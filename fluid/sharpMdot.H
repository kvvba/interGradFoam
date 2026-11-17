/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Jakub Cranmer <jakub@posteo.net>
-------------------------------------------------------------------------------
License
    This file is part of interGradFoam.

    interGradFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    interGradFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with interGradFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// Sharp phase change mass source, as in PSI-BOIL (Sato & Niceno 2013,
// J. Comput. Phys. 249; Lal, Sato & Niceno 2015, Int. J. Heat Fluid Flow 51,
// Eq. 3): mDot = qDot/hlv [kg/m^3/s], which is non-zero only in the cells
// cut by the alpha1 = 0.5 iso-surface (calcQdot.H). There is no smoothing and
// no shifting of the source out of the interface region (this replaces the
// Hardt & Wondra treatment of mulesGradFoam's smoothQdot.H).
//
// The same mDot is used in
//   - pEqn.H:      div(U) = mDot*(1/rho2 - 1/rho1)
//   - alphaSuSp.H: Su = -mDot/rho1, Sp = -mDot*(1/rho2 - 1/rho1)
// (at the start of the next time step).
//
// Cap and redistribution (PSI-BOIL PhaseChange4::mdot_cut / mdot):
// mDot comes from interface heat flux x interface area and does not know how
// much liquid the cell holds. Within one time step a cell may lose at most its
// liquid (evaporation) or fill at most its vapour volume (condensation):
//     mDot > 0:  mDot_c = min(mDot,  rho1*alpha1/dt)
//     mDot < 0:  mDot_c = max(mDot, -rho1*(1 - alpha1)/dt)
// The excess (mDot - mDot_c)*V [kg/s] is passed to face neighbours: towards
// the liquid for evaporation, towards the vapour for condensation. Face f of
// the donor gets weight (d & nf)^2 where d is that direction and nf the
// outward face normal, over faces with (d & nf) > 0 (on a hex mesh these are
// PSI-BOIL's nx^2, ny^2, nz^2). The receivers are capped again and anything
// left over is dropped (and logged).
//
// Time step: this mDot is used in this step's pressure equation and in the
// alpha step at the start of the NEXT time step, which uses the next deltaT.
// With adjustTimeStep that can be up to maxDeltaTGrowth times the current one
// (setMultiRegionDeltaT.H), so the cap uses dtCap = maxDeltaTGrowth*deltaT.
// alphaSuSp.H reports any cell that still exceeds the cap at the actual dt.
{
    const scalarField& V = mesh.V();
    const scalarField& alpha1I = alpha1.primitiveField();
    const scalar rhoL = rho1.value();

    const scalar dtCap =
        runTime.deltaTValue()*(adjustTimeStep ? maxDeltaTGrowth : scalar(1));

    // Limit mDot in cell celli so the source alone keeps alpha1 in [0, 1]
    auto capMdot = [&](const label celli, const scalar m) -> scalar
    {
        const scalar a = min(max(alpha1I[celli], scalar(0)), scalar(1));

        if (m >= 0)
        {
            return min(m, rhoL*a/dtCap);
        }

        return max(m, -rhoL*(1 - a)/dtCap);
    };

    scalarField& mDotI = mDot.primitiveFieldRef();

    // Sharp source: qDot is already zero outside interface cells
    mDotI = qDot.primitiveField()/hlv.value();

    // Uncapped evaporation rate [kg/s]
    const scalar mRaw = gSum(mDotI*V);

    // --- First cap. excess [kg/m^3/s] is what the donor cell could not take.
    scalarField excess(mDotI.size(), Zero);
    label nCapped = 0;

    forAll(mDotI, celli)
    {
        if (mDotI[celli] != 0)
        {
            const scalar mc = capMdot(celli, mDotI[celli]);
            excess[celli] = mDotI[celli] - mc;
            mDotI[celli] = mc;

            if (excess[celli] != 0)
            {
                ++nCapped;
            }
        }
    }

    const scalar mExcess = gSum(excess*V);   // [kg/s], signed
    scalar mNoTarget = 0;  // excess with no eligible neighbour face [kg/s]
    scalar mDropped = 0;   // removed by the second cap [kg/s]
    label nFlipped = 0;    // donors whose iso-surface normal pointed to vapour
    label nRecapped = 0;   // receivers limited by the second cap

    // fvc::grad and the processor exchange are collective, so all processors
    // take this branch if any of them has a capped cell
    if (returnReduce(nCapped, sumOp<label>()) > 0)
    {
        // Used only to orient the iso-surface area vector towards the liquid
        const volVectorField gradAlpha1(fvc::grad(alpha1));

        const vectorField& faceAreas = mesh.faceAreas();
        const labelList& faceOwner = mesh.faceOwner();
        const labelList& faceNeighbour = mesh.faceNeighbour();
        const cellList& cells = mesh.cells();
        const label nInternalFaces = mesh.nInternalFaces();
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        // Coupled patches the excess may cross. Other boundary faces (walls,
        // inlets, outlets, symmetry, empty, wedge, AMI) never receive.
        boolList transferPatch(pbm.size(), false);
        forAll(pbm, patchi)
        {
            transferPatch[patchi] =
                isA<processorPolyPatch>(pbm[patchi])
             || isA<cyclicPolyPatch>(pbm[patchi]);
        }

        // Mass rate [kg/s] received by each local cell
        scalarField received(mDotI.size(), Zero);

        // Mass rate [kg/s] sent out through each coupled boundary face
        scalarField sendBf(mesh.nBoundaryFaces(), Zero);

        // Scratch space for the weights of one donor
        DynamicList<label> donorFaces(12);
        DynamicList<scalar> donorWeights(12);

        forAll(excess, celli)
        {
            if (excess[celli] == 0)
            {
                continue;
            }

            const scalar rate = excess[celli]*V[celli];  // [kg/s]

            if (mag(normal[celli]) < VSMALL)
            {
                mNoTarget += rate;
                continue;
            }

            // Unit normal pointing into the liquid
            vector nL = normal[celli]/mag(normal[celli]);
            if ((nL & gradAlpha1[celli]) < 0)
            {
                nL = -nL;
                ++nFlipped;
            }

            // Evaporation excess goes towards the liquid, condensation
            // excess towards the vapour
            const vector d = (excess[celli] > 0 ? nL : -nL);

            donorFaces.clear();
            donorWeights.clear();
            scalar sumW = 0;

            for (const label facei : cells[celli])
            {
                if (facei >= nInternalFaces)
                {
                    const label patchi = pbm.whichPatch(facei);
                    if (!transferPatch[patchi])
                    {
                        continue;
                    }
                }

                // Outward unit normal of facei seen from celli
                vector nf = faceAreas[facei]/mag(faceAreas[facei]);
                if (faceOwner[facei] != celli)
                {
                    nf = -nf;
                }

                const scalar c = d & nf;

                if (c > 0)
                {
                    donorFaces.append(facei);
                    donorWeights.append(sqr(c));
                    sumW += sqr(c);
                }
            }

            if (sumW < VSMALL)
            {
                mNoTarget += rate;
                continue;
            }

            forAll(donorFaces, k)
            {
                const label facei = donorFaces[k];
                const scalar share = rate*donorWeights[k]/sumW;

                if (facei < nInternalFaces)
                {
                    const label own = faceOwner[facei];
                    const label nei = faceNeighbour[facei];
                    received[own == celli ? nei : own] += share;
                }
                else
                {
                    sendBf[facei - nInternalFaces] += share;
                }
            }
        }

        // Swap the coupled boundary values: each coupled face now holds what
        // the cell on the other side sent through it
        syncTools::swapBoundaryFaceList(mesh, sendBf);

        forAll(sendBf, bFacei)
        {
            if (sendBf[bFacei] != 0)
            {
                const label facei = bFacei + nInternalFaces;
                if (transferPatch[pbm.whichPatch(facei)])
                {
                    received[faceOwner[facei]] += sendBf[bFacei];
                }
            }
        }

        // --- Second cap; the remainder is dropped
        forAll(received, celli)
        {
            if (received[celli] != 0)
            {
                const scalar mNew = mDotI[celli] + received[celli]/V[celli];
                const scalar mc = capMdot(celli, mNew);

                if (mc != mNew)
                {
                    mDropped += (mNew - mc)*V[celli];
                    ++nRecapped;
                }

                mDotI[celli] = mc;
            }
        }
    }

    mDot.correctBoundaryConditions();

    reduce(mNoTarget, sumOp<scalar>());
    reduce(mDropped, sumOp<scalar>());
    reduce(nCapped, sumOp<label>());
    reduce(nFlipped, sumOp<label>());
    reduce(nRecapped, sumOp<label>());

    const scalar mFinal = gSum(mDotI*V);

    // Same line as mulesGradFoam, for existing post-processing scripts
    Info<< "Timeev= " << runTime.timeName() << " "
        << "mEvap= " << mRaw << " " << nl << endl;

    // All mass rates in kg/s (signed, > 0 = evaporation).
    // raw = sum(qDot/hlv*V); final = sum(mDot*V) = raw - noTarget - dropped
    Info<< "sharpMdot: raw = " << mRaw
        << " final = " << mFinal
        << " capped = " << mExcess << " (" << nCapped << " cells)"
        << " redistributed = " << mExcess - mNoTarget
        << " noTarget = " << mNoTarget
        << " dropped = " << mDropped << " (" << nRecapped << " cells)"
        << " flippedNormals = " << nFlipped
        << " dtCap = " << dtCap
        << endl;
}
