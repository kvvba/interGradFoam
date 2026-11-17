// Included in VoF/alphaEqn.H after the high-order alpha1 face flux
// talphaPhi1Un has been computed (used by both the MULESCorr and the
// explicit branch, and limited by MULES afterwards).
//
// Face value for flow leaving a phase change source cell (mDot != 0):
// the alpha1 of the RECEIVING cell instead of the MULES face value.
//
// Why: all of the phase change volume source mDot*(1/rho2 - 1/rho1) sits in
// the cut cell, and the fluid it pushes out through a face is the phase on
// the other side of that face: liquid through the faces towards the liquid,
// vapour through the faces towards the vapour. The vapour produced stays in
// the source cell, so alpha1 of the source cell drops at the rate of vapour
// production and the cell empties completely before the alpha1 = 0.5
// iso-surface (and with it mDot) moves to the next cell.
// With the MULES face value the outflow takes (close to) the upwind value,
// i.e. the alpha1 of the source cell: vapour is pushed into the liquid cell,
// the iso-surface moves on early and the source cell is left with
// alpha1 ~ 0.15 that is never removed (no source, no flux, no compression
// since phic ~ |phi| ~ 0 inside the bubble). Observed on the 7.5 um spherical
// bubble case.
//
// Switch: displacedPhaseFlux in the alpha1 solver dictionary
// (system/<region>/fvSolution, solvers/"alpha.*"), default yes.
//
// Needs extremaCoeff > 0 (0.1 is enough) in the same dictionary: otherwise
// the MULES limiter keeps the source cell above the alpha1 of its vapour-side
// neighbours and the cell stalls at ~0.07 (see readme.org).
if (displacedPhaseFlux)
{
    surfaceScalarField& aPhi = talphaPhi1Un.ref();
    const surfaceScalarField& phiF = phiCN();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const scalarField& mDotI = mDot.primitiveField();
    const scalarField& alpha1I = alpha1.primitiveField();

    scalarField& aPhiI = aPhi.primitiveFieldRef();
    const scalarField& phiI = phiF.primitiveField();

    forAll(own, facei)
    {
        const scalar f = phiI[facei];

        if (f > 0 && mDotI[own[facei]] != 0)
        {
            aPhiI[facei] = f*alpha1I[nei[facei]];
        }
        else if (f < 0 && mDotI[nei[facei]] != 0)
        {
            aPhiI[facei] = f*alpha1I[own[facei]];
        }
    }

    // Coupled (processor, cyclic) patches: neighbour values across the patch
    surfaceScalarField::Boundary& aPhiBf = aPhi.boundaryFieldRef();

    forAll(aPhiBf, patchi)
    {
        fvsPatchScalarField& aPhip = aPhiBf[patchi];

        if (!aPhip.coupled())
        {
            continue;
        }

        const fvsPatchScalarField& phip = phiF.boundaryField()[patchi];
        const labelUList& fc = mesh.boundary()[patchi].faceCells();

        const scalarField alphaNei
        (
            alpha1.boundaryField()[patchi].patchNeighbourField()
        );
        const scalarField mDotNei
        (
            mDot.boundaryField()[patchi].patchNeighbourField()
        );

        forAll(aPhip, i)
        {
            if (phip[i] > 0 && mDotI[fc[i]] != 0)
            {
                aPhip[i] = phip[i]*alphaNei[i];
            }
            else if (phip[i] < 0 && mDotNei[i] != 0)
            {
                aPhip[i] = phip[i]*alpha1I[fc[i]];
            }
        }
    }
}
