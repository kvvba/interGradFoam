/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
		Copyright (C) 2023 Jakub Cranmer <jakub@posteo.net>
-------------------------------------------------------------------------------
License
    This file is part of mulesGradFoam.

    mulesGradFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mulesGradFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mulesGradFoam.  If not, see <http://www.gnu.org/licenses/>.
		
\*---------------------------------------------------------------------------*/
// Included at the start of every fluid region solve. It binds local
// references to the stored fields of fluid region i, so that the
// interFoam-style equation files can use the usual names (U, p_rgh, alpha1,
// ...). It also builds the temporary fields needed for this time step and
// handles the optional bubble seeding.

dynamicFvMesh& mesh = fluidRegions[i];

immiscibleIncompressibleTwoPhaseMixture& mixture(mixtureFluid[i]);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

volScalarField& rho = rhoFluid[i];
volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];
surfaceScalarField& alphaPhi10 = alphaPhi10Fluid[i];

incompressible::turbulenceModel& turbulence = turbulenceFluid[i];
volScalarField& p_rgh = p_rghFluid[i];
volScalarField& p = pFluid[i];

volScalarField& mDot = mDotFluid[i];
volScalarField& qDot = qDotFluid[i];
volScalarField& qL = qLFluid[i];
volScalarField& qV = qVFluid[i];
volVectorField& centre = centreFluid[i];
volVectorField& normal = normalFluid[i];

volScalarField& T = TFluid[i];
volScalarField& TL = TLFluid[i];
volScalarField& TV = TVFluid[i];
volScalarField& kappa = kappaFluid[i];

const bool& alphaRestart = alphaRestartFluid[i];

scalar& cumulativeContErr = cumulativeContErrFluid[i];

volScalarField& gh = ghFluid[i];
surfaceScalarField& ghf = ghfFluid[i];
// const dimensionedScalar& ghRef = ghRefFluid[i];

// surfaceScalarField& brackFS = brackFSFluid[i];

fv::options& fvOptions = fluidFvOptions[i];

// bool frozenFlow = frozenFlowFluid[i];

const label pRefCell = pRefCellFluid[i];
const scalar pRefValue = pRefValueFluid[i];

// Use cached phase properties instead of repeated dictionary lookups
const dimensionedScalar& k1 = k1Fluid[i];
const dimensionedScalar& k2 = k2Fluid[i];
const dimensionedScalar& cp1 = cp1Fluid[i];
const dimensionedScalar& cp2 = cp2Fluid[i];
const dimensionedScalar& Tsat = TsatFluid[i];
const dimensionedScalar& hlv = hlvFluid[i];
const dimensionedScalar& Prt = PrtFluid[i];
const bool turbulentHeat = turbulentHeatFluid[i];

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// Mixture volumetric heat capacity. In calcQdot.H it is overwritten
// cell by cell with the pure-phase values away from the interface; it is
// then only used for the convective boundary fluxes in
// computeFluidHeatFluxes.H.
volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1*cp1 + alpha2*rho2*cp2,
    alpha1.boundaryField().types()
);
rhoCp.oldTime();


//Energy flux

// surfaceScalarField rhoPhiCp
// (
//     IOobject
//     (
//         "rhoPhiCp",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     fvc::interpolate(rhoCp)*phi
// );

surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

pisoControl& piso(pisoFluid[i]);

bool correctPhi
(
    piso.dict().getOrDefault("correctPhi", mesh.dynamic())
);

// bool checkMeshCourantNo
// (
//     piso.dict().getOrDefault("checkMeshCourantNo", false)
// );

// bool moveMeshOuterCorrectors
// (
//     piso.dict().getOrDefault("moveMeshOuterCorrectors", false)
// );

IOMRFZoneList& MRF = MRFfluid[i];

// surfaceScalarField& alphaPhi0(alphaPhi10Fluid[i]);

// MULES Correction, kept per region between time steps so that
// alphaApplyPrevCorr in alphaEqn.H can reuse the previous step's correction
tmp<surfaceScalarField>& talphaPhi1Corr0 = talphaPhi1Corr0Fluid[i];

autoPtr<surfaceVectorField> Uf;

if (fluidRegions[i].dynamic())
{
    Info<< "Constructing face velocity Uf\n" << endl;

    Uf.reset
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(UFluid[i])
        )
    );
}

volScalarField rAU
(
    IOobject
    (
        "rAU",
        runTime.timeName(),
        fluidRegions[i],
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidRegions[i],
    dimensionedScalar("rAU", dimTime/dimDensity, 1)
);

// Uniform fields used by the interface heat flux (calcQdot.H) and by the
// phase temperature equations (TEqn.H):
//  - TsatField: the interface temperature (Dirichlet value at the
//    interface) and the value TL/TV are fixed to in the other phase
//  - kappaL/kappaV: liquid and vapour conductivities
const volScalarField TsatField
(
 IOobject
 (
	"TsatField",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
 ),
 mesh,
 Tsat
);

const volScalarField kappaL
(
 IOobject
 (
	"kappa1Field",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
 ),
 mesh,
 k1
);

const volScalarField kappaV
(
 IOobject
 (
	"kappa2Field",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
 ),
 mesh,
 k2
);

// scalar isoFaceTol = 1.0e-8;
// Interface area density: area of the reconstructed interface in a cell
// divided by the cell volume [1/m]. Set in calcQdot.H, zero elsewhere.
volScalarField SArea
(
 IOobject
 (
	"SArea",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
 mesh,
 dimensionedScalar(dimless/dimLength, Zero)
 );
// --- Bubble seeding: the optional `constant/bubbleSeeding` dictionary is
//     read once at start-up (createFluidFields.H). A spherical vapour region
//     is seeded when the temperature at the seed cell reaches
//     `activationTemperature`. The seeding state is stored per region so
//     that it persists between time steps.
const label& _bubbleSeedCell = bubbleSeedCellFluid[i];
const scalar& _bubbleActivationT = bubbleActivationTFluid[i];
const vector& _bubbleCoordinate = bubbleCoordinateFluid[i];
const vector& _bubbleNormal = bubbleNormalFluid[i];
const scalar& _bubbleContactAngleDeg = bubbleContactAngleDegFluid[i];
const scalar& _bubbleRadius = bubbleRadiusFluid[i];
const bool& _bubbleSeedOnce = bubbleSeedOnceFluid[i];
bool& _bubbleSeeded = bubbleSeededFluid[i];

// Check trigger each timestep (only if configured and not yet seeded or seedOnce==false)
if (_bubbleActivationT < GREAT && (!_bubbleSeeded || !_bubbleSeedOnce))
{
    // Sample temperature at the cached/nearest cell (only on the owning
    // processor).  Then use returnReduce so that ALL processors agree on
    // whether to seed — correctBoundaryConditions() is collective MPI.
    scalar Tcell = -GREAT;  // default: not reached
    if (_bubbleSeedCell >= 0)
    {
        Tcell = T.internalField()[_bubbleSeedCell];
    }

    // Collective decision: take the maximum across all ranks
    scalar TcellGlobal = returnReduce(Tcell, maxOp<scalar>());

    if (TcellGlobal >= _bubbleActivationT)
    {
        Info << "bubbleSeeding: activation temperature reached (T=" << TcellGlobal
             << ") — inserting bubble" << endl;

        // vapour is always represented by alpha1 == 0 (user-specified)
        // compute spherical centre so the interface meets the wall at the
        // supplied contact angle (measured through the liquid):
        // centre = contactPoint + n * (R * cos(theta)). theta = 90 deg gives
        // a hemisphere, theta < 90 deg a cap larger than a hemisphere.
        scalar theta = _bubbleContactAngleDeg * Foam::acos(scalar(-1.0)) / scalar(180.0);
        vector nHat = (_bubbleNormal == vector::zero) ? vector(0,0,1) : _bubbleNormal / mag(_bubbleNormal);
        vector sphereCentre = _bubbleCoordinate + nHat * (_bubbleRadius * Foam::cos(theta));

        // set alpha fields and vapour temperature inside sphere
        // — runs on ALL processors so each modifies its own local cells
        const pointField& Cc = mesh.C();
        forAll(Cc, ci)
        {
            if (mag(Cc[ci] - sphereCentre) <= _bubbleRadius)
            {
                // alpha1 == 0 => vapour, alpha2 == 1 => vapour fraction
                alpha1[ci] = 0.0;
                alpha2[ci] = 1.0;

                // set vapour temperature to saturation for seeded vapour cells
                TV[ci] = TsatField[ci];

                // keep mixture/diagnostic temperature consistent
                T[ci] = TV[ci];
            }
        }

        // mark seeded to avoid repeated insertion (unless seedOnce==false)
        _bubbleSeeded = true;

        // update boundary fields — collective MPI, now called by ALL ranks
        alpha1.correctBoundaryConditions();
        alpha2.correctBoundaryConditions();
        TV.correctBoundaryConditions();
        T.correctBoundaryConditions();
        Info << "bubbleSeeding: bubble inserted (centre=" << sphereCentre
             << ", R=" << _bubbleRadius << ")" << endl;
    }
}
