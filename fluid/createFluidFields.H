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
// Initialise fluid field pointer lists
//
// Phase 1 (alpha1) is the liquid and phase 2 (alpha2 = 1 - alpha1) is the
// vapour. Every list holds one entry per fluid region.

// Two-phase mixture (phase properties, interface curvature, surface tension)
PtrList<immiscibleIncompressibleTwoPhaseMixture> mixtureFluid
(
 fluidRegions.size()
);

PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
// MULES alpha1 face flux, kept between time steps
PtrList<surfaceScalarField> alphaPhi10Fluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<incompressible::turbulenceModel> turbulenceFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
// Mixture temperature, T = alpha1*TL + alpha2*TV
PtrList<volScalarField> TFluid(fluidRegions.size());
// Liquid and vapour phase temperatures, solved separately in TEqn.H
PtrList<volScalarField> TLFluid(fluidRegions.size());
PtrList<volScalarField> TVFluid(fluidRegions.size());
// Conductive heat flux vector -kappa*grad(T), only created if
// writeQwall is present in transportProperties
PtrList<volVectorField> qwFluid(fluidRegions.size());

// Mixture thermal conductivity
PtrList<volScalarField> kappaFluid(fluidRegions.size());
// Volumetric phase change mass source [kg/m^3/s], sharp: non-zero only in
// interface cells and their direct neighbours (sharpMdot.H)
PtrList<volScalarField> mDotFluid(fluidRegions.size());
// Interfacial heat flux per unit volume [W/m^3], qDot = qL + qV
PtrList<volScalarField> qDotFluid(fluidRegions.size());
// Vapour-side and liquid-side contributions to qDot
PtrList<volScalarField> qVFluid(fluidRegions.size());
PtrList<volScalarField> qLFluid(fluidRegions.size());
// Centre and area vector of the reconstructed interface (alpha1 = 0.5
// iso-surface) in each interface cell, zero elsewhere (calcQdot.H)
PtrList<volVectorField> centreFluid(fluidRegions.size());
PtrList<volVectorField> normalFluid(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());
List<scalar> cumulativeContErrFluid(fluidRegions.size());
PtrList<dimensionedScalar> ghRefFluid(fluidRegions.size());
// List<bool> frozenFlowFluid(fluidRegions.size(),false);
List<bool> alphaRestartFluid(fluidRegions.size(),false);

PtrList<pisoControl> pisoFluid(fluidRegions.size());

// Cached phase properties to avoid repeated dictionary lookups
PtrList<dimensionedScalar> k1Fluid(fluidRegions.size());   // liquid conductivity
PtrList<dimensionedScalar> k2Fluid(fluidRegions.size());   // vapour conductivity
PtrList<dimensionedScalar> cp1Fluid(fluidRegions.size());  // liquid heat capacity
PtrList<dimensionedScalar> cp2Fluid(fluidRegions.size());  // vapour heat capacity
PtrList<dimensionedScalar> TsatFluid(fluidRegions.size()); // saturation temperature
PtrList<dimensionedScalar> hlvFluid(fluidRegions.size());  // latent heat
PtrList<dimensionedScalar> PrtFluid(fluidRegions.size());  // turbulent Prandtl number
List<bool> turbulentHeatFluid(fluidRegions.size(), false); // add turbulent conductivity (non-laminar models only)

// PtrList <surfaceScalarField> brackFSFluid(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

// max(T) at the previous time step, used for the presimulation
// convergence check in computeFluidHeatFluxes.H
scalarField maxTFluid(fluidRegions.size());

// MULES correction flux from the previous time step (alphaApplyPrevCorr)
List<tmp<surfaceScalarField>> talphaPhi1Corr0Fluid(fluidRegions.size());

// Bubble seeding settings and state, read from the optional
// constant/<region>/bubbleSeeding dictionary (see below)
List<label> bubbleSeedCellFluid(fluidRegions.size(), -1);
List<scalar> bubbleActivationTFluid(fluidRegions.size(), GREAT); // GREAT = disabled
List<vector> bubbleCoordinateFluid(fluidRegions.size(), vector::zero);
List<vector> bubbleNormalFluid(fluidRegions.size(), vector(0,0,1));
List<scalar> bubbleContactAngleDegFluid(fluidRegions.size(), 90.0);
List<scalar> bubbleRadiusFluid(fluidRegions.size(), 0.0);
List<bool> bubbleSeedOnceFluid(fluidRegions.size(), true);
List<bool> bubbleSeededFluid(fluidRegions.size(), false);

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
	Info<< "*** Reading fluid mesh thermophysical properties for region "
			<< fluidRegions[i].name() << nl << endl;

	Info<< "    Adding to p_rghFluid\n" << endl;
	p_rghFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "p_rgh",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		 )
    );


	Info<< "    Adding to UFluid\n" << endl;
	UFluid.set
    (
		 i,
		 new volVectorField
		 (
			IOobject
			(
			 "U",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		 )
    );

	Info<< "    Adding to phiFluid\n" << endl;
	phiFluid.set
    (
		 i,
		 new surfaceScalarField
		 (
			IOobject
			(
			 "phi",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::READ_IF_PRESENT,
			 IOobject::AUTO_WRITE
			),
			linearInterpolate(UFluid[i])
			& fluidRegions[i].Sf()
		 )
    );

	Info<< "    Adding to hRefFluid\n" << endl;
	hRefFluid.set
    (
		 i,
		 new uniformDimensionedScalarField
		 (
			IOobject
			(
			 "hRef",
			 runTime.constant(),
			 fluidRegions[i],
			 IOobject::READ_IF_PRESENT,
			 IOobject::NO_WRITE
			),
			dimensionedScalar("hRef", dimLength, Zero) // uses name
		 )
    );

	// mDot is read if present: on a restart the first alpha step needs the
	// mDot the stored phi was computed with (alphaSuSp.H). Otherwise it
	// starts at zero. Its boundary values are not used.
	Info << "    Adding to mDot\n" << endl;
	mDotFluid.set
		(
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "mDot",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::READ_IF_PRESENT,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i],
			dimensionedScalar(dimDensity/dimTime, Zero),
			calculatedFvPatchScalarField::typeName
		 )
		);

	Info << "    Adding to qDot\n" << endl;
	qDotFluid.set
		(i,
		 new volScalarField
		 (
			IOobject
			(
			 "qDot",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i],
			dimensionedScalar
			(
			 "qDot",
			 dimensionSet(1,-1,-3,0,0,0,0),
			 0.0
			)
		 )
		);

		qVFluid.set
		(i,
		 new volScalarField
		 (
			IOobject
			(
			 "qV",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i],
			dimensionedScalar
			(
			 "qV",
			 dimensionSet(1,-1,-3,0,0,0,0),
			 0.0
			)
		 )
		);
	qLFluid.set
		(i,
		 new volScalarField
		 (
			IOobject
			(
			 "qL",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i],
			dimensionedScalar
			(
			 "qL",
			 dimensionSet(1,-1,-3,0,0,0,0),
			 0.0
			)
		 )
		);

	centreFluid.set
		(
		 i,
		 new volVectorField
		 (
			IOobject
			(
			 "centre",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i],
			dimensionedVector(dimLength, Zero)
		 )
		);

	normalFluid.set
		(
		 i,
		 new volVectorField
		 (
			IOobject
			(
			 "normal",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i],
			dimensionedVector(dimArea, Zero)
		 )
		);

	dimensionedScalar ghRef
		(
		 mag(g.value()) > SMALL
		 ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
		 : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
		);

	ghRefFluid.set(i, new dimensionedScalar(ghRef));

	Info<< "    Adding to ghFluid\n" << endl;
	ghFluid.set
    (
		 i,
		 new volScalarField
		 (
			"gh",
			(g & fluidRegions[i].C()) - ghRef
		 )
    );

	Info<< "    Adding to ghfFluid\n" << endl;
	ghfFluid.set
    (
		 i,
		 new surfaceScalarField
		 (
			"ghf",
			(g & fluidRegions[i].Cf()) - ghRef
		 )
    );

	mixtureFluid.set
    (
		 i,
		 new immiscibleIncompressibleTwoPhaseMixture
		 (
			UFluid[i],
			phiFluid[i]
		 )
    );


	Info<< "    Adding to turbulenceFluid\n" << endl;
	turbulenceFluid.set
    (
		 i,
		 incompressible::turbulenceModel::New
		 (
			UFluid[i],
			phiFluid[i],
			mixtureFluid[i]
		 ).ptr()
    );


	immiscibleIncompressibleTwoPhaseMixture& mixture(mixtureFluid[i]);
	volScalarField& alpha1(mixture.alpha1());
	volScalarField& alpha2(mixture.alpha2());


	//- This also builds the surface reconstruction scheme

	const dimensionedScalar& rho1 = mixture.rho1();
	const dimensionedScalar& rho2 = mixture.rho2();

	// Thermal properties come from constant/transportProperties: kappa and
	// Cp from each phase sub-dictionary, and Tsat and hlv from the
	// PhaseChangeProperties sub-dictionary (a diff entry, used by
	// mulesGradFoam's Hardt-Wondra smoothing, is ignored)
	dictionary phase1 = mixture.subDict(mixture.phase1Name());
	dictionary phase2 = mixture.subDict(mixture.phase2Name());
	dictionary pcProperties = mixture.subDict("PhaseChangeProperties");

	// Cache thermal and phase-change properties to avoid repeated lookups during time steps
	dimensionedScalar k1
    (
		 "kappa",
		 dimensionSet(1,1,-3,-1,0,0,0),
		 readScalar(phase1.lookup("kappa"))
    );

	dimensionedScalar k2
    (
		 "kappa",
		 dimensionSet(1,1,-3,-1,0,0,0),
		 readScalar(phase2.lookup("kappa"))
    );

	dimensionedScalar cp1
    (
		 "Cp",
		 dimensionSet(0,2,-2,-1,0,0,0),
		 readScalar(phase1.lookup("Cp"))
    );

	dimensionedScalar cp2
    (
		 "Cp",
		 dimensionSet(0,2,-2,-1,0,0,0),
		 readScalar(phase2.lookup("Cp"))
    );

	dimensionedScalar Tsat
    (
		 "Tsat",
		 dimTemperature,
		 readScalar(pcProperties.lookup("Tsat"))
    );

	dimensionedScalar hlv
    (
		 "hlv",
		 dimensionSet(0,2,-2,0,0,0,0),
		 readScalar(pcProperties.lookup("hlv"))
    );

	// Store cached values for use in setRegionFluidFields.H
	k1Fluid.set(i, new dimensionedScalar(k1));
	k2Fluid.set(i, new dimensionedScalar(k2));
	cp1Fluid.set(i, new dimensionedScalar(cp1));
	cp2Fluid.set(i, new dimensionedScalar(cp2));
	TsatFluid.set(i, new dimensionedScalar(Tsat));
	hlvFluid.set(i, new dimensionedScalar(hlv));

	// Turbulent heat transfer: with a RAS/LES model the energy equations use
	// kappa + rho*Cp*nut/Prt. Laminar cases skip this entirely, so they are
	// unaffected. Prt is an optional top-level entry in transportProperties.
	turbulentHeatFluid[i] =
		!isA<incompressible::laminarModel>(turbulenceFluid[i]);

	PrtFluid.set
	(
		i,
		new dimensionedScalar
		(
			"Prt",
			dimless,
			mixture.getOrDefault<scalar>("Prt", 0.85)
		)
	);

	if (turbulentHeatFluid[i])
	{
		Info<< "    Turbulent heat transfer on, Prt = "
			<< PrtFluid[i].value() << nl << endl;
	}

	Info<< "    Adding to kappaFluid\n" << endl;
	kappaFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "kappa",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			 ),
			k1*alpha1 + k2*alpha2
			)
		 );


	Info<< "    Adding to TFluid\n" << endl;
	TFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "T",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		 )
    );

	TLFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "TL",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		 )
    );
	TVFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "TV",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			fluidRegions[i]
		 )
    );

		
	maxTFluid[i] = gMax(TFluid[i]);

	rhoFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "rho",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::NO_READ,
			 IOobject::AUTO_WRITE
			),
			alpha1*rho1 + alpha2*rho2
		 )
    );

	pisoFluid.set
    (
		 i,
		 new pisoControl(fluidRegions[i])
    );

	// pimpleFluid[i].dict().readIfPresent("frozenFlow", frozenFlowFluid[i]);

	pFluid.set
    (
		 i,
		 new volScalarField
		 (
			IOobject
			(
			 "p",
			 runTime.timeName(),
			 fluidRegions[i],
			 IOobject::MUST_READ,
			 IOobject::AUTO_WRITE
			),
			p_rghFluid[i] + rhoFluid[i]*ghFluid[i]
		 )
    );

	pRefCellFluid[i] = -1;
	pRefValueFluid[i] = 0.0;
	setRefCell
    (
		 pFluid[i],
		 p_rghFluid[i],
		 pisoFluid[i].dict(),
		 pRefCellFluid[i],
		 pRefValueFluid[i]
    );



	if (p_rghFluid[i].needReference())
    {
			pFluid[i] += dimensionedScalar
        (
				 "p",
				 pFluid[i].dimensions(),
				 pRefValueFluid[i] - getRefCellValue(pFluid[i], pRefCellFluid[i])
        );
			p_rghFluid[i] = pFluid[i] - rhoFluid[i]*ghFluid[i];
    }

	fluidRegions[i].setFluxRequired(p_rghFluid[i].name());
	fluidRegions[i].setFluxRequired(alpha1.name());


	Info<< "    Adding fvOptions\n" << endl;
	fluidFvOptions.set
    (
		 i,
		 new fv::options(fluidRegions[i])
    );

	turbulenceFluid[i].validate();

	Info<< "    Adding MRF\n" << endl;
	MRFfluid.set
    (
		 i,
		 new IOMRFZoneList(fluidRegions[i])
    );

	// MULES alpha flux from the previous time step. If alphaPhi0 is in the
	// start time directory this is a restart, and the Crank-Nicolson
	// off-centring in alphaEqn.H is applied from the first step.
	IOobject alphaPhi10Header
    (
		 IOobject::groupName("alphaPhi0", alpha1.group()),
		 runTime.timeName(),
		 fluidRegions[i],
		 IOobject::READ_IF_PRESENT,
		 IOobject::AUTO_WRITE
    );

	alphaRestartFluid[i] =
		alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true);

	if (alphaRestartFluid[i])
    {
			Info << "Restarting alpha" << endl;
    }

	alphaPhi10Fluid.set
    (
		 i,
		 new surfaceScalarField
		 (
			alphaPhi10Header,
			phiFluid[i]*fvc::interpolate(alpha1)
		 )
    );

	// brackFSFluid.set
	//   (
	//    i,
	//    new surfaceScalarField
	//    (
	//     fvc::interpolate(mixture.sigmaK())*fvc::snGrad(alpha1)
	//     )
	//    );

	uniformDimensionedScalarField cumulativeContErrIO
    (
		 IOobject
		 (
			"cumulativeContErr",
			runTime.timeName(),
			"uniform",
			fluidRegions[i],
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		 ),
		 dimensionedScalar(dimless, Zero)
    );


	cumulativeContErrFluid[i] = cumulativeContErrIO.value();


	// Optional output of the conductive heat flux vector; it is only updated
	// in computeFluidHeatFluxes.H, i.e. when presimControl/computeFluxes is on
	if(mixture.found("writeQwall"))
    {
			qwFluid.set
        (
				 i,
				 new volVectorField
				 (
					IOobject
					(
					 "qw",
					 runTime.timeName(),
					 fluidRegions[i],
					 IOobject::NO_READ,
					 IOobject::AUTO_WRITE
					),
					-kappaFluid[i]*fvc::grad(TFluid[i])
				 )
        );
    }

	// --- Bubble seeding: read the optional `constant/bubbleSeeding`
	//     dictionary once. The bubble itself is inserted in
	//     setRegionFluidFields.H when the seed cell reaches
	//     `activationTemperature`.
	//
	// Dictionary (optional, place in `constant/<region>/bubbleSeeding`):
	// {
	//     activationTemperature 373.15;    // scalar
	//     coordinate (0.01 0.02 0.0);      // vector (contact point on surface)
	//     surfaceNormal (0 0 1);           // vector (points into fluid)
	//     contactAngle 60;                 // degrees, measured through the liquid
	//     initialRadius 1e-3;              // meters
	//     seedOnce true;                   // optional (default true)
	// }
	{
		IOobject bsHeader
		(
			"bubbleSeeding",
			runTime.constant(),
			fluidRegions[i],
			IOobject::READ_IF_PRESENT,
			IOobject::NO_WRITE
		);

		// Use IOobject header check for presence, then construct dictionary from file
		if (bsHeader.typeHeaderOk<IOdictionary>(true))
		{
			IOdictionary bsDict(bsHeader);

			bubbleActivationTFluid[i] = bsDict.getOrDefault<scalar>("activationTemperature", GREAT);
			bubbleCoordinateFluid[i] = bsDict.getOrDefault<vector>("coordinate", vector::zero);
			bubbleNormalFluid[i] = bsDict.getOrDefault<vector>("surfaceNormal", vector(0,0,1));
			bubbleContactAngleDegFluid[i] = bsDict.getOrDefault<scalar>("contactAngle", 90.0);
			bubbleRadiusFluid[i] = bsDict.getOrDefault<scalar>("initialRadius", 0.0);
			bubbleSeedOnceFluid[i] = bsDict.getOrDefault<bool>("seedOnce", true);

			if (bubbleRadiusFluid[i] > 0 && bubbleActivationTFluid[i] < GREAT)
			{
				// Use meshSearch for efficient spatial lookup of the seed cell.
				// In parallel, only the owning processor returns a valid index.
				meshSearch bubbleMeshSearch(fluidRegions[i]);
				bubbleSeedCellFluid[i] =
					bubbleMeshSearch.findCell(bubbleCoordinateFluid[i]);

				// Check across all ranks whether ANY processor found the cell
				label globalFound = returnReduce
				(
					(bubbleSeedCellFluid[i] >= 0) ? label(1) : label(0),
					maxOp<label>()
				);

				if (globalFound == 0)
				{
					Info<< "bubbleSeeding: could not find a mesh cell for coordinate "
						<< bubbleCoordinateFluid[i]
						<< " on any processor, disabling seeding." << endl;
					bubbleActivationTFluid[i] = GREAT;
				}
				else
				{
					Info<< "bubbleSeeding: configured seed"
						<< " (coord " << bubbleCoordinateFluid[i]
						<< ") with R=" << bubbleRadiusFluid[i]
						<< " contactAngle=" << bubbleContactAngleDegFluid[i] << " deg"
						<< endl;
				}
			}
			else
			{
				Info<< "bubbleSeeding: dictionary present but incomplete (must provide activationTemperature and initialRadius). Ignoring." << endl;
				bubbleActivationTFluid[i] = GREAT;
			}
		}
		else
		{
			// no bubbleSeeding dictionary — that's fine
			Info << "No bubbleSeeding dictionary found" << endl;
		}
	}

}
