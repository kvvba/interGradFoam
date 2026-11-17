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
// Energy equations for the liquid temperature TL and the vapour temperature
// TV. Each phase temperature is solved on the whole mesh, but:
//  - conduction is switched off across faces touching the other phase
//  - the interface acts as a Dirichlet boundary at Tsat, through the
//    implicit sources TsourceL/TsourceV from calcQdot.H
//  - cells of the other phase are fixed to Tsat
// The mixture temperature T is then rebuilt from TL and TV.

// Face conductivities for each phase
surfaceScalarField kappaLf(fvc::interpolate(kappaL));
surfaceScalarField kappaVf(fvc::interpolate(kappaV));

// Turbulent conductivity rho_k*Cp_k*nut/Prt (RAS/LES only), added before
// the conduction switch-off below so it is zeroed on the same faces
if (turbulentHeat)
{
	const surfaceScalarField nutf(fvc::interpolate(turbulence.nut()));
	kappaLf += rho1*cp1/Prt*nutf;
	kappaVf += rho2*cp2/Prt*nutf;
}

const labelList& owner = mesh.owner();
const labelList& neighbour = mesh.neighbour();

// No liquid conduction across faces next to a vapour cell, and no vapour
// conduction across faces next to a liquid cell. Heat transfer between a
// phase and the interface goes only through the interface sources.
// Internal faces first:
forAll(kappaLf, facei)
{
	const label& own = owner[facei];
	const label& nei = neighbour[facei];
	const scalar vapourOwn = vapourCells[own];
	const scalar vapourNei = vapourCells[nei];
	const scalar liquidOwn = liquidCells[own];
	const scalar liquidNei = liquidCells[nei];

	if (vapourOwn || vapourNei)
		// cell or neighbour cell is on the vapour side
		{
			kappaLf[facei] = 0.0;
		}
	if (liquidOwn || liquidNei)
		// cell or neighbour cell is on the liquid side
		{
			kappaVf[facei] = 0.0;
		}
}

// Then coupled (processor, cyclic) boundary faces, where the cell on the
// other side is found through the patch neighbour values of indicator
// fields
volScalarField vapourInd
(
 IOobject
 (
	"vapourInd",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
 ),
 mesh,
 dimensionedScalar(dimless, Zero)
);
volScalarField liquidInd("liquidInd", vapourInd);

forAll(vapourInd, celli)
{
	vapourInd[celli] = vapourCells[celli];
	liquidInd[celli] = liquidCells[celli];
}
vapourInd.correctBoundaryConditions();
liquidInd.correctBoundaryConditions();

forAll(kappaLf.boundaryField(), patchi)
{
	if (vapourInd.boundaryField()[patchi].coupled())
	{
		const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
		const scalarField vapourNei
		(
			vapourInd.boundaryField()[patchi].patchNeighbourField()
		);
		const scalarField liquidNei
		(
			liquidInd.boundaryField()[patchi].patchNeighbourField()
		);
		fvsPatchScalarField& kappaLp = kappaLf.boundaryFieldRef()[patchi];
		fvsPatchScalarField& kappaVp = kappaVf.boundaryFieldRef()[patchi];

		forAll(faceCells, facei)
		{
			if (vapourCells[faceCells[facei]] || vapourNei[facei] > 0.5)
			{
				kappaLp[facei] = 0.0;
			}
			if (liquidCells[faceCells[facei]] || liquidNei[facei] > 0.5)
			{
				kappaVp[facei] = 0.0;
			}
		}
	}
}

volScalarField rhoCpL
(
 IOobject
 (
	"rhoCpL",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
 ),
 mesh,
 rho1 * cp1
);

volScalarField rhoCpV
(
 IOobject
 (
	"rhoCpV",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
 ),
 mesh,
 rho2 * cp2
);

// Phase heat capacity per unit mixture volume, alpha_k*rho_k*cp_k. Strictly
// this is zero in cells of the other phase; there it is set to 1e-7 times
// the pure-phase value instead, to keep it small but nonzero because a
// zero value is numerically unstable. Those cells are fixed to Tsat by
// setValues below anyway.
forAll(rhoCpL, celli)
{
	if (interfaceCells[celli])
		{
			rhoCpL[celli] = alpha1[celli] * rho1.value() * cp1.value();
			rhoCpV[celli] = alpha2[celli] * rho2.value() * cp2.value();
		}
	else if (vapourCells[celli])
		{
			rhoCpL[celli] = rho1.value() * cp1.value() * 1.0e-7;
		}
	else if (liquidCells[celli])
		{
			rhoCpV[celli] = rho2.value() * cp2.value() * 1.0e-7;
		}
}

surfaceScalarField rhoPhiCpL
("rhoPhiCpL",
 fvc::interpolate(rhoCpL)*phi);

surfaceScalarField rhoPhiCpV
("rhoPhiCpV",
 fvc::interpolate(rhoCpV)*phi);

// Phase energy equations in non-conservative form: the Sp term removes
// T*(ddt(rhoCp) + div(rhoPhiCp)) from the conservative form, leaving
// rhoCp*DT/Dt = laplacian(kappa, T) - interface heat sink
fvScalarMatrix TLEqn
(
 fvm::ddt(rhoCpL,TL)
 + fvm::div(rhoPhiCpL,TL)
 - fvm::Sp(fvc::ddt(rhoCpL) + fvc::div(rhoPhiCpL), TL)
 - fvm::laplacian(kappaLf,TL)
 + TsourceL
);

fvScalarMatrix TVEqn
(
 fvm::ddt(rhoCpV,TV)
 + fvm::div(rhoPhiCpV,TV)
 - fvm::Sp(fvc::ddt(rhoCpV) + fvc::div(rhoPhiCpV), TV)
 - fvm::laplacian(kappaVf,TV)
 + TsourceV
);

// Fix each phase temperature to Tsat in cells of the other phase
TLEqn.setValues(VapourCellsLabels, TsatField);
TVEqn.setValues(LiquidCellsLabels, TsatField);

TL.correctBoundaryConditions();
TV.correctBoundaryConditions();

TLEqn.relax();
TVEqn.relax();

TLEqn.solve();
TVEqn.solve();

// Mixture temperature, used for output, the solid coupling boundary
// conditions, bubble seeding and computeFluidHeatFluxes.H
T = TL*alpha1 + TV*alpha2;

T.correctBoundaryConditions();

mixture.correct();

Info << "Min(TL) = " << min(TL).value() << " Max(TL) = " << max(TL).value() << endl;
Info << "Min(TV) = " << min(TV).value() << " Max(TV) = " << max(TV).value() << endl;

