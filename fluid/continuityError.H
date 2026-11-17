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
// Mixture mass continuity of the alpha step: ddt(rho) + div(rhoPhi) with
// rhoPhi from alphaEqn.H. With Su/Sp from alphaSuSp.H consistent with the
// pressure equation this should be ~0 also with phase change.
volScalarField contErr(fvc::ddt(rho) + fvc::div(rhoPhi));

scalar sumLocalContErr = runTime.deltaTValue()*
    mag(contErr)().weightedAverage(mesh.V()).value();

scalar globalContErr = runTime.deltaTValue()*
    contErr.weightedAverage(mesh.V()).value();

cumulativeContErr += globalContErr;

Info<< "time step continuity errors : sum local = " << sumLocalContErr
    << ", global = " << globalContErr
    << ", cumulative = " << cumulativeContErr
    << endl;

// Volume continuity of the pressure equation: div(phi) should equal the
// phase change volume source mDot*(1/rho2 - 1/rho1) (pEqn.H)
{
    volScalarField volContErr
    (
        fvc::div(phi) - mDot*(1.0/rho2 - 1.0/rho1)
    );

    const scalar sumLocalVolErr = runTime.deltaTValue()*
        mag(volContErr)().weightedAverage(mesh.V()).value();

    const scalar globalVolErr = runTime.deltaTValue()*
        volContErr.weightedAverage(mesh.V()).value();

    Info<< "time step volume continuity errors : sum local = "
        << sumLocalVolErr
        << ", global = " << globalVolErr
        << endl;
}
