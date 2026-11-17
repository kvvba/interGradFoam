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

// Phase change sources for the alpha1 (liquid) equation in VoF/alphaEqn.H:
//
//     ddt(alpha1) + div(phi, alpha1) = Su + (Sp + divU)*alpha1
//
// The target equation (rho1 = liquid, rho2 = vapour) is
//
//     ddt(alpha1) + div(alpha1 U) = -mDot/rho1                  (A)
//     div(U) = mDot*(1/rho2 - 1/rho1)                            (B, pEqn.H)
//
// which together conserve mixture mass exactly. (A) has no alpha1*div(U)
// term, so the alpha1*divU term of alphaEqn.H must vanish where the phase
// change volume source sits. It is removed by subtracting the phase change
// part from divU:
//
//     Su   = -mDot/rho1
//     Sp   = 0
//     divU = fvc::div(phi) - mDot*(1/rho2 - 1/rho1)
//
// divU is then only the continuity error of the pressure solution, i.e. the
// role it has in interFoam. mDot is the sharp, capped source from
// sharpMdot.H of the previous time step, i.e. the one the current phi was
// made divergent with.
//
// Why not Sp = -mDot*(1/rho2 - 1/rho1) with divU = fvc::div(phi) (the
// interPhaseChangeFoam pattern)? Sp is implicit (new alpha1) but divU is
// explicit (old alpha1) in MULES::explicitSolve, and only Sp enters the
// MULES::correct step of MULESCorr. The two then do not cancel: the update
// becomes (Su - div(alphaPhi))*dt/(1 + D*dt) with D = mDot*(1/rho2 - 1/rho1),
// an error of about (rho1/rho2)*(change of alpha1 per step) relative to the
// evaporated volume. With rho1/rho2 ~ 1600 and maxCo = 0.1 this was ~60x the
// source on the 7.5 um spherical bubble case (log.sharpGradFoam_splitSp).
// With the form above the phase change part cancels exactly in the
// predictor, the explicit solve and the correction.
//
// Boundedness of the explicit Su is ensured by the cap in sharpMdot.H
// (mDot*dt/rho1 <= alpha1) together with the MULES limiter.
//
// This is the copy that is compiled (-I./fluid comes before -I./fluid/VoF in
// Make/options; the VoF/ copy was removed). It is included in every alpha
// corrector, so it must not print anything.

volScalarField::Internal Su
(
    "alphaSu",
    -mDot()/rho1
);

zeroField Sp;

volScalarField::Internal divU(fvc::div(phi));
divU -= mDot()*(1.0/rho2 - 1.0/rho1);
