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
// One time step of a fluid region:
//   1. advect alpha1 with MULES, with the phase change sources Su/Sp
//      (alphaSuSp.H) from the previous step's mDot
//   2. surface tension force
//   3. interfacial heat flux and sharp phase change mass source mDot
//   4. liquid and vapour temperatures
//   5. momentum predictor and PISO pressure correction, with the
//      phase change volume source mDot*(1/rho2 - 1/rho1) in the pressure
//      equation

mesh.update();

mixture.correct();

// --- VoF advection: ddt(alpha1) + div(alpha1 U) = -mDot/rho1, with mDot
//     from the previous time step (the one the current phi was made
//     divergent with in pEqn.H)
#include "alphaControls.H"

// Outflow from phase change source cells takes the receiving cell's alpha1
// (displacedPhaseFlux.H); "no" gives plain MULES face values for comparison
const bool displacedPhaseFlux
(
	alphaControls.getOrDefault<bool>("displacedPhaseFlux", true)
);

// Alpha sub-cycling (nAlphaSubCycles > 1) is allowed: Su and divU act over
// each sub-cycle, so over the full time step in total. mDot, and with it the
// set of source cells used by displacedPhaseFlux.H, stays fixed over the
// sub-cycles (it is computed once per time step in sharpMdot.H).

#include "alphaSourceCheck.H"
#include "alphaEqnSubCycle.H"
#include "liquidVolumeBalance.H"

// Update curvature and interface normals for the new alpha1
mixture.correct();

// --- Surface tension (CSF) face force brackFS, used in pEqn.H
#include "calcF.H"

// --- Reconstruct the interface and compute the interfacial heat flux
//     qDot from the liquid- and vapour-side temperature gradients
#include "calcQdot.H"

// --- Sharp mDot = qDot/hlv in the interface cells, capped so that the
//     source cannot empty or fill a cell in one time step
#include "sharpMdot.H"

// --- Liquid (TL) and vapour (TV) phase temperatures
#include "TEqn.H"

#include "UEqn.H"

// --- pressure loop
for (int corr=0; corr<nCorr; corr++)
	{
    #include "pEqn.H"
	}

turbulence.correct();

// --- Wall heat flux balance and convergence check for single-phase
//     presimulation
if(computeFluxes)
				{
								#include "computeFluidHeatFluxes.H"
				}
