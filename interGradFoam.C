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

Application
    interGradFoam

Description
    Transient solver for multi-region incompressible two-phase interface
    tracking simulations with phase change.

    Derived from mulesGradFoam. The phase change mass source is sharp: it
    sits in the interface cells only, without the Hardt & Wondra (2008)
    smoothing (fluid/sharpMdot.H), and enters both the pressure equation
    and the MULES alpha equation (fluid/alphaSuSp.H), as in PSI-BOIL
    (Sato & Niceno 2013; Lal, Sato & Niceno 2015).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
// #include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "incompressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
// #include "radiationModel.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "pressureControl.H"
#include "processorFvPatchFields.H"
#include "volPointInterpolation.H"
#include "cutCellIso.H"
#include "isoAlpha.H"
#include "implicitInterfaceDiffFlux.H"
#include "meshSearch.H"
#include "syncTools.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
// #include "explicitInterfaceDiffFlux.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for multi-region two-phase "
				"interface tracking phase change simulations."
    );

    // Support running as "interGradFoam -postProcess"; in that mode the
    // meshes are created by createMeshesPostProcess.H instead of the
    // default single-mesh createMesh.H
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    
    #include "setRootCaseLists.H"
    #include "createTime.H"
    // Fluid (dynamicFvMesh) and solid (fvMesh) regions, as listed in
    // constant/regionProperties
    #include "createMeshes.H"
		// #include "createMeshControls.H"
    // Per-region field lists for all fluid and solid regions
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    // Initial time step from the largest fluid Courant number (CoNum)
    // and the largest solid diffusion number (DiNum) over all regions
    #include "incompressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    // Largest factor by which setMultiRegionDeltaT.H lets deltaT grow from
    // one step to the next. sharpMdot.H caps mDot with this factor because
    // mDot is used by the alpha step of the next (possibly longer) step.
    const scalar maxDeltaTGrowth = 1.2;

    while (runTime.run())
    {
        // Re-read the time step and single-phase presimulation controls
        // every step so that they can be changed while the case is running
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPresimControls.H"
			  // #include "readMeshControls.H"

        // Adjust the time step to satisfy maxCo in the fluid regions and
        // maxDi in the solid regions
        #include "incompressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // if (nOuterCorr != 1)
        // {
        //     forAll(fluidRegions, i)
        //     {
        //         #include "storeOldFluidFields.H"
        //     }
        // }

				// Regions are solved one after another (segregated) and only
				// once per time step: there is no outer loop over regions.
				// The coupling between them happens through the boundary
				// conditions on the region interfaces.
				forAll(fluidRegions, i)
					{
						Info<< "\nSolving for fluid region "
								<< fluidRegions[i].name() << endl;
            #include "setRegionFluidFields.H"
            #include "readFluidMultiRegionPISOControls.H"
            #include "solveFluid.H"
					}

				forAll(solidRegions, i)
					{
						Info<< "\nSolving for solid region "
								<< solidRegions[i].name() << endl;
            #include "setRegionSolidFields.H"
            #include "readSolidMultiRegionPISOControls.H"
            #include "solveSolid.H"
					}
				
        runTime.write();

        runTime.printExecutionTime(Info);

				// Single-phase presimulation: stop once the wall heat balance
				// and max(T) have converged (converged is set in
				// computeFluidHeatFluxes.H)
				if(stopWhenConverged && converged)
        {
            Info<< "All fluxes are converged."
                << "The simulation ends here." << endl;

            runTime.writeAndEnd();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
