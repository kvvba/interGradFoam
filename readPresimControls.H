// Controls for single-phase presimulation, read from the optional
// presimControl sub-dictionary of system/fvSolution:
//
//     presimControl
//     {
//         computeFluxes      true;  // compute and log wall heat fluxes
//         stopWhenConverged  true;  // end the run once converged
//         fluxTolerance      1e-8;  // tolerance on the net wall heat flux [W]
//         TTolerance         1e-8;  // tolerance on the relative change in max(T)
//     }
//
// Every entry defaults to off (or 1e-8) when the dictionary is missing.
// This file is included every time step, so fvSolution is re-read and
// the controls can be changed while the case is running.

fvSolution solutionDict(runTime);

const dictionary* presimControlPtr = solutionDict.findDict("presimControl");

// Compute the wall heat fluxes (computeFluidHeatFluxes.H)
const Switch computeFluxes =
    presimControlPtr ? presimControlPtr->getOrDefault<Switch>("computeFluxes", false) : Switch(false);

// Cleared in computeFluidHeatFluxes.H by any fluid region that does not
// meet both tolerances below, so it is only true if all regions converged.
// Starts false when the fluxes are not computed at all.
bool converged(computeFluxes);

// End the run when converged
const Switch stopWhenConverged =
    presimControlPtr ? presimControlPtr->getOrDefault<Switch>("stopWhenConverged", false) : Switch(false);

// Convergence criteria: magnitude of the net heat flux through all
// non-processor fluid boundaries, and the relative change in max(T)
// between time steps
const scalar fluxTol(presimControlPtr ? presimControlPtr->getOrDefault<double>("fluxTolerance", 1e-8) : 1e-8);
const scalar TTol(presimControlPtr ? presimControlPtr->getOrDefault<double>("TTolerance", 1e-8) : 1e-8);
