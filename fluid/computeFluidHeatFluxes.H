// Single-phase presimulation diagnostics, run when presimControl/
// computeFluxes is on (see readPresimControls.H):
//  - conductive and convective heat flow [W] through every non-processor
//    boundary patch of the fluid region, and their sum over all patches
//  - max(T) and its relative change since the previous time step
// The net heat flow and max(T) are appended to <region>_heatFluxBalance.csv
// and <region>_maxResT.csv. The run counts as converged (thermal steady
// state) when |net heat flow| < fluxTolerance and the max(T) change <
// TTolerance.

scalar sumConduction = 0.0;
scalar sumConvection = 0.0;
scalar sumTot = 0.;

//- Effective conductivity, including the turbulent part as in TEqn.H
tmp<volScalarField> tkappaEff(kappa);
if (turbulentHeat)
{
    tkappaEff = kappa + rhoCp*turbulence.nut()/Prt;
}
const volScalarField& kappaEff = tkappaEff();

//- First, compute the heat fluxes in the fluid region
surfaceScalarField heatFluxD = -fvc::interpolate(kappaEff)*fvc::snGrad(T);

const surfaceScalarField::Boundary& patchHeatFluxD = heatFluxD.boundaryFieldRef();

if(mixture.found("writeQwall"))
{
    qwFluid[i] = -kappaEff*fvc::grad(T);
}


Info<< "\nWall heat fluxes [W] fluid region " << fluidRegions[i].name()
    << endl;
forAll(patchHeatFluxD, patchi)
{
    if(isA<processorFvPatchScalarField>(T.boundaryField()[patchi]))
    {
        continue;
    }

    scalar conduction = gSum
    (
        mesh.magSf().boundaryField()[patchi]
       *heatFluxD.boundaryField()[patchi]
    );

    scalar convection = 0.0;

    if(mesh.foundObject<surfaceScalarField>("phi"))
    {
        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>("phi");

        convection = gSum
        (
            rhoCp.boundaryField()[patchi]
            *T.boundaryField()[patchi]
            *phi.boundaryField()[patchi]
        );
    }


    Info<< mesh.boundary()[patchi].name()
        << " "
        << conduction
        << " "
        << convection
        << " "
        << conduction + convection
        << endl;

    sumConduction += conduction;
    sumConvection += convection;

}

sumTot = sumConduction + sumConvection;
Info<< "sum "
    << sumConduction
    << " "
    << sumConvection
    << " "
    << sumTot
    << nl << endl;

//- Compute maximum temperature
scalar maxTnew(gMax(T));

scalar Tres(mag(maxTFluid[i]-maxTnew)/(maxTFluid[i]+SMALL));

maxTFluid[i] = maxTnew;


Info<< "maxT = " << maxTFluid[i]  << "   residuals on T = "
    << Tres << nl << endl;

if(Pstream::master())
{
    std::ofstream outfile;

    outfile.open
    (
        fluidRegions[i].name() + "_heatFluxBalance.csv", std::ios_base::app
    );

    outfile << runTime.value() << ", "
            << sumConduction + sumConvection<< "\n";

    outfile.close();

    outfile.open
    (
        fluidRegions[i].name() + "_maxResT.csv", std::ios_base::app
    );

    outfile << runTime.value() << ", "
            << maxTFluid[i] << ", "
            << Tres << "\n";


}

// converged starts true each step (readPresimControls.H); any region that
// has not converged clears it
if(!(mag(sumTot) < fluxTol && Tres < TTol))
{
    converged = false;
}
