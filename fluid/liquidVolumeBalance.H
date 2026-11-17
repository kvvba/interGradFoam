// Included in solveFluid.H just after the alpha step (needs liquidVol0 from
// alphaSourceCheck.H).
//
// Global liquid volume balance of the alpha step [m^3]:
//     dV        = sum(alpha1*V) after - before
//     source    = -sum(mDot*V)/rho1*dt        (Su of alphaSuSp.H)
//     outflow   = dt*sum(alphaPhi1Avg) over non-coupled boundary faces
//                 (alpha1 flux averaged over the alpha sub-cycles)
//     imbalance = dV - (source - outflow)
// The imbalance should be ~0. It contains alpha1*divU with divU =
// fvc::div(phi) - mDot*(1/rho2 - 1/rho1), i.e. the pressure solver
// continuity error (alphaSuSp.H), and any clipping of alpha1. With the Euler ddt scheme
// alphaPhi1Avg is the flux actually used; with CrankNicolson it is the
// end-of-step flux and the balance is only approximate.
{
    const scalar dt = runTime.deltaTValue();

    const scalar liquidVol1 = gSum(alpha1.primitiveField()*mesh.V());

    const scalar sourceVol =
       -gSum(mDot.primitiveField()*mesh.V())/rho1.value()*dt;

    scalar outflowVol = 0;
    forAll(alphaPhi1Avg.boundaryField(), patchi)
    {
        const fvsPatchScalarField& apf = alphaPhi1Avg.boundaryField()[patchi];

        if (!apf.coupled())
        {
            outflowVol += sum(apf);
        }
    }
    reduce(outflowVol, sumOp<scalar>());
    outflowVol *= dt;

    const scalar dV = liquidVol1 - liquidVol0;
    const scalar imbalance = dV - (sourceVol - outflowVol);

    Info<< "liquidVolumeBalance: dV = " << dV
        << " source = " << sourceVol
        << " outflow = " << outflowVol
        << " imbalance = " << imbalance
        << " imbalance/source = "
        << imbalance/(mag(sourceVol) + VSMALL)
        << endl;
}
