// Included in solveFluid.H just before the alpha step.
//
// 1. Store the liquid volume for liquidVolumeBalance.H.
// 2. Check the phase change source against the cap at the deltaT the alpha
//    step will actually use. sharpMdot.H capped mDot with
//    maxDeltaTGrowth*deltaT of the previous step, so normally nothing is
//    reported. If cells are reported, the source would by itself push alpha1
//    out of [0, 1]; pEqn.H has already used this mDot, so it is not changed
//    here, only reported (reduce maxCo / maxDeltaT or check the time step
//    control).

const scalar liquidVol0 = gSum(alpha1.primitiveField()*mesh.V());

{
    const scalar dt = runTime.deltaTValue();
    const scalar rhoL = rho1.value();
    const scalarField& V = mesh.V();

    label nOver = 0;
    scalar mOver = 0;

    forAll(mDot, celli)
    {
        const scalar m = mDot[celli];

        if (m != 0)
        {
            const scalar a = min(max(alpha1[celli], scalar(0)), scalar(1));
            const scalar lim = (m > 0 ? rhoL*a/dt : -rhoL*(1 - a)/dt);

            if ((m > 0 && m > lim) || (m < 0 && m < lim))
            {
                ++nOver;
                mOver += (m - lim)*V[celli];
            }
        }
    }

    reduce(nOver, sumOp<label>());
    reduce(mOver, sumOp<scalar>());

    if (nOver > 0)
    {
        Info<< "alphaSourceCheck: mDot exceeds the cap at deltaT = " << dt
            << " in " << nOver << " cells, excess = " << mOver << " kg/s"
            << endl;
    }
}
