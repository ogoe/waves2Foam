//    surfaceScalarField muEff
//    (
//        "muEff",
//        fvc::interpolate(turbulence->mut())
//    );

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
//        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
        pm->ddt(rho, U)
      + 1.0/porosity*fvm::div(rhoPhi/(porosityFace), U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    // mu is needed for the resistance field
    tmp<volScalarField> tmu = mixture.mu();
    const volScalarField& mu = tmu();

    pm->addResistance(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
