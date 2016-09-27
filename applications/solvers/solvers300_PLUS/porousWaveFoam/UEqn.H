    surfaceScalarField muEff
    (
        "muEff",
        mixture.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    // Calculate and cache mu for the porous media
    volScalarField mu(mixture.mu());

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        pm->ddt(rho, U) + 
        1/porosity*fvm::div(rhoPhi/(porosityFace), U)
      + MRF.DDt(rho, U)
//      + turbulence->divDevRhoReff(rho, U)
        - fvm::laplacian(muEff/porosityFace, U)
        - (fvc::grad(U) & (1.0/porosity*fvc::grad(muEff)))
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

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
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
