    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    // Calculate and cache mu for the porous media
    volScalarField mu(twoPhaseProperties.mu());

    /*fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + turbulence->divDevRhoReff(rho, U)
    );*/

    fvVectorMatrix UEqn
    (
          pm->ddt(rho, U)
        + 1.0/porosity*fvm::div(rhoPhi/(porosityFace), U)
        - fvm::laplacian(muEff/porosityFace, U)
        - (fvc::grad(U) & (1.0/porosity*fvc::grad(muEff)))
    );

    UEqn.relax();

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
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );
    }
