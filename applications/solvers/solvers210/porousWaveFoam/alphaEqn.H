{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir(phic*interface.nHatf());

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        // MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);
        MULES::explicitSolve
        (
            porosity,
            alpha1,
            phi,
            phiAlpha,
            zeroField(),
            zeroField(),
            1,
            0
        );

        rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;
    }

    scalar porousVolume = Foam::sum
        (
            alpha1.internalField()*mesh.V()*porosity.internalField()
        );
    porousVolume /= Foam::sum(mesh.V()).value();

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value() << endl
        << " Porosity weighted volume: " << porousVolume
        << endl;
}
