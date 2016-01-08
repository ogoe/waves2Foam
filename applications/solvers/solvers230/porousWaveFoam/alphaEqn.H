{
	if (MULESCorr)
	{
		FatalErrorIn("alphaEqn.H")
			<< "Implicit solution to the alpha field is chosen for porousWaveFoam.\n"
			<< "No implementation is currently available that corrects for the\n"
			<< "porosity in the implicit advection of alpha.water.\n"
			<< "\nSwitch to the explicit advection routine."
			<< endl << endl << exit(FatalError);
	}

    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // Standard face-flux compression coefficient
    surfaceScalarField phic(interface.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (interface.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phic.boundaryField()[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> tphiAlpha;

    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)()
          + 1.0/porosity*fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)()
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;

        tmp<surfaceScalarField> tphiAlphaUD(alpha1Eqn.flux());
        tphiAlpha = tmp<surfaceScalarField>
        (
            new surfaceScalarField(tphiAlphaUD())
        );

        // Cache the upwind-flux
        tphiAlphaCorr0 = tphiAlphaUD;

        alpha2 = 1.0 - alpha1;

        interface.correct();
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phir(phic*interface.nHatf());

        tmp<surfaceScalarField> tphiAlphaUn
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tphiAlphaCorr(tphiAlphaUn() - tphiAlpha());
            volScalarField alpha10(alpha1);

            MULES::correct
            (
//                oneField(),
            	porosity,
            	alpha1,
            	tphiAlphaUn(),
            	tphiAlphaCorr(),
            	zeroField(),
            	zeroField(),
            	1,
            	0
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                tphiAlpha() += tphiAlphaCorr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                tphiAlpha() += 0.5*tphiAlphaCorr();
            }
        }
        else
        {
            tphiAlpha = tphiAlphaUn;

            MULES::explicitSolve
            (
            	porosity,
            	alpha1,
            	phi,
            	tphiAlpha(),
            	zeroField(),
            	zeroField(),
            	1,
            	0
            );
        }

        alpha2 = 1.0 - alpha1;

        interface.correct();
    }

    rhoPhi = tphiAlpha()*(rho1 - rho2) + phi*rho2;

    scalar porousVolume = Foam::sum
        (
            alpha1.internalField()*mesh.V()*porosity.internalField()
        );
    porousVolume /= Foam::sum(mesh.V()).value();

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << " Porosity weighted volume: " << porousVolume
        << endl;
}
