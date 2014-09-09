    Info << "\nReading waveProperties" << endl;

    autoPtr<Foam::waveTheories::externalWaveForcing> externalWave =
        Foam::waveTheories::externalWaveForcing::New
        (
            IOobject
            (
                "externalWaveForcing",
                runTime.constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            runTime,
            mesh
        );
