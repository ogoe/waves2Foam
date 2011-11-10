    Info << "\nReading waveProperties" << endl;

    IOdictionary waveProperties
    (
        IOobject
        (
            "waveProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
