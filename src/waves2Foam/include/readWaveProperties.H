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

    /* Compute a reference point that is placed on the free surface */
    // A zero offset between original and new reference level. Default offset
    dimensionedVector referencePoint("NULL", dimLength, vector::zero);

	// Get the sea level
	scalar sL = readScalar(waveProperties.lookup("seaLevel"));

	if (SMALL < Foam::mag(sL))
	{
	    if (Switch(waveProperties.lookup("seaLevelAsReference")))
        {
        	// Make a positive unit vector along the direction of gravity
            referencePoint.value() = g.value()/Foam::mag(g.value());
            referencePoint.value() = Foam::cmptMag(referencePoint.value());

            // Make the sea level the reference level
        	referencePoint.value() *= sL;
        }
    }
