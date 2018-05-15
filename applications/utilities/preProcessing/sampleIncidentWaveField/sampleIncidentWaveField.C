/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    syntheticWaveField

Description


Author
    Niels Gjoel Jacobsen, Deltares.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "OFstream.H"

#include "uniformDimensionedFields.H"

#include "waveTheory.H"

using namespace Foam;

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    Info << "\nReading g" << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Create the waveProperties dictionary
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

	const dictionary& subDict
	    (waveProperties.subDict("sampleIncidentWaveField"));

    // Make output file
    fileName fn = "syntheticWaveField/0";
    Foam::mkDir(fn);

    autoPtr<OFstream> osPtr;
    osPtr.reset
    (
        new OFstream
        (
            runTime.time().path() + "/" + fn + "/" + "surfaceElevation.dat"
        )
    );

    autoPtr<OFstream> osPtrU;
        osPtrU.reset
        (
            new OFstream
            (
                runTime.time().path() + "/" + fn + "/" + "U.dat"
            )
        );

    autoPtr<OFstream> osPtrV;
        osPtrV.reset
        (
    		new OFstream
       		(
   				runTime.time().path() + "/" + fn + "/" + "V.dat"
       		)
        );

    autoPtr<OFstream> osPtrW;
        osPtrW.reset
        (
       		new OFstream
       		(
    			runTime.time().path() + "/" + fn + "/" + "W.dat"
       		)
        );

    // Get all the relaxation zones
    wordList relaxationNames(subDict.lookup("names"));

    // Read all points to be read
    pointField input(subDict.lookup("points"));

    // Write the header file
    osPtr() << "Time";
    osPtrU() << "Time";
    osPtrV() << "Time";
    osPtrW() << "Time";

    forAll(relaxationNames, relaxi)
    {
    	forAll(input, pointi)
        {
            osPtr() << tab << relaxationNames[relaxi];
            osPtrU() << tab << relaxationNames[relaxi];
            osPtrV() << tab << relaxationNames[relaxi];
            osPtrW() << tab << relaxationNames[relaxi];
        }

    }
    osPtr() << endl;
    osPtrU() << endl;
    osPtrV() << endl;
    osPtrW() << endl;

    osPtr() << "-1";
    osPtrU() << "-1";
    osPtrV() << "-1";
    osPtrW() << "-1";
    forAll(relaxationNames, relaxi)
    {
    	forAll(input, pointi)
    	{
            osPtr() << tab << input[pointi].x();
            osPtrU() << tab << input[pointi].x();
            osPtrV() << tab << input[pointi].x();
            osPtrW() << tab << input[pointi].x();
    	}
    }
    osPtr() << endl;
    osPtrU() << endl;
    osPtrV() << endl;
    osPtrW() << endl;


    osPtr() << "-2";
    osPtrU() << "-2";
    osPtrV() << "-2";
    osPtrW() << "-2";

    forAll(relaxationNames, relaxi)
    {
    	forAll(input, pointi)
    	{
            osPtr() << tab << input[pointi].y();
            osPtrU() << tab << input[pointi].y();
            osPtrV() << tab << input[pointi].y();
            osPtrW() << tab << input[pointi].y();
    	}
    }
    osPtr() << endl;
    osPtrU() << endl;
    osPtrV() << endl;
    osPtrW() << endl;

    osPtr() << "-3";
    osPtrU() << "-3";
    osPtrV() << "-3";
    osPtrW() << "-3";

    forAll(input, pointi)
    {
        forAll(relaxationNames, relaxi)
    	{
            osPtr() << tab << input[pointi].z();
            osPtrU() << tab << input[pointi].z();
            osPtrV() << tab << input[pointi].z();
            osPtrW() << tab << input[pointi].z();
    	}
    }
    osPtr() << endl;
    osPtrU() << endl;
    osPtrV() << endl;
    osPtrW() << endl;

    // Prepare all the wave theories
    List<autoPtr<Foam::waveTheories::waveTheory> > theories(relaxationNames.size());

    forAll(relaxationNames, relaxi)
    {
        theories[relaxi] = Foam::waveTheories::waveTheory::New(relaxationNames[relaxi], mesh);
    }

    // Make the time axis
    scalar t = 0;
    scalar dt = readScalar(subDict.lookup("deltaT"));
    scalar T = readScalar(subDict.lookup("endTime"));

    label count = 0;

    // Loop over all times until finished
    while (t <= T)
    {

        osPtr() << t;
        osPtrU() << t;
        osPtrV() << t;
        osPtrW() << t;

        forAll(relaxationNames, relaxi)
        {
            scalarField eta = theories[relaxi]->eta(input, t);

            forAll (eta, etai)
            {
                osPtr() << tab << eta[etai];

                if (-(theories[relaxi]->returnDir() & input[etai]) < eta[etai])
                {
                	vector U = theories[relaxi]->U(input[etai], t);
                	osPtrU() << tab << U.x();
                	osPtrV() << tab << U.y();
                	osPtrW() << tab << U.z();
                }
                else
                {
                	vector U = vector::zero;
                	osPtrU() << tab << U.x();
                	osPtrV() << tab << U.y();
                	osPtrW() << tab << U.z();
                }
            }
        }

        osPtr() << endl;
        osPtrU() << endl;
        osPtrV() << endl;
        osPtrW() << endl;

        t += dt;
    }

    Info << nl << nl << "End" << nl << endl;

    return 0;
}
