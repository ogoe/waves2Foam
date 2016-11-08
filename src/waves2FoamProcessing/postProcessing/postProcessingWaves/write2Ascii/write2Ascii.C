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

\*---------------------------------------------------------------------------*/

#include "write2Ascii.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(write2Ascii, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    write2Ascii,
    postProcessingWaves
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


write2Ascii::write2Ascii
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action )
{
}


write2Ascii::~write2Ascii()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void write2Ascii::evaluate()
{
    // Writing the time instances to ascii file
    IOField<scalar> time
    (
        IOobject
        (
            callName_ + "_time",
            rT_.constant(),
            addDir_,
            rT_,
            IOobject::MUST_READ
        )
    );

    autoPtr<OFstream> asciiPtr_;
    asciiPtr_.reset
    (
        new OFstream(directDir_ + "/" + callName_ + "_time.txt")
    );

    forAll (time, timei)
    {
        asciiPtr_() << time[timei] << endl;
    }

    // Writing the data fields
    label count(0);

    while (true)
    {
        std::stringstream ss;
        ss << callName_ << "_" << count;

        IOobject fileHeader
        (
            ss.str(),
            rT_.constant(),
            addDir_,
            rT_,
            IOobject::MUST_READ
        );

#if OFPLUSBRANCH == 1
    #if OFVERSION<1606
        if (!fileHeader.headerOk())
        {
            break;
        }
    #else
        // Step required to update the header
        fileHeader.typeHeaderOk<volScalarField>(false);

        if (!(fileHeader.headerClassName() == "scalarField" 
              || fileHeader.headerClassName() == "vectorField"))
        {
            // The check for the header is discus
            break;
        }
    #endif
#else
        if (!fileHeader.headerOk())
        {
            break;
        }
#endif

        asciiPtr_.reset(new OFstream( directDir_ + "/" + ss.str() + ".txt" ) );

        if (fileHeader.headerClassName() == "scalarField")
        {
            IOField<scalar> field(fileHeader);

            forAll (field, datai)
            {
                asciiPtr_() << field[datai] << endl;
            }
        }
        else if (fileHeader.headerClassName() == "vectorField")
        {
            IOField<vector> field(fileHeader);

            forAll (field, datai)
            {
                asciiPtr_() << field[datai].x() << tab << field[datai].y()
                            << tab << field[datai].z() << endl;
            }
        }
        else
        {
            FatalErrorIn("void write2Ascii::evaluate()" )
                    << "The field type" << fileHeader.headerClassName() << endl
                    << "is not supported." << endl << exit(FatalError);
        }

        count++;
    }

    // Writing the content of the dictionary
    IOdictionary dict
    (
        IOobject
        (
            callName_ + "_dict",
            rT_.constant(),
            addDir_,
            rT_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Writing time step
    scalar dt = readScalar( dict.lookup("deltaT") );

    asciiPtr_.reset
    (
        new OFstream(directDir_ + "/" + callName_ + "_deltaT.txt")
    );

    asciiPtr_() << dt << endl;

    // Writing indices and spatial coordinates
    labelList indices( dict.lookup("index") );

    // Either locations or names are stated in the dictionary
    if (dict.found("x"))
    {
        scalarField x( dict.lookup("x") );
        scalarField y( dict.lookup("y") );
        scalarField z( dict.lookup("z") );

        asciiPtr_.reset
        (
            new OFstream(directDir_ + "/" + callName_ + "_indexXYZ.txt")
        );

        forAll (indices, indexi)
        {
            asciiPtr_() << indices[indexi] << tab << x[indexi] << tab
                        << y[indexi] << tab << z[indexi] << endl;
        }
    }
    else
    {
        wordList names( dict.lookup("names") );

        asciiPtr_.reset
        (
            new OFstream(directDir_ + "/" + callName_ + "_indexNames.txt")
        );

        forAll (indices, indexi)
        {
            asciiPtr_() << indices[indexi] << tab << names[indexi] << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
