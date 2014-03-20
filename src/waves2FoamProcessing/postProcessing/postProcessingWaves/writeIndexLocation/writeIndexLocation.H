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

Class
    Foam::writeIndexLocation

Description
    Reads the OF-formatted files and write the content in ASCII format without
    any headers, etc.

SourceFiles
    writeIndexLocation.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef writeIndexLocation_H
#define writeIndexLocation_H

#include "postProcessingWaves.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class writeIndexLocation
:
    public postProcessingWaves
{

protected:

    // Protected data
        word inputDir_;

    // Protected member functions

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        writeIndexLocation(const writeIndexLocation&);

        //- Disallow default bitwise assignment
        void operator=(const writeIndexLocation&);

    // Private member data

public:

    //- Runtime type information
    TypeName("writeIndexLocation");

    // Constructors

        //- Construct from components
        writeIndexLocation
        (
            const Time&,
            const dictionary&,
            const word&
        );


    // Destructor

        virtual ~writeIndexLocation();

    // Member Functions
        virtual void evaluate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
