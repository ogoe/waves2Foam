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
    Foam::write2Ascii

Description
    Reads the OF-formatted files and write the content in ASCII format without
    any headers, etc.

SourceFiles
    write2Ascii.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef write2Ascii_H
#define write2Ascii_H

#include "postProcessingWaves.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class write2Ascii
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
        write2Ascii(const write2Ascii&);

        //- Disallow default bitwise assignment
        void operator=(const write2Ascii&);

    // Private member data

public:

    //- Runtime type information
    TypeName("write2Ascii");

    // Constructors

        //- Construct from components
        write2Ascii
        (
            const Time&,
            const dictionary&,
            const word&
        );


    // Destructor

        virtual ~write2Ascii();

    // Member Functions
        virtual void evaluate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
