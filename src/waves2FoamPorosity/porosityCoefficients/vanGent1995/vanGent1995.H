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
    Foam::vanGent1995

Description
    Porosity coefficients following the work by van Gent, M (1995).

SourceFiles
    vanGent1995.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef vanGent1995_H
#define vanGent1995_H

#include "porosityCoefficient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class vanGent1995
:
    public porosityCoefficient
{

protected:

    // Protected data

    // Protected member functions

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        vanGent1995(const vanGent1995&);

        //- Disallow default bitwise assignment
        void operator=(const vanGent1995&);

    // Private member data

public:

    //- Runtime type information
    TypeName("vanGent1995");

    // Constructors

        //- Construct from components
        vanGent1995
        (
            const dictionary &
        );


    // Destructor

        virtual ~vanGent1995();

    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
