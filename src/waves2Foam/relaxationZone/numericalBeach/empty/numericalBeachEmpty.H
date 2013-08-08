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
    Foam::numericalBeaches::numericalBeachEmpty

Description
    A dummy numerical beach to be used, if a relaxation zone should not include
    artificial viscosity contributions in the momentum equation

SourceFiles
    numericalBeachEmpty.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef numericalBeachEmpty_H
#define numericalBeachEmpty_H

#include "numericalBeach.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace numericalBeaches
{

/*---------------------------------------------------------------------------*\
                     Class numericalBeachEmpty Declaration
\*---------------------------------------------------------------------------*/

class numericalBeachEmpty
:
    public numericalBeach
{

protected:


public:

    //- Runtime type information
    TypeName("numericalBeachEmpty");

    // Constructors

        //- from components
        numericalBeachEmpty
        (
            const word& ,
            const fvMesh& mesh_
        );


    // Destructor

        ~numericalBeachEmpty()
        {}


    // Member Functions

        virtual void correct
        (
            const labelList&,
            const scalarField&,
            volScalarField&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numericalBeaches
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
