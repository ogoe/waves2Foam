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
    Foam::relaxationSchemes::relaxationSchemeSpatial

Description
    An object which results in relaxation, which is only dependent on the
    spatial location within a relaxationShape.

SourceFiles
    relaxationSchemeSpatial.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef relaxationSchemeSpatial_H
#define relaxationSchemeSpatial_H

#include "relaxationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationSchemes
{

/*---------------------------------------------------------------------------*\
                           Class relaxationSchemeSpatial Declaration
\*---------------------------------------------------------------------------*/

class relaxationSchemeSpatial
:
    public relaxationScheme
{

protected:

    // Protected data
        scalar exponent_;
        scalarField weight_;

    // Protected member functions


public:

    //- Runtime type information
    TypeName("relaxationSchemeSpatial");

    // Constructors

        //- from components
        relaxationSchemeSpatial
        (
            const word& ,
            const fvMesh& mesh_,
            vectorField& U,
            scalarField& alpha
        );

    // Destructor

        virtual ~relaxationSchemeSpatial()
        {}

    // Member Functions
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
