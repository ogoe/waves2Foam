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
    Foam::relaxationSchemes::relaxationSchemeSpatialInterpolation

Description
    This relaxation scheme utilises interpolation over a pre-defined interval
    for the intersection. This speed of the execution of the intersections 
    with a huge factor. Especially when using irregular waves with a lot of 
    wave components.
    Should be combined with a relaxation shape that implements the interpolation
    method.

SourceFiles
    relaxationSchemeSpatialInterpolation.C

Author
    Niels Gjøl Jacobsen, Deltares.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef relaxationSchemeSpatialInterpolation_H
#define relaxationSchemeSpatialInterpolation_H

#include "relaxationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationSchemes
{

/*---------------------------------------------------------------------------*\
            Class relaxationSchemeSpatialInterpolation Declaration
\*---------------------------------------------------------------------------*/

class relaxationSchemeSpatialInterpolation
:
    public relaxationScheme
{

protected:

    // Protected data
        scalar exponent_;
        scalarField weight_;

        scalarField surfaceElevation_;

    // Protected member functions

        //- Returns a scalarField of the signed distance to an arbitrary
        //  surface
        virtual void signedPointToSurfaceDistance
        (
            const pointField&,
            scalarField&
        );


        //- Returns a scalar of the signed distance to an arbitrary surface
        virtual scalar signedPointToSurfaceDistance
        (
            const point&
        ) const;

public:

    //- Runtime type information
    TypeName("relaxationSchemeSpatialInterpolation");

    // Constructors

        //- from components
        relaxationSchemeSpatialInterpolation
        (
            const word& ,
            const fvMesh& mesh_,
            vectorField& U,
            scalarField& alpha
        );

    // Destructor

        virtual ~relaxationSchemeSpatialInterpolation()
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
