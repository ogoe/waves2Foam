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
    Foam::relaxationWeights::relaxationWeightCourantCorrection

Description
    This class performs an additional weighting of the relaxation zone weights,
    which are based on the local cell Courant number.

    This is based on a presentation by Sopheak Seng, Technical University of
    Denmark.

SourceFiles
    relaxationWeightCourantCorrection.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef relaxationWeightCourantCorrection_H
#define relaxationWeightCourantCorrection_H

#include "fvMesh.H"
#include "fvCFD.H"
//#include "volFields.H"

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationWeights
{

/*---------------------------------------------------------------------------*\
                      Class relaxationShape Declaration
\*---------------------------------------------------------------------------*/

class relaxationWeightCourantCorrection
{

protected:

private:
    // Private data
        const fvMesh& mesh_;

        Switch courantCorrection_;

    // Private Member Functions

        scalar courantNumber
        (
               const surfaceScalarField& phi,
               const surfaceScalarField& magSf,
               const surfaceScalarField& deltaCoeffs,
               const cell& cc
        ) const;

        //- Disallow default bitwise copy construct
        relaxationWeightCourantCorrection
        (
            const relaxationWeightCourantCorrection&
        );

        //- Disallow default bitwise assignment
        void operator=(const relaxationWeightCourantCorrection&);

public:

    //- Runtime type information
    TypeName("relaxationWeightCourantCorrection");

    // Constructors

        //- Construct from components
        relaxationWeightCourantCorrection
        (
            const fvMesh&,
            const dictionary&
        );

    // Destructor

        virtual ~relaxationWeightCourantCorrection();

    // Member Functions
        void courantCorrection
        (
            const labelList&,
            scalarField&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationWeights
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
