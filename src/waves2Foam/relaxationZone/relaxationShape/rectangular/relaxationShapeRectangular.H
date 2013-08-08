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
    Foam::relaxationShapes::relaxationShapeRectangular

Description
    An object used to define RECTANGULAR relaxation zones, which can be
    arbitrarily orientated in space. See e.g.

    @article { jacobsenFuhrmanFredsoe2011,
        Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
        title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
        Journal = {{Int. J. for Numer. Meth. Fluids}},
        Year = {2012},
        Volume = {70},
        Number = {9},
        Pages = {1073-1088},
        DOI = {{10.1002/fld.2726}},
    }

SourceFiles
    relaxationShapeRectangular.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef relaxationShapeRectangular_H
#define relaxationShapeRectangular_H

#include "relaxationShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

/*---------------------------------------------------------------------------*\
                           Class relaxationShapeRectangular Declaration
\*---------------------------------------------------------------------------*/

class relaxationShapeRectangular
:
    public relaxationShape
{

protected:

    // Protected data
        vectorField cornerNodes_;
        vector      orient_;
        vector      crossOrient_;
        vector      xAxis_; // The 'x-axis in the rotated coordinate system'
        vector      yAxis_; // The 'y-axis in the rotated coordinate system'

        scalar      width_;


    // Protected member functions
        bool insideZone
        (
            const label&
        ) const;

        virtual void computeSigmaCoordinate();

        virtual void findComputationalCells();

public:

    //- Runtime type information
    TypeName("relaxationShapeRectangular");

    // Constructors

        //- from components
        relaxationShapeRectangular
        (
            const word& ,
            const fvMesh& mesh_
        );


    // Destructor

        ~relaxationShapeRectangular()
        {}


    // Member Functions
        virtual const pointField& pointSet();

        virtual scalar interpolation(const scalarField&, const point&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
