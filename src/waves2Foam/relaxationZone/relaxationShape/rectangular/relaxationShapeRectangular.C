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

#include "relaxationShapeRectangular.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationShapeRectangular, 0);
addToRunTimeSelectionTable
(
    relaxationShape,
    relaxationShapeRectangular,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


relaxationShapeRectangular::relaxationShapeRectangular
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    relaxationShape(subDictName, mesh_)
{
    orient_ = vector(coeffDict_.lookup("orientation"));
    cornerNodes_.setSize(4);

    // Geometric properties needed to uniquely define relaxation-shape
    orient_ /= Foam::mag(orient_);

    cornerNodes_[0] = vector(coeffDict_.lookup("startX"));
    cornerNodes_[2] = vector(coeffDict_.lookup("endX"));
    cornerNodes_[3] = cornerNodes_[0]
        + (orient_ & (cornerNodes_[2] - cornerNodes_[0]))*orient_;
    cornerNodes_[1] = cornerNodes_[0] + (cornerNodes_[2] - cornerNodes_[3]);

    crossOrient_ = (cornerNodes_[1] - cornerNodes_[0]);
    crossOrient_ /= Foam::mag(crossOrient_);

    width_ = Foam::mag(cornerNodes_[3] - cornerNodes_[0]);

    xAxis_ = Foam::cmptMultiply( orient_, orient_ );
    yAxis_ = Foam::cmptMultiply( crossOrient_, crossOrient_ );

    coeffDict_.lookup("relaxType") >> relaxType_;

    // Find computational cells inside the relaxation-shape
    findComputationalCells();

    // Computate the sigma coordinate
    computeSigmaCoordinate();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationShapeRectangular::findComputationalCells()
{
    const vectorField& cc = mesh_.C();

    cells_.setSize(5000);
    label count(0);

    forAll (cc, celli)
    {
        if (insideZone( celli ))
        {
            cells_[count++] = celli;

            if (count == cells_.size())
            {
                cells_.setSize( static_cast<label>( count*1.1 ) );
            }
        }
    }

    cells_.setSize(count);
}


void relaxationShapeRectangular::computeSigmaCoordinate()
{
    const vectorField& C = mesh_.C();
    sigma_.setSize(cells_.size(), 0);

    forAll (cells_, celli)
    {
        vector cc( C[cells_[celli]] );
        cc -= ( ( cc & direction_ )*direction_ );

        sigma_[celli] = Foam::mag(((xAxis_ & cc) - (xAxis_ & cornerNodes_[0])))
            /( Foam::mag(xAxis_)*width_ );

        if (relaxType_ == "INLET")
        {
            sigma_[celli] = Foam::mag( sigma_[celli] - 1.0 );
        }
    }
}


bool relaxationShapeRectangular::insideZone
(
    const label& celli
) const
{
    bool inside( false );

    label negative(0);
    label positive(0);

    scalar pX(mesh_.C()[celli] & xAxis_);
    scalar pY(mesh_.C()[celli] & yAxis_);

    scalarField nX(cornerNodes_ & xAxis_);
    scalarField nY(cornerNodes_ & yAxis_);

    // Uses method 3 on
    // http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/
    // to consider whether or not the point is inside the relaxation zone
    forAll (nX, pointi)
    {
        scalar temp(0);
        label pointj((pointi == nX.size()-1) ? 0 : pointi + 1);
        temp = (pY - nY[pointi])*(nX[pointj] - nX[pointi])
            - (pX - nX[pointi])*(nY[pointj] - nY[pointi]);

        if (temp > 0)
        {
            positive++;
        }
        if (temp < 0)
        {
            negative++;
        }
    }

    if (positive == 0 || negative == 0)
    {
        inside = true;
    }

    return inside;
}


const pointField& relaxationShapeRectangular::pointSet()
{
    notImplemented("pointSet is not implemented for this shape");
}


scalar relaxationShapeRectangular::interpolation
(
    const scalarField& source,
    const point& p0
) const
{
    notImplemented("interpolation is not implemented for this shape");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// ************************************************************************* //
