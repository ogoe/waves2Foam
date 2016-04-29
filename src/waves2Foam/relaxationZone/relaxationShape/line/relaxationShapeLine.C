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

#include "relaxationShapeLine.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationShapeLine, 0);
addToRunTimeSelectionTable
(
    relaxationShape,
    relaxationShapeLine,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


relaxationShapeLine::relaxationShapeLine
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    relaxationShapeRectangular(subDictName, mesh_),

    N_(readLabel(coeffDict_.lookup("N"))),

    dataRange_(2, -GREAT),

    localCoordinates_(N_, 0),

    globalCoordinates_(N_, point::zero),

    updated_(-1)
{
    dataRange_[0] = GREAT;
    dataRange();

    localCoordinates();

    globalCoordinates();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationShapeLine::dataRange()
{
    const labelListList& cellPoints = mesh_.cellPoints();
    const pointField& points = mesh_.points();

    scalar l, Min(dataRange_[0]), Max(dataRange_[1]);

    forAll (cells_, celli)
    {
        const labelList& cp = cellPoints[cells_[celli]];

        forAll (cp, pointi)
        {
            l = localCoordinate(points[cp[pointi]]);
            Min = Foam::min(l,Min);
            Max = Foam::max(l,Max);
        }
    }

    reduce(Min, minOp<scalar>());
    reduce(Max, maxOp<scalar>());

    scalar dist = Max - Min;
    Max += 0.005 * dist;
    Min -= 0.005 * dist;

    dataRange_[0] = Min;
    dataRange_[1] = Max;
}


void relaxationShapeLine::localCoordinates()
{
    scalar dlc = (dataRange_[1] - dataRange_[0])/static_cast<scalar>(N_ - 1);

    forAll (localCoordinates_, pointi)
    {
        localCoordinates_[pointi] = dataRange_[0]
            + static_cast<scalar>(pointi)*dlc;
    }
}


void relaxationShapeLine::globalCoordinates()
{
    forAll (globalCoordinates_, pointi)
    {
        globalCoordinates_[pointi] =
            globalCoordinate(localCoordinates_[pointi]);
    }
}


scalar relaxationShapeLine::localCoordinate(const point& pp) const
{
    return ((pp - cornerNodes_[0]) & orient_);
}


point relaxationShapeLine::globalCoordinate(const scalar& lc) const
{
    return (lc * orient_ + cornerNodes_[0]);
}


const pointField& relaxationShapeLine::pointSet()
{
    globalCoordinates();

    return globalCoordinates_;
}


scalar relaxationShapeLine::interpolation
(
    const scalarField& source,
    const point& p0
) const
{
    scalar lc = localCoordinate(p0);
    scalar res(0);

//  OLD IMPLEMENTATION - BRUTE FORCE SEARCH
//    for (label i = 0; i < source.size() - 1; i++)
//    {
//        if (localCoordinates_[i] <= lc && lc < localCoordinates_[i + 1])
//        {
//            scalar dlc = localCoordinates_[i + 1] - localCoordinates_[i];
//            res = source[i] + (source[i + 1] - source[i])/dlc
//                *(lc - localCoordinates_[i]);
//            break;
//        }
//    }

    // Utilise the knowledge that the interpolation axis is equidistant when
    // location the index-range of 'lc' on this axis
    scalar length1 = localCoordinates_[localCoordinates_.size() - 1]
         - localCoordinates_[0];
    scalar length2 = lc - localCoordinates_[0];
    label index = static_cast<label>(length2/length1*localCoordinates_.size());

    // Perform the interpolation
    scalar dlc = localCoordinates_[index + 1] - localCoordinates_[index];
    res = source[index] + (source[index + 1] - source[index])/dlc
    		*(lc - localCoordinates_[index]);

    return res;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// ************************************************************************* //
