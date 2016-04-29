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

#include "relaxationSchemeSpatialInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationSchemes
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationSchemeSpatialInterpolation, 0);
addToRunTimeSelectionTable
(
    relaxationScheme,
    relaxationSchemeSpatialInterpolation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

relaxationSchemeSpatialInterpolation::relaxationSchemeSpatialInterpolation
(
    const word& subDictName,
    const fvMesh& mesh,
    vectorField& U,
    scalarField& alpha
)
:
    relaxationScheme(subDictName, mesh, U, alpha),

    exponent_(coeffDict_.lookupOrDefault<scalar>("alphaExponent",3.5)),

    surfaceElevation_(0)
{
    const scalarField sigma(relaxShape_->sigma());

    weight_.setSize(sigma.size());

    forAll (weight_, celli)
    {
        weight_[celli] = 1.0 -
            (Foam::exp(Foam::pow(sigma[celli],exponent_)) - 1.0)
            /(Foam::exp(1.0) - 1.0);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationSchemeSpatialInterpolation::signedPointToSurfaceDistance
(
    const pointField& pp,
    scalarField& sd
)
{
    forAll (pp, pointi)
    {
        sd[pointi] = signedPointToSurfaceDistance(pp[pointi]);
    }
}


scalar relaxationSchemeSpatialInterpolation::signedPointToSurfaceDistance
(
    const point& pp
) const
{
    scalar temp = relaxShape_->interpolation(surfaceElevation_, pp);
    temp += (waveProps_->returnDir() & pp);
    temp *= -1.0;

    return temp;
}


void relaxationSchemeSpatialInterpolation::correct()
{
    // Obtain relaxation zone cells and local sigma coordinate
    // The number of cells and the sigma coordinate can have changed
    // for dynamic meshes
    const labelList& cells = relaxShape_->cells();
    const scalarField& sigma = relaxShape_->sigma();

    const pointField& interpolationPoints = relaxShape_->pointSet();
    surfaceElevation_.setSize(interpolationPoints.size(), 0.0);

    forAll (surfaceElevation_, pointi)
    {
        surfaceElevation_[pointi] = waveProps_->eta
             (interpolationPoints[pointi], db().time().value());
    }

    // Compute the relaxation weights - only changes for moving/changing meshes
    if (weight_.size() != sigma.size())
    {
    	weight_.clear();
        weight_.setSize(sigma.size(), 0.0);
    }

    relaxWeight_->weights(cells, sigma, weight_);

    // Perform the correction
    const scalarField& V (mesh_.V());
    const vectorField& C (mesh_.C());
    const cellList& cc(mesh_.cells());
    const pointField& pp(mesh_.points());
    const faceList& fL(mesh_.faces());

    forAll (cells, celli)
    {
        const label cellNo = cells[celli];
        const pointField& p = cc[cellNo].points(fL, pp);

        // Evaluate the cell height and the signedDistance to the surface from
        // the cell centre
        scalar cellHeight
            (
                Foam::max(p & waveProps_->returnDir())
              - Foam::min(p & waveProps_->returnDir())
            );
        scalar sdc(signedPointToSurfaceDistance(C[cellNo]));

        // Target variables
        scalar alphaTarget(0.0);
        vector UTarget(waveProps_->windVelocity(mesh_.time().value()));

        // Only do cutting, if surface is close by
        if (Foam::mag(sdc) <= 2*cellHeight)
        {
            localCellNeg lc(mesh_, cellNo);

            dividePolyhedral(point::zero, vector::one, lc);

            if (lc.ccNeg().size() >= 4)
            {
                UTarget = waveProps_->U(lc.centreNeg(), mesh_.time().value());
                alphaTarget = lc.magNeg()/V[cellNo];
            }
        }
        else if (sdc < 0)
        {
            alphaTarget = 1.0;
            UTarget     = waveProps_->U(C[cellNo], mesh_.time().value());
        }
        else
        {
            alphaTarget = 0.0;
            UTarget     = waveProps_->windVelocity(mesh_.time().value());
        }

        U_[cellNo] = (1 - weight_[celli])*UTarget + weight_[celli]*U_[cellNo];

        alpha_[cellNo] = (1 - weight_[celli])*alphaTarget
            + weight_[celli]*alpha_[cellNo];
    }

    // REMEMBER NOT TO PUT correctBoundaryConditions() HERE BUT ONE LEVEL UP
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationSchemes
} // End namespace Foam

// ************************************************************************* //
