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

#include "concatenatedDistribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(concatenatedDistribution, 0);
addToRunTimeSelectionTable
(
    pointDistributions,
    concatenatedDistribution,
    pointDistributions
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


concatenatedDistribution::concatenatedDistribution
(
//    const fvMesh& mesh,
    const dictionary& dict
)
:
    pointDistributions(dict),

    ppS_(0),
    ppE_(0)
{
}


concatenatedDistribution::~concatenatedDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


pointField concatenatedDistribution::evaluateStart()
{
    ppS_.setSize(0);

    wordList toc(pointDict_.toc());

    forAll (toc, item)
    {
        word name(toc[item]);

        if (pointDict_.isDict(name))
        {
            // Get the point distribution
            const dictionary subDict(pointDict_.subDict(name));

            autoPtr<Foam::pointDistributions> pd
                (
                    Foam::pointDistributions::New(subDict)
                );

            pointField pp(pd->evaluateStart());

            // Update the size of the ppS_
            label N = ppS_.size();
            ppS_.setSize(N + pp.size());

            // Insert all points
            forAll (pp, pointi)
            {
                ppS_[pointi + N] = pp[pointi];
            }
        }
    }

    return ppS_;
}


pointField concatenatedDistribution::evaluateEnd()
{
    ppE_.setSize(0);

    wordList toc(pointDict_.toc());

    forAll (toc, item)
    {
        word name(toc[item]);

        if (pointDict_.isDict(name))
        {
            // Get the point distribution
            const dictionary subDict(pointDict_.subDict(name));

            autoPtr<Foam::pointDistributions> pd
            (
                    Foam::pointDistributions::New(subDict)
            );

            pointField pp(pd->evaluateEnd());

            // Update the size of the ppE_
            label N = ppE_.size();
            ppE_.setSize(N + pp.size());

            // Insert all points
            forAll (pp, pointi)
            {
                ppE_[pointi + N] = pp[pointi];
            }
        }
    }

    return ppE_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
