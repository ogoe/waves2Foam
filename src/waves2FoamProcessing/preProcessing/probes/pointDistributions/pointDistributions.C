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

#include "pointDistributions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointDistributions, 0);
defineRunTimeSelectionTable(pointDistributions, pointDistributions);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


pointDistributions::pointDistributions
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),

    pointDict_( dict )
{
}


autoPtr<pointDistributions> pointDistributions::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    word pd( dict.lookup("pointDistribution") );

    pointDistributionsConstructorTable::iterator cstrIter =
            pointDistributionsConstructorTablePtr_->find( pd );

    if (cstrIter == pointDistributionsConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "pointDistributions::New(const dictionary&)"
        )   << "Unknown point distribution: " << pd
            << endl << endl
            << "Valid methods are :" << endl
            << pointDistributionsConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<pointDistributions>(cstrIter()( mesh, dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
