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

#include "relaxationWeight.H"

#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationWeights
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationWeight, 0);
defineRunTimeSelectionTable(relaxationWeight, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


relaxationWeight::relaxationWeight
(
    const word& subDictName,
    const fvMesh& mesh
)
:
    IOdictionary
    (
        mesh.thisDb().lookupObject<IOobject>("waveProperties")
    ),

    mesh_(mesh),
    coeffDict_(subDict(subDictName + "Coeffs").subDict("relaxationZone")),
    PI_( M_PI ),

    rwcc_( mesh, coeffDict_ )
{
}


relaxationWeight::~relaxationWeight()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


autoPtr<relaxationWeight> relaxationWeight::New
(
    const word& subDictName,
    const fvMesh& mesh_
)
{
    word relaxationWeightTypeName;

    // Enclose the creation of the waveProp to ensure it is
    // deleted before the relaxationWeight is created otherwise the dictionary
    // is entered in the database twice
    {
        const dictionary coeffDict_
        (
            (mesh_.thisDb().lookupObject<IOdictionary>("waveProperties"))
             .subDict(subDictName + "Coeffs")
             .subDict("relaxationZone")
        );

        relaxationWeightTypeName = coeffDict_.lookupOrDefault<word>
            (
                "relaxationWeight",
                "Exponential"
            );

    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find
        (
            "relaxationWeight"+relaxationWeightTypeName
        );

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "relaxationWeight::New(const word&, const fvMesh&)"
        )   << "Unknown relaxation weight type 'relaxationWeight"
            << relaxationWeightTypeName << "'"<< endl << endl
            << "Valid relaxation weight types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<relaxationWeight>(cstrIter()(subDictName,mesh_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationWeight::weights
(
    const labelList& cells,
    const scalarField& sigma,
    scalarField& weights
)
{
    // First obtain the weights from the specific weight specification
    computeWeights(cells, sigma, weights);

    // Perform the courant number dependent weight correction
    rwcc_.courantCorrection(cells, weights);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationWeights
} // End namespace Foam

// ************************************************************************* //
