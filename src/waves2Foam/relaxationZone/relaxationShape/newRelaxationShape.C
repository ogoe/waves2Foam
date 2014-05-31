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

#include "error.H"
#include "relaxationShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<relaxationShape> relaxationShape::New
(
    const word& subDictName,
    const fvMesh& mesh_
)
{
    word relaxationShapeTypeName;

    // Enclose the creation of the waveProp to ensure it is
    // deleted before the relaxationShape is created otherwise the dictionary
    // is entered in the database twice
    {
        const dictionary coeffDict_
        (
            (mesh_.thisDb().lookupObject<IOdictionary>("waveProperties"))
             .subDict(subDictName + "Coeffs")
             .subDict("relaxationZone")
        );

        coeffDict_.lookup("relaxationShape") >> relaxationShapeTypeName;

    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find
        (
            "relaxationShape"+relaxationShapeTypeName
        );

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "relaxationShape::New(const word&, const fvMesh&)"
        )   << "Unknown relaxation shape type 'relaxationShape"
            << relaxationShapeTypeName << "'" << endl << endl
            << "Valid relaxation shape types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<relaxationShape>(cstrIter()(subDictName,mesh_));
}


autoPtr<relaxationShape> relaxationShape::New
(
    const word& subDictName,
    const word relaxationShapeTypeName,
    const fvMesh& mesh_
)
{
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find
        (
            "relaxationShape"+relaxationShapeTypeName
        );

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "relaxationShape::New(const word&, const fvMesh&)"
        )   << "Unknown relaxation shape type 'relaxationShape"
            << relaxationShapeTypeName << "'" << endl << endl
            << "Valid relaxation shape types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<relaxationShape>(cstrIter()(subDictName,mesh_));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// ************************************************************************* //
