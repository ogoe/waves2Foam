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
#include "relaxationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationSchemes
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<relaxationScheme> relaxationScheme::New
(
    const word& subDictName,
    const fvMesh& mesh,
    vectorField& U,
    scalarField& alpha
)
{
    word relaxationSchemeTypeName;

    // Enclose the creation of the waveProp to ensure it is
    // deleted before the relaxationScheme is created otherwise the dictionary
    // is entered in the database twice
    {
        const dictionary coeffDict_
        (
            (mesh.thisDb().lookupObject<IOdictionary>("waveProperties"))
            .subDict(subDictName + "Coeffs")
            .subDict("relaxationZone")
        );

        coeffDict_.lookup("relaxationScheme") >> relaxationSchemeTypeName;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find
        (
            "relaxationScheme"+relaxationSchemeTypeName
        );

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "relaxationScheme::New"
            "("
            "    const word&, "
            "    const fvMesh&,"
            "    vectorField& U,"
            "    scalarField& alpha"
            ")"
        )   << "Unknown relaxationScheme '" << relaxationSchemeTypeName << "'"
            << endl << endl
            << "Valid relaxation schemes are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<relaxationScheme>(cstrIter()(subDictName,mesh,U,alpha));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationSchemes
} // End namespace Foam

// ************************************************************************* //
