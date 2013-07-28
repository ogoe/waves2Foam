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
#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


autoPtr<waveTheory> waveTheory::New
(
    const word& subDictName,
    const fvMesh& mesh_
)
{
    word waveTheoryTypeName;

    // Enclose the creation of the waveProp to ensure it is
    // deleted before the waveTheory is created otherwise the dictionary
    // is entered in the database twice
    {
        const dictionary coeffDict_
        (
            (mesh_.thisDb().lookupObject<IOdictionary>("waveProperties"))
           .subDict(subDictName + "Coeffs")
        );

        coeffDict_.lookup("waveType") >> waveTheoryTypeName;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(waveTheoryTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "waveTheory::New(const word&, const fvMesh&)"
        )   << "Unknown wave theory type " << waveTheoryTypeName
            << endl << endl
            << "Valid wave theories types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<waveTheory>(cstrIter()(subDictName,mesh_));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
