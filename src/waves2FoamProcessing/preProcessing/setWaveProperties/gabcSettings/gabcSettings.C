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

#include "gabcSettings.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(gabcSettings, 0);
defineRunTimeSelectionTable(gabcSettings, gabcSettings);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void gabcSettings::lineFormatting( Ostream& os, const word& key)
{
    os << indent << key << token::SPACE;

    for (int i=key.size(); i<Nspaces_-1; i++)
    {
        os << token::SPACE;
    }
}


void gabcSettings::addITstream
(
    Ostream& os,
    const word& key,
    const ITstream& it
)
{
    lineFormatting(os, key);

    forAll (it, ii)
    {
        os << it[ii];

        if (ii < it.size() - 1)
        {
            os << token::SPACE;
        }
    }

    os << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


void gabcSettings::writeGiven( Ostream& os, word name )
{
    ITstream it( dict_.lookup(name) );

    addITstream( os, name, it );
}
//
//
void gabcSettings::writeDerived( Ostream& os, word name, scalar val)
{
    lineFormatting(os, name);

    os << val << token::END_STATEMENT << nl;
}


void gabcSettings::writeDerived( Ostream& os, word name, word type)
{
    lineFormatting(os, name);

    os << type << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


gabcSettings::gabcSettings
(
    dictionary& dict,
    scalar depth
)
:
    dict_(dict),
    depth_(depth)
{
    Nspaces_ = 20;
}


gabcSettings::~gabcSettings()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


autoPtr<gabcSettings> gabcSettings::New
(
    dictionary& dict,
    scalar depth
)
{
    word gabcTypeName;
    dict.lookup("preProcessMethod") >> gabcTypeName;

#if OFVERSION < 2206
    gabcSettingsConstructorTable::iterator cstrIter =
        gabcSettingsConstructorTablePtr_->find
        (
            gabcTypeName
        );

    if (cstrIter == gabcSettingsConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "gabcSettings::New(dictionary&)"
        )
        << "Unknown GABC method " << gabcTypeName
        << endl << endl
        << "Valid GABC methods are :" << endl
        << gabcSettingsConstructorTablePtr_->toc()
        << exit(FatalError);
    }

    return autoPtr<gabcSettings>(cstrIter()(dict, depth));
#else
    auto* cstrIter = gabcSettingsConstructorTable
        (
            gabcTypeName
        );

    if (!cstrIter)
    {
        FatalErrorIn
        (
            "gabcSettings::New(dictionary&)"
        )
        << "Unknown GABC method " << gabcTypeName
        << endl << endl
        << "Valid GABC methods are :" << endl
        << gabcSettingsConstructorTablePtr_->toc()
        << exit(FatalError);
    }

    return autoPtr<gabcSettings>(cstrIter(dict, depth));
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
