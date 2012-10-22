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

#include "setWaveProperties.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(setWaveProperties, 0);
defineRunTimeSelectionTable(setWaveProperties, setWaveProperties);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void setWaveProperties::writeBeginning( Ostream & os)
{
	wordList names( dict_.name().components('/') );

	fileName ends( names[names.size() -1]);

	word dictName( ends.components(':')[1] );

	os << dictName << nl << token::BEGIN_BLOCK << incrIndent << nl;
}

void setWaveProperties::writeGiven( Ostream & os )
{
	wordList toc( dict_.toc() );

	forAll( toc, item )
	{
		if ( !dict_.isDict( toc[item]) )
		{
			ITstream it( dict_.lookup( toc[item] ) );
			addITstream( os, toc[item], it);
		}
	}
}

void setWaveProperties::writeGiven( Ostream & os, word name, label N)
{
	ITstream it( dict_.lookup(name) );

	addITstream( os, name, it, N);
}

void setWaveProperties::writeRelaxationZone( Ostream & os )
{
	word rl("relaxationZone");

	if ( dict_.found( rl ) )
	{
		os << nl << indent << rl << nl << indent << token::BEGIN_BLOCK << incrIndent << nl;

		dictionary sd( dict_.subDict(rl) );

		wordList toc( sd.toc() );

		forAll( toc, item )
		{
			ITstream it( sd.lookup(toc[item]) );

			addITstream( os, toc[item], it );
		}

		os << decrIndent << indent << token::END_BLOCK << endl;
	}


}

void setWaveProperties::writeEnding( Ostream & os )
{
	os  << decrIndent << token::END_BLOCK << nl << endl;
}

void setWaveProperties::addITstream( Ostream & os, const word & key, const ITstream & it, label N )
{
	os << indent << key;

	for( int i=0; i<N; i++)
		os << tab;

	forAll( it, ii)
	{
		os << it[ii];

		if ( ii < it.size() - 1)
			os << " ";
	}

	os << ";" << nl;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

setWaveProperties::setWaveProperties
(
	const fvMesh & mesh,
	dictionary & dict,
	bool write
)
:
	write_(write),
	dict_(dict),

// Takes care of the fact that the gravity vector is defined differently between OF1.5 and OF1.6+
	g_( uniformDimensionedVectorField( mesh.thisDb().lookupObject<uniformDimensionedVectorField>("g")).value() )
{
	G_  = Foam::mag(g_);
	PI_ = M_PI;
}


setWaveProperties::~setWaveProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

autoPtr<setWaveProperties> setWaveProperties::New
(
	const fvMesh & mesh,
	dictionary & dict,
	bool write
)
{
    word waveTheoryTypeName;
    dict.lookup("waveType") >> waveTheoryTypeName;

    setWavePropertiesConstructorTable::iterator cstrIter =
    		setWavePropertiesConstructorTablePtr_->find(waveTheoryTypeName+"Properties");

    if (cstrIter == setWavePropertiesConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "setWaveProperties::New(const fvMesh &, dictionary &, bool)"
        )   << "Unknown wave property type " << waveTheoryTypeName << "Properties"
            << endl << endl
            << "Valid wave property types are :" << endl
            << setWavePropertiesConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<setWaveProperties>(cstrIter()(mesh, dict, write));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
