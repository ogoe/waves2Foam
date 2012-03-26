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

#include "irregularProperties.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(irregularProperties, 0);
addToRunTimeSelectionTable(setWaveProperties, irregularProperties, setWaveProperties);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

irregularProperties::irregularProperties
(
	const fvMesh & mesh,
	dictionary & dict,
	bool write
)
:
	setWaveProperties(mesh, dict, write),
	mesh_(mesh)
{
	Info << "\nConstructing: " << this->type() << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void irregularProperties::ws
(
	const string & fileName,
	const string & propName,
	const scalarField & field
) const
{
	IOField<scalar> output
	(
		IOobject
		(
			fileName+propName,
			"constant",
			"additionalWaveProperties",
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		field
	);
	output.write();
}

void irregularProperties::wv
(
	const string & fileName,
	const string & propName,
	const vectorField & field
) const
{
	IOField<vector> output
	(
		IOobject
		(
			fileName+propName,
			"constant",
			"additionalWaveProperties",
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		field
	);
	output.write();
}



void irregularProperties::set()
{
	scalarField amp(0);
	scalarField frequency(0);
	scalarField phaselag(0);
	vectorField waveNumber(0);

	autoPtr<waveSpectra> spectra( waveSpectra::New(mesh_, dict_, amp, frequency, phaselag, waveNumber) );

	spectra->set();

	if ( write_ )
	{
		mkDir("constant/additionalWaveProperties");

		string fileName = (dict_.name()).name().replace("waveProperties::","");

		ws(fileName, "Amplitude", amp);

		ws(fileName, "Frequency", frequency);

		ws(fileName, "Phaselag", phaselag);

		wv(fileName, "WaveNumber", waveNumber);

//		IOField<scalar> amplitude
//		(
//			IOobject
//			(
//				fileName+"Amplitude",
//				"constant",
//				"additionalWaveProperties",
//				mesh_,
//				IOobject::NO_READ,
//				IOobject::NO_WRITE
//			),
//			amp
//		);
//		amplitude.write();
//
//		IOField<scalar> freq
//		(
//			IOobject
//			(
//				fileName+"Frequency",
//				"constant",
//				"additionalWaveProperties",
//				mesh_,
//				IOobject::NO_READ,
//				IOobject::NO_WRITE
//			),
//			frequency
//		);
//		freq.write();
//
//		IOField<scalar> pl
//		(
//			IOobject
//			(
//				fileName+"Phaselag",
//				"constant",
//				"additionalWaveProperties",
//				mesh_,
//				IOobject::NO_READ,
//				IOobject::NO_WRITE
//			),
//			phaselag
//		);
//		pl.write();
//
//		IOField<vector> wn
//		(
//			IOobject
//			(
//				fileName+"WaveNumber",
//				"constant",
//				"additionalWaveProperties",
//				mesh_,
//				IOobject::NO_READ,
//				IOobject::NO_WRITE
//			),
//			waveNumber
//		);
//		wn.write();
	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
