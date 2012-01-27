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
		// Amplitude string
		std::stringstream samp;

		samp << "nonuniform List<scalar> " << amp.size() << "(";
		forAll(amp, index)
		{
			samp << amp[index] << " ";
		}
		samp << ")";

		// Frequency string
		std::stringstream sfreq;

		sfreq << "nonuniform List<scalar> " << frequency.size() << "(";
		forAll(frequency, index)
		{
			sfreq << 2.0 * PI_ * frequency[index] << " ";
		}
		sfreq << ")";

		// Phaselag string
		std::stringstream sphi;

		sphi << "nonuniform List<scalar> " << phaselag.size() << "(";
		forAll(phaselag, index)
		{
			sphi << phaselag[index] << " ";
		}
		sphi << ")";

		// waveNumber string
		std::stringstream sk;

		sk << "nonuniform List<vector> " << waveNumber.size() << "(";
		forAll(waveNumber, index)
		{
			sk << "(" << waveNumber[index].x() << " " << waveNumber[index].y() << " " << waveNumber[index].z() << ") ";
		}
		sk << ")";

		// Write strings to dictionary
		dict_.add("amplitude" , samp.str() , true);
		dict_.add("omega"     , sfreq.str(), true);
		dict_.add("phi"       , sphi.str() , true);
		dict_.add("waveNumber", sk.str()   , true);
	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
