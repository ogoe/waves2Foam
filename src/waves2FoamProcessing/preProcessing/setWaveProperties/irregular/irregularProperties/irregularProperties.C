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
addToRunTimeSelectionTable
(
    setWaveProperties,
    irregularProperties,
    setWaveProperties
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


irregularProperties::irregularProperties
(
    const Time& rT,
    dictionary& dict,
    bool write
)
:
    setWaveProperties(rT, dict, write),
    rT_(rT)
{
    Info << "\nConstructing: " << this->type() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void irregularProperties::set( Ostream& os )
{
    // Write the beginning of the sub-dictionary
    writeBeginning( os );

    // Write the already given parameters
    writeGiven( os, "waveType" );

    writeGiven( os, "spectrum");

    writeGiven( os, "N" );

    writeGiven( os, "Tsoft");

    if (dict_.found("writeSpectrum" ))
    {
        writeGiven( os, "writeSpectrum");
    }

    // Make a pointer to the spectral theory
    scalarField amp(0);
    scalarField frequency(0);
    scalarField phaselag(0);
    vectorField waveNumber(0);

    autoPtr<waveSpectra> spectra
    (
        waveSpectra::New(rT_, dict_, amp, frequency, phaselag, waveNumber)
    );

    // Write properties specific to chosen spectral theory
    wordList specificInput( spectra->list() );

    forAll (specificInput, speci)
    {
        writeGiven( os, specificInput[speci] );
    }

    // Computing the spectral quantities
    spectra->set( os );

    if (write_)
    {
        writeDerived( os, "amplitude", amp);
        writeDerived( os, "frequency", frequency);
        writeDerived( os, "phaselag", phaselag);
        writeDerived( os, "waveNumber", waveNumber);
    }

    // Write the frequency axis information
    if (dict_.found("frequencyAxis"))
    {
        word fa("frequencyAxis");

        os << nl << indent << fa << nl << indent << token::BEGIN_BLOCK
           << incrIndent << nl;

        dictionary sd(dict_.subDict(fa));

        wordList toc( sd.toc() );

        forAll (toc, item)
        {
            ITstream it(sd.lookup(toc[item]));

            addITstream(os, toc[item], it);
        }

        os << decrIndent << indent << token::END_BLOCK << endl;
    }

    // Write the relaxation zone
    writeRelaxationZone( os );

    // Write the closing bracket
    writeEnding( os );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
