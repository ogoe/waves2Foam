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

#include "powerSpectraLS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(powerSpectraLS, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    powerSpectraLS,
    postProcessingWaves
);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


void powerSpectraLS::evaluateScalar()
{
    Info << "        - Power spectra computed for scalar quantities" << endl;

    spectralMethodsLeastSquaresBased smls( rT_, actionProperties_ );

    List<scalarField> input = readScalarFields( indices_ );

    scalarField time = readIOScalarField( callName_ + "_time" );

    List<scalarField> spectra =
        smls.powerSpectra(time, input, N_, 1.0 /period_);

    scalarField frequencies = smls.frequencies( N_ );

    writeScalar( frequencies, spectra);
}


void powerSpectraLS::writeScalar
(
    const scalarField& frequencies,
    const List<scalarField>& spectra
)
{
    Info << "        - Writing computed spectra to: " << directDir_.c_str()
         << this->type() << endl;

    mkDir( directDir_ + this->type() );

    autoPtr<OFstream> spectrumPtr_;

    forAll (indices_, indexi)
    {
        std::stringstream ss;
        ss << callName_ << "_" << indices_[indexi];

        spectrumPtr_.reset
        (
            new OFstream
            (
                directDir_ + "/" + this->type() + "/"
              + ss.str() + "_spectrum.dat"
            )
        );

        const scalarField& spectrum( spectra[indexi] );

        forAll (frequencies, freqi)
        {
            spectrumPtr_() << frequencies[freqi] << tab << spectrum[freqi]
                           << endl;
        }
    }
}


void powerSpectraLS::evaluateVector()
{
    Info << "        - Power spectra computed for vector quantities" << endl;

    spectralMethodsLeastSquaresBased smls( rT_, actionProperties_ );

    List<vectorField> input = readVectorFields( indices_ );

    scalarField time = readIOScalarField( callName_ + "_time" );

    List<vectorField> spectra =
        smls.powerSpectra(time, input, N_, 1.0/period_);

    scalarField frequencies = smls.frequencies( N_ );

    writeVector( frequencies, spectra);
}


void powerSpectraLS::writeVector
(
    const scalarField& frequencies,
    const List<vectorField>& spectra
)
{
    Info << "        - Writing computed spectra to: " << directDir_.c_str()
         << this->type() << endl;

    mkDir( directDir_ + this->type() );

    autoPtr<OFstream> spectrumPtr_;

    forAll (indices_, indexi)
    {
        std::stringstream ss;
        ss << callName_ << "_" << indices_[indexi];

        spectrumPtr_.reset
        (
            new OFstream
            (
                directDir_ + "/" + this->type() + "/"
              + ss.str() + "_spectrum.dat"
            )
        );

        const vectorField& spectrum( spectra[indexi] );

        forAll (frequencies, freqi)
        {
            spectrumPtr_() << frequencies[freqi]  << tab
                           << spectrum[freqi].x() << tab
                           << spectrum[freqi].y() << tab
                           << spectrum[freqi].z() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


powerSpectraLS::powerSpectraLS
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action ),

#   include "../../../dataDict.H"
    ,

    N_( readLabel( actionProperties_.lookup("nFreq") ) ),

    period_( readScalar( actionProperties_.lookup("period") ))
{
    readIndices( dataDict_, indices_ );
}


powerSpectraLS::~powerSpectraLS()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void powerSpectraLS::evaluate()
{
    if (dataType() == "scalar")
    {
        evaluateScalar();
    }
    else if (dataType() == "vector")
    {
        evaluateVector();
    }
    else
    {
        notImplemented("Only scalars and vectors are supported.");
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
