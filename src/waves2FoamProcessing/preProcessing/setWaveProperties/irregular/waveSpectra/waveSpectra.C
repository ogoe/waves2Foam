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

#include "waveSpectra.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(waveSpectra, 0);
defineRunTimeSelectionTable(waveSpectra, waveSpectra);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


waveSpectra::waveSpectra
(
    const Time& rT,
    dictionary& dict,
    scalarField& amp,
    scalarField& freq,
    scalarField& phi,
    vectorField& k
)
:
    rT_(rT),
    dict_(dict),
    amp_(amp),
    freq_(freq),
    phi_(phi),
    k_(k),

#if OFPLUSBRANCH==1
    #if OFVERSION<1812
        G_(Foam::mag(uniformDimensionedVectorField
            (
                rT_.db().lookupObject<uniformDimensionedVectorField>("g")
            ).value())),
    #else
        G_(Foam::mag(Foam::meshObjects::gravity::New(rT_).value())),
    #endif
#else
        G_(Foam::mag(uniformDimensionedVectorField
            (
                rT_.db().lookupObject<uniformDimensionedVectorField>("g")
            ).value())),
#endif

    PI_( M_PI ),

    phases_(Foam::phases::New(rT_, dict_))
{
}


waveSpectra::~waveSpectra()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void waveSpectra::writeSpectrum
(
    Ostream& os,
    const scalarField& freq,
    const scalarField& S
) const
{
    if (dict_.subDict("frequencyAxis").lookupOrDefault<Switch>("writeSpectrum",false))
    {
        S.writeEntry("spectramValues", os);
        os << nl;

        freq.writeEntry("fspectrum", os);
        os << nl;
    }
}


autoPtr<waveSpectra> waveSpectra::New
(
    const Time& rT,
    dictionary& dict,
    scalarField& amp,
    scalarField& freq,
    scalarField& phi,
    vectorField& k
)
{
    word spectrumName;
    dict.lookup("spectrum") >> spectrumName;

    waveSpectraConstructorTable::iterator cstrIter =
            waveSpectraConstructorTablePtr_->find(spectrumName);

    if (cstrIter == waveSpectraConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "waveSpectra::New(const fvMesh&, dictionary&, bool)"
        )   << "Unknown wave spectrum '" << spectrumName << "'"
            << endl << endl
            << "Valid wave spectra are:" << endl
            << waveSpectraConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<waveSpectra>(cstrIter()(rT, dict, amp, freq, phi, k));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
