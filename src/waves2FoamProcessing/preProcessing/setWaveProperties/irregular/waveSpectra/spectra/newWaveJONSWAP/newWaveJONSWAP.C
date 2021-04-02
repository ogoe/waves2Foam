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

#include "newWaveJONSWAP.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(newWaveJONSWAP, 0);
addToRunTimeSelectionTable(waveSpectra, newWaveJONSWAP, waveSpectra);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

scalarField newWaveJONSWAP::spectralValue
(
    const scalar& Hs,
    const scalar& Tp,
    const scalar& gamma,
    const scalarField& freq
) const
{
    // Additional parameters
    scalar fp = (1.0/Tp);
    scalar alpha = 0.0624/(0.230 + 0.0336*gamma - 0.185/(1.9 + gamma));

    scalarField sigma(freq.size(), 0.07), beta(freq.size(), 0.0);

    forAll (sigma, ii)
    {
        if (freq[ii] >= fp)
        {
            sigma[ii] = 0.09;
        }
    }

    beta = Foam::exp(- Foam::pow(freq - fp, 2.0)
        /(2*Foam::pow(sigma, 2.0)*Foam::pow(fp, 2.0)));

    // Compute spectrum
    scalarField S(alpha*Foam::pow(Hs,2.0)*Foam::pow(fp,4.0)*Foam::pow(freq,-5.0)
        *Foam::pow(gamma,beta)*Foam::exp(- 5.0/4.0*Foam::pow(fp/freq, 4.0)));

    return S;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


newWaveJONSWAP::newWaveJONSWAP
(
    const Time& rT,
    dictionary& dict,
    vector g,
    scalarField& amp,
    scalarField& freq,
    scalarField& phi,
    vectorField& k
)
:
    waveSpectra(rT, dict, g, amp, freq, phi, k)
{
    Info << "\nConstructing: " << this->type() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


wordList newWaveJONSWAP::list()
{
    wordList res(5);

    res[0] = "etaMax";
    res[1] = "Tp";
    res[2] = "gamma";
    res[3] = "depth";
    res[4] = "direction";

    return res;
}


void newWaveJONSWAP::set(Ostream& os)
{
    // Get the input parameters
    scalar etaMax(readScalar(dict_.lookup("etaMax")));

    scalar Tp(readScalar(dict_.lookup("Tp")));

    scalar gamma(readScalar(dict_.lookup("gamma")));

    scalar depth(readScalar(dict_.lookup("depth")));

    vector direction(vector(dict_.lookup("direction")));

    label N = readLabel(dict_.lookup("N"));

    // Calculate the frequency axis
    autoPtr<Foam::frequencyAxis> fA = Foam::frequencyAxis::New(rT_, dict_);
    scalarField nodeFrequency(N + 1, 0);

    // An intermediate step needed for certain discretisation types
    // Placed in scopes such that the temporary variables to not 'survive'
    {
        equidistantFrequencyAxis equiFA(rT_, dict_);

        scalarField tempFreqAxis = equiFA.freqAxis(10000);
        scalarField tempSpectrum
            = this->spectralValue(1.0, Tp, gamma, tempFreqAxis);

        nodeFrequency = fA->freqAxis(tempFreqAxis, tempSpectrum, N);
    }

    // Prepare variables
    freq_.setSize(N);
    amp_.setSize(N);
    phi_.setSize(N);
    k_.setSize(N);

    // Calculate the spectrum (unit have height)
    scalarField S = this->spectralValue(1.0, Tp, gamma, nodeFrequency);

    // Prepare stokesFirst to compute wave numbers
    Foam::stokesFirstProperties stp(rT_, dict_, g_);

    scalar Stotal(0);

    for (int i = 1; i < S.size(); i++)
    {
        Stotal += (nodeFrequency[i - 1] + nodeFrequency[i])*
                0.5*(S[i - 1] + S[i]);
    }

    // Compute return variables
    for (int i = 1; i < N + 1; i++)
    {
        // The frequency is the mid-point between two nodes
        freq_[i - 1] = 0.5*(nodeFrequency[i - 1] + nodeFrequency[i]);

        // Amplitude is the square root of the trapezoidal integral
        amp_[i - 1] =
            etaMax*(nodeFrequency[i - 1] + nodeFrequency[i])*
                0.5*(S[i - 1] + S[i])/Stotal;

        // Wave number based on linear wave theory
        k_[i - 1] = direction*stp.linearWaveNumber(depth, freq_[i-1]);

        // The phase is computed based on the phase-function
        phi_[i - 1] = phases_->phase(freq_[i - 1], k_[i - 1]);
    }

    writeSpectrum(os, nodeFrequency, S);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
