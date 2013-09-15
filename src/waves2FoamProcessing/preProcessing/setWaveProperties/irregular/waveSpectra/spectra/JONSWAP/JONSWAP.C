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

#include "JONSWAP.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(JONSWAP, 0);
addToRunTimeSelectionTable(waveSpectra, JONSWAP, waveSpectra);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


JONSWAP::JONSWAP
(
    const Time& rT,
    dictionary& dict,
    scalarField& amp,
    scalarField& freq,
    scalarField& phi,
    vectorField& k
)
:
    waveSpectra(rT, dict, amp, freq, phi, k)
{
    Info << "\nConstructing: " << this->type() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


wordList JONSWAP::list()
{
    wordList res(5);

    res[0] = "Hs";
    res[1] = "Tp";
    res[2] = "gamma";
    res[3] = "depth";
    res[4] = "direction";

    return res;
}


void JONSWAP::set(Ostream& os)
{
    // Get the input parameters
    scalar Hs( readScalar(dict_.lookup("Hs")) );
    scalar Tp( readScalar(dict_.lookup("Tp")) );
    scalar gamma( readScalar(dict_.lookup("gamma")) );
    scalar depth( readScalar(dict_.lookup("depth")) );
    vector direction( vector( dict_.lookup("direction")));
    Switch freqEqui(Switch(dict_.lookup("equidistantFrequencyAxis")));

    label  N( static_cast<label>(readScalar(dict_.lookup("N"))) );

    freq_.setSize(N);
    amp_.setSize(N);
    phi_.setSize(N);
    k_.setSize(N);

    N++;

    // Additional parameters
    scalar fp    = ( 1.0/Tp );
    scalar alpha = 0.0624/(0.230 + 0.0336*gamma - 0.185/( 1.9 + gamma) );

    scalar flow( 0.3*fp ), fhigh( 3.0*fp );

    if (dict_.found("lowerFrequencyCutoff"))
    {
        flow = readScalar(dict_.lookup("lowerFrequencyCutoff"));
    }

    if (dict_.found("upperFrequencyCutoff"))
    {
        fhigh = readScalar(dict_.lookup("upperFrequencyCutoff"));
    }

    scalarField f(N, 0.0), sigma(N, 0.07), beta(N, 0.0);

    if (!freqEqui)
    {
    	// Stretched frequency axis
        label Nlow ( ceil( (fp - flow)/( fhigh - fp)*N ) );
        label Nhigh( N - Nlow );

        for (int i=0; i < Nlow; i++)
        {
            f[i] = ( fp - flow )*Foam::sin( 2*PI_/( 4.0*Nlow )*i ) + flow;
        }

        for (int i=0; i<=Nhigh; i++)
        {
            f[Nlow - 1 + i] =
                (fhigh - fp)*(- Foam::cos(2*PI_/(4*Nhigh)*i) + 1) + fp;
        }
    }
    else
    {
    	// Equidistant frequency axis
    	scalar df = (fhigh - flow)/static_cast<scalar>(N-1);
    	for (int i = 0; i < N; i++)
    	{
    		f[i] = static_cast<scalar>(i)*df + flow;
    	}
    }

    forAll (sigma, ii)
    {
        if (f[ii] >= fp)
        {
            sigma[ii] = 0.09;
        }
    }

    beta = Foam::exp(- Foam::pow(f - fp, 2.0)
        /(2*Foam::pow(sigma, 2.0)*Foam::pow(fp, 2.0)));

    // Compute spectrum
    scalarField S = alpha*Foam::pow(Hs,2.0)*Foam::pow(fp,4.0)*Foam::pow(f,-5.0)
        *Foam::pow(gamma,beta)*Foam::exp( - 5.0/4.0*Foam::pow( fp/f, 4.0 ) );

    Foam::stokesFirstProperties stp( rT_, dict_ );

    // Compute return variables
    for (int i = 1; i < N; i++)
    {
        freq_[i - 1] = 0.5*( f[i - 1] + f[i] );
        amp_[i - 1] = Foam::sqrt( ( S[i-1] + S[i] )*( f[i] - f[i - 1] ) );
        k_[i - 1] = direction*stp.linearWaveNumber(depth, freq_[i-1]);

        phi_[i - 1] = phases_->phase(freq_[i - 1], k_[i - 1]);
    }

    if (dict_.lookupOrDefault<Switch>("writeSpectrum",false))
    {
        S.writeEntry("spectramValues", os);
        os << nl;

        f.writeEntry("fspectrum", os);
        os << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
