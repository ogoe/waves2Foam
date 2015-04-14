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

#include "cosineStretchedFrequencyAxis.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(cosineStretchedFrequencyAxis, 0);
addToRunTimeSelectionTable(frequencyAxis, cosineStretchedFrequencyAxis, frequencyAxis);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


cosineStretchedFrequencyAxis::cosineStretchedFrequencyAxis
(
    const Time& rT,
    dictionary& dict
)
:
    frequencyAxis(rT, dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalarField cosineStretchedFrequencyAxis::freqAxis
(
    const scalarField&,
    const scalarField&,
    const label& N
) const
{
    scalarField freq(N + 1, 0.0);

    scalar fp = 1.0/readScalar(dict_.lookup("Tp"));

    // Calculate the number of upper and lower frequencies
    label Nlow(ceil((fp - fl_)/(fu_ - fp)*(N + 1)));
    label Nhigh(N + 1 - Nlow);

    for (int i = 0; i < Nlow; i++)
    {
        freq[i] = (fp - fl_)*Foam::sin(2*M_PI/(4.0*Nlow)*i) + fl_;
    }

    for (int i = 0; i <= Nhigh; i++)
    {
        freq[Nlow - 1 + i] =
            (fu_ - fp)*(- Foam::cos(2*M_PI/(4*Nhigh)*i) + 1) + fp;
    }

    return freq;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
