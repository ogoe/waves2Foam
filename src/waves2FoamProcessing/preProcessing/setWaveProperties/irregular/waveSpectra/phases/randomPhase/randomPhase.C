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

#include "randomPhase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(randomPhase, 0);
addToRunTimeSelectionTable(phases, randomPhase, phases);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


randomPhase::randomPhase
(
    const Time& rT,
    dictionary& dict
)
:
    phases(rT, dict)
{
    Info << "\nConstructing: " << this->type() << endl;
    srand(time(NULL));

    if (dict.found("seedForRandomPhase"))
    {
        label seed = readLabel(dict.lookup("seedForRandomPhase"));
        srand(seed);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar randomPhase::phase(const scalar& freq, const vector& k)
{
    return
    (
        2.0*M_PI*static_cast<scalar>(rand())
       /static_cast<scalar>(RAND_MAX)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
