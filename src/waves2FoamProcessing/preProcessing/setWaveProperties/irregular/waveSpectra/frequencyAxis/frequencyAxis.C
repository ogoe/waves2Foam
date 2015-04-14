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

#include "frequencyAxis.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(frequencyAxis, 0);
defineRunTimeSelectionTable(frequencyAxis, frequencyAxis);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


frequencyAxis::frequencyAxis
(
    const Time& rT,
    dictionary& dict
)
:
    rT_(rT),
    dict_(dict)
{
    if (dict_.subDict("frequencyAxis").found("lowerFrequencyCutoff"))
    {
        fl_ = readScalar
            (
                dict_.subDict("frequencyAxis").lookup("lowerFrequencyCutoff")
            );
    }
    else
    {
        scalar Tp = readScalar(dict_.lookup("Tp"));
        fl_ = 0.3/Tp;
    }

    if (dict_.subDict("frequencyAxis").found("upperFrequencyCutoff"))
    {
        fu_ = readScalar
            (
                dict_.subDict("frequencyAxis").lookup("upperFrequencyCutoff")
            );
    }
    else
    {
        scalar Tp = readScalar(dict_.lookup("Tp"));
        fu_ = 3.0/Tp;
    }
}


frequencyAxis::~frequencyAxis()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


autoPtr<frequencyAxis> frequencyAxis::New
(
    const Time& rT,
    dictionary& dict
)
{
    word discretisation = dict.subDict("frequencyAxis").lookup("discretisation");

    frequencyAxisConstructorTable::iterator cstrIter =
            frequencyAxisConstructorTablePtr_->find(discretisation);

    if (cstrIter == frequencyAxisConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "frequencyAxis::New(const Time&, dictionary&)"
        )   << "Unknown discretisation method '" << discretisation << "'"
            << endl << endl
            << "Valid discretisation methods are:" << endl
            << frequencyAxisConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<frequencyAxis>(cstrIter()(rT, dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
