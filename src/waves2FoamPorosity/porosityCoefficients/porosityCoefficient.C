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

#include "porosityCoefficient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(porosityCoefficient, 0);
defineRunTimeSelectionTable(porosityCoefficient, porosityCoefficient);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

porosityCoefficient::porosityCoefficient
(
    const dictionary & poroProp
)
:
    poroProperties_( poroProp ),

    linearCoefficient_( dimensionedVector("null", dimless, vector::zero ) ),

    quadraticCoefficient_( dimensionedVector("null", dimless, vector::zero ) )
{

}


porosityCoefficient::~porosityCoefficient()
{}

autoPtr<porosityCoefficient> porosityCoefficient::New
(
    const dictionary & poroProp
)
{
    word poroForm = poroProp.lookup("resistanceFormulation");

    porosityCoefficientConstructorTable::iterator cstrIter =
            porosityCoefficientConstructorTablePtr_->find(poroForm);

    if (cstrIter == porosityCoefficientConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "porosityCoefficient::New(const dictionary &)"
        )   << "Unknown resistance formulation: " << poroForm
            << endl << endl
            << "Valid methods are :" << endl
            << porosityCoefficientConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<porosityCoefficient>(cstrIter()( poroProp ));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
