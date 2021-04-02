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

#include "polynomialDefault.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(polynomialDefault, 0);
addToRunTimeSelectionTable
(
    gabcSettings,
    polynomialDefault,
    gabcSettings
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


polynomialDefault::polynomialDefault
(
    dictionary& dict,
    scalar depth
)
:
    gabcSettings(dict, depth)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polynomialDefault::set( Ostream& os )
{
    // Write the preProcessMethod for future reference
    this->writeGiven(os, "preProcessMethod");

    // Write the range for future reference
    this->writeGiven(os, "defaultRange");
    this->writeDerived(os, "shapeFunction", "P03Celerity");

    // Write the default coefficients
    word range(dict_.lookup("defaultRange"));

    if (range == "kh000_030")
    {
        writeDerived(os, "a0",  1.6207);
        writeDerived(os, "a1",  0.0000);
        writeDerived(os, "a2", -2.7209);
        writeDerived(os, "a3",  1.5429);
    }
    else if (range == "kh000_050")
    {
        writeDerived(os, "a0",  1.3797);
        writeDerived(os, "a1",  0.0000);
        writeDerived(os, "a2", -0.3806);
        writeDerived(os, "a3", -0.7182);
    }
    else if (range == "kh000_100")
    {
        writeDerived(os, "a0",  1.0613);
        writeDerived(os, "a1",  0.0000);
        writeDerived(os, "a2",  2.5231);
        writeDerived(os, "a3", -3.4319);
    }
    else if (range == "kh030_080")
    {
        writeDerived(os, "a0",  0.7071);
        writeDerived(os, "a1",  0.0000);
        writeDerived(os, "a2",  3.0863);
        writeDerived(os, "a3", -3.6110);
    }
    else if (range == "kh050_100")
    {
        writeDerived(os, "a0",  0.4437);
        writeDerived(os, "a1",  0.0000);
        writeDerived(os, "a2",  3.6696);
        writeDerived(os, "a3", -3.9378);
    }
    else
    {
        FatalErrorIn("void polynomialDefault::set( Ostream& os ") << endl
            << "Default range not found. Valid default ranges are:\n"
            << "(\n"
            << "    kh000_030\n"
            << "    kh000_050\n"
            << "    kh000_100\n"
            << "    kh030_080\n"
            << "    kh050_100\n"
            << ")"
            << exit(FatalError);


    }

    this->writeDepth(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
