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

#include "rationalDefault.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rationalDefault, 0);
addToRunTimeSelectionTable
(
    gabcSettings,
    rationalDefault,
    gabcSettings
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


rationalDefault::rationalDefault
(
    dictionary& dict,
    scalar depth
)
:
    gabcSettings(dict, depth)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void rationalDefault::set( Ostream& os )
{
    // Write the preProcessMethod for future reference
    this->writeGiven(os, "preProcessMethod");

    // Write the range for future reference
    this->writeGiven(os, "defaultRange");
    this->writeDerived(os, "shapeFunction", "rationalCelerity");

    // Write the default coefficients
    word range(dict_.lookup("defaultRange"));

    if (range == "kh000_030")
    {
        writeDerived(os, "a0",  1.5867);
        writeDerived(os, "a1",  0.8559);
        writeDerived(os, "a2", -1.4986);
        writeDerived(os, "a3",  0.5394);
        writeDerived(os, "a4",  0.7049);
    }
    else if (range == "kh000_050")
    {
        writeDerived(os, "a0",  1.5805);
        writeDerived(os, "a1", -0.6535);
        writeDerived(os, "a2", -0.5843);
        writeDerived(os, "a3", -0.4135);
        writeDerived(os, "a4",  0.4502);
    }
    else if (range == "kh000_100")
    {
        writeDerived(os, "a0",  1.3143);
        writeDerived(os, "a1", -1.2931);
        writeDerived(os, "a2",  0.0064);
        writeDerived(os, "a3", -0.9839);
        writeDerived(os, "a4",  0.1758);

    }
    else
    {
        FatalErrorIn("void rationalDefault::set( Ostream& os ") << endl
            << "Default range not found. Valid default ranges are:\n"
            << "(\n"
            << "    kh000_030\n"
            << "    kh000_050\n"
            << "    kh000_100\n"
            << ")"
            << exit(FatalError);


    }

    this->writeDepth(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
