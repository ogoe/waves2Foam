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

#include "oceanWave3DProperties.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(oceanWave3DProperties, 0);
addToRunTimeSelectionTable
(
    setWaveProperties,
    oceanWave3DProperties,
    setWaveProperties
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


oceanWave3DProperties::oceanWave3DProperties
(
    const Time& rT,
    dictionary& dict,
    bool write
)
:
    setWaveProperties(rT, dict, write)
{
    Info << "\nConstructing: " << this->type() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void oceanWave3DProperties::set( Ostream& os )
{
    // Write the beginning of the sub-dictionary
    writeBeginning( os );

    // Write the already given parameters
    writeGiven( os, "waveType" );

    os << endl;

    // Write the intervals
    writeGiven(os, "nIntervals");

    scalarField startTimes("startTimes", dict_, readLabel(dict_.lookup("nIntervals")));
    writeDerived(os, "startTimes", startTimes);

    scalarField endTimes("endTimes", dict_, readLabel(dict_.lookup("nIntervals")));
    writeDerived(os, "endTimes", endTimes);

    os << endl;

    // Write the ramp information
    writeGiven(os, "rampInterval");

    if (Switch(dict_.lookup("rampInterval")))
    {
    	writeGiven(os, "Tsoft");
    }

    os << endl;

    // Write the information for the mapping
    writeGiven(os, "mappingZone");

    // Set the translation of the OF-mesh
    writeGiven(os, "translateOpenFoamMesh");

    // Write the closing bracket
    writeEnding( os );
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
