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

#include "relaxationShape.H"

#if OFVERSION != 15
    #include "uniformDimensionedFields.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationShape, 0);
defineRunTimeSelectionTable(relaxationShape, dictionary);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

relaxationShape::relaxationShape
(
	const word & subDictName,
	const fvMesh & mesh
)
:
    IOdictionary
    (
#if OFVERSION == 15
        mesh.db().lookupObject<IOobject>("waveProperties")
#else
        mesh.thisDb().lookupObject<IOobject>("waveProperties")
#endif
    ),

    mesh_(mesh),
    coeffDict_(subDict(subDictName + "Coeffs").subDict("relaxationZone"))
{
	// Takes care of the fact that the gravity vector is defined differently between OF1.5 and OF1.6+
#if OFVERSION == 15
    vector g( dimensionedVector( mesh_.db().lookupObject<IOdictionary>("environmentalProperties").lookup("g") ).value() );
#else
	vector g( uniformDimensionedVectorField( mesh_.thisDb().lookupObject<uniformDimensionedVectorField>("g")).value() );
#endif
    direction_ = g / mag(g);
}


relaxationShape::~relaxationShape()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// ************************************************************************* //
