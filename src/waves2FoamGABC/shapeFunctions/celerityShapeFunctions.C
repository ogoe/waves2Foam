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

#include "celerityShapeFunctions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(celerityShapeFunctions, 0);
defineRunTimeSelectionTable(celerityShapeFunctions, celerityShapeFunctions);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


word celerityShapeFunctions::dictName() const
{
    return "GABCSettings";
}


scalarField celerityShapeFunctions::zCoordinateFace() const
{
	// Get the face centre field on the boundary and calculate the vertical
	// coordinate
	label patchID = mesh_.boundaryMesh().findPatchID(patchName_);
	const vectorField cf = mesh_.Cf().boundaryField()[patchID];

	scalarField zCoor = (cf & vertVector_);

	return zCoor;
}


scalar celerityShapeFunctions::returnDepth() const
{
    return readScalar(dict_.subDict("GABCSettings").lookup("depth"));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

celerityShapeFunctions::celerityShapeFunctions
(
    const fvMesh& mesh,
    const word& patchName,
    const dictionary& dict,
    const vector& vertVector,
    const scalar& seaLevel
)
:
    mesh_(mesh),

    patchName_(patchName),

    dict_(dict),

    seaLevel_(seaLevel)
{
	vector temp = Foam::cmptMultiply(vertVector, vertVector);
	vertVector_ = temp/Foam::mag(temp);
}


celerityShapeFunctions::~celerityShapeFunctions()
{}


autoPtr<celerityShapeFunctions> celerityShapeFunctions::New
(
    const fvMesh& mesh,
    const word& patchName,
    const dictionary& dict,
    const vector& vert,
    const scalar& seaLevel
)
{
	const dictionary& sd = dict.subDict("GABCSettings");

    word formulation;
    sd.lookup("shapeFunction") >> formulation;

    celerityShapeFunctionsConstructorTable::iterator cstrIter =
            celerityShapeFunctionsConstructorTablePtr_->find(formulation);

    if (cstrIter == celerityShapeFunctionsConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "celerityShapeFunctions::New(const fvMesh &, ...)"
        )   << "Unknown formulation for celerity shape function: " << formulation
            << endl << endl
            << "Valid methods are :" << endl
            << celerityShapeFunctionsConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<celerityShapeFunctions>
        (
            cstrIter()(mesh, patchName, dict, vert, seaLevel)
        );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
