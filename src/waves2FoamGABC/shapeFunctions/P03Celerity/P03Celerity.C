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

#include "P03Celerity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(P03Celerity, 0);
addToRunTimeSelectionTable
(
    celerityShapeFunctions,
    P03Celerity,
    celerityShapeFunctions
);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

scalar P03Celerity::func(const scalar& z) const
{
    return a0_ + a1_*z + a2_*Foam::pow(z, 2.0) + a3_*Foam::pow(z, 3.0);
}


scalarField P03Celerity::func(const scalarField& z) const
{
    scalarField res(z.size(), 0.0);

    forAll (res, facei)
    {
    	res[facei] = func(z[facei]);
    }

    return res;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

P03Celerity::P03Celerity
(
    const fvMesh& mesh,
    const word& patchName,
    const dictionary& dict,
    const vector& vertVec,
    const scalar& seaLevel
)
:
    celerityShapeFunctions(mesh, patchName, dict, vertVec, seaLevel),

    depth_(returnDepth()),

    a0_(readScalar(dict.subDict(dictName()).lookup("a0"))),

    a1_(readScalar(dict.subDict(dictName()).lookup("a1"))),

    a2_(readScalar(dict.subDict(dictName()).lookup("a2"))),

    a3_(readScalar(dict.subDict(dictName()).lookup("a3")))
{
}


P03Celerity::~P03Celerity()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> P03Celerity::rhoCelerity
(
	const fvPatch& p,
    const scalarField& rho,
    const scalarField& alpha,
    const scalar& g
)
{
	// Prepare the return field
    tmp<scalarField> tres(new scalarField);
    scalarField& res = tres.ref();

    // Initialise with the shape function calculated at z/h + 1 = 1
    // (still water level)
    res.setSize(rho.size(), func(1));

    // Get the vertical axis with z/h = -1.0 at the bottom
    scalarField z = zCoordinateFace();
    z = (z - seaLevel_)/depth_;

    // Evaluate shape function over entire depth and truncate at z/h=0
    scalarField S = func(z + 1);

    res = Foam::pos(z)*res + Foam::neg(z)*S;

    // Produce the dimensionally correct value
    res *= rho*Foam::sqrt(g*depth_);

    return tres;
}


void P03Celerity::rhoCelerity
(
    const fvPatch& p,
    const scalarField& rho,
    const scalarField& alpha,
    const scalar& g,
    scalarField& rhoCSqr,
    scalarField& rhoC,
    scalarField& C
)
{
    // Get the 2DV version of the propagation speed
    tmp<scalarField> ttemp = this->rhoCelerity(p, rho, alpha, g);
    const scalarField& temp(ttemp());

    // Set size of all return fields
    C.setSize(temp.size());
    rhoC.setSize(temp.size());
    rhoCSqr.setSize(temp.size());

    // Calculate the return fields
    C = temp/rho;
    rhoC = temp;
    rhoCSqr = rhoC*C;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
