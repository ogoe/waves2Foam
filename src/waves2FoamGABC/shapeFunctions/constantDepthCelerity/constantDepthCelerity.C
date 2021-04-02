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

#include "constantDepthCelerity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(constantDepthCelerity, 0);
addToRunTimeSelectionTable
(
    celerityShapeFunctions,
    constantDepthCelerity,
    celerityShapeFunctions
);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantDepthCelerity::constantDepthCelerity
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

    fact_(readScalar(dict.subDict(dictName()).lookup("celerityFactor")))
{
}


constantDepthCelerity::~constantDepthCelerity()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> constantDepthCelerity::rhoCelerity
(
	const fvPatch& p,
    const scalarField& rho,
    const scalarField& alpha,
    const scalar& g
)
{
    tmp<scalarField> tres(new scalarField);
    scalarField& res = tres.ref();

    res.setSize(rho.size(), 0.0);

    res = fact_*rho*Foam::sqrt(g*depth_);

    return tres;
}


void constantDepthCelerity::rhoCelerity
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
