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

#include "porosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    defineTypeNameAndDebug(porosityModel, 0);
    defineRunTimeSelectionTable(porosityModel, porosityModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::porosityModel::porosityModel
(
	const fvMesh& mesh
)
{
}


Foam::porosityModel::~porosityModel()
{}


autoPtr<porosityModel> porosityModel::New
(
    const fvMesh& mesh
)
{
    word porosityModelTypeName;

    // Enclose the creation of the dictionary to ensure it is deleted before
    // the actual porosity model is created
    {
    	// Not successful with porosityZones, because the reading of the
    	// porous zones properties was made impossible, if an additional
    	// keyword was added.
//        IOdictionary dict
//        (
//        	IOobject
//        	(
//        	    "porosityZones",
//        	    mesh.time().constant(),
//        	    mesh,
//        	    IOobject::MUST_READ,
//        	    IOobject::NO_WRITE
//        	)
//        );
//
//        dict.lookup("porosityModel") >> porosityModelTypeName;

    	mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
    		.lookup("porosityModel") >> porosityModelTypeName;
    }

    porosityModelConstructorTable::iterator cstrIter =
    		porosityModelConstructorTablePtr_->find(porosityModelTypeName);

    if (cstrIter == porosityModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "porosityModel::New(const fvMesh&)"
        )   << "Unknown porosity model of type " << porosityModelTypeName
            << endl << endl
            << "Valid porosity models are :" << endl
            << porosityModelConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<porosityModel>(cstrIter()(mesh));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
