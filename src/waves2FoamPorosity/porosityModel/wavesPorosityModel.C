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

#include "wavesPorosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    defineTypeNameAndDebug(wavesPorosityModel, 0);
    defineRunTimeSelectionTable(wavesPorosityModel, wavesPorosityModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::wavesPorosityModel::wavesPorosityModel
(
	const fvMesh& mesh
)
:
    porosity_
    (
        IOobject
        (
            "porosity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("NULL", dimless, 1.0),
        "zeroGradient"
    )
{
}


Foam::wavesPorosityModel::~wavesPorosityModel()
{}


autoPtr<wavesPorosityModel> wavesPorosityModel::New
(
    const fvMesh& mesh
)
{
    word wavesPorosityModelTypeName;

    // Enclose the creation of the dictionary to ensure it is deleted before
    // the actual porosity model is created
    {
        if (mesh.thisDb().foundObject<IOdictionary>("waveProperties"))
        {
        	mesh.thisDb().lookupObject<IOdictionary>("waveProperties")
        		.lookup("porosityModel") >> wavesPorosityModelTypeName;
        }
        else
        {
        	IOdictionary wp
        	(
        	    IOobject
        	    (
        	        "waveProperties",
        	        mesh.time().constant(),
        	        mesh,
        	        IOobject::MUST_READ,
        	        IOobject::NO_WRITE
        	    )
        	);

        	wp.lookup("porosityModel") >> wavesPorosityModelTypeName;
        }
    }

    wavesPorosityModelConstructorTable::iterator cstrIter =
    		wavesPorosityModelConstructorTablePtr_->find(wavesPorosityModelTypeName);

    if (cstrIter == wavesPorosityModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "wavesPorosityModel::New(const fvMesh&)"
        )   << "Unknown porosity model of type " << wavesPorosityModelTypeName
            << endl << endl
            << "Valid porosity models are :" << endl
            << wavesPorosityModelConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<wavesPorosityModel>(cstrIter()(mesh));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
