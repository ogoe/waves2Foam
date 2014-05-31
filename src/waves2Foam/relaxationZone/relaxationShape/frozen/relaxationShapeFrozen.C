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

#include "relaxationShapeFrozen.H"
#include "addToRunTimeSelectionTable.H"

#include "relaxationShapeRectangular.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationShapeFrozen, 0);
addToRunTimeSelectionTable
(
    relaxationShape,
    relaxationShapeFrozen,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


relaxationShapeFrozen::relaxationShapeFrozen
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    relaxationShape(subDictName, mesh_)
{
    // Everything is placed in the constructor in order to avoid any updating
    // the relaxation zone coefficients during mesh motion.
    if (mesh_.time().timeName() == "0")
    {
        word actualType(coeffDict_.lookup("actualRelaxationShape"));

        if ("relaxationShape"+actualType == typeName)
        {
            FatalErrorIn("relaxationShapeFrozen::relaxationShapeFrozen")
                << "The chosen type of the relaxation shape type "
                << "'actualRelaxationShape' is the calling type.\n"
                << "This is not allowed.\n\n"
                << endl << exit(FatalError);
        }

        autoPtr<Foam::relaxationShapes::relaxationShape> rsr
            = Foam::relaxationShapes::relaxationShape
              ::New(subDictName, actualType, mesh_);

        const labelList& recCells = rsr->cells();
        const scalarField& recSigma = rsr->sigma();

        cells_.setSize(recCells.size(), -1);
        cells_ = recCells;

        sigma_.setSize(recSigma.size(), -1);
        sigma_ = recSigma;

        // Writing the cell indices to prepare for restart of computation
        IOField<label> cellsOutput
        (
            IOobject
            (
                "cells_" + subDictName,
                mesh_.time().constant(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );

        cellsOutput.setSize(cells_.size());
        forAll (cells_, celli)
        {
            cellsOutput[celli] = cells_[celli];
        }

        cellsOutput.write();

        IOField<scalar> sigmaOutput
        (
            IOobject
            (
                "sigma_" + subDictName,
                mesh_.time().constant(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );

        sigmaOutput.setSize(sigma_.size());
        sigmaOutput = sigma_;

        sigmaOutput.write();
    }
    else
    {
        IOField<label> cellsInput
        (
            IOobject
            (
                "cells_" + subDictName,
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        cells_.setSize(cellsInput.size());
        cells_ = cellsInput;

        IOField<scalar> sigmaInput
        (
            IOobject
            (
                "sigma_" + subDictName,
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        sigma_.setSize(sigmaInput.size());
        sigma_ = sigmaInput;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationShapeFrozen::findComputationalCells()
{
    // Intentionally left blank
}


void relaxationShapeFrozen::computeSigmaCoordinate()
{
    // Intentionally left blank
}


const pointField& relaxationShapeFrozen::pointSet()
{
    notImplemented("pointSet is not implemented for this shape");
}


scalar relaxationShapeFrozen::interpolation
(
    const scalarField& source,
    const point& p0
) const
{
    notImplemented("interpolation is not implemented for this shape");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// ************************************************************************* //
