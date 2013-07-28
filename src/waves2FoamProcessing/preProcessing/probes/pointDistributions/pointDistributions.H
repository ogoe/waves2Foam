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

Class
    Foam::pointDistributions

Description
    Abstract base class for distribution of points on/inside geometric objects.
    This is used for both the definition of wave gauges and point probes.

SourceFiles
    pointDistributions.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef pointDistributions_H
#define pointDistributions_H

#include "fvMesh.H"
#include "fvCFD.H"
#include "volFields.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "OFstream.H"

#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class porosityCoefficient Declaration
\*---------------------------------------------------------------------------*/

class pointDistributions
{

protected:

    // Protected data
        const fvMesh& mesh_;
        const dictionary& pointDict_;

    // Protected member functions

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointDistributions(const pointDistributions&);

        //- Disallow default bitwise assignment
        void operator=(const pointDistributions&);

public:

    //- Runtime type information
    TypeName("pointDistributions");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointDistributions,
            pointDistributions,
            (
                   const fvMesh& mesh,
                const dictionary& dict
            ),
            ( mesh, dict )
        );

    // Constructors

        //- Construct from components
        pointDistributions
        (
            const fvMesh&,
            const dictionary&
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<pointDistributions> New
        (
            const fvMesh&,
            const dictionary&
        );

    // Destructor

        virtual ~pointDistributions() {};

    // Member Functions

        // Access
        virtual pointField evaluate() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
