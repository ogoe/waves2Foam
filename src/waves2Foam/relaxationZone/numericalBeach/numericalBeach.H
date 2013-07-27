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
    Foam::numericalBeaches::numericalBeach

Description
    An abstract base class for different types of numerical beaches

SourceFiles
    numericalBeach.C
    newNumericalBeach.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef numericalBeach_H
#define numericalBeach_H

#include "volFields.H"

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

//#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace numericalBeaches
{

/*---------------------------------------------------------------------------*\
                      Class numericalBeach Declaration
\*---------------------------------------------------------------------------*/

class numericalBeach
:
    public IOdictionary
{

protected:

    // Protected data
        const fvMesh& mesh_;

        const dictionary coeffDict_;

        word beachType_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        numericalBeach(const numericalBeach&);

        //- Disallow default bitwise assignment
        void operator=(const numericalBeach&);

public:

    //- Runtime type information
    TypeName("numericalBeach");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            numericalBeach,
            dictionary,
            (
                const word& subDictName, const fvMesh& mesh_
            ),
            (subDictName, mesh_)
        );

    // Constructors

        //- Construct from components
        numericalBeach
        (
            const word& type,
            const fvMesh& mesh_
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<numericalBeach> New
        (
            const word& subDictName,
            const fvMesh& mesh_
        );

    // Destructor

        virtual ~numericalBeach();


    // Member Functions

        virtual void correct
        (
            const labelList&,
            const scalarField&,
            volScalarField&
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
