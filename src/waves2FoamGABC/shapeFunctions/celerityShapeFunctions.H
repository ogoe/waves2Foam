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
    Foam::celerityShapeFunctions

Description


SourceFiles
    celerityShapeFunctions.C

Author
    Niels Gjoel Jacobsen, Deltares
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef celerityShapeFunctions_H
#define celerityShapeFunctions_H

#include "fvMesh.H"
#include "fvCFD.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class celerityShapeFunctions
{

protected:

    // Protected data
        const fvMesh& mesh_;

        const word patchName_;

        const dictionary& dict_;

        vector vertVector_;
        scalar seaLevel_;

    // Protected member functions
        word dictName() const;

        scalarField zCoordinateFace() const;

        scalar returnDepth() const;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        celerityShapeFunctions(const celerityShapeFunctions&);

        //- Disallow default bitwise assignment
        void operator=(const celerityShapeFunctions&);

public:

    //- Runtime type information
    TypeName("celerityShapeFunctions");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            celerityShapeFunctions,
            celerityShapeFunctions,
            (
                const fvMesh& mesh,
                const word& patchName,
                const dictionary& dict,
                const vector& vertical,
                const scalar& seaLevel
            ),
            (mesh, patchName, dict, vertical, seaLevel)
        );

    // Constructors

        //- Construct from components
        celerityShapeFunctions
        (
            const fvMesh&,
            const word&,
            const dictionary&,
            const vector&,
            const scalar&
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<celerityShapeFunctions> New
        (
            const fvMesh&,
            const word&,
            const dictionary&,
            const vector&,
            const scalar&
        );

    // Destructor

        virtual ~celerityShapeFunctions();

    // Member Functions
        bool pressureMethod() const {return false;};

        virtual tmp<scalarField> rhoCelerity
        (
        	const fvPatch&,
            const scalarField&,
            const scalarField&,
            const scalar&
        ) = 0;

        virtual void rhoCelerity
        (
            const fvPatch&,
            const scalarField&,
            const scalarField&,
            const scalar&,
            scalarField&,
            scalarField&,
            scalarField&
        ) = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
