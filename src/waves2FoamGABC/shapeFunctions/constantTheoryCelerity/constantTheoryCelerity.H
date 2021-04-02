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
    Foam::constantTheoryCelerity

Author
    Niels Gjoel Jacobsen, Deltares
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef constantTheoryCelerity_H
#define constantTheoryCelerity_H

#include "celerityShapeFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class constantTheoryCelerity
:
    public celerityShapeFunctions
{

protected:

    // Protected data
    	word sourceName_;

    	scalar celerity_;

    // Protected member functions

private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantTheoryCelerity(const constantTheoryCelerity&);

        //- Disallow default bitwise assignment
        void operator=(const constantTheoryCelerity&);

public:

    //- Runtime type information
    TypeName("constantTheoryCelerity");

    // Constructors

        //- Construct from components
        constantTheoryCelerity
        (
            const fvMesh&,
            const word&,
            const dictionary&,
            const vector&,
            const scalar&
        );

    // Destructor

        virtual ~constantTheoryCelerity();

    // Member Functions
        virtual tmp<scalarField> rhoCelerity
        (
        	const fvPatch&,
            const scalarField&,
            const scalarField&,
            const scalar&
        );

        virtual void rhoCelerity
        (
            const fvPatch&,
            const scalarField&,
            const scalarField&,
            const scalar&,
            scalarField&,
            scalarField&,
            scalarField&
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
