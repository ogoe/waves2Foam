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
    Foam::constantDepthCelerity

Author
    Niels Gjoel Jacobsen, Deltares
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef constantDepthCelerity_H
#define constantDepthCelerity_H

#include "celerityShapeFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class constantDepthCelerity
:
    public celerityShapeFunctions
{

protected:

    // Protected data
    	scalar depth_;
    	scalar fact_;

    // Protected member functions

private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantDepthCelerity(const constantDepthCelerity&);

        //- Disallow default bitwise assignment
        void operator=(const constantDepthCelerity&);

public:

    //- Runtime type information
    TypeName("constantDepthCelerity");

    // Constructors

        //- Construct from components
        constantDepthCelerity
        (
            const fvMesh&,
            const word&,
            const dictionary&,
            const vector&,
            const scalar&
        );

    // Destructor

        virtual ~constantDepthCelerity();

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
