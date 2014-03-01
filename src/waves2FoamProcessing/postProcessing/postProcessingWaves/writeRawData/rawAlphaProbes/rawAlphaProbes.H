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
    Foam::rawAlphaProbes

Description
    Reads, concatenates and output raw alpha1 probe data.

SourceFiles
    rawAlphaProbes.C

Author
    Bjarne Jensen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef rawAlphaProbes_H
#define rawAlphaProbes_H

#include "postProcessingWaves.H"
#include "crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class rawAlphaProbes
:
    public postProcessingWaves
{

protected:

    // Protected data
        word inputDir_;

        wordList timeDirs_;

        Switch removeDuplicate_;

    // Protected member functions

private:

    // Private Member Functions
        void resizeFields
        (
            List<std::pair<scalar, label> >&,
            List<scalarField>&,
            label
        );

        void writeRawData
        (
            const List<std::pair<scalar, label> >&,
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const List<scalarField>&
        );

        //- Disallow default bitwise copy construct
        rawAlphaProbes(const rawAlphaProbes&);

        //- Disallow default bitwise assignment
        void operator=(const rawAlphaProbes&);

    // Private member data

public:

    //- Runtime type information
    TypeName("rawAlphaProbes");

    // Constructors

        //- Construct from components
        rawAlphaProbes
        (
            const Time&,
            const dictionary&,
            const word&
        );


    // Destructor

        virtual ~rawAlphaProbes();

    // Member Functions
        virtual void evaluate();

        void readAlphaProbesData
        (
            List<std::pair<scalar, label> >& timeLabel,
            scalarField& x,
            scalarField& y,
            scalarField& z,
            List<scalarField>& alphas
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
