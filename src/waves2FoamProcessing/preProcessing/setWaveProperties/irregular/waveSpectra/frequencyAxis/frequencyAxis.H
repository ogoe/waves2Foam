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
    Foam::frequencyAxis

Description
    An abstract base class for setting the frequencyAxis for irregular wave trains.

SourceFiles
    frequencyAxis.C

Author
    Niels Gjøl Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef frequencyAxis_H
#define frequencyAxis_H

#include "fvMesh.H"
#include "volFields.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "setWaveProperties.H"

#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class frequencyAxis Declaration
\*---------------------------------------------------------------------------*/

class frequencyAxis
{

protected:

    // Protected data
        const Time& rT_;

        dictionary& dict_;

        scalar fp_;
        scalar fl_;
        scalar fu_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        frequencyAxis(const frequencyAxis&);

        //- Disallow default bitwise assignment
        void operator=(const frequencyAxis&);

public:

    //- Runtime type information
    TypeName("frequencyAxis");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            frequencyAxis,
            frequencyAxis,
            (
                const Time& rT,
                dictionary& dict
            ),
            (rT, dict)
        );

    // Constructors

        //- Construct from components
        frequencyAxis
        (
            const Time&,
            dictionary&
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<frequencyAxis> New
        (
            const Time&,
            dictionary&
        );

    // Destructor

        virtual ~frequencyAxis();

    // Member functions

        virtual scalarField freqAxis
        (
            const scalarField&,
            const scalarField&,
            const label&
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
