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
    Foam::removeData

Description
    Remove the OF and ASCII formatted data.

    OF data is removed by default and ASCII formatted data is by default
    retained in the folder.

SourceFiles
    removeData.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef removeData_H
#define removeData_H

#include "postProcessingWaves.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class removeData
:
    public postProcessingWaves
{

protected:

    // Protected data
        word inputDir_;

        Switch removeOF_;

        Switch removeAscii_;

    // Protected member functions

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        removeData(const removeData&);

        //- Disallow default bitwise assignment
        void operator=(const removeData&);

    // Private member data

public:

    //- Runtime type information
    TypeName("removeData");

    // Constructors

        //- Construct from components
        removeData
        (
            const Time&,
            const dictionary&,
            const word&
        );


    // Destructor

        virtual ~removeData();

    // Member Functions
        virtual void evaluate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
