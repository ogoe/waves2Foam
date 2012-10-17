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
    Foam::setWaveProperties

Description
	An abstract base class for a pre-processing tool for setting the wave
	properties in the dictionary waveProperties.

    A description of the general wave framework is given in

    @article { jacobsenFuhrmanFredsoe2011,
		Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
		title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
		Journal = {{Int. J. for Numer. Meth. Fluids}},
		Year = {2011},
		Volume = {In print},
		Pages = {},
	}

SourceFiles
    setWaveProperties.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef setWaveProperties_H
#define setWaveProperties_H

#include "fvMesh.H"
#include "volFields.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "uniformDimensionedFields.H"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class setWaveProperties Declaration
\*---------------------------------------------------------------------------*/

class setWaveProperties
{

protected:

    // Protected data
		bool write_;
		dictionary & dict_;
		vector g_;
		scalar G_;

		scalar PI_;

    // Protected member functions

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        setWaveProperties(const setWaveProperties&);

        //- Disallow default bitwise assignment
        void operator=(const setWaveProperties&);

public:

    //- Runtime type information
    TypeName("setWaveProperties");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            setWaveProperties,
            setWaveProperties,
            (
            	const fvMesh & mesh, dictionary & dict, bool write
            ),
            (mesh, dict, write)
        );

    // Constructors

        //- Construct from components
        setWaveProperties
        (
        	const fvMesh &,
        	dictionary &,
        	bool
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<setWaveProperties> New
        (
        	const fvMesh &,
        	dictionary &,
        	bool
        );

    // Destructor

        virtual ~setWaveProperties();


    // Member Functions

        virtual void set() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
