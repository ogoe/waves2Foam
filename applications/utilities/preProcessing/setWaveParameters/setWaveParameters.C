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

Application
    setWaveParameters

Description
    This function loops through all of the sub-dictionaries in waveProperties
    dictionary in <root case>/constant. The needed wave parameters are computed
    based on a set of required input data.

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

Additional information
    Implementation published and validated in the following journal article:

    @article { jacobsenFuhrmanFredsoe2011,
		Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
		title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
		Journal = {{Int. J. for Numer. Meth. Fluids}},
		Year = {2011},
		Volume = {In print},
		Pages = {},
	}

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvMesh.H"

#include "setWaveProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    Info << "\nReading waveProperties\n" << endl;

    IOdictionary waveProperties
    (
        IOobject
        (
            "waveProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

#if OFVERSION == 15
#    include "readEnvironmentalProperties.H"
#else
#   include "readGravitationalAcceleration.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    wordList toc = waveProperties.toc();
    
    /* Loop over all subdicts in waveProperties. For each of them compute the
       wave parameters relevant for that particular wave theory. */
    forAll(toc, item )
    {
    	if ( waveProperties.isDict(toc[item]) )
    	{
    		dictionary & sd = waveProperties.subDict(toc[item]);

    		autoPtr<setWaveProperties> props( setWaveProperties::New(mesh, sd, true) );

    		props->set();
    	}
    }

    // Write waveProperties with the above computed changes
    waveProperties.regIOobject::write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
