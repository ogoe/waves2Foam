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
        Year = {2012},
        Volume = {70},
        Number = {9},
        Pages = {1073-1088},
        DOI = {{10.1002/fld.2726}},
    }

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvMesh.H"
#include "IOstream.H"

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

#   include "readGravitationalAcceleration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    IOobject wOut
    (
            "wavePropertiesTEMP",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
    );

    // Write waveProperties with the above computed changes
    OFstream os
    (
        wOut.objectPath(),
        ios_base::out|ios_base::trunc,
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED
    );

    // Write the OF banner
    wOut.writeBanner( os );

    // Write the file information. Class name is not correct when
    // using wOut.writeHeader( os ); hence manual entries
    os << "FoamFile" << nl;
    os << token::BEGIN_BLOCK << incrIndent << nl;
    os << indent << "version" << tab << IOstream::currentVersion << token::END_STATEMENT << nl;
    os << indent << "format" << tab << "ascii;" << nl;
    os << indent << "class" << tab << "dictionary;" << nl;
    os << indent << "object" << tab << "waveProperties;" << nl;
    os << decrIndent << indent << token::END_BLOCK << nl;

    // Write the divider
    wOut.writeDivider( os );
    os << nl;
    
    /* Loop over all subdicts in waveProperties. For each of them compute the
       wave parameters relevant for that particular wave theory. */
    wordList toc = waveProperties.toc();

    forAll(toc, item )
    {
    	// If a sub-dictionary, then compute parameters and write the subdict
    	if ( waveProperties.isDict(toc[item]) )
    	{
    		dictionary & sd = waveProperties.subDict(toc[item]);

    		autoPtr<setWaveProperties> props( setWaveProperties::New(mesh, sd, true) );

    		props->set( os );
    	}
    	else
    	{
    		label Nspaces = 20;

    		// Read the entry and write to the dummy output file
    	    ITstream read = waveProperties.lookup(toc[item]);
    	    os << toc[item] << token::SPACE;

    		for( int i=toc[item].size(); i<Nspaces-1; i++)
    			os << token::SPACE;
    	    
    	    forAll(read, ri )
    	    {
    	        if ( ri < read.size() - 1)
        	        os << read[ri] << token::SPACE;
        	    else
        	        os << read[ri];
        	}
    	        
    	    os << token::END_STATEMENT << nl << endl;
    	}
    }

//    // TO BE REMOVE EVENTUALLY
//    waveProperties.regIOobject::write();

    // Write end divider
    wOut.writeEndDivider( os );

    // Move the dummy output file to the waveProperties file
    fileName src( wOut.objectPath() );
    fileName dst( wOut.path() + "/waveProperties" );

    Foam::mv( src, dst );

    // End
    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
