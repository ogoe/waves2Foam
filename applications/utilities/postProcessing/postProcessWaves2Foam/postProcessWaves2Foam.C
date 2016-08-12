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
    postProcessWaves2Foam

Description
    Utility to carry out a wide range of post-processing operations, which are
    typically related to coastal and offshore engineering topics. The post-
    processing utility can, however, also be used on other data sets, e.g.
    general probing of velocities.

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "argList.H"

#if EXTBRANCH==1
    #if 310<OFVERSION
        #include "foamTime.H"
    #else
        #include "Time.H"
    #endif
#elif OFPLUSBRANCH==1
    #include "Time.H"
#else
    #include "Time.H"
#endif

//#if EXTBRANCH==1 && OFVERSION>310
//    #include "foamTime.H"
//#else
//    #include "Time.H"
//#endif

#include "fvMesh.H"

#include "postProcessingWaves.H"
#include "uniformDimensionedFields.H"

using namespace Foam;

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"

    // Needed by e.g. the reflection analysis
    Info << "\nReading g" << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read the dataProcessing dictionary
    IOdictionary postProcProperties
    (
        IOobject
        (
            "postProcessingProperties",
            "constant",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Potentially clean out prior directory postProcessedWaves2Foam
    Switch deleteOutput
        (
            postProcProperties.lookup("deleteParentOutputDirectory")
        );

    if (deleteOutput)
    {
        fileName dirName( "postProcessedWaves2Foam" );

        if (isDir(dirName))
        {
            Info << "Deleting 'postProcessedWaves2Foam'.\n" << endl;

            Foam::rmDir( dirName );
        }
    }

    // Get a list of sub-dictionaries in dictionary DP
    wordList toc( postProcProperties.toc() );

    // Loop over all items in TOC
    forAll (toc, itemi)
    {
        word tocName( toc[itemi] );

        if (postProcProperties.isDict(tocName))
        {
            Info << "Processing: " << tocName << endl;

            const dictionary& subDict( postProcProperties.subDict( tocName ) );

            wordList actionList( subDict.lookup("actionList") );

            // Loop over the action list for each process
            forAll (actionList, actionItem)
            {
                Info << "    Processing sub-action: "
                     << actionList[actionItem] << endl;
                autoPtr<postProcessingWaves> action
                    (
                        postProcessingWaves::New
                        (
                            runTime,
                            subDict,
                            actionList[actionItem]
                        )
                    );

                action->evaluate();
            }
            Info << endl;
        }
    }

    Info << nl << "End" << endl;

    return 0;
}

