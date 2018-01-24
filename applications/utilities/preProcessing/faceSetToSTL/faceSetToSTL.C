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
    faceSetToSTL

Description
    A small utility to create simple STL-surfaces from a pointField and a
    faceList

    It can also based on a faceList with one entry perform a simple trans-
    lation of this face and create a closed STL. See the example file
    "stlDefinitions".

Author
    Niels Gjoel Jacobsen, Technical University of Denmark.  All rights reserved.

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

#include "fvMesh.H"
#include "faceTriangulation.H"
#include "triSurface.H"

using namespace Foam;

void extrudeFacesAndPoints
(
    const dictionary&,
    faceList&,
    pointField&
);

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    IOdictionary stlDefs
    (
        IOobject
        (
            "stlDefinitions",
            runTime.constant(),
            "triSurface",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    wordList toc = stlDefs.toc();

    forAll (toc, item)
    {
        if (stlDefs.isDict(toc[item]))
        {
            Info << "\nCreates the STL surface for " << toc[item] << endl;

            pointField pp(stlDefs.subDict(toc[item]).lookup("points"));
            faceList faces(stlDefs.subDict(toc[item]).lookup("faces"));

            triFaceList tfl(0);
            label count(0);

            if (
                   stlDefs.subDict(toc[item])
                   .lookupOrDefault<Switch>("extrude", false)
               )
            {
                if (faces.size() <= 1)
                {
                    extrudeFacesAndPoints
                        (
                            stlDefs.subDict(toc[item]),
                            faces,
                            pp
                        );
                }
                else
                {
                    Info << "\nWARNING: Using extrude, but"
                         << " multiple faces are defined\n" << endl;
                }
            }

            forAll (faces, facei)
            {
                faceTriangulation triangulation(pp, faces[facei], true);

                tfl.setSize(count + triangulation.size());

                forAll (triangulation, triI)
                {
                    tfl[count++] = triangulation[triI];
                }
            }

            triSurface ts(tfl, pp);

            Info << "Writes the STL surface for " << toc[item] << endl;

            ts.write( "constant/triSurface/"+toc[item]+".stl" );
        }
    }

    Info << nl << "End" << endl;

    return 0;
}

void extrudeFacesAndPoints
(
    const dictionary& dict, faceList& fL, pointField& pp
)
{
	// Get the extrude vector
	vector extrude( dict.lookup("extrudeVector") );

	// Check for correct (positive) orientation of the final shape
    vector faceCentre(fL[0].centre(pp));
    vector faceNormal(fL[0].normal(pp));

    scalar projection = (extrude & faceNormal)
    	/(Foam::mag(extrude)*Foam::mag(faceNormal));

    // If the projection is positive (!), the face is ordered correctly.
    // Note that this is the case, since the original face is swapped in
    // orientation at the end of this method.
    if (projection < -SMALL)
    {
        Info << "Orientation of the STL-surface is swapped to\n"
        	 << "yield outward pointing normals." << endl;

        // Swap the direction of the original face
        face dummy = fL[0];
        face& target = fL[0];

        forAll (target, pointi)
        {
          	target[target.size() - 1 - pointi] = dummy[pointi];
        }
    }


	// Get the number of points and extrude the points
    label N = pp.size();
    pp.setSize(2*N);

    for (int i=0; i < N; i++)
    {
        pp[N + i] = pp[i] + extrude;
    }

    label M = fL[0].size();

    fL.setSize(2 + M);

    fL[1].setSize(M);

    face& fOrg(fL[0]);
    face& fExt(fL[1]);

    forAll (fOrg, pointi)
    {
//        fExt[N - 1 - pointi] = fOrg[pointi] + N;
        fExt[pointi] = fOrg[pointi] + N;
    }

    for (int i = 0; i < M ; i++)
    {
        face& f(fL[i + 2]);
        f.setSize(4);

        f[0] = fOrg[i];
        f[1] = fOrg[(i + 1)%M];
        f[2] = fExt[(i + 1)%M];
        f[3] = fExt[i];
    }

    // Swap the direction of the original face
    face dummy = fOrg;

    forAll (fOrg, pointi)
    {
//        fExt[N - 1 - pointi] = fOrg[pointi] + N;
    	fOrg[N - 1 - pointi] = dummy[pointi];
//        fExt[pointi] = fOrg[pointi] + N;
    }
}
