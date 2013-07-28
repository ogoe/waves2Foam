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
    localCellNeg.C

Description
    See localCellNeg.H

Author
    Niels Gjoel Jacobsen, Technical University of Denmark

\*---------------------------------------------------------------------------*/

#include "localCellNeg.H"
#include <map>

namespace Foam
{
// * * * * * * * * * * * * * * * Constructurs  * * * * * * * * * * * * * * * //


localCellNeg::localCellNeg
(
    const fvMesh& mesh,
    const label& cellI
)
:
    ccNeg_(0),
    negCount_(0),
    nextFace_(0)
{
    localizeCell(mesh, cellI);
}


localCellNeg::localCellNeg
(
    const cell cc,
    const faceList fL,
    const pointField pp,
    bool checkCell
)
:
    cc_(cc),
    fL_(fL),
    pp_(pp),
    ccNeg_(0),
    negCount_(0),
    nextFace_(0)
{
    cellConnectivity();

    if (checkCell)
    {
        notImplemented("Check cell in localCell");
    }
}


localCellNeg::localCellNeg()
:
    ccNeg_(0),
    negCount_(0),
    nextFace_(0)
{}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


void localCellNeg::localizeCell
(
    const fvMesh& mesh,
    const label& cellI
)
{
    // Reference to some global mesh data
    const labelList& own = mesh.owner();
    const faceList& ff = mesh.faces();
    const pointField& pp = mesh.points();
    const label nFaces      = mesh.nInternalFaces();

    // Get local cell properties
    cc_ = mesh.cells()[cellI];
    fL_.setSize(cc_.nFaces());

    // Initialise the faceList
    forAll (fL_, facei)
    {
        fL_[facei] = ff[cc_[facei]];
    }

    // Orient all face, so the normal vector is pointing outward

    forAll (fL_, facei)
    {
        label nFace( cc_[facei] );

        // If the other holds, the normal points outward per definition
        if (nFace < nFaces && own[nFace] != cellI)
        {
            face& f( fL_[facei] );
            f = f.reverseFace();
        }
    }

    // Make cellFaces local
    forAll (cc_, celli)
    {
        cc_[celli] = celli;
    }

    // Make local pointField and global-local pointMap
    const labelList pLabels(cc_.labels(fL_));
    pp_.setSize(pLabels.size());
    std::map <label, label> pointMap;

    forAll (pLabels, pointi)
    {
        pp_[pointi] = pp[pLabels[pointi]];
        pointMap[pLabels[pointi]] = pointi;
    }

    // Change point labels in the faces to be consistent with local
    // point numbering
    forAll (fL_, facei)
    {
        face& f(fL_[facei]);

        forAll (f, pointi)
        {
            f[pointi] = pointMap.find(f[pointi])->second;
        }
    }

    eL_ = cc_.edges(fL_);

    // Generate faceEdges - better to use original data structure?
    faceEdges_.setSize(cc_.size());

    edgeFaces_.setSize( eL_.size() );
    labelList counter(eL_.size(), 0);

    forAll (edgeFaces_, edgei)
    {
        labelList& eF = edgeFaces_[edgei];
        eF.setSize(2);
    }

    forAll (fL_, facei)
    {
        const face& f(fL_[facei]);
        const edgeList& eLf = f.edges();

        labelList& fE( faceEdges_[facei] );
        fE.setSize(eLf.size());

        forAll (eLf, edgei)
        {
            label count(0);

            while (true)
            {
                if (eLf[edgei] == eL_[count])
                {
                    fE[edgei] = count;

                    labelList& eF( edgeFaces_[count]);
                    eF[counter[count]++] = facei;

                    break;
                }
                count++;
            }
        }
    }

    cellConnectivity();
}


void localCellNeg::cellConnectivity()
{
    eL_ = cc_.edges(fL_);

    // Generate faceEdges - better to use original data structure?
    faceEdges_.setSize(cc_.size());

    edgeFaces_.setSize( eL_.size() );
    labelList counter(eL_.size(), 0);

    forAll (edgeFaces_, edgei)
    {
        labelList& eF = edgeFaces_[edgei];
        eF.setSize(2);
    }

    forAll (fL_, facei)
    {
        const face& f(fL_[facei]);
        const edgeList& eLf = f.edges();

        labelList& fE( faceEdges_[facei] );
        fE.setSize(eLf.size());

        forAll (eLf, edgei)
        {
            label count(0);

            while (true)
            {
                if (eLf[edgei] == eL_[count])
                {
                    fE[edgei] = count;

                    labelList& eF( edgeFaces_[count]);
                    eF[counter[count]++] = facei;

                    break;
                }
                count++;
            }
        }
    }
}

// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //


void localCellNeg::localizeCell
(
    const word cellSide
)
{
    // Set the cell to the negative or positive intersected cell
    if (cellSide == "neg")
    {
        this->cc_ = this->ccNeg_;
    }
    else
    {
        Info << "FATAL ERROR!!!" << endl;
    }

    // Clean out the face list to contain only the necessary faces
    faceList fL( this->cc_.size());

    forAll (this->cc_, facei)
    {
        fL[facei] = fL_[this->cc_[facei]];
    }

    this->fL_ = fL;

    // Renumber the cell
    forAll (this->cc_, celli)
    {
        this->cc_[celli] = celli;
    }

    // Make local pointField and global-local pointMap
    labelList pLabels(this->cc_.labels(this->fL_));
    pointField pp(pLabels.size());
    std::map <label, label> pointMap;

    forAll (pLabels, pointi)
    {
        pp[pointi] = pp_[pLabels[pointi]];
        pointMap[pLabels[pointi]] = pointi;
    }

    this->pp_ = pp;

    // Change point labels in the faces to be consistent with local
    // point numbering
    forAll (this->fL_, facei)
    {
        face& f(this->fL_[facei]);

        forAll (f, pointi)
        {
            f[pointi] = pointMap.find(f[pointi])->second;
        }
    }

    cellConnectivity();

    ccNeg_.setSize(0);
    negCount_ = 0;
    nextFace_ = cc_.size();
}


void localCellNeg::initCell( const fvMesh& mesh, const label& cellI)
{
    localizeCell(mesh, cellI);
}


void localCellNeg::emptyCell()
{
    cc_.setSize(0);
    fL_.setSize(0);
    pp_.setSize(0);

    eL_.setSize(0);
    faceEdges_.setSize(0);
    edgeFaces_.setSize(0);

    ccNeg_.setSize(0);
    negCount_ = 0;
    nextFace_ = 0;
}

/*void localCellNeg::write(const fvMesh& mesh, const List<localCellNeg>& lcs)
{
    label pCount(0);

    forAll (lcs, celli)
    {
        const localCellNeg& lc(lcs[celli]);
        pCount += lc.points().size();
    }

    faceList fL( lcs.size() );

    pointField pp(pCount);
    pCount = 0;

    forAll (lcs, celli)
    {
        const localCellNeg& lc(lcs[celli]);

        // Put points into a common pointField
        pointField p( lc.points() );

        forAll (p, pointi)
        {
            pp[pointi + pCount] = p[pointi];
        }

        // Add interface-face to the face list. Correct for the new point
        // numbering
        face& f( fL[celli]);

        f = lc.iface();

        forAll (f, pointi)
        {
            f[pointi] += pCount;
        }

        // Update point offset
        pCount += p.size();
    }

    pointIOField points
    (
        IOobject
        (
            "points",
            mesh.time().timeName(),
            "surfaceReconstruction",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    faceIOList faceLists
    (
        IOobject
        (
            "faces",
            mesh.time().timeName(),
            "surfaceReconstruction",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    points = pp;
    points.write();

    faceLists = fL;
    faceLists.write();

}*/


} // End namespace




