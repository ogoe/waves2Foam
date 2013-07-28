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
    localCellNeg.H

Description

SourceFiles
    localCellNeg.C

Author
    Niels Gjoel Jacobsen, Technical University of Denmark

\*---------------------------------------------------------------------------*/

#ifndef localCellNeg_H
#define localCellNeg_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class localCellNeg Declaration
\*---------------------------------------------------------------------------*/

class localCellNeg
{
private:

    // Private member data
        cell cc_;
        faceList fL_;
        pointField pp_;
        edgeList eL_;
        labelListList faceEdges_;
        labelListList edgeFaces_;

        cell ccNeg_;
        label negCount_, nextFace_;

    // Private member functions
        void localizeCell( const fvMesh&, const label&);

        void cellConnectivity();

protected:

public:

    // Constructor
        localCellNeg();

        localCellNeg(const fvMesh&, const label&);

        localCellNeg
        (
            const cell,
            const faceList,
            const pointField,
            bool checkCell = false
        );

        void initCell( const fvMesh&, const label&);


    // Cut-control
        void clearCut()
        {
            ccNeg_.setSize(0);
            fL_.setSize(cc_.size());
            pp_.setSize(cc_.labels(fL_).size());
            negCount_ = 0;
            nextFace_ = cc_.size();
        };

        void initCut()
        {
            ccNeg_.setSize(cc_.size() + 2);
            fL_.setSize( 3*cc_.size() );
            nextFace_ = cc_.size();
            negCount_ = 0;
        }

        void doneCut()
        {
            ccNeg_.setSize( negCount_ );
            fL_.setSize( nextFace_ );
        }

        void addNeg( const face f )
        {
            ccNeg_[negCount_++] = nextFace_;
            fL_[ nextFace_++] = f;

            if (nextFace_ == fL_.size())
            {
                fL_.setSize( nextFace_ + cc_.size() );
            }
        }

        void addNeg( const label faceI )
        {
            ccNeg_[negCount_++] = faceI;
        }

        void fullNeg()
        {
            ccNeg_ = cc_;
            negCount_ = cc_.size();
        }

        void emptyCell();

        void localizeCell( const word type = "neg");

    // Access functions
        scalar mag()    const { return cc_.mag(pp_, fL_);    };
        scalar magNeg() const { return ccNeg_.mag(pp_, fL_); };

        vector centre()    const { return cc_.centre(pp_, fL_);    };
        vector centreNeg() const { return ccNeg_.centre(pp_, fL_); };

        const pointField& points() const { return pp_; };
        pointField&       points()       { return pp_; };

        const edgeList& edges()     const { return eL_;        };
        const labelListList& faceEdges() const { return faceEdges_; };
        const labelListList& edgeFaces() const { return edgeFaces_; };

        const cell& cc()    const { return cc_;    };
        const cell& ccNeg() const { return ccNeg_; };
        cell&       ccNeg()       { return ccNeg_; };

        const faceList& faces() const { return fL_; };
        faceList&       faces()       { return fL_; };
};


} // End namespace

#endif
