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
    localCell.H

Description

SourceFiles
    localCell.C

Author
    Niels Gjoel Jacobsen, Technical University of Denmark

\*---------------------------------------------------------------------------*/

#ifndef localCell_H
#define localCell_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class localCell Declaration
\*---------------------------------------------------------------------------*/

class localCell
{
private:

    // Private member data
        cell cc_;
        faceList fL_;
        pointField pp_;
        edgeList eL_;
        labelListList faceEdges_;
        labelListList edgeFaces_;

        cell ccNeg_, ccPos_;
        label negCount_, posCount_, nextFace_;

    // Private member functions
        void localizeCell( const fvMesh&, const label&);

        void cellConnectivity();

protected:

public:

    // Constructor
        localCell();

        localCell(const fvMesh&, const label&);

        localCell
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
            ccPos_.setSize(0);
            fL_.setSize(cc_.size());
            pp_.setSize(cc_.labels(fL_).size());
            negCount_ = 0;
            posCount_ = 0;
            nextFace_ = cc_.size();
        };

        void initCut()
        {
            ccNeg_.setSize(cc_.size() + 2);
            ccPos_.setSize(cc_.size() + 2);
            fL_.setSize( 3*cc_.size() );
            nextFace_ = cc_.size();
            negCount_ = 0;
            posCount_ = 0;
        }

        void doneCut()
        {
            ccNeg_.setSize( negCount_ );
            ccPos_.setSize( posCount_ );
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

        void addPos( const face f )
        {
            ccPos_[posCount_++] = nextFace_;
            fL_[ nextFace_++ ] = f;

            if (nextFace_ == fL_.size())
            {
                fL_.setSize( nextFace_ + cc_.size() );
            }
        }

        void addPos( const label faceI )
        {
            ccPos_[posCount_++] = faceI;
        }

        void fullNeg()
        {
            ccNeg_ = cc_;
            negCount_ = cc_.size();
        }

        void fullPos()
        {
            ccPos_ = cc_;
            posCount_ = cc_.size();
        }

        void emptyCell();

        void localizeCell( const word );

    // Access functions
        scalar mag() const { return cc_.mag(pp_, fL_); };
        scalar magNeg() const { return ccNeg_.mag(pp_, fL_); };
        scalar magPos() const { return ccPos_.mag(pp_, fL_); };

        vector centre()    const { return cc_.centre( pp_, fL_);  };
        vector centreNeg() const { return ccNeg_.centre(pp_, fL_); };
        vector centrePos() const { return ccPos_.centre(pp_, fL_); };

        const pointField& points() const {return pp_; };
        pointField& points() { return pp_; };

        const edgeList& edges() const { return eL_; };
        const labelListList& faceEdges() const { return faceEdges_; };
        const labelListList& edgeFaces() const { return edgeFaces_; };

        const cell& cc() const { return cc_;};
        const cell& ccNeg() const { return ccNeg_; };
        const cell& ccPos() const { return ccPos_; };

        cell& ccNeg() { return ccNeg_; };

        const faceList& faces() const { return fL_;};
        faceList& faces() { return fL_; };

        const face& faces( const label faceI ) const { return fL_[faceI]; };

        const face iface() const { return fL_[ccNeg_[ ccNeg_.size() - 1]]; };

    // Access functions used for testing

        scalar fmag(const label faceI) const {return fL_[faceI].mag(pp_); };
        vector fcentre(const label faceI) const
        {
            return fL_[faceI].centre(pp_);
        };
        vector fnormal(const label faceI) const
        {
            return fL_[faceI].normal(pp_);
        };

    // Write
        void write( const fvMesh& , const List<localCell>& );
};


} // End namespace

#endif
