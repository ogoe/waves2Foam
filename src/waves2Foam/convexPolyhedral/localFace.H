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
    localFace

Description

SourceFiles
    localFace.C

Author
    Niels Gjoel Jacobsen, Technical University of Denmark

\*---------------------------------------------------------------------------*/

#ifndef localFace_H
#define localFace_H

#include "fvMesh.H"


namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class localFace Declaration
\*---------------------------------------------------------------------------*/

class localFace
{
private:

    pointField pp_;
    face negF_;
    face posF_;
    label numberOfIntersections_;

public:

    // Constructors
        localFace( )
        :
            pp_(0),
            negF_(0),
            posF_(0),
            numberOfIntersections_(0)
        {};

    // Clearing private variables
        void clearCut()
        {
            negF_.setSize(0);
            posF_.setSize(0);
            pp_.setSize(0);
            numberOfIntersections_ = 0;
        };

    // Set local face properties
        void setAttributes
        (
            const pointField& pp,
            const face& negF,
            const face& posF,
            const label& nOI
        )
        {
            pp_    = pp;
            negF_  = negF;
            posF_  = posF;
            numberOfIntersections_ = nOI;
        };

    // Demand driven data
        inline bool isPosFace() const
        {
            return posF_.size() > 2 ? true : false;
        };

        inline bool isNegFace() const
        {
            return negF_.size() > 2 ? true : false;
        };

        const face negFace() const {return negF_; };

        const face posFace() const {return posF_; };

        face& negFace() { return negF_; };
        face& posFace() { return posF_; };

        const pointField& points() const { return pp_; };
        pointField& points() { return pp_; };
        void points( const pointField& pp) { pp_ = pp; };

        inline point posCentre() const { return posF_.centre(pp_); };

        inline point negCentre() const
        {
            return ( isNegFace() ? negF_.centre(pp_) : vector::zero );
        };

        inline scalar posMag()   const { return posF_.mag(pp_);    };

        inline scalar negMag() const
        {
            return ( isNegFace() ? negF_.mag(pp_) : 0.0 );
        };

        inline label noi() const { return numberOfIntersections_; };
        void noi(label noi) { numberOfIntersections_ = noi;};

};

}


#endif
