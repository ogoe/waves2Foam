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
    Foam::gabcHelperFunctions

Description


SourceFiles
    gabcHelperFunctions.C

Author
    Niels Gjoel Jacobsen, Deltares
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef gabcHelperFunctions_H
#define gabcHelperFunctions_H

#include "volFields.H"
#include "uniformDimensionedFields.H"

#include "crossVersionCompatibility.H"


namespace Foam
{

class gabcHelperFunctions
{
private:

    // Private member variables
        const fvMesh& mesh_;
        vector vertDir_;
        const label patchID_;

        const scalar alphaLimit_;
        word gDirName_;

        // Patch neighbours
        void makeNeighbourList() const;


        label findPatchID(const label&) const;

        mutable labelListList faceFaces_;
        mutable Field<vectorField> adjacentEdges_;

        mutable scalar yMin_, yMax_;

        scalarField solveLS
        (
            const labelList& active,
            const scalarField& coorT,
            const scalarField& coorZ,
            const scalarField& vel,
            const label facei
        ) const;

        scalarField solveLS
        (
            const labelList& active,
            const vectorField& adjacentEdges,
            const scalarField& coorT,
            const scalarField& coorZ,
            const scalarField& vel,
            const vector& proj,
            const label facei
        ) const;

        void solveLS
        (
            const List<scalarField>& A,
            scalarField& b
        ) const;

    // Private member functions
        void checkTangential(const vectorField&) const;

        //- Disallow default bitwise copy construct
        gabcHelperFunctions(const gabcHelperFunctions &);

        //- Disallow default bitwise assignment
        void operator=(const gabcHelperFunctions &);

protected:



public:

        gabcHelperFunctions
        (
            const fvMesh&,
            const vector,
            const label
        );


        // Member functions
        inline scalar alphaLimit() const {return alphaLimit_;};

        tmp<vectorField> normalGradient
        (
            const scalarField&,
            const scalarField&
        ) const;

        tmp<vectorField> tangentialProjection() const;

        void zeroVorticity
        (
            const vectorField&,
            const scalarField&,
            vectorField&
        ) const;

        tmp<scalarField> localCoordinate() const;

        const labelListList& neighbourList() const;
};

} // End namespace Foam

#endif
