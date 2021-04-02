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
    gabcHelperFunctions

Description

Author
    Niels Gjoel Jacobsen, Deltares

\*---------------------------------------------------------------------------*/

#include "gabcHelperFunctions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gabcHelperFunctions::gabcHelperFunctions
(
    const fvMesh& mesh,
    const vector g,
    const label patchID
)
:
    mesh_(mesh),
    vertDir_(g),
    patchID_(patchID),

    alphaLimit_(Foam::waves2Foam::alphaInterface())
{
    vertDir_ /= Foam::mag(g);
    vertDir_ = Foam::cmptMultiply(vertDir_, vertDir_);

    if (vertDir_.y() < vertDir_.z())
    {
        gDirName_ = "z";
    }
    else if (vertDir_.z() < vertDir_.y())
    {
        gDirName_ = "y";
    }
    else
    {
        notImplemented("The code does not support gravity along the x-axis");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void gabcHelperFunctions::makeNeighbourList() const
{
    const fvPatch& patch(mesh_.boundary()[patchID_]);

    // Force the size to the number of faces of the patch
    faceFaces_.setSize(0);
    faceFaces_.setSize(patch.size());

    adjacentEdges_.setSize(0);
    adjacentEdges_.setSize(patch.size());

    label index = patch.index();
    label size = patch.size();
    label start = patch.patch().start();

    // Get the overall mesh properties
    const labelListList& ef = mesh_.edgeFaces();

    // Loop over all faces on the patch and find the connectivity
    for (long facei = start; facei < (start + size); facei++)
    {
        const labelList& fe = mesh_.faceEdges(facei);

        labelHashSet set(0);

        forAll (fe, edgei)
        {
            const labelList& ff(ef[fe[edgei]]);

            forAll (ff, facej)
            {
                label faceJ = ff[facej];

                if (start <= faceJ && faceJ < start + size && !set.found(faceJ))
                {
                    set.insert(faceJ);
                }
            }
        }

        labelList& fF = faceFaces_[facei - start];
        fF.setSize(set.toc().size());
        fF = set.toc() - start;
    }

    // Loop over all faces on the patch and find those, where faceEdges
    // is only of size 2 and not belonging to the same boundary
    const edgeList& edges = mesh_.edges();
    const pointField& points = mesh_.points();

    for (long facei = start; facei < (start + size); facei++)
    {
        const labelList& fe = mesh_.faceEdges(facei);

//        labelList tmpFaces(fe.size(), -1);
//        labelList tmpPatch(fe.size(), -1);
        vectorField tmpEdges(fe.size(), vector::zero);
        label count(0);

        forAll (fe, edgei)
        {
            const labelList& ff(ef[fe[edgei]]);

            if (ff.size() == 2)
            {
                forAll (ff, jj)
                {
                    if (start <= ff[jj] && ff[jj] < start + size)
                    {
                        // Nothing (but easier construct ;)
                    }
                    else
                    {
                        tmpEdges[count++] = edges[fe[edgei]].centre(points);
                    }
                }
            }
        }

        tmpEdges.setSize(count);

        adjacentEdges_[facei - start].setSize(tmpEdges.size());
        adjacentEdges_[facei - start] = tmpEdges;
    }
}


label gabcHelperFunctions::findPatchID(const label& faceID) const
{
    forAll (mesh_.boundary(), patchID)
    {
        const fvPatch& patch(mesh_.boundary()[patchID]);

        label size = patch.size();
        label start = patch.patch().start();

        if (start <= faceID && faceID < start + size)
        {
            return patchID;
        }
    }
}


const labelListList& gabcHelperFunctions::neighbourList() const
{
    const fvPatch& patch(mesh_.boundary()[patchID_]);

    if (faceFaces_.size() != patch.size() || mesh_.moving())
    {
        makeNeighbourList();
    }

    return faceFaces_;
}


tmp<vectorField> gabcHelperFunctions::normalGradient
(
    const scalarField& vel,
    const scalarField& alpha
) const
{
    tmp<vectorField> tproj = this->tangentialProjection();
    const vectorField& proj = tproj();

    // Coordinates
    const vectorField& cf = mesh_.Cf().boundaryField()[patchID_];
    scalarField coorZ = (cf & vertDir_);
    scalarField coorT = (cf & proj);

    // Get reference to the neighbour list
    const labelListList& nei(this->neighbourList());
    const Field<vectorField>& neiEdges(adjacentEdges_);

    // Prepare the return field
    tmp<vectorField> tres(new vectorField);
    vectorField& res(tres.ref());
    res.setSize(vel.size(), vector::zero);

    forAll(nei, facei)
    {
        if (alphaLimit() <= alpha[facei])
        {
            // Find all wetted cells
            const labelList& n = nei[facei];
            labelHashSet activeSet(0);

            forAll (n, facej)
            {
                if (alphaLimit() <= alpha[n[facej]])
                {
                    activeSet.insert(n[facej]);
                }
            }

            labelList active(activeSet.toc().size());
            active = activeSet.toc();

            if (neiEdges[facei].size() == 0)
            {
                scalarField sol(solveLS(active, coorT, coorZ, vel, facei));
                res[facei] = sol[0]*proj[facei] + sol[1]*vertDir_;
            }
            else
            {
                scalarField sol
                (
                    solveLS
                    (
                        active,
                        neiEdges[facei],
                        coorT,
                        coorZ,
                        vel,
                        proj[facei],
                        facei
                    )
                );

                res[facei] = sol[0]*proj[facei] + sol[1]*vertDir_;
            }
        }
        else
        {
            res[facei] = vector::zero;
        }
    }

    return tres;
}


tmp<vectorField> gabcHelperFunctions::tangentialProjection() const
{
    // Prepare the return field
    tmp<vectorField> tres(new vectorField);
    vectorField& res(tres.ref());

    vectorField nf = mesh_.Sf().boundaryField()[patchID_];
    nf /= Foam::mag(nf);

    res.setSize(nf.size());

    if (gDirName_ == "y")
    {
        res = (nf ^ vertDir_);
    }
    else
    {
        res = (vertDir_ ^ nf);
    }

    return tres;
}


void gabcHelperFunctions::checkTangential(const vectorField& proj) const
{
    scalarField sum(proj.component(0) + proj.component(1) + proj.component(2));

    forAll (proj, facei)
    {
        scalar a = Foam::cmptMax(proj[facei]);

        if (SMALL < a - sum[facei])
        {

            FatalErrorIn("void gabcHelperFunctions::checkTangential(const vectorField& proj) const")
                << "The patch is not aligned with the Cartesian coordinate system."
                << endl << endl << exit(FatalError);
        }
    }
}


scalarField gabcHelperFunctions::solveLS
(
    const labelList& active,
    const vectorField& adjacentEdges,
    const scalarField& coorT,
    const scalarField& coorZ,
    const scalarField& vel,
    const vector& proj,
    const label faceI
) const
{
    // Generate a subset of local coordinates and include the distance
    // for the boundaries
    // Currently only implemented for slip-type velocity boundaries!
    label N = active.size();
    labelList activeL(N + adjacentEdges.size(), 0);
    scalarField coorTL(N + adjacentEdges.size(), 0.0);
    scalarField coorZL(N + adjacentEdges.size(), 0.0);
    scalarField velL(N + adjacentEdges.size(), 0.0);

    label faceIL = -1;

    forAll (active, facei)
    {
        coorTL[facei] = coorT[active[facei]];
        coorZL[facei] = coorZ[active[facei]];
        velL[facei] = vel[active[facei]];
        activeL[facei] = facei;

        if (active[facei] == faceI)
        {
            faceIL = facei;
        }
    }

    // Loop over all the adjacent boundaries to create the additional
    // distances for same-order gradient approximation
    forAll (adjacentEdges, edgei)
    {
        vector e(adjacentEdges[edgei]);

        // Create a temporary coordinate for the ghost point
        vector tmpCoor = mesh_.Cf().boundaryField()[patchID_][faceI];
        tmpCoor += 2*((e - tmpCoor) & proj)*proj + 2*((e - tmpCoor) & vertDir_)*vertDir_;

        // Calculate the projected coordinates
        coorTL[N + edgei] = (tmpCoor & proj);
        coorZL[N + edgei] = (tmpCoor & vertDir_);

        // Insert as active
        activeL[N + edgei] = N + edgei;

        // Insert velocities such that they reflect a slip boundary condition
        scalar distT = Foam::mag(coorTL[N + edgei] - coorTL[faceIL]);
        scalar distZ = Foam::mag(coorZL[N + edgei] - coorZL[faceIL]);

        // The boundary is either (top or) bottom, so a slip means that the
        // velocity is simply copied.
        if (distT < distZ)
        {
            velL[N + edgei] = velL[faceIL];
        }
        // The boundary is on the sides, so the velocity flips sign
        else
        {
            velL[N + edgei] = -velL[faceIL];
        }
    }

    // Call the least-square system with the local data
    return solveLS(activeL, coorTL, coorZL, velL, faceIL);
}





scalarField gabcHelperFunctions::solveLS
(
    const labelList& active,
    const scalarField& coorT,
    const scalarField& coorZ,
    const scalarField& vel,
    const label facei
) const
{
    // Find the distance of all wetted cells to facei
    scalarField distT(active.size(), 0.0);
    scalarField distZ(active.size(), 0.0);

    forAll (active, ii)
    {
        distT[ii] = coorT[active[ii]] - coorT[facei];
        distZ[ii] = coorZ[active[ii]] - coorZ[facei];
    }

    scalarField res(2, 0);

    // There are no lateral cells available to calculate the gradient
    if (Foam::max(Foam::mag(distT)) < 1e-5)
    {
        List<scalarField> A(2);
        scalarField b(active.size());

        for (int ii = 0; ii < 2; ii++)
        {
            scalarField& atmp(A[ii]);
            atmp.setSize(active.size());
        }

        scalarField& A0 = A[0];
        scalarField& A1 = A[1];

        // Construct and solve
        forAll (active, nI)
        {
            A0[nI] = 1;
            A1[nI] = coorZ[active[nI]] - coorZ[facei];

            b[nI]  = vel[active[nI]];
        }

        dictionary tmpDict;
        solveLS(A, b);

        // b[0] is the average value
        res[1] = b[1];

        return res;
    }
    // Can only solve laterally
    else if (Foam::max(Foam::mag(distZ)) < 1e-5)
    {
        List<scalarField> A(2);
        scalarField b(active.size());

        for (int ii = 0; ii < 2; ii++)
        {
            scalarField& atmp(A[ii]);
            atmp.setSize(active.size());
        }

        scalarField& A0 = A[0];
        scalarField& A1 = A[1];

        // Construct and solve
        forAll (active, nI)
        {
            A0[nI] = 1;
            A1[nI] = coorT[active[nI]] - coorT[facei];

            b[nI]  = vel[active[nI]];
        }

        dictionary tmpDict;
        solveLS(A, b);

        // b[0] is the average value
        res[0] = b[1];

        return res;
    }
    else
    {
        List<scalarField> A(3);
        scalarField b(active.size());

        for (int ii = 0; ii < 3; ii++)
        {
            scalarField& atmp(A[ii]);
            atmp.setSize(active.size());
        }

        scalarField& A0 = A[0];
        scalarField& A1 = A[1];
        scalarField& A2 = A[2];

        // Construct and solve
        forAll (active, nI)
        {
            A0[nI] = 1;
            A1[nI] = coorT[active[nI]] - coorT[facei];
            A2[nI] = coorZ[active[nI]] - coorZ[facei];

            b[nI]  = vel[active[nI]];
        }

        dictionary tmpDict;

        solveLS(A, b);

        // b[0] is the average value
        res[0] = b[1]; res[1] = b[2];

        return res;
    }
}


void gabcHelperFunctions::solveLS
(
    const List<scalarField>& A,
    scalarField& b
) const
{
    // Create the least-squares right and left hand sides
    label N( A.size() );

#if EXTBRANCH==1
    scalarSquareMatrix AtA(N, 0.0);
#elif OFPLUSBRANCH==1
    SquareMatrix<scalar> AtA(N, N);
#else
    SquareMatrix<scalar> AtA(N, N);
#endif
    scalarField Atb( N, 0.0 );

    // Fill the matrix elements
    for (int i=0; i<N; i++)
    {
        const scalarField& ai( A[i] );

        for (int j=0; j<N; j++)
        {
            const scalarField& aj( A[j] );
            AtA[i][j] = Foam::sum( ai*aj );
        }

        Atb[i] = Foam::sum(ai*b);
    }

    // Solve the square system
#if EXTBRANCH==1
    Foam::scalarSquareMatrix::LUsolve(AtA, Atb);
#elif OFPLUSBRANCH==1
    Foam::LUsolve(AtA, Atb);
#else
    Foam::LUsolve(AtA, Atb);
#endif

    // Return solution in input field
    b = Atb;
}


void gabcHelperFunctions::zeroVorticity
(
    const vectorField& U,
    const scalarField& alpha,
    vectorField& gradient
) const
{
    gradient.setSize(U.size(), vector::zero);

    // Get a unit, positive vector
    vectorField nf = mesh_.Sf().boundaryField()[patchID_];
    nf /= Foam::mag(nf);
    nf = Foam::cmptMultiply(nf, nf);

    // Calculate the positive normal velocity
    scalarField vel(U & nf);

    // Define the coordinate system
    tmp<vectorField> tproj = this->tangentialProjection();
    const vectorField& proj = tproj();

    const vectorField& cf = mesh_.Cf().boundaryField()[patchID_];
    scalarField coorZ = (cf & vertDir_);
    scalarField coorT = (cf & proj);

    // Get reference to the neighbour list
    const labelListList& nei(this->neighbourList());

    forAll(nei, facei)
    {
        if (alphaLimit() <= alpha[facei])
        {
            // Find all wetted cells
            const labelList& n = nei[facei];
            labelHashSet activeSet(0);

            forAll (n, facej)
            {
                if (alphaLimit() <= alpha[n[facej]])
                {
                    activeSet.insert(n[facej]);
                }
            }

            labelList active(activeSet.toc().size());
            active = activeSet.toc();

            scalarField sol(solveLS(active, coorT, coorZ, vel, facei));

            gradient[facei] += (proj[facei]*sol[0] + vertDir_*sol[1]);
        }
    }
}


tmp<scalarField> gabcHelperFunctions::localCoordinate() const
{
    tmp<scalarField> tres(new scalarField);
    scalarField& res(tres.ref());

    // Get the points and face centres
    const fvPatch& patch = mesh_.boundary()[patchID_];
    const pointField& pp(patch.patch().localPoints());
    const vectorField cf = patch.patch().faceCentres();

    // Get projection vector and tangential coordinates
    tmp<vectorField> tproj(this->tangentialProjection());
    const vectorField& proj(tproj.ref());

    // Here, it is important to recognise that only flat boundaries are allowed!
    scalarField ypp(pp & proj[0]);
    res.setSize(cf.size());
    res = (cf & proj);

    // Determine bounding box
    // Maximum value
    bool zeroSize(false);
    if (ypp.size() == 0)
    {
        ypp.setSize(1, -GREAT);
        zeroSize = true;
    }

    yMax_ = Foam::gMax(ypp);

    // Minimum value
    if (zeroSize)
    {
        ypp.setSize(1, +GREAT);
    }

    yMin_ = Foam::gMin(ypp);

    // Centre the local coordinate around the middle of the patch
    res -= (0.5*(yMax_ + yMin_));

    // Return the field
    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
