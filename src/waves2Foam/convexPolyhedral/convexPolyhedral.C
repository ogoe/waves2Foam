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
    convexPolyhedral

Description
    See convexPolyhedral.H

Author
    Niels Gjoel Jacobsen, Technical University of Denmark

\*---------------------------------------------------------------------------*/

#include "convexPolyhedral.H"
#include <map>

namespace Foam
{

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


labelList convexPolyhedral::edgeCutLabel
(
    const edgeList& eL,
    const labelList& pType,
    const scalarField& sD,
    pointField& pf
)
{
    labelList edgeCut( eL.size(), -1);

    scalarField d( Foam::mag(sD) );

    label pCount(pf.size());
    pf.setSize( pf.size() + eL.size() );

    forAll (eL, edgei)
    {
        edge e(eL[edgei]);

        if (pType[e.start()]*pType[e.end()] == -1)
        {
            point newP;

            if (!iterateDistance_)
            {
                label start(e.start()), end(e.end());

                newP = pf[start] + d[start]/( d[start] + d[end] )
                    *(pf[end] - pf[start]);
            }
            else
            {
                point p0( pf[e.start()]), p1( pf[e.end()] );
                scalar d0( sD[e.start()]), d1( sD[e.end()]);

                newP = p0 + Foam::mag(d0)/(Foam::mag(d0) + Foam::mag(d1))
                    *( p1 - p0 );

                scalar dist( signedPointToSurfaceDistance(newP));
                label count(0);

                while (Foam::mag( dist ) > 5.0e-15)
                {
                    if (Foam::sign( dist ) == Foam::sign( d0 ))
                    {
                        d0 = dist;
                        p0 = newP;
                    }
                    else
                    {
                        d1 = dist;
                        p1 = newP;
                    }

                    newP = p0 + Foam::mag(d0)/(Foam::mag(d0) + Foam::mag(d1))
                        *( p1 - p0 );
                    dist = signedPointToSurfaceDistance( newP );

                    if (count++ > 100)
                    {
                        break;
                    }
                }

            }

            edgeCut[edgei] = pCount;
            pf[pCount++]   = newP;
        }
    }
    pf.setSize(pCount);

    return edgeCut;
}


void convexPolyhedral::faceCut
(
    const labelList& pType,
    const edgeList& eL,
    const labelList& edgeCut,
    localCell& lc
)
{
    lc.initCut();

    localFace lf;

    edgeList cuttedEdges(lc.cc().size());
    label cutCount(0);

    bool noiProblem( false );

    // Loop over all external faces and find neg/pos cuts.
    // Since they are created with the same orientation as the
    // original cell, all cell normals points outward.

    const labelListList& faceEdges(lc.faceEdges());
    const faceList& fL = lc.faces();

    forAll (lc.cc(), facei)
    {
        const face& f( fL[facei] );
        const edgeList& eLt  (f.edges());

        labelList edgeCutt(f.size());

        const labelList& fE( faceEdges[facei] );

        forAll (eLt, edgei)
        {
            edgeCutt[edgei] = edgeCut[ fE[edgei] ];
        }

        edge& cutted( cuttedEdges[cutCount]);

        faceCut(pType, f, eLt, edgeCutt, lc.points(), lf, cutted);

        if (lf.isNegFace() && lf.isPosFace())
        {
            lc.addNeg( lf.negFace() );
            lc.addPos( lf.posFace() );

            cutCount++;

        }
        else if (lf.isNegFace())
        {
            lc.addNeg( facei );
        }
        else if (lf.isPosFace())
        {
            lc.addPos( facei );
        }
        if (lf.noi() > 2)
        {
            noiProblem = true;
        }
    }

    if (noiProblem)
    {
        lc.initCut();

        if (signedPointToSurfaceDistance(lc.centre() ) < 0.0)
        {
            lc.fullNeg();
        }
        else
        {
            lc.fullPos();
        }
    }
    else
    {

        // Add special edges, where both end are on the plane
        forAll (eL, edgei)
        {
            if (pType[eL[edgei].start()] == 0 && pType[eL[edgei].end()] == 0)
            {
                cuttedEdges[cutCount++] = eL[edgei];
            }
        }

        // Take care of the interface face
        cuttedEdges.setSize(cutCount);

        if (cuttedEdges.size() >= 3)
        {
            face f = combineEdgeList(cuttedEdges);

            if (Foam::sign( f.normal( lc.points()) & normalToPlane_ ) == 1)
            {
                lc.addPos(f.reverseFace());
                lc.addNeg(f);
            }
            else
            {
                lc.addPos(f);
                lc.addNeg(f.reverseFace());
            }
        }
    }
}


void convexPolyhedral::faceCut
(
    const labelList& pType,
    const edgeList& eL,
    const labelList& edgeCut,
    localCellNeg& lc
)
{
    lc.initCut();

    localFace lf;

    edgeList cuttedEdges(lc.cc().size());
    label cutCount(0);

    bool noiProblem( false );

    // Loop over all external faces and find neg/pos cuts.
    // Since they are created with the same orientation as the
    // original cell, all cell normals points outward.

    const labelListList& faceEdges(lc.faceEdges());
    const faceList& fL = lc.faces();

    forAll (lc.cc(), facei)
    {
        const face& f( fL[facei] );
        const edgeList& eLt  (f.edges());

        labelList edgeCutt(f.size());

        const labelList& fE( faceEdges[facei] );

        forAll (eLt, edgei)
        {
            edgeCutt[edgei] = edgeCut[ fE[edgei] ];
        }

        edge& cutted( cuttedEdges[cutCount]);

        faceCut(pType, f, eLt, edgeCutt, lc.points(), lf, cutted);

        if (lf.isNegFace() && lf.isPosFace())
        {
            lc.addNeg( lf.negFace() );

            cutCount++;
        }
        else if (lf.isNegFace())
        {
            lc.addNeg( facei );
        }
        else
        {
        }

        if (lf.noi() > 2)
        {
            noiProblem = true;
        }

    }

    if (noiProblem)
    {
        lc.initCut();

        if (signedPointToSurfaceDistance(lc.centre() ) < 0.0)
        {
            lc.fullNeg();
        }
        //else lc is a negative cell, hence further work is redundant
    }
    else
    {
        // Add special edges, where both ends are on the plane
        forAll (eL, edgei)
        {
            if (pType[eL[edgei].start()] == 0 && pType[eL[edgei].end()] == 0)
            {
                cuttedEdges[cutCount++] = eL[edgei];
            }
        }

        // Take care of the interface face
        cuttedEdges.setSize(cutCount);

        if (cuttedEdges.size() >= 3)
        {
            face f = combineEdgeList(cuttedEdges);

            if (Foam::sign( f.normal( lc.points()) & normalToPlane_ ) == 1)
            {
                lc.addNeg(f);
            }
            else
            {
                lc.addNeg(f.reverseFace());
            }
        }
    }
}


void convexPolyhedral::faceCut
(
    const labelList& pType,
    const face& f,
    const edgeList& eL,
    const labelList& edgeCut,
    const pointField& pf,
    localFace& lf,
    edge& cutted
)
{
    lf.points(pf);
    lf.noi(0);

    // Used for checking the type of intersection. This part is needed
    // because pType.size() != f.size() when polyhedral intersections
    // are considered.
    label pTypeSqr(0);
    label pTypeSum(0);

    forAll (f, pointi)
    {
        pTypeSqr += ( pType[f[pointi]]*pType[f[pointi]] );
        pTypeSum += pType[f[pointi]];
    }

    // Faces completely on the surface is added to both positive and
    // negative side
    if (pTypeSqr == 0)
    {
        lf.posFace() = f;
        lf.negFace() = f;
    }
    // Only the negative side has non-zero face properties
    else if (pTypeSum == - pTypeSqr)
    {
        lf.negFace() = f;
        lf.posFace().setSize(0);
    }
    // Only the positive side has non-zero face properties
    else if (pTypeSum == pTypeSqr)
    {
        lf.negFace().setSize(0);
        lf.posFace() = f;
    }
    else
    {
        lf.posFace().setSize( f.size() + 2);
        lf.negFace().setSize( f.size() + 2);
        face& posF( lf.posFace() );
        face& negF( lf.negFace() );
        label posCount(0), negCount(0), noi(0);

        forAll (eL, edgei)
        {
            edge e(eL[edgei]);
            label start(e.start());

            if (pType[start] == -1)
            {
                negF[negCount++] = start;
            }
            else if (pType[start] == 1)
            {
                posF[posCount++] = start;
            }
            else
            {
                negF[negCount++] = start;
                posF[posCount++] = start;

                if (noi < 2)
                {
                    noi == 0 ? cutted.start() = start : cutted.end() = start;
                }

                noi++;
            }

            if (edgeCut[edgei] != -1)
            {
                negF[negCount++] = edgeCut[edgei];
                posF[posCount++] = edgeCut[edgei];

                if (noi < 2)
                {
                    noi == 0 ? cutted.start() = edgeCut[edgei] :
                         cutted.end() = edgeCut[edgei];
                }

                noi++;
            }
        }

        negF.setSize(negCount);
        posF.setSize(posCount);

        lf.noi(noi);
    }
}


face convexPolyhedral::combineEdgeList
(
    const edgeList& eL
)
{
    face fInter(eL.size());

    fInter[0] = eL[0].start();
    fInter[1] = eL[0].end();

    label count(1);
    label prevEdge(-1);

    while (true)
    {
        for (label i = 1; i < eL.size(); i++)
        {
            if (i != prevEdge)
            {
                if (eL[i].start() == fInter[count])
                {
                    fInter[++count] = eL[i].end();
                    prevEdge = i;
                    break;
                }
                else if (eL[i].end() == fInter[count])
                {
                    fInter[++count] = eL[i].start();
                    prevEdge = i;
                    break;
                }
            }
        }

        if (count == eL.size() - 1)
        {
            break;
        }
    }

    return fInter;

}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //


void convexPolyhedral::signedPointToSurfaceDistance
//void convexPolyhedral::signedPointToSurfaceDistance
(
    const pointField& pp,
    scalarField& signedDistance
)
{
    if (pp.size() != signedDistance.size())
    {
        signedDistance.setSize( pp.size() );
    }

    forAll (pp, pointi)
    {
        signedDistance[pointi] = signedPointToSurfaceDistance( pp[pointi] );
    }
}


scalar convexPolyhedral::signedPointToSurfaceDistance
(
    const point& p
) const
{
    return (p - pointOnPlane_) & normalToPlane_;
}


// This avoids having to work with very small numbers, hence if a given
// distance is less than 5.0e-14, it is assumed that the point is lying
// on the surface and its labelList attribute is 0. Else the attribute
// is the sign of the distance.
void convexPolyhedral::floatingPointToLabel
(
    const scalarField& s,
    labelList& l
)
{
    forAll (s, pointi)
    {
        l[pointi] = floatingPointToLabel( s[pointi] );
    }
}


label convexPolyhedral::floatingPointToLabel
(
    const scalar& s
)
{
    return Foam::mag(s) < 5.0e-14 ? 0 : Foam::sign(s);
}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //


localFace convexPolyhedral::divideFace
(
    const label& faceLabel,
    const point& pointOnPlane,
    const vector& normalToPlane
)
{
    // Update private member functions
    setPoint ( pointOnPlane  );
    setNormal( normalToPlane );

    const pointField& pp(mesh_.points());
    const faceList& fL(mesh_.faces());

    // Make a local face representation
    face f( fL[faceLabel] );
    pointField pf( f.size() );

    forAll (f, pointi)
    {
        pf[pointi] = pp[f[pointi]];
        f[pointi]  = pointi;
    }

    // Get the signed distance and point type
    scalarField sD   ( f.size(), 0.0);
    labelList   pType( f.size(), 0.0);

    signedPointToSurfaceDistance(pf, sD);
    floatingPointToLabel( sD, pType );

    // Add potential new points to pf and return their labels in edgeCut
    const edgeList& eL(f.edges());
    labelList edgeCut = edgeCutLabel(eL, pType, sD, pf);

    // Return variable. To be populated below - constructed to empty
    localFace lf;
    edge cutted;

    faceCut(pType, f, eL, edgeCut, pf, lf, cutted);

    return lf;
}


localFace convexPolyhedral::divideFace
(
    const label& faceLabel
)
{
    const pointField& pp(mesh_.points());
    const faceList& fL(mesh_.faces());

    // Make a local face representation
    face f( fL[faceLabel] );
    pointField pf( f.size() );

    forAll (f, pointi)
    {
        pf[pointi] = pp[f[pointi]];
        f[pointi]  = pointi;
    }

    // Get the signed distance and point type
    scalarField sD   ( f.size(), 0.0);
    labelList   pType( f.size(), 0.0);

    signedPointToSurfaceDistance(pf, sD);
    floatingPointToLabel( sD, pType );

    // Add potential new points to pf and return their labels in edgeCut
    const edgeList& eL(f.edges());
    labelList edgeCut = edgeCutLabel(eL, pType, sD, pf);

    // Return variable. To be populated below - constructed to empty
    localFace lf;
    edge cutted;

    faceCut(pType, f, eL, edgeCut, pf, lf, cutted);

    return lf;
}


localCell convexPolyhedral::dividePolyhedral
(
    const label& cellLabel,
    const point& pointOnPlane,
    const vector& normalToPlane
)
{
    // Create a localised cell
    localCell lc(mesh_, cellLabel);

    dividePolyhedral(pointOnPlane, normalToPlane, lc);

    return lc;
}


void convexPolyhedral::dividePolyhedral
(
    const point& pointOnPlane,
    const vector& normalToPlane,
    localCell& lc
)
{
    lc.clearCut();

    // Update private member functions
    setPoint ( pointOnPlane  );
    setNormal( normalToPlane );

    // Get data from localised cell
    const edgeList& eL( lc.edges() );
    pointField& pp( lc.points() );

    // Get the signed distance and point type
    scalarField sD   ( pp.size(), 0.0);
    labelList   pType( pp.size(), 0.0);

    signedPointToSurfaceDistance(pp, sD);
    floatingPointToLabel( sD, pType );

    if (Foam::sum( pType) == - Foam::sum( pType*pType ))
    {
        lc.fullNeg();
    }
    else if (Foam::sum( pType ) == Foam::sum( pType*pType ))
    {
        lc.fullPos();
    }
    else
    {
        // Add potential new points to pf and return their labels in edgeCut
        labelList edgeCut = edgeCutLabel(eL, pType, sD, pp);

        faceCut(pType, eL, edgeCut, lc);
    }

    lc.doneCut();
}


void convexPolyhedral::dividePolyhedral
(
    const point& pointOnPlane,
    const vector& normalToPlane,
    localCellNeg& lc
)
{
    lc.clearCut();

    // Update private member functions
    setPoint ( pointOnPlane  );
    setNormal( normalToPlane );

    // Get data from localised cell
    const edgeList& eL( lc.edges() );
    pointField& pp( lc.points() );

    // Get the signed distance and point type
    scalarField sD   ( pp.size(), 0.0);
    labelList   pType( pp.size(), 0.0);

    signedPointToSurfaceDistance(pp, sD);
    floatingPointToLabel( sD, pType );

    if (Foam::sum( pType) == - Foam::sum( pType*pType ))
    {
        lc.fullNeg();
    }
    else if (Foam::sum( pType ) == Foam::sum( pType*pType ))
    {
        // Nothing to be done
        // The cell is fully positive, and since localCellNeg only contains the
        // negative part of the intersection, further work is redundant
    }
    else
    {
        // Add potential new points to pf and return their labels in edgeCut
        labelList edgeCut = edgeCutLabel(eL, pType, sD, pp);

        faceCut(pType, eL, edgeCut, lc);
    }


    lc.doneCut();
}


void convexPolyhedral::unionSet
(
    const localCell& cell0,
    localCell& cell1
)
{
    const cell& c( cell0.ccNeg() );
    const faceList& fL( cell0.faces() );
    const pointField& pp( cell0.points() );

    forAll (c, facei)
    {
        const vector n( fL[c[facei]].normal(pp) );
        const point Cf( fL[c[facei]].centre(pp) );

        dividePolyhedral(Cf, n, cell1);

        if (cell1.ccNeg().size())
        {
            cell1.localizeCell("neg");
        }
        else
        {
            cell1.emptyCell();
            break;
        }

    }
}


void convexPolyhedral::unionSet
(
    const localCellNeg& cell0,
    localCellNeg& cell1
)
{
    const cell& c( cell0.ccNeg() );
    const faceList& fL( cell0.faces() );
    const pointField& pp( cell0.points() );

    forAll (c, facei)
    {
        const vector n( fL[c[facei]].normal(pp) );
        const point Cf( fL[c[facei]].centre(pp) );

        dividePolyhedral(Cf, n, cell1);

        if (cell1.ccNeg().size())
        {
            cell1.localizeCell("neg");
        }
        else
        {
            cell1.emptyCell();
            break;
        }
    }
}

} // End namespace



