/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sampledSurfaceElevation.H"
#include "globalIndex.H"
#include "interpolation.H"
#include "volFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
Foam::tmp<GeoField>
Foam::sampledSurfaceElevation::getOrLoadField(const word& fieldName) const
{
    tmp<GeoField> tfield;

    if (loadFromFiles_)
    {
        tfield.reset
        (
            new GeoField
            (
                IOobject
                (
                    fieldName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_
            )
        );
    }
    else
    {
        // Slightly paranoid here
        tfield.cref(mesh_.cfindObject<GeoField>(fieldName));
    }

    return tfield;
}


template<class Type>
void Foam::sampledSurfaceElevation::performAction
(
    const VolumeField<Type>& fld,
    unsigned request
)
{
    const scalar timeValue = fld.time().timeOutputValue();

    // The interpolator for this field
    autoPtr<interpolation<Type>> interpPtr;

    if (!samplePointScheme_.empty() && samplePointScheme_ != "cell")
    {
        interpPtr.reset(interpolation<Type>::New(samplePointScheme_, fld));
    }

    const unsigned int width(IOstream::defaultPrecision() + 7);

    if (Pstream::master())
    {
        surfaceElevationFilePtr_() << timeValue;
    }

    // Find the values for the cells
    forAll(*this, seti)
    {
        const sampledSet& s = (*this)[seti];
        const globalIndex& globIdx = globalIndices_[seti];
        const labelList& globOrder = gatheredSorting_[seti];

        Field<Type> values(s.size());
        Field<scalar> coords(s.size(), 0);

        if (interpPtr)
        {
            forAll(s, samplei)
            {
                const point& p = s[samplei];
                const label celli = s.cells()[samplei];
                const label facei = s.faces()[samplei];

                coords[samplei] = s.scalarCoord(samplei);

                if (celli == -1 && facei == -1)
                {
                    // Special condition for illegal sampling points
                    values[samplei] = pTraits<Type>::max;
                }
                else
                {
                    values[samplei] = interpPtr().interpolate(p, celli, facei);
                }
            }
        }
        else
        {
            forAll(s, samplei)
            {
                const label celli = s.cells()[samplei];

                coords[samplei] = s.scalarCoord(samplei);

                if (celli == -1)
                {
                    values[samplei] = pTraits<Type>::max;
                }
                else
                {
                    values[samplei] = fld[celli];
                }
            }
        }

        // Collect data from all processors
        globIdx.gatherInplace(values);
        globIdx.gatherInplace(coords);

        // Perform the surface integration of the alpha field on master
        if (Pstream::master())
        {
            // Use sorted order
            values = UIndirectList<Type>(values, globOrder)();
            coords = UIndirectList<scalar>(coords, globOrder)();

            // Trapezoidal integration
            Type eta = pTraits<Type>::zero;

            scalar tolerance(0.0001);

            // Write "-GREAT" is the integration set is less than 2 points
            if (values.size() < 2)
            {
                eta = -GREAT*pTraits<Type>::one;
            }
            // Write "-GREAT" if both points are above or below the water
            // surface
            else if
            (
               (values[0] < tolerance && values[values.size()-1] < tolerance)
               ||
               (
                   values[0] > 1.0 - tolerance &&
                   values[values.size()-1] > 1.0 - tolerance
               )
            )
            {
                eta = -GREAT*pTraits<Type>::one;
            }
            else
            {
                for (int ii = 0; ii < values.size() - 2; ii++)
                {
                    scalar dist(coords[ii + 1] - coords[ii]);

                    eta += 0.5*(values[ii + 1] + values[ii])*dist;
                }

                // Add the bottom point
                eta += pTraits<Type>::one*coords[0];
            }
            // Write the field to the output file
            surfaceElevationFilePtr_() << tab << eta;
        }
    }

    // Finish line
    if (Pstream::master())
    {
        surfaceElevationFilePtr_() << endl;
    }
}


template<class GeoField>
void Foam::sampledSurfaceElevation::performAction
(
    const IOobjectList& objects,
    unsigned request
)
{
    wordList fieldNames;
    if (loadFromFiles_)
    {
        fieldNames = objects.sortedNames<GeoField>(fieldSelection_);
    }
    else
    {
        fieldNames = mesh_.thisDb().sortedNames<GeoField>(fieldSelection_);
    }

    for (const word& fieldName : fieldNames)
    {
        tmp<GeoField> tfield = getOrLoadField<GeoField>(fieldName);

        if (tfield)
        {
            performAction<typename GeoField::value_type>(tfield(), request);
        }
    }
}


// ************************************************************************* //
