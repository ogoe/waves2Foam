/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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
#include "volFields.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::sampledSurfaceElevation::volFieldSampler<Type>::volFieldSampler
(
    const word& interpolationScheme,
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const PtrList<sampledSet>& samplers
)
:
    List<Field<Type>>(samplers.size()),
    name_(field.name())
{
    autoPtr<interpolation<Type>> interpolator
    (
        interpolation<Type>::New(interpolationScheme, field)
    );

    forAll(samplers, setI)
    {
        Field<Type>& values = this->operator[](setI);
        const sampledSet& samples = samplers[setI];

        values.setSize(samples.size());
        forAll(samples, sampleI)
        {
            const point& samplePt = samples[sampleI];
            label celli = samples.cells()[sampleI];
            label facei = samples.faces()[sampleI];

            if (celli == -1 && facei == -1)
            {
                // Special condition for illegal sampling points
                values[sampleI] = pTraits<Type>::max;
            }
            else
            {
                values[sampleI] = interpolator().interpolate
                (
                    samplePt,
                    celli,
                    facei
                );
            }
        }
    }
}


template<class Type>
Foam::sampledSurfaceElevation::volFieldSampler<Type>::volFieldSampler
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const PtrList<sampledSet>& samplers
)
:
    List<Field<Type>>(samplers.size()),
    name_(field.name())
{
    forAll(samplers, setI)
    {
        Field<Type>& values = this->operator[](setI);
        const sampledSet& samples = samplers[setI];

        values.setSize(samples.size());
        forAll(samples, sampleI)
        {
            label celli = samples.cells()[sampleI];

            if (celli ==-1)
            {
                values[sampleI] = pTraits<Type>::max;
            }
            else
            {
                values[sampleI] = field[celli];
            }
        }
    }
}


template<class Type>
Foam::sampledSurfaceElevation::volFieldSampler<Type>::volFieldSampler
(
    const List<Field<Type>>& values,
    const word& name
)
:
    List<Field<Type>>(values),
    name_(name)
{}


template<class Type>
Foam::fileName Foam::sampledSurfaceElevation::writeSampleFile
(
    const coordSet& masterSampleSet,
    const PtrList<volFieldSampler<Type>>& masterFields,
    const label setI,
    const fileName& timeDir,
    const writer<Type>& formatter
)
{
    wordList valueSetNames(masterFields.size());
    List<const Field<Type>*> valueSets(masterFields.size());

    forAll(masterFields, fieldi)
    {
        valueSetNames[fieldi] = masterFields[fieldi].name();
        valueSets[fieldi] = &masterFields[fieldi][setI];
    }

    fileName fName
    (
        timeDir/formatter.getFileName(masterSampleSet, valueSetNames)
    );

    OFstream ofs(fName);
    if (ofs.opened())
    {
        formatter.write
        (
            masterSampleSet,
            valueSetNames,
            valueSets,
            ofs
        );
        return fName;
    }
    else
    {
        WarningInFunction
            << "File " << ofs.name() << " could not be opened. "
            << "No data will be written" << endl;
        return fileName::null;
    }
}


template<class T>
void Foam::sampledSurfaceElevation::combineSampledValues
(
    const PtrList<volFieldSampler<T>>& sampledFields,
    const labelListList& indexSets,
    PtrList<volFieldSampler<T>>& masterFields
)
{
    forAll(sampledFields, fieldi)
    {
        List<Field<T>> masterValues(indexSets.size());

        forAll(indexSets, setI)
        {
            // Collect data from all processors
            List<Field<T>> gatheredData(Pstream::nProcs());
            gatheredData[Pstream::myProcNo()] = sampledFields[fieldi][setI];
            Pstream::gatherList(gatheredData);

            if (Pstream::master())
            {
                Field<T> allData
                (
                    ListListOps::combine<Field<T>>
                    (
                        gatheredData,
                        Foam::accessOp<Field<T>>()
                    )
                );

                masterValues[setI] = UIndirectList<T>
                (
                    allData,
                    indexSets[setI]
                )();
            }
        }

        masterFields.set
        (
            fieldi,
            new volFieldSampler<T>
            (
                masterValues,
                sampledFields[fieldi].name()
            )
        );
    }
}


template<class Type>
void Foam::sampledSurfaceElevation::sampleAndWrite(fieldGroup<Type>& fields)
{
	bool interpolate = interpolationScheme_ != "cell";

	// Create or use existing writer
	if (fields.formatter.empty())
	{
	    fields = writeFormat_;
	}

	// Storage for interpolated values - brute force size to one!
	PtrList<volFieldSampler<Type>> sampledFields(1);

	if (Pstream::master() && verbose_)
	{
		Pout<< "surfaceElevation::sampleAndWrite: "
				<< Foam::waves2Foam::aName() << endl;
	}

	if (loadFromFiles_)
	{
		GeometricField<Type, fvPatchField, volMesh> vf
		(
			IOobject
			(
				Foam::waves2Foam::aName(),
				mesh_.time().timeName(),
				mesh_,
				IOobject::MUST_READ,
				IOobject::NO_WRITE,
				false
			),
			mesh_
		);

		if (interpolate)
		{
			sampledFields.set
			(
			    0,
			    new volFieldSampler<Type>
			(
			    interpolationScheme_,
			    vf,
			    *this
			)
			);
		}
		else
		{
			sampledFields.set
			(
				0,
				new volFieldSampler<Type>(vf, *this)
			);
		}
	}
	else
	{
		if (interpolate)
		{
			sampledFields.set
			(
				0,
				new volFieldSampler<Type>
			(
				interpolationScheme_,
				mesh_.lookupObject
				<GeometricField<Type, fvPatchField, volMesh>>
				(Foam::waves2Foam::aName()),
				*this
			)
			);
		}
		else
		{
			sampledFields.set
			(
				0,
				new volFieldSampler<Type>
			(
				mesh_.lookupObject
				<GeometricField<Type, fvPatchField, volMesh>>
				(Foam::waves2Foam::aName()),
				*this
			)
			);
		}
	}

	// Combine sampled fields from processors.
	// Note: only master results are valid

	PtrList<volFieldSampler<Type>> masterFields(sampledFields.size());
	combineSampledValues(sampledFields, indexSets_, masterFields);

    if (Pstream::master())
	{
        // Make a field to collect all results
		scalarField result(masterSampledSets_.size(), 0.0);

		// Run over all sampled sets and perform integration
    	forAll (masterSampledSets_, seti)
		{
    		const coordSet & cs(masterSampledSets_[seti]);

    		List<const Field<scalar>*> valueSets(masterFields.size());
    		valueSets[0] = &masterFields[0][seti];

    		List<const Field<scalar>*> columns(valueSets.size());
    		columns[0] = valueSets[0];

    		const Field<scalar>& alpha = *columns[0];
    		scalar tolerance(0.0001);

    		// Write "-GREAT" is the integration set is less than 2 points
    		if (alpha.size() < 2)
    		{
    			result[seti] = -GREAT;
    		}
    		// Write "-GREAT" if both points are above or below the water
    		// surface
    		else if
    		(
    			(alpha[0] < tolerance && alpha[alpha.size()-1] < tolerance)
    			||
    			(
    			    alpha[0] > 1.0 - tolerance &&
    			    alpha[alpha.size()-1] > 1.0 - tolerance
    			)
    		)
    		{
    			result[seti] = -GREAT;
    		}
    		// Perform the integration
    		else
    		{
    			scalar value(0);
    			scalar minScalarCoord(cs.scalarCoord(0));

    			for (int pointi=0; pointi < alpha.size() - 1; pointi++)
    			{
    				value +=
    				(
    				    cs.scalarCoord(pointi + 1)
    				  - cs.scalarCoord(pointi)
    				)
    				*(alpha[pointi + 1] + alpha[pointi]);

    				minScalarCoord =
    			        Foam::min
    					(
    					    minScalarCoord,
    					    cs.scalarCoord(pointi + 1)
    					);
    			}
    			value *= 0.5;

    			result[seti] = value + minScalarCoord;
    		}
		}

    	// Write time - Enforce 10 digit time precision for long time
    	// series
    	label precision = surfaceElevationFilePtr_().precision(10);
        surfaceElevationFilePtr_() << mesh_.time().time().value();
        surfaceElevationFilePtr_().precision(precision);

        // Write the surface elevation
        forAll (result, seti)
        {
            surfaceElevationFilePtr_() << "\t" << result[seti];
        }

        surfaceElevationFilePtr_() << "\n";
	}

}


// ************************************************************************* //
