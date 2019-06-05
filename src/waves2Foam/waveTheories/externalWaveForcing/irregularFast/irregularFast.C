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

\*---------------------------------------------------------------------------*/

#include "irregularFast.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(irregularFast, 0);

addToRunTimeSelectionTable
(
    externalWaveForcing,
    irregularFast,
    externalWaveForcing
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


irregularFast::irregularFast
(
    IOobject io,
    Time& rT,
    const fvMesh& mesh
)
:
    externalWaveForcing(io, rT, mesh),

    waveProps_(io.db().lookupObject<IOdictionary>("waveProperties")),
    coeffDict_(waveProps_.subDict("externalForcingCoeffs")),

    ignoreMeshMotion_(Switch(waveProps_.lookup("ignoreMeshMotion"))),

    relaxationNames_(waveProps_.lookup("relaxationNames")),

    N_(readLabel(coeffDict_.lookup("N"))),

    h_(readScalar(coeffDict_.lookup("depth"))),

    amp_("amplitude", coeffDict_, N_),

    omega_("frequency", coeffDict_, N_),

    phi_("phaselag", coeffDict_, N_),

    k_("waveNumber", coeffDict_, N_),

    Tsoft_( readScalar(coeffDict_.lookup("Tsoft"))),

    seaLevel_(readScalar(waveProps_.lookup("seaLevel"))),

    rampFactor_(factor(rT_.time().value())),

    tolerance_(readScalar(waveProps_.lookup("searchTolerance")))
{
	// Find the relaxation zone names, which uses the externalSource as forcing
    externalSourceZones_.setSize(relaxationNames_.size(), -1);
    label count(0);

    forAll (relaxationNames_, zonei)
    {
        const dictionary& sd
            (
                waveProps_.subDict(relaxationNames_[zonei] + "Coeffs")
            );

        if (word(sd.lookup("waveType")) == "externalSource")
        {
             externalSourceZones_[count++] = zonei;
        }
    }

    externalSourceZones_.setSize(count);

    // Set additional wave properties
    K_.setSize(N_, 0.0);
    compDir_.setSize(N_, vector::zero);
    period_.setSize(N_, 0);
    velAmp_.setSize(N_, 0);

    // Get the cyclic frequency
    omega_ *= (2.0*M_PI);

    // Compute the length of k_
    K_ = Foam::mag(k_);

    // Get the direction of each wave component
    compDir_ = k_ / K_;

    // Compute the period
    forAll (period_, index)
    {
        period_[index] = 2.0*M_PI/omega_[index];
    }

    // Compute the velocity amplitude
    forAll (velAmp_, index)
    {
        velAmp_[index] = M_PI*2.0*amp_[index]/period_[index]
            /Foam::sinh(K_[index]*h_);
    }

    scalar timeBefore = mesh_.time().elapsedCpuTime();

    // Make unit vectors
    calcUnitVectors();

    // Calculate the unique coordinates, which makes the execution of the
    // irregular wave summation faster
    calcUniqueCoordinates();

    // Initialise all components
    calcAllFields();

    // Info on the initialisation procedure
    Info << "Constructing irregularFast in " << mesh_.time().elapsedCpuTime()
    	 - timeBefore << " s.\n" << endl;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar irregularFast::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (0.0 < Tsoft_ && time < Tsoft_)
    {
        factor = Foam::sin(2*M_PI/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}


void irregularFast::calcUnitVectors()
{
	// Get the gravity field
#if OFPLUSBRANCH==1
    #if OFVERSION<1812
        vector g(uniformDimensionedVectorField
            (
                mesh_.thisDb().lookupObject<uniformDimensionedVectorField>("g")
            ).value());
    #else
        vector g(Foam::meshObjects::gravity::New(mesh_.thisDb().time()).value());
    #endif
#else
        vector g(uniformDimensionedVectorField
            (
                mesh_.thisDb().lookupObject<uniformDimensionedVectorField>("g")
            ).value());
#endif
//    vector g =
//    	uniformDimensionedVectorField
//        (
//            mesh_.thisDb().lookupObject<uniformDimensionedVectorField>("g")
//        ).value();

    // The vertical coordinate is based on the gravity field
    eZ_ = Foam::cmptMultiply(g, g);
    eZ_ /= Foam::mag(eZ_);

    // Assume that unit in x-direction is actually x and correct for any
    // projection on the vertical axis
    eX_ = vector(1, 0, 0);

    if (Foam::mag(eZ_ & eX_) > SMALL)
    {
    	eX_ -= (eZ_ & eX_)*eZ_;
    	eX_ /= Foam::mag(eX_);
    }

    // Define the unit-y as the cross product between these two unit vectors
    eY_ = eX_^eZ_;
    eY_ /= Foam::mag(eY_);
}


void irregularFast::calcUniqueCoordinates()
{
	// Set the size of the horizontal coordinate list
    horizontalCoordinates_.setSize(2);

    // Get references to each primary coordinate
    scalarField& x = horizontalCoordinates_[0];
    scalarField& y = horizontalCoordinates_[1];
    scalarField& z = verticalCoordinate_;

    // Prepare some counters
    label count(0), countZ(0);

    // Re-initialise the lists
    x.setSize(10000, 0);
    y.setSize(10000, 0);
    z.setSize(10000, 0);

    // Get refernece to the coordinates
    const vectorField& C = mesh_.C().internalField();

    // Reference to add the points in the horizontal plane
    const labelListList& cellPoints = mesh_.cellPoints();
    const pointField& pp = mesh_.points();

    // Get the coordinates of all cells in the relaxation zones
    forAll (externalSourceZones_, zonei)
    {
    	autoPtr<relaxationShapes::relaxationShape> rs =
            Foam::relaxationShapes::relaxationShape::New
            (
                relaxationNames_[externalSourceZones_[zonei]], mesh_
            );

    	const labelList cells = rs->cells();

    	forAll(cells, celli)
    	{
            x[count] = (C[cells[celli]] & eX_);
            y[count++] = (C[cells[celli]] & eY_);
            z[countZ++] = (C[cells[celli]] & eZ_);

            // Resize if necessary
            if (count == x.size())
            {
            	x.setSize(2*x.size());
            	y.setSize(2*y.size());
            }

            if (countZ == z.size())
            {
            	z.setSize(2*z.size());
            }

            // Add the points in the horizontal plane
            const labelList& cp = cellPoints[cells[celli]];

            forAll (cp, pointi)
            {
            	x[count] = (pp[cp[pointi]] & eX_);
            	y[count++] = (pp[cp[pointi]] & eY_);

            	// Resize if necessary
            	if (count == x.size())
            	{
            		x.setSize(2*x.size());
            		y.setSize(2*y.size());
            	}
            }

    	}
    }

    // Finally, truncate the list of coordinates
    x.setSize(count);
    y.setSize(count);
    z.setSize(countZ);

    // Unique-sorting of the lists
    this->uniqueSorting(x);
    this->uniqueSorting(y);
    this->uniqueSorting(z);
}


void irregularFast::uniqueSorting(scalarField& in)
{
	// No sorting needed (possible) if the size of the input field is zero.
	if (in.size() == 0)
	{
		return;
	}

	// Make a copy of the in-field
    List<scalar> final(in);

    // Sort the values
    std::sort(final.begin(), final.end());

    // Find unique values (within a search tolerance)
    labelList indices(final.size(), 0);
    label count(1);

    for (label i = 1; i < final.size(); i++)
    {
    	if (Foam::mag(final[i] - final[i - 1]) > tolerance_)
    	{
            indices[count++] = i;
    	}
    }

    // Get the indices of the unique coordinates
    indices.setSize(count);

    // Correct the field to return sorted and unique values
    in.setSize(indices.size());

    forAll (indices, i)
    {
    	in[i] = final[indices[i]];
    }
}


void irregularFast::calcAllFields()
{
    calcHyperbolicFunctions();

    calcWaveNumberTrigonometricFunctions();

    calcWaveFrequencyTrigonometricFunctions();

    calcSurfaceElevation();

    calcVelocityComponents();
}


void irregularFast::calcHyperbolicFunctions()
{
	// Set the size of the hyperbolic functions. One for each vertical
	// coordinate
    Cosh_.setSize(verticalCoordinate_.size());
    Sinh_.setSize(verticalCoordinate_.size());

    // Make a local coordinate system, which is 0 at zero level
    scalarField Z = verticalCoordinate_ - seaLevel_;

    // Loop over all vertical coordinates
    forAll (Z, zi)
    {
    	// Set the size of the fields
    	scalarField& fieldC = Cosh_[zi];
    	scalarField& fieldS = Sinh_[zi];

    	fieldC.setSize(K_.size(), 0.0);
    	fieldS.setSize(K_.size(), 0.0);

    	// Compute the cosh and sinh contributions
        fieldC = Foam::cosh(K_*(Z[zi] + h_));
        fieldS = Foam::sinh(K_*(Z[zi] + h_));
    }
}


void irregularFast::calcWaveNumberTrigonometricFunctions()
{
	// Compute the two horizontal components of the wave number
    scalarField kx = (k_ & eX_);
    scalarField ky = (k_ & eY_);

    const scalarField& x = horizontalCoordinates_[0];
    const scalarField& y = horizontalCoordinates_[1];

    // Set the size of the fields.
    // Level 1: X
    // Level 2: Y
    // Level 3: k-field
    CosK_.setSize(x.size());
    SinK_.setSize(x.size());

    forAll (x, xi)
    {
    	List<scalarField>& fieldC = CosK_[xi];
    	List<scalarField>& fieldS = SinK_[xi];

    	fieldC.setSize(y.size());
    	fieldS.setSize(y.size());

    	forAll (y, yi)
    	{
    		scalarField& fC = fieldC[yi];
    		scalarField& fS = fieldS[yi];

    		fC.setSize(kx.size(), 0.0);
    		fS.setSize(kx.size(), 0.0);

    		fC = Foam::cos(kx*x[xi] + ky*y[yi]);
    		fS = Foam::sin(kx*x[xi] + ky*y[yi]);
    	}
    }
}


void irregularFast::calcWaveFrequencyTrigonometricFunctions()
{
    CosOmega_.setSize(omega_.size(), 0.0);
    SinOmega_.setSize(omega_.size(), 0.0);

    forAll (omega_, indexi)
    {
        CosOmega_[indexi] = Foam::cos(omega_[indexi]*rT_.time().value() + phi_[indexi]);
        SinOmega_[indexi] = Foam::sin(omega_[indexi]*rT_.time().value() + phi_[indexi]);
    }
}


void irregularFast::calcSurfaceElevation()
{
    const scalarField& x = horizontalCoordinates_[0];
    const scalarField& y = horizontalCoordinates_[1];

	Eta_.setSize(x.size());

    forAll (Eta_, nX)
    {
    	scalarField& eta = Eta_[nX];
        eta.setSize(y.size());

    	forAll (eta, nY)
    	{
    		// Get the reference to the trigonometric functions
    		const scalarField& Sin = SinK_[nX][nY];
    		const scalarField& Cos = CosK_[nX][nY];

    		// Loop over all components
    		eta[nY] = Foam::sum(amp_*(CosOmega_*Cos + SinOmega_*Sin));

    	    eta[nY] *= rampFactor_;
    	    eta[nY] += seaLevel_;
    	}
    }
}


void irregularFast::calcVelocityComponents()
{
    const scalarField& x = horizontalCoordinates_[0];
    const scalarField& y = horizontalCoordinates_[1];

    Ucomponent_.setSize(x.size());
    Wcomponent_.setSize(x.size());

    forAll (Ucomponent_, nX)
    {
     	List<vectorField>& Ucmp = Ucomponent_[nX];
     	List<scalarField>& Wcmp = Wcomponent_[nX];

     	Ucmp.setSize(y.size());
     	Wcmp.setSize(y.size());

        forAll (Ucmp, nY)
     	{
        	// Get the reference to the trigonometric functions
        	const scalarField& Sin = SinK_[nX][nY];
        	const scalarField& Cos = CosK_[nX][nY];

        	vectorField& U = Ucmp[nY];
        	scalarField& W = Wcmp[nY];

            U.setSize(Sin.size());
            W.setSize(Sin.size());

            U = compDir_*velAmp_*(CosOmega_*Cos + SinOmega_*Sin);
            W = velAmp_*(SinOmega_*Cos - CosOmega_*Sin);
     	}

    }

//	// Get the reference to the trigonometric functions
//	const scalarField& Sin = SinK_[nX][nY];
//	const scalarField& Cos = CosK_[nX][nY];
//
//	U = Foam::sum(compDir_*velAmp_*Cosh*(CosOmega_*Cos + SinOmega_*Sin))
//	    - eZ_*(Foam::sum(velAmp_*Sinh*(SinOmega_*Cos - CosOmega_*Sin)));

}


label irregularFast::bisectionFindIndex
(
    const scalar& val,
    const scalarField& field
) const
{
    label res = -1;

    label Nf = field.size() - 1;

    label N0 = 0, NN = Nf;

    // Check the one end point
    if (Foam::mag(val - field[N0]) < tolerance_)
    {
    	res = N0;
    	return res;
    }

    // Check the other end point
    if (Foam::mag(val - field[NN]) < tolerance_)
    {
    	res = NN;
    	return res;
    }

    // Only perform the search, if the point is in the search interval
    if (val < field[0] || field[NN] < val)
    {
    	return res;
    }

    // Perform a bi-section search in terms on the indices in the list
    while (true)
    {
    	// Find the average index
    	label avI = N0 + (NN - N0)/2;

    	// If this is the solution, return the index
    	if (Foam::mag(val - field[avI]) < tolerance_)
    	{
    		res = avI;
    		break;
    	}

    	// Update the upper/lower bounds
    	if (val < field[avI])
    	{
    		NN = avI;
    	}
    	else
    	{
    		N0 = avI;
    	}

    	// If the gap between the indices is less than two, stop the search
    	if (Foam::mag(NN - N0) < 2)
    	{
    		break;
    	}
    }

    // Return the index
    return res;
}


void irregularFast::findIndexing
(
    const point& p,
    label& nX,
    label& nY,
    label& nZ,
    bool checkVertical
) const
{
	// Get reference to the primary coordinate fields
	const scalarField& x = horizontalCoordinates_[0];
	const scalarField& y = horizontalCoordinates_[1];
	const scalarField& z = verticalCoordinate_;

	// Compute the primary coordinates of the point 'p'
	scalar X = (p & eX_);
	scalar Y = (p & eY_);
	scalar Z = (p & eZ_);

	// Find the x-label
	nX = bisectionFindIndex(X, x);

	// Find the y-label
    if (y.size() == 3)
	{
		nY = 1;
	}
    else
    {
    	nY = bisectionFindIndex(Y, y);
    }

    // Find the z-level
    if (checkVertical)
    {
    	nZ = bisectionFindIndex(Z, z);
    }
    else
    {
    	nZ = -1;
    }
}


scalar irregularFast::interpolateSurfaceElevation
(
    const point& p,
    const scalar& t,
    const label& nX,
    const label& nY
) const
{
	scalar eta(0);

	scalar x = (p & eX_);
	scalar y = (p & eY_);

	label Nx = horizontalCoordinates_[0].size();
	label Ny = horizontalCoordinates_[1].size();

    if (
    	   horizontalCoordinates_[0][0] < x &&
    	   x < horizontalCoordinates_[0][Nx - 1] &&
    	   Ny == 3
       )
    {
        label nLower = 0;

        forAll (horizontalCoordinates_[0], xi)
        {
        	if (horizontalCoordinates_[0][xi] < x && x < horizontalCoordinates_[0][xi + 1])
        	{
        		nLower = xi;
        		break;
        	}
        }

        scalar eta0 = surfaceElevation(nLower, 0);
        scalar eta1 = surfaceElevation(nLower + 1, 0);

        scalar x0 = horizontalCoordinates_[0][nLower];
        scalar x1 = horizontalCoordinates_[0][nLower + 1];

        eta = eta0*(x - x0)/(x1 - x0) + eta1*(x1 - x)/(x1 - x0);
    }
    else if (
    	   horizontalCoordinates_[0][0] < x &&
    	   x < horizontalCoordinates_[0][Nx - 1] &&
    	   Ny > 3
       )
    {
    	notImplemented("Irregular interpolation 3D fields");
    }
    else
    {
    	forAll (amp_, index)
    	{
    		scalar arg = omega_[index]*t - (k_[index] & p) + phi_[index];
    		eta += amp_[index]*Foam::cos(arg);
    	}

    	eta *= rampFactor_;
    	eta += seaLevel_;
    }

    return eta;
}

scalar irregularFast::surfaceElevation(const label& nX, const label& nY) const
{
    // Values are pre-computed in irregularFast::step()
	return Eta_[nX][nY];
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void irregularFast::step()
{
	scalar beforeTime = rT_.elapsedCpuTime();

	rampFactor_ = factor(rT_.time().value());

	if (!mesh_.changing() || ignoreMeshMotion_)
	{
        // No need to re-evaluate the part, which depends on coordinates. Only
		// necessary to re-evaluate the part, which depends on time.
        calcWaveFrequencyTrigonometricFunctions();

        // Re-evalute the pre-computed surface elevation
        calcSurfaceElevation();

        // Re-evaluate the pre-computed velocity components
        calcVelocityComponents();
	}
	else // Static mesh
	{
		notImplemented
		(
	        "void irregularFast::step() for changing/moving meshes"
	    );

        calcAllFields();
	}

	Info << "External step: " << rT_.elapsedCpuTime() - beforeTime << " s."
	     << endl;
}


scalar irregularFast::eta
(
    const point& x,
    const scalar& time
) const
{
	scalar eta = 0;

	label nX, nY, nZ;

	findIndexing(x, nX, nY, nZ, false);

	if (nX == -1 || nY == -1)
	{
		eta = interpolateSurfaceElevation(x, time, nX, nY);
	}
	else
	{
		eta = surfaceElevation(nX, nY);
	}

    return eta;
}


//scalar irregularFast::ddxPd
//(
//    const point& x,
//    const scalar& time,
//    const vector& unitVector
//) const
//{
//    return 0.0;
//}


scalar irregularFast::pExcess
(
    const point& x,
    const scalar& time
) const
{
    return 0.0;
}


vector irregularFast::U
(
    const point& x,
    const scalar& time
) const
{
	// Find the labels in the search table
	label nX, nY, nZ;

	findIndexing(x, nX, nY, nZ, true);

    // Define the return variable
    vector U(vector::zero);

    if (nX == -1 || nY == -1 || nZ == -1)// || true)
    {
    	scalar Z = (x & eZ_) - seaLevel_;

    	// Loop over all wave components
    	forAll (amp_, index)
    	{
    		//      OLD IMPLEMENTATION
    		scalar arg0 = omega_[index]*time - (k_[index] & x) + phi_[index];
    		scalar arg1 = K_[index]*(Z + h_);

    		scalar Uhorz = velAmp_[index]*Foam::cosh(arg1)*Foam::cos(arg0);
    		scalar Uvert = - velAmp_[index]*Foam::sinh(arg1)*Foam::sin(arg0);

    		// Note "-" because of "g" working in the opposite direction
    		U += Uhorz*compDir_[index] + Uvert*eZ_;
    	}

    }
    else
    {
    	// Get the reference to the hyperbolic functions
    	const scalarField& Sinh = Sinh_[nZ];
    	const scalarField& Cosh = Cosh_[nZ];

    	const vectorField& Ucmp = Ucomponent_[nX][nY];
    	const scalarField& Wcmp = Wcomponent_[nX][nY];

    	// New implementation with fewer operations
        U = Foam::sum(Ucmp*Cosh) - eZ_*Foam::sum(Wcmp*Sinh);
    }

    // Multiply by the ramping factor
    U *= rampFactor_;

    return U;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
