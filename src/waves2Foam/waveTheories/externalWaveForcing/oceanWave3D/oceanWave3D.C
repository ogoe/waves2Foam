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

#include "oceanWave3D.H"
#include "addToRunTimeSelectionTable.H"

#include "fvCFD.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

extern "C" double __globalvariables_MOD_time;
extern "C" double __globalvariables_MOD_dt;
extern "C" int __globalvariables_MOD_tstep;
extern "C" int __globalvariables_MOD_ic;
extern "C" int __globalvariables_MOD_nsteps;
extern "C" void oceanwave3dt0setup_();
extern "C" void oceanwave3dtakeatimestep_();
extern "C" void closeiofiles_();
extern "C" void closevariables_();
extern "C" void interpolationinitialize_();
extern "C" void calculatekinematics_();
extern "C" void openfoaminterface_eta_(double(*)[3] , double *);
extern "C" void openfoaminterface_u_(double(*)[3] , double *, double *, double *);
extern "C" void writeoceanwave3d_(int *);


namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oceanWave3D, 0);

addToRunTimeSelectionTable
(
    externalWaveForcing,
    oceanWave3D,
    externalWaveForcing
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


oceanWave3D::oceanWave3D
(
    IOobject io,
    Time& rT,
    const fvMesh& mesh
)
:
    externalWaveForcing(io, rT, mesh),

    waveProps_(io.db().lookupObject<IOdictionary>("waveProperties")),
    coeffDict_(waveProps_.subDict("externalForcingCoeffs")),

    seaLevel_(readScalar(waveProps_.lookup("seaLevel"))),

    N_(readLabel(coeffDict_.lookup("nIntervals"))),

    startTimes_("startTimes", coeffDict_, N_),

    endTimes_("endTimes", coeffDict_, N_),

    ramp_(Switch(coeffDict_.lookup("rampInterval"))),

    Tsoft_(0),

    translateOFMesh_
    (
    	coeffDict_.lookup("translateOpenFoamMesh")
//    	coeffDict_.lookupOrDefault<vector>
//        (
//        	"translateOpenFoamMesh", vector::zero
//        )
    ),

    OFtoOCW_(tensor::zero),

    OCWtoOF_(tensor::zero)
{
	if (N_ > 1)
	{
		FatalErrorIn("oceanWave3D::oceanWave3D(IOobject io, Time& rT, const fvMesh& mesh)")
			<< "Unwanted behaviour has been observed in the re-starting, when\n"
			<< "multiple intervals are used. Consequently, make one simulation"
			<< "per interval instead.\n"
			<< endl << exit(FatalError);
	}

	// Initialise to the first set of intervals
    startTime_ = -GREAT;
    endTime_ = GREAT;

    forAll (startTimes_, timei)
    {
    	if (rT_.time().value() < endTimes_[timei])
    	{
    		startTime_ = startTimes_[timei];
    		endTime_ = endTimes_[timei];
    		N_ = timei;
    		break;
    	}
    }

    if (rT_.endTime().value() < startTime_)
    {
		FatalErrorIn("oceanWave3D::oceanWave3D(IOobject io, Time& rT, const fvMesh& mesh)")
			<< "The simulation will terminate before the selected interval starts.\n"
			<< "Correct the endTime in system/controlDict.\n"
			<< endl << exit(FatalError);
    }

    if (rT_.endTime().value() < endTime_)
    {
		FatalErrorIn("oceanWave3D::oceanWave3D(IOobject io, Time& rT, const fvMesh& mesh)")
			<< "The simulation will terminate before the selected interval ends.\n"
			<< "Correct the endTime in system/controlDict.\n"
			<< endl << exit(FatalError);
    }

	if (ramp_)
	{
		Tsoft_ = readScalar(coeffDict_.lookup("Tsoft"));
	}

	// Make the mapping tensors
	mappingTensors();

	// Start OceanWave3D
	oceanwave3dt0setup_();

	// Initialise the interpolation routine
	interpolationinitialize_();

	// Get the uniform time step specified in OceanWave3D input file
	maxDT_ = __globalvariables_MOD_dt;

	label ocwDuration = __globalvariables_MOD_nsteps;

	if (maxDT_*ocwDuration < rT_.endTime().value())
	{
		FatalErrorIn("oceanWave3D::oceanWave3D(IOobject io, Time& rT, const fvMesh& mesh)")
			<< "The duration of the OpenFoam simulation ("
			<< rT_.endTime().value() << " s) exceeds the duration of the \n"
			<< "OceanWave3D simulation (" << maxDT_*ocwDuration << " s).\n"
			<< exit(FatalError) << endl << endl;
	}

	// Update the OceanWave3D to current time (restart)
	alignTimes();
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar oceanWave3D::factor(const scalar& time) const
{
    scalar factor(1.0);

    scalar rampTime = time;

    if (ramp_)
    {
        rampTime -= startTime_;
    }

    if (0.0 < Tsoft_)
    {
     	factor = Foam::sin(2*M_PI/(4.0*Tsoft_)*Foam::min(Tsoft_, rampTime));
    }

    return factor;
}


dimensionedScalar oceanWave3D::OCWTimeStep() const
{
    return dimensionedScalar
    	(
    		"null",
    		dimTime,
    		Foam::min(1.2*rT_.deltaT().value(), maxDT_)
        );
}


void oceanWave3D::mappingTensors()
{
	// Get a unit vector along the direction of gravity.
	vector g = uniformDimensionedVectorField
        (
	        mesh_.thisDb().lookupObject<uniformDimensionedVectorField>("g")
	    ).value();

	g = cmptMultiply(g, g);
	g /= Foam::mag(g);

	// Discard a gravity vector pointing partyly in the X-direction
	if ((g & vector(1, 0, 0)) < SMALL)
	{
		// M_11 = 1;
        OFtoOCW_.xx() = 1;

        // Set the remaining coefficients of M depending on the direction of
        // gravity.
        if (SMALL < (g & vector(0, 1, 0)))
        {
        	OFtoOCW_.zy() = 1;
        	OFtoOCW_.yz() = 1;
        }
        else
        {
        	OFtoOCW_.yy() = 1;
        	OFtoOCW_.zz() = 1;
        }
	}
	else
	{
		FatalErrorIn("void oceanWave3D::mappingTensors()")
			<< "The gravity points (partly) in the x-direction. This is not \n"
			<< "supported together with the oceanWave3D external forcing class."
			<< endl << exit(FatalError);
	}

	// Create the inverse map
	OCWtoOF_ = Foam::inv(OFtoOCW_);

	// Check that the translation vector is horizontal
    if (SMALL < Foam::mag(translateOFMesh_ & g))
    {
    	FatalErrorIn("void oceanWave3D::mappingTensors()")
    	    << "The translation vector of the computational mesh for OpenFoam\n"
    		<< "is not horizontal. "
    		<< endl << exit(FatalError);
    }
}


void oceanWave3D::updateIntervals()
{
	// Update the start-end times, if there are any intervals left. Otherwise,
	// set the start time and end times to a very, very big number.
	if (endTime_ < rT_.time().value())
	{
        N_++;

		if (N_ < startTimes_.size())
		{
		    startTime_ = startTimes_[N_];
		    endTime_ = endTimes_[N_];
		}
		else
		{
			startTime_ = GREAT;
			endTime_ = GREAT;
		}
	}
}


void oceanWave3D::alignTimes()
{
	Info << "OF-time: " << rT_.time().value() << endl;
	Info << "OCW-time: " << __globalvariables_MOD_time << endl;
    Info << "Interval: [" << startTime_ << ", " << endTime_ << "]" << endl;

    // If the starting times are identical, simply do nothing
    if (Foam::mag(rT_.time().value() - __globalvariables_MOD_time) < 1.0e-9)
    {
    	Info << "Identical starting times for OceanWave3D and OpenFoam" << endl;

    	return;
    }

    // If the OpenFoam starting time exceeds that of OceanWave3D, take a number
    // of OceanWave3D steps.
    if (__globalvariables_MOD_time < rT_.time().value())
    {
    	// Evaluate the mismatch in start time
    	scalar dT = rT_.time().value() - __globalvariables_MOD_time;

    	// Calculate the number of OceanWave3D time steps and the time step
        label NtimeSteps = std::ceil(dT/maxDT_);
        scalar alignDt = dT/NtimeSteps;

        // Loop OceanWave3D
        for (long n=0; n < NtimeSteps; n++)
        {
            takeTimeStep(alignDt, false);
        }

        // Calculate the kinematics, which could be needed for setWaveField at
        // at a mapping time larger than 0
        calculateKinematics();

        // Chech again for the alignment of the times.
        alignTimes();
    }
    else
    {
    	FatalErrorIn("void oceanWave3D::alignTimes()")
    			<< "The starting time in OceanWave3D exceeds the starting time\n"
    			<< "in OpenFoam: \n\n"
    			<< "Start time in OceanWave3D: " << __globalvariables_MOD_time << " s.\n"
    			<< "Start time in OpenFoam:    " << rT_.time().value() << " s."
    			<< endl << exit(FatalError);
    }
}


void oceanWave3D::timeStepOceanWave3D()
{
    // If the current time is less than the next startTime_ for OpenFoam,
	// perform a lot of time steps with OceanWave3D, map the solution and
	// go back to OpenFoam
	if (rT_.time().value() < startTime_)
	{
		// Make sure to take the one time step that the time in OpenFoam has
		// already been increased with:
		takeTimeStep(rT_.deltaT().value(), false);

	    while (SMALL < startTime_ - rT_.time().value())
	    {
	    	// Set the size of the time step
	    	dimensionedScalar dt = OCWTimeStep();

	    	// Go forward in time
		    rT_.setDeltaT(dt);
            rT_++;

            Info << this->type() << ": Time = " << rT_.time().value() << " s.\n"
                 << this->type() << ": deltaT = " << rT_.deltaT().value() << " s."
                 << endl << endl;

            // Take a time step with OpenFoam 3D
            takeTimeStep(false);

            // Make sure that the model ends, if the endTime of the OpenFoam
            // simulation is exceeded. Instead of going back to OpenFoam and do
            // one more time step, a simple 'FatalErrorIn' approach has been
            // used. It is not pretty, but it works.
            if (rT_.endTime().value() < rT_.time().value())
            {
            	FatalErrorIn("void oceanWave3D::step()")
            		<< "This is not a FatalError, but it is a simple way to stop"
            		<< " the execution. " << exit(FatalError);
            }
	    }

	    // Calculate the kinematics from the last time step
	    calculateKinematics();

	    // Map the solution to the OF-domain
        mapSolution();

        updatePhi();

        updateTimeAndTimeStep();
	}
}


void oceanWave3D::mapSolution()
{
	// Obtain references to the various field. Note the const_case!
    volScalarField& rho = const_cast<volScalarField&>
        (
        	mesh_.thisDb().lookupObject<volScalarField>("rho")
        );

    volScalarField& alpha = const_cast<volScalarField&>
        (
           	mesh_.thisDb().lookupObject<volScalarField>
            (
                Foam::waves2Foam::aName()
            )
        );

    volScalarField& pd = const_cast<volScalarField&>
        (
           	mesh_.thisDb().lookupObject<volScalarField>
            (
                Foam::waves2Foam::pName()
            )
        );

    volVectorField& U = const_cast<volVectorField&>
        (
           	mesh_.thisDb().lookupObject<volVectorField>
            (
                "U"
            )
        );

    // Map the solution using the setWaveField functionality
    Info << "Mapping the OceanWave3D solution to the OpenFoam domain" << endl;
    Info << "Applied ramping factor: " << factor(rT_.time().value()) << endl;
    Info << endl;

    word name(coeffDict_.lookup("mappingZone"));
    setWaveField mapSolution(mesh_, name, U, alpha, pd);

    mapSolution.correct();

    U.correctBoundaryConditions();
    alpha.correctBoundaryConditions();
    pd.correctBoundaryConditions();

    // Obtain the densities
    const dictionary& dict
        = mesh_.thisDb().lookupObject<IOdictionary>("transportProperties");

    dimensionedScalar rho0(dict.subDict(waves2Foam::airPhase()).lookup("rho"));
    dimensionedScalar rho1(dict.subDict(waves2Foam::waterPhase()).lookup("rho"));

    // Update the density based on the mapped alpha-field
    rho = (alpha*rho1 + (1 - alpha)*rho0);
    rho.correctBoundaryConditions();
}


void oceanWave3D::updatePhi()
{
	// Obtain references to the phi field. Note the const_case!
	surfaceScalarField& phi = const_cast<surfaceScalarField&>
	    (
	      	mesh_.thisDb().lookupObject<surfaceScalarField>("phi")
	    );

    const volScalarField& pd =
    	mesh_.thisDb().lookupObject<volScalarField>
        (
            Foam::waves2Foam::pName()
        );

    const volScalarField& rho =
    	mesh_.thisDb().lookupObject<volScalarField>
        (
            "rho"
        );

    const volVectorField& U = mesh_.thisDb().lookupObject<volVectorField>("U");

    surfaceScalarField phiTemp
    (
        IOobject
        (
            "phiTemp",
            rT_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(U) & mesh_.Sf()
    );


    // Map the phiTemp -> phi
#if EXTBRANCH==1
    phi.internalField() = phiTemp.internalField();
#elif OFPLUSBRANCH==1
    #if OFVERSION<1706
        phi.internalField() = phiTemp.internalField();
    #else
        phi.ref() = phiTemp.internalField();
    #endif
#else
    #if OFVERSION<400
        phi.internalField() = phiTemp.internalField();
    #else
        phi.ref() = phiTemp.internalField();
    #endif
#endif
    

    forAll (phi.boundaryField(), patchi)
    {
#if OFPLUSBRANCH==1
    #if OFVERSION<1706
    	phi.boundaryField()[patchi] == phiTemp.boundaryField()[patchi];
    #else
    	phi.boundaryFieldRef()[patchi] == phiTemp.boundaryField()[patchi];
    #endif
#else
  	phi.boundaryField()[patchi] == phiTemp.boundaryField()[patchi];
#endif
    }

    // Perform the correction to phi utilising the correctPhi code (foam-extend-3.1)
    {
    	// Get the controls
    	dictionary pimple = mesh_.solutionDict().subDict("PIMPLE");

 	    int nNonOrthCorr =
   	        pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

 	    // Get information on the reference pressure
 	    label pdRefCell = 0;
 	    scalar pdRefValue = 0.0;

 	    if (pd.needReference())
 	    {
            pdRefCell = readLabel(pimple.lookup("pdRefCell"));
            pdRefValue = readScalar(pimple.lookup("pdRefValue"));
 	    }

 	    // Set up the correction of phi
        wordList pcorrTypes
        (
            pd.boundaryField().size(),
            zeroGradientFvPatchScalarField::typeName
        );

        for (label i=0; i<pd.boundaryField().size(); i++)
        {
            if (pd.boundaryField()[i].fixesValue())
            {
                pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        volScalarField pcorr
        (
            IOobject
            (
                "pcorr",
                rT_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("pcorr", pd.dimensions(), 0.0),
            pcorrTypes
        );

        dimensionedScalar rUAf("(1|A(U))", dimTime/rho.dimensions(), 1.0);

        adjustPhi(phi, U, pcorr);

        for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pcorrEqn
            (
                fvm::laplacian(rUAf, pcorr) == fvc::div(phi)
            );

            pcorrEqn.setReference(pdRefCell, pdRefValue);
            pcorrEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
            	phi -= pcorrEqn.flux();
            }
        }
    }

    // Write the mapped field
    U.write();
    pd.write();
    rho.write();
    mesh_.thisDb().lookupObject<volScalarField>(waves2Foam::aName()).write();
}


void oceanWave3D::updateTimeAndTimeStep()
{
	// Make constant references to change the names of local variables to fit
	// with the standard approaches
	const fvMesh& mesh = mesh_;
	const surfaceScalarField& phi
	    = mesh_.thisDb().lookupObject<surfaceScalarField>("phi");

	// Make name-change to rT_. Has to be a non-const, because setDeltaT.H
	// modifies the time step
	Time& runTime = rT_;

	// Calculate the Courant number
    #include "CourantNo.H"

	// Set the time step according to the new phi-field
#if EXTBRANCH==1
    #if OFVERSION<400
	#include "readTimeControls.H"
    #else
        #include "createTimeControls.H"   
    #endif
#elif OFPLUSBRANCH==1
    #include "createTimeControls.H"
#else
    #if OFVERSION<300
        #include "readTimeControls.H"
    #else
        #include "createTimeControls.H"
    #endif
#endif
    #include "setDeltaT.H"

	// Update time:
	runTime++;
}


void oceanWave3D::takeTimeStep(bool calcKinematics)
{
	takeTimeStep(rT_.deltaT().value(), calcKinematics);
}


void oceanWave3D::takeTimeStep(const scalar dt, bool calcKinematics)
{

	// Set the time step for OceanWave3D and increase its time step counter
    __globalvariables_MOD_dt = dt;
    __globalvariables_MOD_tstep = __globalvariables_MOD_tstep + 1;

    // Take one time step in OceanWave3D
    oceanwave3dtakeatimestep_();

    // If required, calculate the wave kinematics in OceanWave3D
    if (calcKinematics)
    {
    	calculateKinematics();
    }
}


void oceanWave3D::calculateKinematics()
{
	calculatekinematics_();
}


void oceanWave3D::writeExternal() const
{
	// Create the output directory, if it does not exist
	fileName outputDir = rT_.path();

	if (Pstream::parRun())
	{
		outputDir = rT_.path()/"../OCW3Dhotstart/";
	}
	else
	{
		outputDir = rT_.path()/"OCW3Dhotstart/";
	}

	if (!Foam::isDir(outputDir))
	{
		Foam::mkDir(outputDir);
	}

    // Write the OCW3D data, if the current time is an output time and the
	// running process is the master node
	if (Pstream::master() && rT_.outputTime())
    {
        // The index of the written field is based on the timeIndex in OpenFoam
		// Consequently, it is easy to match a hot-start file for OCW3D with
		// OpenFoam, because the time index is written in <time>/uniform/time
    	int writeCounter = rT_.timeIndex();

    	// Write the data
    	writeoceanwave3d_(&writeCounter);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void oceanWave3D::step()
{
	scalar beforeTime = rT_.elapsedCpuTime();

	// Update intervals, if the old one is finished
    updateIntervals();

    // Perform time steps in OceanWave3D until OpenFoam has to be used
    timeStepOceanWave3D();

    // Take a single time step in OceanWave3D and return to OpenFoam
    takeTimeStep(true);

    writeExternal();

	Info << "External step: " << rT_.elapsedCpuTime() - beforeTime << " s."
	     << endl;
}


void oceanWave3D::close()
{
	// If used for setWaveField at t = 0, there is nothing to close
	if (0 < __globalvariables_MOD_time)
	{
        closevariables_();
    }
}


scalar oceanWave3D::eta
(
    const point& x,
    const scalar& time
) const
{
	// Rotate the point x according to the predefined rotation matrix
	vector xx = OFtoOCW_ & (x + translateOFMesh_);

	// Create output
    double eta(0);

    // Create location in Fortran-format based on rotated coordinates
    double x0[3];
    x0[0] = xx[0];
    x0[1] = xx[1];
    x0[2] = xx[2] - seaLevel_;  // Notice displacement of coordinate system in OceanWave3D

    // Evaluate the surface elevation
    openfoaminterface_eta_(&x0,&eta);

    // Multiply by the ramping factor
    eta *= factor(time);
    eta += seaLevel_;

    // Return data
    return eta;
}


//scalar oceanWave3D::ddxPd
//(
//    const point& x,
//    const scalar& time,
//    const vector& unitVector
//) const
//{
//    return 0.0;
//}


scalar oceanWave3D::pExcess
(
    const point& x,
    const scalar& time
) const
{
    return 0.0;
}


vector oceanWave3D::U
(
    const point& x,
    const scalar& time
) const
{
	// Rotate the point x according to the predefined rotation matrix
	vector xx = OFtoOCW_ & (x + translateOFMesh_);

    // Map the coordinates to fortran format
    double x0[3];

    x0[0] = xx[0];
    x0[1] = xx[1];
    x0[2] = xx[2] - seaLevel_;

    // Make the return variables
    double utemp(0);
    double vtemp(0);
    double wtemp(0);

    // Evaluate the velocity in the given point
    openfoaminterface_u_(&x0,&utemp,&vtemp,&wtemp);

    // Map the solution of OF-format but still in rotated form
    vector U(vector::zero);
    U.x() = utemp;
    U.y() = vtemp;
    U.z() = wtemp;

    // Rotate the solution back to the defined coordinate system in OpenFoam
    U = OCWtoOF_ & U;

    // Multiply by the ramping factor
    U *= factor(time);

    // Return the solution
    return U;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
