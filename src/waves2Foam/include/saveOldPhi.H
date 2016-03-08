    // Save the face flux at the boundaries; required for dedicated boundary conditions
    PtrList<IOField<scalar> > phiOldBoundary(mesh.boundary().size());

	forAll (phiOldBoundary, patchi)
	{
		phiOldBoundary.set
		(
		    patchi,
		    new IOField<scalar>
		    (
		         IOobject
		         (
		             "phiBoundary_" + mesh.boundary()[patchi].name(),
		             runTime.timeName(),
		             mesh,
		             IOobject::NO_READ,
		             IOobject::NO_WRITE
		         ),
		         phi.boundaryField()[patchi]
		    )
		);
	}


