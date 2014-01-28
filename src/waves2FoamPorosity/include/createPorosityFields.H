porosityZones pZones(mesh);

// Obtain the porosity field
tmp<volScalarField> tporosity(pZones.porosity());
volScalarField& porosity = tporosity();

// Writes the porosity field to starting time folder
porosity.write();

// Create a field with the porosity on the faces
surfaceScalarField porosityFace(fvc::interpolate(pZones.porosity()()));
