// Obtain the porosity field
tmp<volScalarField> tporosity(pm->porosity());
const volScalarField& porosity = tporosity();

// Create a field with the porosity on the faces
const surfaceScalarField porosityFace(fvc::interpolate(porosity));
