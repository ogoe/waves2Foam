// Update the porosity field - accounting for moving meshes
pm->updatePorosity();

// Get a reference to the porosity field
const volScalarField& porosity = pm->porosity();

// Create a field with the porosity on the faces
const surfaceScalarField porosityFace(fvc::interpolate(porosity));

