// Get an autoPtr to a porosity model
autoPtr<wavesPorosityModel> pm = Foam::wavesPorosityModel::New(mesh);

// Inside a scope to get the field deleted again.
{
	Info << "HERE" << endl;
#   include "calcPorosity.H"

	porosity.write();
}




