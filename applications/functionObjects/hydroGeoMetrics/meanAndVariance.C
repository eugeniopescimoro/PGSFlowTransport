template<class Type>
void Foam::functionObjects::hydroGeoMetrics::meanAndVariance
(
    const GeometricField<Type, fvPatchField, volMesh>& K
)

{
  typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

	Info<< "Computing meanAndVariance on "
		<< K.name()
		<< endl;


	const Field<Type>& k = K.primitiveField();

	const scalar vol(gSum(mesh_.V()));

  fileName filePath =
              mesh_.time().path()
              /
              "postProcessing/fieldMetrics"
              /
              mesh_.time().timeName();

  if (Pstream::parRun())
  {
      filePath =
          mesh_.time().path()
          /
          "../postProcessing/fieldMetrics"
          /
          mesh_.time().timeName();
  }

  mkDir(filePath);

  OFstream fout
              (
              filePath
              /
              "meanAndVariance_"
              +
              K.name()
              );

  for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
  {
    scalar meank(gSum(k.component(cmpt)*mesh_.V())/vol);
    scalar vark
        (
        gSum( pow(k.component(cmpt),2)*mesh_.V() ) / vol
        -
        pow(meank,2)
        );

    fout << meank << "\t" << vark << endl;
  }

  Info << "meanAndVariance written to " << filePath << endl;
}
