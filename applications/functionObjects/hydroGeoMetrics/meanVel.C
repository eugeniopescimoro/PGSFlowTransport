template<class Type>
void Foam::functionObjects::hydroGeoMetrics::meanVel
(
    const GeometricField<Type, fvPatchField, volMesh>& U
)

{
  typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

	Info<< "Computing meanVel on "
		<< U.name()
		<< endl;


	const Field<Type>& u = U.primitiveField();

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
              "meanVel_"
              +
              U.name()
              );

  for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
  {
    scalar meanu(gSum(u.component(cmpt)*mesh_.V())/vol);
    scalar varu
        (
        gSum( pow(u.component(cmpt),2)*mesh_.V() ) / vol
        -
        pow(meanu,2)
        );

    fout << meanu << "\t" << varu << endl;
  }

  Info << "meanVel written to " << filePath << endl;
}
