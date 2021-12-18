template<class Type>
bool Foam::functionObjects::hydroGeoMetrics::checkField
(
    autoPtr<GeometricField<Type, fvPatchField, volMesh>>& fPtr
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if (mesh_.foundObject<VolFieldType>(fieldName_))
    {
        DebugInfo
            << "Field " << fieldName_ << " already in database"
            << endl;

        return true;
    }
    else
    {
        IOobject fieldHeader
        (
            fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if
        (
            fieldHeader.typeHeaderOk<VolFieldType>(false)
         && fieldHeader.headerClassName() == VolFieldType::typeName
        )
        {
			Info<< "Reading field " << fieldName_
				<< endl;

			fPtr.reset
		    (
		        new GeometricField<Type, fvPatchField, volMesh>
				(
					fieldHeader,
					mesh_
				)
		    );

            return true;
        }
    }

	return false;
}

template<class Type>
bool Foam::functionObjects::hydroGeoMetrics::computeMetrics
(
    autoPtr<GeometricField<Type, fvPatchField, volMesh>>& fPtr
)
{
    if(!checkField(fPtr))
    {
        return false;
    }

    forAll(operations_,opi)
    {
      if  (
          operations_[opi]=="twoPointCorrelation"
          ||
          operations_[opi]=="autocorrelation"
          )
      {
        twoPointCorrelation
        (
            mesh_.lookupObject<GeometricField<Type, fvPatchField, volMesh>>
            (
                fieldName_
            )
        );
      }
      // else if (operations_[opi]=="entropy")
      // {
      //   entropy
      //   (
      //       mesh_.lookupObject<GeometricField<Type, fvPatchField, volMesh>>
      //       (
      //           fieldName_
      //       )
      //   );
      // }
      else if (
              operations_[opi]=="meanAndVariance"
              ||
              operations_[opi]=="mean"
              ||
              operations_[opi]=="variance"
              )
      {
        meanAndVariance
        (
            mesh_.lookupObject<GeometricField<Type, fvPatchField, volMesh>>
            (
                fieldName_
            )
        );
      }
      else
      {
        Info << "Operation " << operations_[opi] << " not supported!"
             << endl << "Available operations: mean, autocorrelation";
      }
    }


    return true;
}
