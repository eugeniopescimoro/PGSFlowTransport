/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         |
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    macroScalarTransport

Description
    Solves the steady or transient transport equation for a passive scalar using
    an adaptive time step while setting a threshold as a lower limit for the
    passive scalar to be computed. It also condiders variable mechanical
    dispersion according to the model: MD = alpha*V

Developers
    Eugenio Pescimoro
    Matteo Icardi

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    //#include "createTimeControls.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nCalculating scalar transport\n" << endl;

    #include "myCourantN.H"
    scalar maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    scalarField sumQ
    (
        fvc::surfaceSum(mag(q))().primitiveField()
    );

    runTime.setDeltaT (maxCo/(0.5*gMax(sumQ/mesh.V().field())));

    label patchId = mesh.boundaryMesh().findPatchID(outletPatch);
    const scalarField& TOut(c.boundaryField()[patchId]);
    const vectorField& Sf( mesh.Sf().boundaryField()[patchId]);
    const scalarField& qout(q.boundaryField()[patchId]);
//    const vectorField& Uout(U.boundaryField()[patchId]);

//    scalarField fluxOut((Uout&Sf));
//    scalarField fluxOut(qout);
    //Info << "Flux out = " << gSum(fluxOut) << nl << endl;

    vector meanVel( (fvc::domainIntegrate(U)).value() / gSum(mesh.V()) );
    Info << "Mean vel = " << meanVel << nl << endl;

    //scalar massFlag(0.0);
    scalar totalMassOld(0.0);
    // surfaceScalarField D(phi*alpha*mag(q)/mesh.magSf());
    //- Update cell-based dispersion
    surfaceSymmTensorField D
      (
       "D",
        phi
        *
        (
          (alphaT * mag(q) ) * symmTensor::I
          +
          (alpha-alphaT) * (mag(q) * sqr(mesh.Sf()) / sqr(mesh.magSf())  )
        ) / mesh.magSf()
      );
    //D.write();


    while (simple.loop(runTime))
    {
        Info << "Adaptive time = " << runTime.timeName() << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(phi, c)
              + fvm::div(q, c)
              - fvm::laplacian(phi*Dm, c)
              - fvm::laplacian(D, c)
             ==
                fvOptions(c)
            );

            TEqn.relax();
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(c);
        }

        scalar totalMass( (fvc::domainIntegrate(c)).value() );
        scalar cFluxOut(gSum(TOut*qout)/(gSum(qout)+SMALL) );

        Info << "Total mass = " << totalMass << endl;
        Info << "Flux out = " << cFluxOut << nl <<  endl;
        //Info << "maximum " << gMax(c) << endl;
        //Info << "minimum "<< gMin(c) << nl << endl;

        if ( cFluxOut  > THRS )
        {
            c.write();
            FatalErrorIn("macroScalarTransport")
                << "Concentration under the threshold" << endl
                << exit(FatalError);
        }

        //scalar massFlagNew(fmod(totalMass,mTHRS));
        //if ( massFlagNew < massFlag )
        if ( (totalMass-totalMassOld) > mTHRS )
        {
          c.write();
          totalMassOld = totalMass;
        }
        //massFlag = massFlagNew;

        runTime.write();
    }

    Info << "End\n" << nl << endl;

    return 0;
}


// ************************************************************************* //
