    // Relative momentum predictor
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevSigma(U)
      + REL->Su()
     ==
        fvModels.source(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    solve(UEqn == -fvc::grad(p));

    fvConstraints.constrain(U);
