# sediBuoyantBoussinesqSimpleFoam
Steady-state solver for sedimentation, buoyant, turbulent flow of incompressible fluids.

Developed based on the native OpenFOAM solver of buoyantBoussinesqSimpleFoam and equilibrium-Euler appraoch

# Notes

For complex problem, this solver is struggling with convergence. (Centrally, the numerical Reynolds number is large enough that inertial kicks in.) In this case, it is recommanded to run the transient solver to steady state (AKA pseudo time stepping, in which the inclusion of the time deriviatie will introduce the diagnal contribution to the coefficient matrix (LHS) and generate a better conditioned system).
