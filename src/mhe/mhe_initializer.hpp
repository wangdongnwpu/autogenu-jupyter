#ifndef MHE_INITIALIZER_H
#define MHE_INITIALIZER_H

#include "newton_gmres_for_disturbance_estimation.hpp"
#include "zero_horizon_disturbance_estimation_problem.hpp"
#include "matrixfree_gmres.hpp"
#include "linear_algebra.hpp"

class MHEInitializer {
public:
  // Sets parameters and allocates vectors. 
  // Arguments:
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  kmax: A parameter for the GMRES method. This parameter represents the
  //     dimension of the Krylov subspace and maximum iteration number of the
  //     GMRES method.
  // Argments:
  //   newton_residual_tolerance: A convergence criteria for the Newton iteration. 
  //     Newton iteration terminates when the error is less than this value.
  //   max_newton_iteration: Maximum number of the Newton iteration. Newton 
  //     iteration for the initialization terminates when the number of the 
  //     iteration is equal to this value.
  MHEInitializer(const double finite_difference_increment, const int kmax, 
                 const double newton_residual_tolerance, 
                 const int max_newton_iteration);

  // Sets parameters and allocates vectors. 
  // Arguments:
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  kmax: A parameter for the GMRES method. This parameter represents the
  //     dimension of the Krylov subspace and maximum iteration number of the
  //     GMRES method.
  MHEInitializer(const double finite_difference_increment, const int kmax);

  // Free vectors.
  ~MHEInitializer();

  // Sets parameters for Newton iteration. 
  //   newton_residual_tolerance: A convergence criteria for the Newton iteration. 
  //     Newton iteration terminates when the error is less than this value.
  //   max_newton_iteration: Maximum number of the Newton iteration. Newton 
  //     iteration for the initialization terminates when the number of the 
  //     iteration is equal to this value.
  void setCriterionsOfNewtonTermination(const double newton_residual_tolerance, 
                                        const int max_newton_iteration);

  // Sets the initial guess soluion, which is composed by the control input
  // vector and the Lgrange multiplier with respect the equality constraints.
  void setInitialGuessSolution(const double* initial_guess_solution_vec);

  // Solves the optimal control problem with horizon whose length is zero
  // under initial_time and initial_state_vec.
  void computeInitialSolution(const double initial_time, 
                              const double* initial_measured_control_input, 
                              const double* initial_measured_output, 
                              const double* initial_estimated_state_vec, 
                              double* initial_solution_vec);

  // Returns the dimenstion of the solution, which is equivalent to the 
  // dim_control_input and the dim_constraints.
  int dim_solution() const;

private:
  NewtonGMRESForDisturbanceEstimation newton_;
  MatrixFreeGMRES<NewtonGMRESForDisturbanceEstimation, const double, 
                  const double*, const double*, const double*, 
                  const double*> mfgmres_;
  const int dim_constraints_, dim_disturbance_, dim_solution_;
  int max_newton_iteration_;
  double newton_residual_tolerance_;
  double *initial_guess_solution_vec_, *solution_update_vec_;
};

#endif // MHE_INITIALIZER_H