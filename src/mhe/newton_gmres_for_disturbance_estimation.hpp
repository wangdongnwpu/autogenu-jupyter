#ifndef NEWTON_GMRES_FOR_DISTURBANCE_ESTIMATION_H
#define NEWTON_GMRES_FOR_DISTURBANCE_ESTIMATION_H

#include <cmath>
#include "linear_algebra.hpp"
#include "zero_horizon_disturbance_estimation_problem.hpp"

class NewtonGMRESForDisturbanceEstimation {
public:
  // Constructs NewtonGMRES with setting parameters and allocates vectors and 
  // matrices.
  // Arguments:
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  NewtonGMRESForDisturbanceEstimation(const double finite_difference_increment);

  // Free vectors and matrices.
  ~NewtonGMRESForDisturbanceEstimation();

  // Returns the squared norm of the optimality residual under time, state_vec, 
  // and the current solution.
  double errorNorm(const double time, 
                   const double* measured_control_input_vec,
                   const double* measured_output_vec,
                   const double* estimated_state_vec, 
                   const double* solution_vec);

  // Computes a vector correspongin to b in Ax=b. This function is called in
  // MatrixfreeGMRES.
  void bFunc(const double time, 
             const double* measured_control_input_vec,
             const double* measured_output_vec,
             const double* estimated_state_vec, 
             const double* current_solution_vec, 
             const double* current_solution_update_vec, double* b_vec);

  // Computes a vector correspongin to Ax in Ax=b. This function is called in
  // MatrixfreeGMRES.
  void AxFunc(const double time, 
              const double* measured_control_input_vec,
              const double* measured_output_vec,
              const double* estimated_state_vec, 
              const double* current_solution_vec, const double* direction_vec,
              double* ax_vec);

  // Returns dimension of the state.
  int dim_state() const;

  // Returns dimension of the control input.
  int dim_control_input() const;

  // Returns dimension of the constraints.
  int dim_constraints() const;

  // Returns dimension of the output.
  int dim_output() const;

  // Returns dimension of the disturbance.
  int dim_disturbance() const;

  // Returns dimension of the solution of the optimal control problem.
  int dim_solution() const;

private:
  ZeroHorizonDisturbanceEstimationProblem problem_;
  const int dim_solution_;
  double finite_difference_increment_;
  double *incremented_solution_vec_, *optimality_residual_, 
      *optimality_residual_1_;
};

#endif // NEWTON_GMRES_FOR_DISTURBANCE_ESTIMATION_H