#ifndef ZERO_HORIZON_DISTURBANCE_ESTIMATION_PROBLEM_H
#define ZERO_HORIZON_DISTURBANCE_ESTIMATION_PROBLEM_H

#include "linear_algebra.hpp"
#include "state_estimation_problem.hpp"

class ZeroHorizonDisturbanceEstimationProblem final 
    : public StateEstimationProblem {
public:
  // Allocate a vector.
  ZeroHorizonDisturbanceEstimationProblem();

  // Free a vector.
  ~ZeroHorizonDisturbanceEstimationProblem();

  // Computes the optimaliy residual under time, state_vec, and solution_vec 
  // that represents the control input and Lgrange multiplier with respect to
  // equality constraints. The result is set in optimality_residual.
  void computeOptimalityResidual(const double time, 
                                 const double* measured_control_input_vec, 
                                 const double* measured_output_vec, 
                                 const double* estimated_state_vec, 
                                 const double* solution_vec,
                                 double* optimality_residual);

  // Return the dimension of the solution, 
  // i.e., dim_constraints+dim_disturbance.
  int dim_solution() const override; 

private:
  int dim_solution_;
  double *lambda_vec_;
};

#endif // ZERO_HORIZON_DISTURBANCE_ESTIMATION_PROBLEM_H