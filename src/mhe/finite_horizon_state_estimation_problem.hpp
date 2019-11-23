#ifndef FINITE_HORIZON_STATE_ESTIMATION_PROBLEM_H
#define FINITE_HORIZON_STATE_ESTIMATION_PROBLEM_H

#include "state_estimation_problem.hpp"
#include "time_varying_smooth_horizon.hpp"
#include "linear_algebra.hpp"

class FiniteHorizonStateEstimationProblem final : public StateEstimationProblem {
public:
  // Constructs FiniteHorizonStateEstimationProblem with setting parameters and 
  // allocates vectors and matrices.
  // Arguments:
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  FiniteHorizonStateEstimationProblem(const double T_f, const double alpha, 
                                      const int N);
  // Constructs FiniteHorizonStateEstimationProblem with setting parameters and 
  // allocates vectors and matrices.
  // Arguments:
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  //  initial_time: Initial time for the length of the horizon.
  FiniteHorizonStateEstimationProblem(const double T_f, const double alpha, 
                                      const int N, const double initial_time);

  // Free vectors and matrices.
  ~FiniteHorizonStateEstimationProblem();

  // Computes the optimaliy residual under time, control_input_mat, and 
  // output_mat, and solution_vec that represents the disturbance sequence. 
  // The result is set in optimality_residual.
  void computeOptimalityResidual(const double time, 
                                 double const* const* control_input_mat,
                                 double const* const* output_mat,
                                 const double* solution_vec,
                                 double* optimality_residual);

  // Reset the length of the horizon by resetting parameters related to the 
  // horizon.
  void resetHorizonLength(const double T_f, const double alpha, 
                          const double initial_time);

  // Reset the length of the horizon by resetting parameters related to the 
  // horizon.
  void resetHorizonLength(const double initial_time);

  // Returns the dimension of the solution, which is equivalent to 
  // N*(dim_control_input+dim_constraints).
  int dim_solution() const override;

  // Returns the grid number of the horizon.
  int N() const;

private:
  TimeVaryingSmoothHorizon horizon_;
  int dim_constraints_and_disturbance_, dim_solution_, N_;
  double *dx_vec_, **state_mat_, **lambda_mat_;
};

#endif // FINITE_HORIZON_STATE_ESTIMATION_PROBLEM_H