#ifndef STATE_ESTIMATE_CONTINUATION_H
#define STATE_ESTIMATE_CONTINUATION_H

#include <cmath>
#include "linear_algebra.hpp"
#include "finite_horizon_state_estimation_problem.hpp"

class StateEstimateContinuation {
public:
  // Constructs SingleShootingContinuation with setting parameters and allocates 
  // vectors and matrices.
  // Arguments:
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  //  initial_time: Initial time for the length of the horizon.
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  zeta: A parameter for stabilization of the C/GMRES method. It may work
  //    well to set this parameters as the reciprocal of the sampling period.
  StateEstimateContinuation(const double T_f, const double alpha, const int N,
                            const double finite_difference_increment,
                            const double zeta);

  // Constructs SingleShootingContinuation with setting parameters and allocates 
  // vectors and matrices.
  // Arguments:
  //  T_f, alpha: Parameters for the length of the horizon. The length horizon
  //    at time t is given by T_f * (1-exp(-alpha*t)).
  //  N: The number of the discretization of the horizon.
  //  initial_time: Initial time for the length of the horizon.
  //  finite_difference_increment: Step length of the finite difference 
  //     approximation of the OCP for the initialization.
  //  zeta: A parameter for stabilization of the C/GMRES method. It may work
  //    well to set this parameters as the reciprocal of the sampling period.
  StateEstimateContinuation(const double T_f, const double alpha, const int N,
                            const double initial_time, 
                            const double finite_difference_increment,
                            const double zeta);

  ~StateEstimateContinuation();

  // Integrates the solution for given optimal update vector of the solution 
  // and the integration length.
  void integrateSolution(double* solution_vec, 
                         const double* solution_update_vec, 
                         const double integration_length);

  // Computes and returns the squared norm of the errors in optimality under 
  // the state_vec and the current solution.
  double computeErrorNorm(const double time, 
                          double const* const* control_input_mat,
                          double const* const* outpu_mat,
                          const double* solution_vec);

  // Reset the length of the horizon by resetting parameters related to the 
  // horizon.
  void resetHorizonLength(const double T_f, const double alpha, 
                          const double initial_time);

  // Reset the length of the horizon by resetting parameters related to the 
  // horizon.
  void resetHorizonLength(const double initial_time);

  // Computes a vector correspongin to b in Ax=b. This function is called in
  // MatrixfreeGMRES.
  void bFunc(const double time, double const* const* control_input_mat,
             double const* const* output_mat, 
             double const* const* incrementd_control_input_mat,
             double const* const* incremented_output_mat, 
             const double* current_solution_vec, 
             const double* current_solution_update_vec, 
             double* b_vec);

  // Computes a vector correspongin to Ax in Ax=b. This function is called in
  // MatrixfreeGMRES.
  void AxFunc(const double time, double const* const* control_input_mat,
              double const* const* output_mat, 
              double const* const* incrementd_control_input_mat,
              double const* const* incremented_output_mat, 
              const double* current_solution_vec, const double* direction_vec,
              double* ax_vec);

  // Returns the dimension of the state.
  int dim_state() const;

  // Returns the dimension of the control input.
  int dim_control_input() const;

  // Returns the dimension of the equality constraints.
  int dim_constraints() const;

  // Returns the dimension of the disturbance.
  int dim_disturbance() const;

  // Returns the dimension of the equality constraints.
  int dim_output() const;

  // Returns the dimension of the solution, which is equivalent to 
  // N*(dim_control_input+dim_constraints).
  int dim_solution() const;

  // Returns the grid number of the horizon.
  int N() const;

private:
  FiniteHorizonStateEstimationProblem estimation_problem_;
  const int dim_state_, dim_control_input_, dim_constraints_, dim_disturbance_, 
      dim_output_, dim_solution_;
  double finite_difference_increment_, zeta_, incremented_time_; 
  double *incremented_solution_vec_, *optimality_residual_, 
      *optimality_residual_1_, *optimality_residual_2_;
};

#endif // STATE_ESTIMATE_CONTINUATION_H