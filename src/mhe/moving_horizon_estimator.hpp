#ifndef MOVING_HORIZON_ESTIMATOR_H
#define MOVING_HORIZON_ESTIMATOR_H

#include "matrixfree_gmres.hpp"
#include "state_estimate_continuation.hpp"
#include "time_series_data.hpp"
#include "linear_algebra.hpp"
#include "time_varying_smooth_horizon.hpp"
#include "mhe_initializer.hpp"

class MovingHorizonEstimator {
public:
  MovingHorizonEstimator(const double T_f, const double alpha, const int N,
                         const double finite_difference_increment,
                         const double zeta, const int kmax, 
                         const double min_sampling_period);
  // Free vectors and matrices.
  ~MovingHorizonEstimator();

  void updateEstimation(const double time, double* estimated_state_vec,
                        const double sampling_period);

  void getStateEstimation(double* estimated_state_vec) const;

  void setParametersForInitialization(const double* initial_guess_solution, 
                                      const double newton_residual_tolerance,
                                      const int max_newton_iteration);

  // Initializes the solution of the C/GMRES method by solving the optimal
  // control problem with the horizon whose length is zero. soltuion_vec_ and 
  // errors_in_optimality_ is fullfilled with the solution of this OCP. The 
  // control input to be applied to the actual system is assigned in 
  // optimal_control_input_vec.
  void initializeSolution(const double initial_time,  
                          const double* initial_control_input_vec,
                          const double* initial_output_vec,
                          const double* initial_estimated_state_vec);

  // Returns the squared norm of the optimality residual under time, state_vec, 
  // and the current solution.
  double getErrorNorm(const double time);

  void appendMeasuredControlInput(const double measured_time, 
                                  const double* measured_control_input);

  void appendMeasuredOutput(const double measured_time, 
                            const double* measured_output);

  // Prohibits copy due to memory allocation.
  MovingHorizonEstimator(const MovingHorizonEstimator&) = delete;
  MovingHorizonEstimator& operator=(const MovingHorizonEstimator&) = delete;

private:
  StateEstimateContinuation continuation_problem_;
  MatrixFreeGMRES<StateEstimateContinuation, const double, 
                  double const* const*, double const* const*,
                  double const* const*, double const* const*> mfgmres_;
  TimeSeriesData control_input_data_, output_data_;
  TimeVaryingSmoothHorizon horizon_;
  MHEInitializer solution_initializer_;
  const int dim_state_, dim_control_input_, dim_constraints_, dim_disturbance_, 
      dim_output_, N_;
  const double finite_difference_increment_;
  double *solution_vec_, *solution_update_vec_, *initial_solution_vec_;
  double **control_input_mat_, **output_mat_, **incremented_control_input_mat_, 
      **incremented_output_mat_;
};

#endif // MOVING_HORIZON_ESTIMATOR_H