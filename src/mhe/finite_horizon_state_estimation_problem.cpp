#include "finite_horizon_state_estimation_problem.hpp"

FiniteHorizonStateEstimationProblem::FiniteHorizonStateEstimationProblem(
    const double T_f, const double alpha, const int N)
  : StateEstimationProblem(),
    horizon_(T_f, alpha),
    dim_constraints_and_disturbance_(
        model_.dimDisturbance()+model_.dimConstraints()),
    dim_solution_(N*(model_.dimDisturbance()+model_.dimConstraints())),
    N_(N),
    dx_vec_(linearalgebra::NewVector(model_.dimState())),
    state_mat_(linearalgebra::NewMatrix(N+1, model_.dimState())),
    lambda_mat_(linearalgebra::NewMatrix(N+1, model_.dimState())) {
}

FiniteHorizonStateEstimationProblem::FiniteHorizonStateEstimationProblem(
    const double T_f, const double alpha, const int N, 
    const double initial_time)
  : StateEstimationProblem(),
    horizon_(T_f, alpha, initial_time),
    dim_constraints_and_disturbance_(
        model_.dimDisturbance()+model_.dimConstraints()),
    dim_solution_(N*(model_.dimDisturbance()+model_.dimConstraints())),
    N_(N),
    dx_vec_(linearalgebra::NewVector(model_.dimState())),
    state_mat_(linearalgebra::NewMatrix(N+1, model_.dimState())),
    lambda_mat_(linearalgebra::NewMatrix(N+1, model_.dimState())) {
}

FiniteHorizonStateEstimationProblem::~FiniteHorizonStateEstimationProblem() {
  linearalgebra::DeleteVector(dx_vec_);
  linearalgebra::DeleteMatrix(state_mat_);
  linearalgebra::DeleteMatrix(lambda_mat_);
}

void FiniteHorizonStateEstimationProblem::computeOptimalityResidual(
    const double time, double const* const* control_input_mat,
    double const* const* output_mat, const double* solution_vec,
    double* optimality_residual) {
  double horizon_length = horizon_.getLength(time);
  double delta_tau = horizon_length / N_;
  // Compute the state trajectory over the horizon on the basis of the 
  // time, solution_vec and the state_vec.
  for (int i=0; i<dim_state_; ++i) {
    state_mat_[N_][i] = solution_vec[N_*dim_constraints_and_disturbance_+i];
  }
  double tau = time;
  for (int i=N_-1; i>=0; --i, tau-=delta_tau) {
    model_.stateFunc(tau, state_mat_[i+1], 
                     control_input_mat[i+1],
                     &(solution_vec[i*dim_constraints_and_disturbance_]), 
                     dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      state_mat_[i][j] = state_mat_[i+1][j] - delta_tau * dx_vec_[j];
    }
  }
  // Compute the Lagrange multiplier over the horizon on the basis of 
  // time, solution_vec and the state_vec.
  model_.phixFunc(time-horizon_length, state_mat_[0], lambda_mat_[0]);
  for (int i=0; i<dim_state_; ++i) {
    lambda_mat_[0][i] *= -1;
  }
  tau = time - horizon_length + delta_tau;
  for (int i=1; i<=N_; ++i, tau+=delta_tau) {
    model_.hxFunc(
        tau, state_mat_[i], control_input_mat[i], output_mat[i],
        &(solution_vec[(i-1)*dim_constraints_and_disturbance_]), 
        &(solution_vec[
            (i-1)*dim_constraints_and_disturbance_+dim_disturbance_]), 
        lambda_mat_[i-1], dx_vec_);
    for (int j=0; j<dim_state_; ++j) {
      lambda_mat_[i][j] = lambda_mat_[i-1][j] - delta_tau * dx_vec_[j];
    }
  }
  // Compute the erros in optimality over the horizon on the basis of the 
  // control_input_vec and the state_vec.
  tau = time;
  for (int i=N_; i>0; --i, tau-=delta_tau) {
    model_.hwFunc(
        tau, state_mat_[i], control_input_mat[i], output_mat[i],
        &(solution_vec[(i-1)*dim_constraints_and_disturbance_]), 
        &(solution_vec[
            (i-1)*dim_constraints_and_disturbance_+dim_disturbance_]), 
        lambda_mat_[i-1], 
        &(optimality_residual[(i-1)*dim_constraints_and_disturbance_]));
    model_.constraintFunc(
        tau, state_mat_[i], control_input_mat[i], output_mat[i], 
        &(solution_vec[(i-1)*dim_constraints_and_disturbance_]), 
        &(optimality_residual[
            (i-1)*dim_constraints_and_disturbance_+dim_disturbance_]));
  }
  model_.etaxFunc(time, state_mat_[N_], control_input_mat[N_], output_mat[N_],
                  dx_vec_);
  for (int i=0; i<dim_state_; ++i) {
    optimality_residual[N_*dim_constraints_and_disturbance_+i]
        = lambda_mat_[N_][i] - dx_vec_[i];
  }
}

void FiniteHorizonStateEstimationProblem::resetHorizonLength(
    const double T_f, const double alpha, const double initial_time) {
  horizon_.resetLength(T_f, alpha, initial_time);
}

void FiniteHorizonStateEstimationProblem::resetHorizonLength(
    const double initial_time) {
  horizon_.resetLength(initial_time);
}

int FiniteHorizonStateEstimationProblem::dim_solution() const {
  return dim_solution_;
}

int FiniteHorizonStateEstimationProblem::N() const {
  return N_;
}