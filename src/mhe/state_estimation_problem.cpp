#include "state_estimation_problem.hpp"

StateEstimationProblem::StateEstimationProblem()
  : model_(),
    dim_state_(model_.dimState()),
    dim_control_input_(model_.dimControlInput()),
    dim_constraints_(model_.dimConstraints()),
    dim_output_(model_.dimOutput()),
    dim_disturbance_(model_.dimDisturbance()) {
}

int StateEstimationProblem::dim_state() const {
  return dim_state_;
}

int StateEstimationProblem::dim_control_input() const {
  return dim_control_input_;
}

int StateEstimationProblem::dim_constraints() const {
  return dim_constraints_;
}

int StateEstimationProblem::dim_output() const {
  return dim_output_;
}

int StateEstimationProblem::dim_disturbance() const {
  return dim_disturbance_;
}