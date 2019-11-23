#include "time_series_data.hpp"

TimeSeriesData::TimeSeriesData(const int dim_data, const int max_num_data)
  : dim_data_(dim_data),
    max_num_data_(max_num_data),
    current_num_data_(0),
    time_series_data_(linearalgebra::NewMatrix(max_num_data, dim_data)),
    time_series_(linearalgebra::NewVector(max_num_data)) {
}

TimeSeriesData::~TimeSeriesData() {
  linearalgebra::DeleteMatrix(time_series_data_);
  linearalgebra::DeleteVector(time_series_);
}

void TimeSeriesData::appendData(const double time, const double* data) {
  if (time<time_series_[0]) {
    std::cout << "error! time must be larger than time_series_[0]!" << std::endl;
    return;
  }
  if (current_num_data_>=max_num_data_) {
    removeOldestData();
  }
  for (int i=0; i<dim_data_; ++i) {
    time_series_data_[current_num_data_][i] = data[i];
  }
  time_series_[current_num_data_] = time;
  ++current_num_data_;
}

void TimeSeriesData::generateUniformlySpacedTimeSeriesData(
    const double time, const double T, const int N, 
    double** generated_data) const {
  if (current_num_data_==0) {
    for (int i=0; i<N; ++i) {
      for (int j=0; j<dim_data_; ++j) {
        generated_data[i][j] = 0.0;
      }
    }
    return;
  } 
  else if (current_num_data_==1) {
    for (int i=0; i<N; ++i) {
      for (int j=0; j<dim_data_; ++j) {
        generated_data[i][j] = time_series_data_[0][j];
      }
    }
    return;
  }
  double delta_tau = T / N;
  int max_interpolation_range = 0;
  int min_interpolation_range = current_num_data_;
  if (time-T<time_series_[0]) {
    max_interpolation_range 
        = static_cast<int>((time_series_[0]-time+T)/delta_tau) + 1;
  }
  if (time>time_series_[current_num_data_-1]) {
    min_interpolation_range
      = static_cast<int>((time-time_series_[current_num_data_-1])/delta_tau) + 1;
  }
  double tau = time - T;
  for (int i=0; i<min_interpolation_range; ++i, tau+=delta_tau) {
    linearExtrapolateLowerData(time_series_[0], time_series_data_[0], 
                               time_series_[1], time_series_data_[1], tau,
                               generated_data[i]);
  }
  for (int i=min_interpolation_range; i<max_interpolation_range; 
          ++i, tau+=delta_tau) {
    int lower_idx = lowerTimeIndex(tau);
    int upper_idx = upperTimeIndex(tau);
    linearInterpolate(time_series_[lower_idx], time_series_data_[lower_idx], 
                      time_series_[upper_idx], time_series_data_[upper_idx], 
                      tau, generated_data[i]);
  }
  for (int i=max_interpolation_range; i<N; ++i) {
    linearExtrapolateLowerData(time_series_[current_num_data_-2], 
                               time_series_data_[current_num_data_-2], 
                               time_series_[current_num_data_-1], 
                               time_series_data_[current_num_data_-1], tau,
                               generated_data[i]);
  }
}

int TimeSeriesData::size() const {
  return current_num_data_;
}

inline void TimeSeriesData::linearInterpolate(const double lower_time, 
                                              const double* lower_data,
                                              const double upper_time,
                                              const double* upper_data, 
                                              const double inter_time, 
                                              double* result_data) const {
  for (int i=0; i<dim_data_; ++i) {
    result_data[i] 
        = lower_data[i] 
            + (upper_data[i]-lower_data[i]) * (inter_time-lower_time)
              / (upper_time-lower_time);
  }
}

inline void TimeSeriesData::linearExtrapolateLowerData(
    const double lower_time, const double* lower_data,const double upper_time,
    const double* upper_data, const double extra_time, 
    double* result_data) const {
  for (int i=0; i<dim_data_; ++i) {
    result_data[i] 
        = lower_data[i] 
            - (upper_data[i]-lower_data[i]) * (lower_time-extra_time)
              / (upper_time-lower_time);
  }
}

inline void TimeSeriesData::linearExtrapolateUpperData(
    const double lower_time, const double* lower_data,const double upper_time,
    const double* upper_data, const double extra_time, 
    double* result_data) const {
  for (int i=0; i<dim_data_; ++i) {
    result_data[i] 
        = upper_data[i] 
            + (upper_data[i]-lower_data[i]) * (extra_time-upper_time)
              / (upper_time-lower_time);
  }
}

inline void TimeSeriesData::removeOldestData() {
  for (int i=1; i<current_num_data_; ++i) {
    for (int j=0; j<dim_data_; ++j) {
      time_series_data_[i-1][j] = time_series_data_[i][j];
    }
  }
  for (int i=1; i<current_num_data_; ++i) {
    time_series_[i-1] = time_series_[i];
  }
  --current_num_data_;
}

int TimeSeriesData::lowerTimeIndex(const double time) const {
  for (int i=0; i<current_num_data_; ++i) {
    if (time_series_[i] >= time) {
      return i - 1;
    }
  }
  return current_num_data_ - 1;
}

int TimeSeriesData::upperTimeIndex(const double time) const {
  for (int i=current_num_data_-1; i>=0; --i) {
    if (time_series_[i] < time) {
      return i + 1;
    }
  }
  return 0;
}