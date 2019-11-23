#ifndef TIME_SERIES_DATA_H
#define TIME_SERIES_DATA_H

#include <iostream>
#include "linear_algebra.hpp"

class TimeSeriesData {
public:
  TimeSeriesData(const int dim_data, const int max_num_data);
  ~TimeSeriesData();

  void appendData(const double time, const double* data);

  void generateUniformlySpacedTimeSeriesData(const double time, const double T, 
                                             const int N, 
                                             double** generated_data) const;

  int size() const;

  // Prohibits copy due to memory allocation.
  TimeSeriesData(const TimeSeriesData&) = delete;
  TimeSeriesData& operator=(const TimeSeriesData&) = delete;

private:
  int dim_data_, max_num_data_, current_num_data_;
  double **time_series_data_;
  double *time_series_;

  inline void linearInterpolate(const double lower_time, 
                                const double* lower_data,
                                const double upper_time,
                                const double* upper_data, 
                                const double inter_time, 
                                double* result_data) const;

  inline void linearExtrapolateLowerData(const double lower_time, 
                                         const double* lower_data,
                                         const double upper_time,
                                         const double* upper_data, 
                                         const double extra_time, 
                                         double* result_data) const;

  inline void linearExtrapolateUpperData(const double lower_time, 
                                         const double* lower_data,
                                         const double upper_time,
                                         const double* upper_data, 
                                         const double extra_time, 
                                         double* result_data) const;

  inline void removeOldestData();

  int lowerTimeIndex(const double time) const;
  int upperTimeIndex(const double time) const;
};

#endif // TIME_SERIES_DATA_H 