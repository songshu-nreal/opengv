#ifndef OPENGV_SAC_SPRT_H
#define OPENGV_SAC_SPRT_H

#include <cmath>
#include <cstddef>
#include <vector>

/**
 * \brief The namespace of this library.
 */
namespace opengv {
/**
 * \brief The namespace for the sample consensus methods.
 */
namespace sac {
// Sequential Probability Ratio Test as proposed in
//
//   "Randomized RANSAC with Sequential Probability Ratio Test",
//   Matas et al., 2005
class SPRT {
public:
  struct Options {
    bool is_sprt = false;
    // Probability of rejecting a good model.
    double delta = 0.01;

    // A priori assumed minimum inlier ratio
    double epsilon = 0.1;

    // The ratio of the time it takes to estimate a model from a random sample
    // over the time it takes to decide whether one data sample is an
    // inlier or not. Matas et al. propose 200 for the 7-point algorithm.
    double eval_time_ratio = 100;

    // Number of models per random sample, that have to be verified. E.g. 1-3
    // for the 7-point fundamental matrix algorithm, or 1-10 for the 5-point
    // essential matrix algorithm.
    int num_models_per_sample = 1;
  };

  explicit SPRT(const Options &options);

  void Update(const Options &options);

  bool Evaluate(const std::vector<double> &residuals, const double max_residual,
                size_t *num_inliers, size_t *num_eval_samples) const;

private:
  void UpdateDecisionThreshold();

  Options options_;
  double delta_epsilon_;
  double delta_1_epsilon_1_;
  double decision_threshold_;
};
} // namespace sac
} // namespace opengv

#include "implementation/sprt.hpp"

#endif //OPENGV_SAC_SPRT_H
