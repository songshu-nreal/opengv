// Copyright (c) 2019, Torsten Sattler
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of the copyright holder nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// author: Torsten Sattler, torsten.sattler.de@googlemail.com

// Computes the number of RANSAC iterations required for a given inlier
// ratio, the probability of missing the best model, and sample size.
// Assumes that min_iterations <= max_iterations.
inline uint32_t NumRequiredIterations(const double inlier_ratio,
                               const double prob_missing_best_model,
                               const int sample_size,
                               const uint32_t min_iterations,
                               const uint32_t max_iterations) {
  if (inlier_ratio <= 0.0) {
    return max_iterations;
  }
  if (inlier_ratio >= 1.0) {
    return min_iterations;
  }

  const double kProbNonInlierSample =
      1.0 - std::pow(inlier_ratio, static_cast<double>(sample_size));
  // If the probability of sampling a non-all-inlier sample is at least
  // 0.99999999999999, RANSAC will take at least 1e+13 iterations for
  // realistic values for prob_missing_best_model (0.5 or smaller).
  // In practice, max_iterations will be smaller.
  if (kProbNonInlierSample >= 0.99999999999999) {
    return max_iterations;
  }

  const double kLogNumerator = std::log(prob_missing_best_model);
  const double kLogDenominator = std::log(kProbNonInlierSample);

  double num_iters = std::ceil(kLogNumerator / kLogDenominator + 0.5);
  uint32_t num_req_iterations =
      std::min(static_cast<uint32_t>(num_iters), max_iterations);
  num_req_iterations = std::max(min_iterations, num_req_iterations);
  return num_req_iterations;
}

// This function implements Fisher-Yates shuffling, implemented "manually"
// here following: https://lemire.me/blog/2016/10/10/a-case-study-in-the-
// performance-cost-of-abstraction-cs-stdshuffle/
inline void RandomShuffle(std::mt19937* rng, std::vector<int>* random_sample) {
  std::vector<int>& sample = *random_sample;
  const int kNumElements = static_cast<int>(sample.size());
  for (int i = 0; i < (kNumElements - 1); ++i) {
    std::uniform_int_distribution<int> dist(i, kNumElements - 1);
    int idx = dist(*rng);
    std::swap(sample[i], sample[idx]);
  }
}

inline void RandomShuffleAndResize(const int target_size, std::mt19937* rng,
                            std::vector<int>* random_sample) {
  RandomShuffle(rng, random_sample);
  random_sample->resize(target_size);
}

inline void RandomShuffleAndResize(
    const std::vector<int>& sample_sizes, std::mt19937* rng,
    std::vector<std::vector<int>>* random_sample) {
  const int kNumDataTypes = static_cast<int>(random_sample->size());
  for (int i = 0; i < kNumDataTypes; ++i) {
    const int kNumData = static_cast<int>((*random_sample)[i].size());
    const int kSampleSize = std::min(kNumData, sample_sizes[i]);
    RandomShuffleAndResize(kSampleSize, rng, &((*random_sample)[i]));
  }
}

template<typename P>
opengv::sac::LoRansac<P>::LoRansac(
    int maxIterations, int minIterations, double threshold, double probability) :
    SampleConsensus<P>(maxIterations, minIterations, threshold, probability) {}

template<typename P>
opengv::sac::LoRansac<P>::LoRansac(
    const LORansacOptions& lo_options, int maxIterations, int minIterations,
    double threshold, double probability):
    SampleConsensus<P>(maxIterations, minIterations, threshold, probability),
        lo_options_(lo_options) {
  assert(lo_options_.non_min_sample_multiplier_ > 0);
  assert(lo_options_.min_sample_multiplicator_ > 0);
  assert(lo_options_.threshold_multiplier_ > 0);
}

template<typename P>
opengv::sac::LoRansac<P>::~LoRansac()= default;

template<typename PROBLEM_T>
bool
opengv::sac::LoRansac<PROBLEM_T>::computeModel(
    int debug_verbosity_level)
{
  typedef PROBLEM_T problem_t;
  typedef typename problem_t::model_t model_t;
  // Initialize parameters
  const int kMinSampleSize = sac_model_->getSampleSize();
  current_iterations_ = 0; // num_iteration
  double k = 1.0; // max_num_iterations compute by probability
  uint32_t max_num_iterations = std::max(max_iterations_, min_iterations_);

  // Best model found so far
  model_t best_model_coefficients; // best minimal model
  double best_model_score = std::numeric_limits<double>::max(); // best minimal score
  std::vector<int> best_model_inliers;
  std::vector<double> best_model_inlier_distances_to_model;

  // Best mininal model found so far
  model_t best_min_model_coefficients; // best minimal model
  double best_min_model_score = std::numeric_limits<double>::max(); // best minimal score
  std::vector<int> best_min_model_inliers;

  // The number of iterations required to achieve the desired probability of
  // local optimization
  int number_lo_iterations = 0;

  unsigned skipped_count = 0;
  // supress infinite loops by just allowing 10 x maximum allowed iterations for
  // invalid model parameters!
  const unsigned max_skip = max_num_iterations * 10;

  // Iterate
  while( current_iterations_ < k && skipped_count < max_skip ) {
    if (current_iterations_ == lo_options_.lo_starting_iterations_ &&
        best_min_model_score < std::numeric_limits<double>::max()) {
      ++number_lo_iterations;
      LocalOptimization(&best_model_coefficients, &best_model_score);

      // Updates the number of RANSAC iterations.
      sac_model_->selectWithinDistance(
          best_model_coefficients, threshold_, best_model_inliers, best_model_inlier_distances_to_model);
      double inlier_ratio = static_cast<double>(best_model_inliers.size()) /
                            static_cast<double>(sac_model_->getIndices()->size());
      k = NumRequiredIterations(
          inlier_ratio, 1.0 - probability_, kMinSampleSize, min_iterations_,
          max_iterations_);
    }

    std::vector<int> current_selection;
    model_t current_model_coefficients;
    double current_model_score = std::numeric_limits<double>::max();

    // Get X samples which satisfy the model criteria
    sac_model_->getSamples(current_iterations_, current_selection);

    if (current_selection.empty()) {
      fprintf(stderr,
              "[sm::RandomSampleConsensus::computeModel] No samples could be selected!\n");
      break;
    }

    // Search for inliers in the point cloud for the current plane model M
    if (!sac_model_->computeModelCoefficients(current_selection, current_model_coefficients)) {
      //++iterations_;
      ++skipped_count;
      continue;
    }

    // Select the inliers that are within threshold_ from the model
    //sac_model_->selectWithinDistance( model_coefficients, threshold_, inliers );
    //if(inliers.empty() && k > 1.0)
    //  continue;

    // compute score
    ScoreModel(current_model_coefficients, threshold_, &current_model_score);

    // Updates the best model found so far.
    if (current_model_score < best_min_model_score ||
        current_iterations_ == lo_options_.lo_starting_iterations_) {
      const bool kBestMinModel = current_model_score < best_min_model_score;

      // Better match ?
      if (kBestMinModel) {
        // Save the current model/inlier/coefficients selection as being the best so far
        best_min_model_score = current_model_score;
        best_min_model_coefficients = current_model_coefficients;

        // Updates the best model.
        UpdateBestModel(best_min_model_score, best_min_model_coefficients,
                        &best_model_score, &best_model_coefficients);
      }

      const bool kRunLO =
          (current_iterations_ >= lo_options_.lo_starting_iterations_ &&
           best_min_model_score < std::numeric_limits<double>::max());

      if ((!kBestMinModel) && (!kRunLO)) {
        continue;
      }

      // Performs local optimization. By construction, the local optimization
      // method returns the best model between all models found by local
      // optimization and the input model, i.e., score_refined_model <=
      // best_min_model_score holds.
      if (kRunLO) {
        ++number_lo_iterations;
        double score = best_min_model_score;

        LocalOptimization(&best_min_model_coefficients, &score);

        // Updates the best model.
        UpdateBestModel(score, best_min_model_coefficients, &best_model_score,
                        &best_model_coefficients);
      }

      // Updates the number of RANSAC iterations.
      sac_model_->selectWithinDistance(
          best_model_coefficients, threshold_, best_model_inliers, best_model_inlier_distances_to_model);

      double inlier_ratio = static_cast<double>(best_model_inliers.size()) /
                            static_cast<double>(sac_model_->getIndices()->size());
      k = NumRequiredIterations(
          inlier_ratio, 1.0 - probability_, kMinSampleSize, min_iterations_,
          max_iterations_);
    }

    if(debug_verbosity_level > 1)
      fprintf(stdout,
              "[sm::RandomSampleConsensus::computeModel] Trial %d out of %f (best is: %d so far).\n",
              current_iterations_, k, best_model_inliers.size() );
    if(current_iterations_ > k)
    {
      if(debug_verbosity_level > 0)
        fprintf(stdout,
                "[sm::RandomSampleConsensus::computeModel] RANSAC reached the maximum number of trials.\n");
      break;
    }
    ++current_iterations_;
  }

  // As proposed by Lebeda et al., Local Optimization is not executed in
  // the first lo_starting_iterations_ iterations. If LO-MSAC needs less than
  // lo_starting_iterations_ iterations, we run LO now.
  if (current_iterations_ <= lo_options_.lo_starting_iterations_ &&
      best_model_score < std::numeric_limits<double>::max()) {
    ++number_lo_iterations;
    LocalOptimization(&best_model_coefficients, &best_model_score);
  }

//  fprintf(stdout,
//          "[sm::RandomSampleConsensus::computeModel] Number of local optimization trials: %d\n",
//          number_lo_iterations);

  if (lo_options_.final_least_squares_) {
    model_t refined_model = best_model_coefficients;
    sac_model_->optimizeModelCoefficients(best_model_inliers, best_model_coefficients, refined_model);

    double score = std::numeric_limits<double>::max();
    ScoreModel(refined_model, threshold_, &score);
    if (score < best_model_score) {
      best_model_score = score;
      best_model_coefficients = refined_model;

      sac_model_->selectWithinDistance(
          best_model_coefficients, threshold_, best_model_inliers, best_model_inlier_distances_to_model);

    }
  }
//  fprintf(stdout,
//          "[sm::RandomSampleConsensus::computeModel] final_least_squares.\n");

  // Save the best model found
  model_coefficients_ = best_model_coefficients;
  inliers_ = best_model_inliers;
  inlier_distances_to_model_ = best_model_inlier_distances_to_model;
  score_ = best_model_score;

  // Compute the final set of inliers that correspond to the best model found
  return (true);
}

// See algorithms 2 and 3 in Lebeda et al.
// The input model is overwritten with the refined model if the latter is
// better, i.e., has a lower score.
template<typename PROBLEM_T>
void opengv::sac::LoRansac<PROBLEM_T>::LocalOptimization(
    model_t* best_minimal_model,
    double* score_best_minimal_model) const {
  // Check non minimal solver sample size.
  const int kNumData = sac_model_->getIndices()->size();
  const int kMinNonMinSampleSize = sac_model_->getLoSampleSize();
  if (kMinNonMinSampleSize > kNumData) {
    return;
  }

  const int kMinSampleSize = sac_model_->getSampleSize();
  const double kThresh = threshold_;
  const double kThreshMult = lo_options_.threshold_multiplier_;

  // Performs an initial least squares fit of the best model found by the
  // minimal solver so far and then determines the inliers to that model
  // under a (slightly) relaxed inlier threshold.
  model_t m_init = *best_minimal_model;
  LeastSquaresFit(kThresh * kThreshMult, &m_init);

  double score = std::numeric_limits<double>::max();
  ScoreModel(m_init, kThresh, &score);
  UpdateBestModel(score, m_init, score_best_minimal_model,
                  best_minimal_model);

  std::vector<int> inliers_base;
  std::vector<double> inlier_distances_to_model;
  sac_model_->selectWithinDistance(
      m_init, kThresh * kThreshMult, inliers_base, inlier_distances_to_model);

  // Determines the size of the non-miminal samples drawn in each LO step.
  const int kNonMinSampleSize =
      std::max(kMinSampleSize,
               std::min(kMinSampleSize * lo_options_.non_min_sample_multiplier_,
                        static_cast<int>(inliers_base.size()) / 2));

  // Performs the actual local optimization (LO).
  std::vector<int> sample;
  for (int r = 0; r < lo_options_.num_lo_steps_; ++r) {
    sample = inliers_base;
    RandomShuffleAndResize(kNonMinSampleSize, &sac_model_->rng_alg_, &sample);

    model_t m_non_min;
    // Search for inliers in the point cloud for the current plane model M
    const bool local_optimizing = true;
    if (!sac_model_->computeLoModelCoefficients(sample, m_non_min)) {
      continue;
    }

    ScoreModel(m_non_min, kThresh, &score);
    UpdateBestModel(score, m_non_min, score_best_minimal_model,
                    best_minimal_model);

    // Iterative least squares refinement.
    LeastSquaresFit(kThresh, &m_non_min);

    // The current threshold multiplier and its update.
    double thresh = kThreshMult * kThresh;
    double thresh_mult_update =
        (kThreshMult - 1.0) * kThresh /
        static_cast<int>(lo_options_.num_lsq_iterations_ - 1);
    for (int i = 0; i < lo_options_.num_lsq_iterations_; ++i) {
      LeastSquaresFit(thresh, &m_non_min);

      ScoreModel(m_non_min, kThresh, &score);
      UpdateBestModel(score, m_non_min, score_best_minimal_model,
                      best_minimal_model);
      thresh -= thresh_mult_update;
    }
  }
}

template<typename PROBLEM_T>
void opengv::sac::LoRansac<PROBLEM_T>::ScoreModel(
    const model_t& model, const double inlier_threshold, double* score) const {
  const int kNumData = sac_model_->getIndices()->size();
  std::vector<double> scores;
  sac_model_->getDistancesToModel(model, scores);
  for (int i = 0; i < kNumData; ++i) {
    scores[i] = std::min(scores[i], inlier_threshold);
  }
  *score = std::accumulate(scores.begin(), scores.end(), 0.0);
}

template<typename PROBLEM_T>
void opengv::sac::LoRansac<PROBLEM_T>::LeastSquaresFit(
    const double thresh, model_t* model) const {
  const int kLSqSampleSize =
      lo_options_.min_sample_multiplicator_ * sac_model_->getSampleSize();
  std::vector<int> inliers;
  std::vector<double> inlier_distances_to_model;
  sac_model_->selectWithinDistance(
      *model, thresh, inliers, inlier_distances_to_model);
  if (inliers.size() < sac_model_->getSampleSize()) {
    return;
  }
  int lsq_data_size = std::min(kLSqSampleSize, static_cast<int>(inliers.size()));
  RandomShuffleAndResize(lsq_data_size, &sac_model_->rng_alg_, &inliers);

  model_t optimized_model;
  sac_model_->optimizeModelCoefficients(inliers, *model, optimized_model);
  *model = optimized_model;
}

template<typename PROBLEM_T>
void opengv::sac::LoRansac<PROBLEM_T>::UpdateBestModel(
    const double score_curr, const model_t& m_curr, double* score_best, model_t* m_best) const {
  if (score_curr < *score_best) {
    *score_best = score_curr;
    *m_best = m_curr;
  }
}