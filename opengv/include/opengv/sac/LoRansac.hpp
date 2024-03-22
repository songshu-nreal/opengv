/******************************************************************************
 * Authors:  Laurent Kneip & Paul Furgale                                     *
 * Contact:  kneip.laurent@gmail.com                                          *
 * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 * * Redistributions of source code must retain the above copyright           *
 *   notice, this list of conditions and the following disclaimer.            *
 * * Redistributions in binary form must reproduce the above copyright        *
 *   notice, this list of conditions and the following disclaimer in the      *
 *   documentation and/or other materials provided with the distribution.     *
 * * Neither the name of ANU nor the names of its contributors may be         *
 *   used to endorse or promote products derived from this software without   *
 *   specific prior written permission.                                       *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
 * SUCH DAMAGE.                                                               *
 ******************************************************************************/

//Note: has been derived from ROS

/**
 * \file LoRansac.hpp
 * \brief Implementation of the Ransac algorithm as outlined in [15].
 */

#ifndef OPENGV_SAC_LORANSAC_HPP_
#define OPENGV_SAC_LORANSAC_HPP_

#include <vector>
#include <opengv/sac/SampleConsensus.hpp>
#include <opengv/sac/sprt.hpp>
#include <cstdio>
#include <numeric>
#include <random>

/**
 * \brief The namespace of this library.
 */
namespace opengv
{
/**
 * \brief The namespace for the sample consensus methods.
 */
namespace sac
{
    // See Lebeda et al., Fixing the Locally Optimized RANSAC, BMVC, Table 1 for
// details on the variables.
class LORansacOptions {
public:
    LORansacOptions()
        : num_lo_steps_(10),
          threshold_multiplier_(std::sqrt(2.0)),
          num_lsq_iterations_(4),
          min_sample_multiplicator_(7),
          non_min_sample_multiplier_(3),
          lo_starting_iterations_(50u),
          final_least_squares_(false) {}
    int num_lo_steps_;
    double threshold_multiplier_;
    int num_lsq_iterations_;
    // The maximum number of data points used for least squares refinement is
    // min_sample_multiplicator_ * min_sample_size. Lebeda et al. recommend
    // setting min_sample_multiplicator_ to 7 (empirically determined for
    // epipolar geometry estimation.
    int min_sample_multiplicator_;
    // The solver needs to report the minimal size of the non-minimal sample
    // needed for its non-minimal solver. In practice, we draw a sample of size
    // min(non_min_sample_size * non_min_sample_multiplier_, N / 2), where N is
    // the number of data points.
    int non_min_sample_multiplier_;
    // As suggested in Sec. 4.4 in Lebeda et al., Local Optimization is only
    // performed after the first K_start iterations (set to 50 by Lebeda et al.)
    // to reduce overhead.
    uint32_t lo_starting_iterations_;
    bool final_least_squares_;
};
    
/**
 * The Ransac sample consensus method, as outlined in [15].
 */
template<typename PROBLEM_T>
class LoRansac : public SampleConsensus<PROBLEM_T>
{
public:
  /** A child of SampleConsensusProblem */
  typedef PROBLEM_T problem_t;
  /** The model we trying to fit */
  typedef typename problem_t::model_t model_t;
  
  using SampleConsensus<problem_t>::max_iterations_;
  using SampleConsensus<problem_t>::min_iterations_;
  using SampleConsensus<problem_t>::threshold_;
  using SampleConsensus<problem_t>::current_iterations_;
  using SampleConsensus<problem_t>::sac_model_;
  using SampleConsensus<problem_t>::model_;
  using SampleConsensus<problem_t>::model_coefficients_;
  using SampleConsensus<problem_t>::inliers_;
  using SampleConsensus<problem_t>::inlier_distances_to_model_;
  using SampleConsensus<problem_t>::score_;
  using SampleConsensus<problem_t>::probability_;

  /**
   * \brief Constructor.
   */
  LoRansac(
      int maxIterations = 1000,
      int minIterations = 100,
      double threshold = 1.0,
      double probability = 0.99);

  LoRansac(
      const LORansacOptions& lo_options,
      int maxIterations = 1000,
      int minIterations = 100,
      double threshold = 1.0,
      double probability = 0.99);

  LoRansac(
      const LORansacOptions& lo_options,
      const SPRT::Options& sprt_options,
      int maxIterations = 1000,
      int minIterations = 100,
      double threshold = 1.0,
      double probability = 0.99);

  /**
   * \brief Destructor.
   */
  virtual ~LoRansac();

  /**
   * \brief Fit the model.
   */
  bool computeModel( int debug_verbosity_level = 0 );

  void LocalOptimization(
        model_t* best_minimal_model, double* score_best_minimal_model) const;

  void ScoreModel(
        const model_t& model, double squared_inlier_threshold, double* score) const;

  void LeastSquaresFit(double thresh, model_t* model) const;

  void UpdateBestModel(
      double score_curr, const model_t& m_curr, double* score_best, model_t* m_best) const;

  void getModelScoreVec(
      const model_t& best_model_coefficients, std::vector<double>& best_model_scores,
      std::vector<double>& best_model_inlier_distances_to_model,
      std::vector<int>& best_model_inliers) const;

  bool UpdateWithSPRT(
      const SPRT& sprt, std::vector<double>& best_model_scores) const;

  LORansacOptions lo_options_;
  // SPRT
  SPRT::Options sprt_options_;
};

} // namespace sac
} // namespace opengv

#include "implementation/LoRansac.hpp"
#endif /* OPENGV_SAC_LORANSAC_HPP_ */
