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
 * \file MultiSampleConsensus.hpp
 * \brief This is a base class for sample consensus methods such as Ransac.
 *        Derivatives call the three basic functions of a sample-consensus
 *        problem (sample drawing, computation of a hypothesis, and verification
 *        of a hypothesis). This version is intended for use with the
 *        RelativeMultiAdapterBase, and attempts to do sampling in multiple
 *        camera-pairs in each hypothesis instantiation.
 */

#ifndef OPENGV_SAC_MULTISAMPLECONSENSUS_HPP_
#define OPENGV_SAC_MULTISAMPLECONSENSUS_HPP_

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

/**
 * Super-class for sample consensus methods, such as Ransac. This one is using
 * multi-indices for homogeneous sampling over groups of samples.
 */
template<typename PROBLEM_T>
class MultiSampleConsensus
{
public:
  /** A child of MultiSampleConsensusProblem */
  typedef PROBLEM_T problem_t;
  /** The model we trying to fit */
  typedef typename problem_t::model_t model_t;

  /**
   * \brief Constructor.
   * \param[in] maxIterations The maximum number of hypothesis generations
   * \param[in] threshold Some threshold value for classifying samples as
   *                      an inlier or an outlier.
   * \param[in] probability The probability of being able to draw at least one
   *                        sample that is free of outliers (see [15])
   */
  MultiSampleConsensus(
      int maxIterations = 1000,
      double threshold = 1.0,
      double probability = 0.99 );
  /**
   * \brief Destructor
   */
  virtual ~MultiSampleConsensus();

  /**
   * \brief Fit the model to the data.
   * \param[in] debug_verbosity_level Sets the verbosity level.
   * \return bool True if success.
   */
  virtual bool computeModel( int debug_verbosity_level = 0 ) = 0;

  // \todo accessors
  //private:
  
  /** the maximum number of iterations */
  size_t max_iterations_;
  /** the current number of iterations */
  size_t current_iterations_;
  /** the threshold for classifying inliers */
  double threshold_;
  /** the current probability (defines remaining iterations) */
  double probability_;
  /** the currently best model coefficients */
  model_t model_coefficients_;
  /** the group-wise indices for the currently best hypothesis */
  std::vector< std::vector<int> > model_;
  /** the group-wise indices of the samples that have been clasified as inliers */
  std::vector< std::vector<int> > inliers_;
  /** the multi-sample-consensus problem we are trying to solve */
  std::shared_ptr<PROBLEM_T> sac_model_;
};

} // namespace sac
} // namespace opengv

#include "implementation/MultiSampleConsensus.hpp"

#endif /* OPENGV_MULTISAMPLECONSENSUS_HPP_ */
