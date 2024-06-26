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
 * \file MultiRansac.hpp
 * \brief Implementation of the Ransac algorithm as outlined in [15]. This
 *        version is intended for use with the RelativeMultiAdapterBase, and
 *        attempts to do sampling in multiple camera-pairs in each hypothesis
 *        instantiation.
 */

#ifndef OPENGV_SAC_MULTIRANSAC_HPP_
#define OPENGV_SAC_MULTIRANSAC_HPP_

#include <vector>
#include <opengv/sac/MultiSampleConsensus.hpp>
#include <cstdio>

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
 * The Ransac sample consensus method, as outlined in [15]. This one is using
 * multi-indices for homogeneous sampling over groups of samples.
 */
template<typename PROBLEM_T>
class MultiRansac : public MultiSampleConsensus<PROBLEM_T>
{
public:
  /** A child of MultiSampleConsensusProblem */
  typedef PROBLEM_T problem_t;
  /** The model we trying to fit */
  typedef typename problem_t::model_t model_t;
  
  using MultiSampleConsensus<problem_t>::max_iterations_;
  using MultiSampleConsensus<problem_t>::threshold_;
  using MultiSampleConsensus<problem_t>::current_iterations_;
  using MultiSampleConsensus<problem_t>::sac_model_;
  using MultiSampleConsensus<problem_t>::model_;
  using MultiSampleConsensus<problem_t>::model_coefficients_;
  using MultiSampleConsensus<problem_t>::inliers_;
  using MultiSampleConsensus<problem_t>::probability_;

  /**
   * \brief Constructor.
   */
  MultiRansac(
      int maxIterations = 1000,
      double threshold = 1.0,
      double probability = 0.99 );
  /**
   * \brief Destructor.
   */
  virtual ~MultiRansac();

  /**
   * \brief Fit the model.
   */
  bool computeModel( int debug_verbosity_level = 0 );
};

} // namespace sac
} // namespace opengv

#include "implementation/MultiRansac.hpp"

#endif /* OPENGV_SAC_MULTIRANSAC_HPP_ */
