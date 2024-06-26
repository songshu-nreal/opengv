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

template<typename P>
opengv::sac::Ransac<P>::Ransac(
    int maxIterations, double threshold, double probability) :
    SampleConsensus<P>(maxIterations, threshold, probability)
{}

template<typename P>
opengv::sac::Ransac<P>::~Ransac(){}


template<typename PROBLEM_T>
bool
opengv::sac::Ransac<PROBLEM_T>::computeModel(
    int debug_verbosity_level)
{
  typedef PROBLEM_T problem_t;
  typedef typename problem_t::model_t model_t;

  current_iterations_ = 0;
  int n_best_inliers_count = -INT_MAX;
  double k = 1.0;
  size_t max_num_iterations = std::max(max_iterations_, min_iterations_);
  const int kMinSampleSize = sac_model_->getSampleSize();

  std::vector<int> selection;
  model_t model_coefficients;

  if (sac_model_->getInitialModel(model_coefficients_)) {
    sac_model_->selectWithinDistance(
        model_coefficients_, threshold_, inliers_, inlier_distances_to_model_);
    double inlier_ratio = static_cast<double>(inliers_.size()) /
                          static_cast<double>(sac_model_->getIndices()->size());

    double p_no_outliers = 1.0 - pow(inlier_ratio, static_cast<double> (kMinSampleSize));
    p_no_outliers =
        (std::max) (std::numeric_limits<double>::epsilon(), p_no_outliers);
    // Avoid division by -Inf
    p_no_outliers =
        (std::min) (1.0 - std::numeric_limits<double>::epsilon(), p_no_outliers);
    // Avoid division by 0.
    k = log(1.0 - probability_) / log(p_no_outliers);
//    std::cout << "Inlier ratio: " << inliers_.size() << " / " << sac_model_->getIndices()->size()
//              << " = " << inlier_ratio << ", k = " << k << std::endl;
    ++current_iterations_;
  }

  int n_inliers_count = 0;
  unsigned skipped_count = 0;
  // supress infinite loops by just allowing 10 x maximum allowed iterations for
  // invalid model parameters!
  const unsigned max_skip = max_num_iterations * 10;

  // Iterate
  while( current_iterations_ < k && skipped_count < max_skip ) {
    // Get X samples which satisfy the model criteria
    sac_model_->getSamples( current_iterations_, selection );

    if(selection.empty()) {
      break;
    }

    // Search for inliers in the point cloud for the current plane model M
    if(!sac_model_->computeModelCoefficients( selection, model_coefficients )) {
      //++current_iterations;
      ++ skipped_count;
      continue;
    }

    // Select the inliers that are within threshold_ from the model
    //sac_model_->selectWithinDistance( model_coefficients, threshold_, inliers );
    //if(inliers.empty() && k > 1.0)
    //  continue;

    n_inliers_count = sac_model_->countWithinDistance(
        model_coefficients, threshold_ );

    // Better match ?
    if(n_inliers_count > n_best_inliers_count) {
      n_best_inliers_count = n_inliers_count;

      // Save the current model/inlier/coefficients selection as being the best so far
      model_              = selection;
      model_coefficients_ = model_coefficients;

      // Compute the k parameter (k=log(z)/log(1-w^n))
      double w = static_cast<double> (n_best_inliers_count) /
          static_cast<double> (sac_model_->getIndices()->size());
      double p_no_outliers = 1.0 - pow(w, static_cast<double> (selection.size()));
      p_no_outliers =
          (std::max) (std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by -Inf
      p_no_outliers =
          (std::min) (1.0 - std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by 0.
      k = log(1.0 - probability_) / log(p_no_outliers);
    }

    ++current_iterations_;

    if(current_iterations_ > max_num_iterations) {
      break;
    }
  }

  if(model_.empty()) {
    inliers_.clear();
    return (false);
  }

  // Get the set of inliers that correspond to the best model found so far
  sac_model_->selectWithinDistance(
      model_coefficients_, threshold_, inliers_, inlier_distances_to_model_);

  return (true);
}
