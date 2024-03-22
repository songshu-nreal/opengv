/******************************************************************************
 * Author:   Laurent Kneip                                                    *
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


#include <opengv/absolute_pose/modules/gp3p_lee/modules.hpp>
#include <opengv/math/polynomial.h>
#include <array>

// Check whether the rays are close to parallel.
bool opengv::absolute_pose::modules::gp3p_lee::CheckParallelRays(
    const Eigen::Vector3d& ray1,
    const Eigen::Vector3d& ray2,
    const Eigen::Vector3d& ray3) {
  const double kParallelThreshold = 1e-5;
  return ray1.cross(ray2).isApproxToConstant(0, kParallelThreshold) &&
         ray1.cross(ray3).isApproxToConstant(0, kParallelThreshold);
}

// Check whether the points are close to collinear.
bool opengv::absolute_pose::modules::gp3p_lee::CheckCollinearPoints(
    const Eigen::Vector3d& X1,
    const Eigen::Vector3d& X2,
    const Eigen::Vector3d& X3) {
  const double kMinNonCollinearity = 1e-5;
  const Eigen::Vector3d X12 = X2 - X1;
  const double non_collinearity_measure =
      X12.cross(X1 - X3).squaredNorm() / X12.squaredNorm();
  return non_collinearity_measure < kMinNonCollinearity;
}

void opengv::absolute_pose::modules::gp3p_lee::ComposePlueckerLine(
    const Eigen::Matrix3d& f,
    const Eigen::Matrix3d& v,
    opengv::plueckers_t& pluecker_vec) {
  pluecker_vec.clear();
  for (int i = 0; i < 3; ++i) {
    const Eigen::Vector3d proj_center = v.col(i);
    const Eigen::Vector3d bearing_normalized = f.col(i);
    opengv::pluecker_t pluecker;
    pluecker << bearing_normalized, proj_center.cross(bearing_normalized);
    pluecker_vec.push_back(pluecker);
  }

  return;
}


Eigen::Vector3d opengv::absolute_pose::modules::gp3p_lee::PointFromPlueckerLineAndDepth(
    const opengv::pluecker_t& pluecker, const double depth) {
  return pluecker.head<3>().cross(pluecker.tail<3>()) +
         depth * pluecker.head<3>();
}

// Compute the coefficients from the system of 3 equations, nonlinear in the
// depth of the points. Inputs are three Pluecker lines and the locations of
// their corresponding points in 3D. The system of equations comes from the
// distance constraints between 3D points:
//
//    || f_i - f_j ||^2 = || (q_i x q_i' + lambda_i * q_i) -
//                           (q_j x q_j' + lambda_j * q_j) ||^2
//
// where [q_i; q_i'] is the Pluecker coordinate of bearing i and f_i is the
// coordinate of the corresponding 3D point in the global coordinate system. A
// 3D point in the local camera coordinate system along this line is
// parameterized through the depth scalar lambda_i as:
//
//    B_fi = q_i x q_i' + lambda_i * q_i.
//
void opengv::absolute_pose::modules::gp3p_lee::ComputePolynomialCoefficients(
    const opengv::plueckers_t& plueckers,
    const opengv::points_t& points3D,
    Eigen::Matrix<double, 3, 6>& K) {
  assert(plueckers.size()==3);
  assert(points3D.size()==3);

  const std::array<int, 3> is = {{0, 0, 1}};
  const std::array<int, 3> js = {{1, 2, 2}};

  for (int k = 0; k < 3; ++k) {
    const int i = is[k];
    const int j = js[k];
    const Eigen::Vector3d moment_difference =
        plueckers[i].head<3>().cross(plueckers[i].tail<3>()) -
        plueckers[j].head<3>().cross(plueckers[j].tail<3>());
    K(k, 0) = 1;
    K(k, 1) = -2 * plueckers[i].head<3>().dot(plueckers[j].head<3>());
    K(k, 2) = 2 * moment_difference.dot(plueckers[i].head<3>());
    K(k, 3) = 1;
    K(k, 4) = -2 * moment_difference.dot(plueckers[j].head<3>());
    K(k, 5) = moment_difference.squaredNorm() -
        (points3D[i] - points3D[j]).squaredNorm();
  }

  return;
}

// Solve quadratics of the form: x^2 + bx + c = 0.
int SolveQuadratic(const double b, const double c, double* roots) {
  const double delta = b * b - 4 * c;
  // Do not allow complex solutions.
  if (delta >= 0) {
    const double sqrt_delta = std::sqrt(delta);
    roots[0] = -0.5 * (b + sqrt_delta);
    roots[1] = -0.5 * (b - sqrt_delta);
    return 2;
  } else {
    return 0;
  }
}

// Given lambda_j, return the values for lambda_i, where:
//     k1 lambda_i^2 + (k2 lambda_j + k3) lambda_i
//      + k4 lambda_j^2 + k5 lambda_j + k6          = 0.
void ComputeLambdaValues(const Eigen::Matrix<double, 3, 6>::ConstRowXpr& k,
                         const double lambda_j,
                         std::vector<double>* lambdas_i) {
  // Note that we solve x^2 + bx + c = 0, since k(0) is one.
  double roots[2];
  const int num_solutions =
      SolveQuadratic(k(1) * lambda_j + k(2),
                     lambda_j * (k(3) * lambda_j + k(4)) + k(5), roots);
  for (int i = 0; i < num_solutions; ++i) {
    if (roots[i] > 0) {
      lambdas_i->push_back(roots[i]);
    }
  }
}

// Given the coefficients of the polynomial system return the depths of the
// points along the Pluecker lines. Use Sylvester resultant to get and 8th
void opengv::absolute_pose::modules::gp3p_lee::ComputeDepthsSylvester(
    const Eigen::Matrix<double, 3, 6>& K,
    opengv::depths_t& depths) {
  const Eigen::Matrix<double, 9, 1> coeffs =
      opengv::absolute_pose::modules::gp3p_lee::ComputeDepthsSylvesterCoeffs(K);

  Eigen::VectorXd roots_real;
  Eigen::VectorXd roots_imag;
  if (!opengv::math::FindPolynomialRootsCompanionMatrix(coeffs, &roots_real, &roots_imag)) {
    return;
  }

  depths.reserve(roots_real.size());
  for (Eigen::VectorXd::Index i = 0; i < roots_real.size(); ++i) {
    const double kMaxRootImagRatio = 1e-3;
    if (std::abs(roots_imag(i)) > kMaxRootImagRatio * std::abs(roots_real(i))) {
      continue;
    }

    const double lambda_3 = roots_real(i);
    if (lambda_3 <= 0) {
      continue;
    }

    std::vector<double> lambdas_2;
    ComputeLambdaValues(K.row(2), lambda_3, &lambdas_2);

    // Now we have two depths, lambda_2 and lambda_3. From the two remaining
    // equations, we must get the same lambda_1, otherwise the solution is
    // invalid.
    for (const double lambda_2 : lambdas_2) {
      std::vector<double> lambdas_1_1;
      ComputeLambdaValues(K.row(0), lambda_2, &lambdas_1_1);
      std::vector<double> lambdas_1_2;
      ComputeLambdaValues(K.row(1), lambda_3, &lambdas_1_2);
      for (const double lambda_1_1 : lambdas_1_1) {
        for (const double lambda_1_2 : lambdas_1_2) {
          const double kMaxLambdaRatio = 1e-2;
          if (std::abs(lambda_1_1 - lambda_1_2) <
              kMaxLambdaRatio * std::max(lambda_1_1, lambda_1_2)) {
            const double lambda_1 = (lambda_1_1 + lambda_1_2) / 2;
            depths.emplace_back(lambda_1, lambda_2, lambda_3);
          }
        }
      }
    }
  }

  return;
}
