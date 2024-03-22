// Copyright (c) 2023, ETH Zurich and UNC Chapel Hill.
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
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#ifndef OPENGV_ABSOLUTE_POSE_MODULES_GP3P_LEE_MODULES_HPP_
#define OPENGV_ABSOLUTE_POSE_MODULES_GP3P_LEE_MODULES_HPP_

#include <stdlib.h>
#include <Eigen/Eigen>
#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include <opengv/types.hpp>


namespace opengv
{
namespace absolute_pose
{
namespace modules
{
namespace gp3p_lee
{

Eigen::Matrix<double, 9, 1> ComputeDepthsSylvesterCoeffs(
    const Eigen::Matrix<double, 3, 6>& K);

bool CheckParallelRays(const Eigen::Vector3d& ray1, const Eigen::Vector3d& ray2,
                       const Eigen::Vector3d& ray3);

bool CheckCollinearPoints(const Eigen::Vector3d& X1, const Eigen::Vector3d& X2,
                      const Eigen::Vector3d& X3);

void ComposePlueckerLine(const Eigen::Matrix3d& f,
                         const Eigen::Matrix3d& v,
                         plueckers_t& pluecker_vec);
Eigen::Vector3d PointFromPlueckerLineAndDepth(const opengv::pluecker_t& pluecker,
                         const double depth);

void  ComputePolynomialCoefficients(
    const opengv::plueckers_t& plueckers,
    const opengv::points_t& points3D,
    Eigen::Matrix<double, 3, 6>& K);

void ComputeDepthsSylvester(
    const Eigen::Matrix<double, 3, 6>& K,
    opengv::depths_t& depths);
}
}
}
}

#endif /* OPENGV_ABSOLUTE_POSE_MODULES_GP3P_LEE_MODULES_HPP_ */
