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


#include <opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp>
#include <opengv/absolute_pose/methods.hpp>

bool opengv::sac_problems::
absolute_pose::AbsolutePoseSacProblem::getInitialModel(model_t & outModel) const {
  if (_adapter.hasT()) {
    rotation_t rotation = _adapter.getR();
    translation_t translation = _adapter.gett();
    outModel.col(3) = translation;
    outModel.block<3, 3>(0, 0) = rotation;
  }
  return _adapter.hasT();
}

bool opengv::sac_problems::
absolute_pose::AbsolutePoseSacProblem::computeModelCoefficients(
    const std::vector<int> &indices,
    algorithm_t algorithm,
    model_t & outModel) const {
  transformations_t solutions;

  switch(algorithm) {
    case TWOPT:
    {
      rotation_t rotation = _adapter.getR();
      translation_t translation = opengv::absolute_pose::p2p(_adapter,indices);

      transformation_t solution;
      translation_t t_bc = _adapter.getCamOffset(indices[0]);
      rotation_t R_bc = _adapter.getCamRotation(indices[0]);
      solution.col(3) = translation - rotation * R_bc.transpose() * t_bc;
      solution.block<3,3>(0,0) = rotation * R_bc.transpose();

      solutions.push_back(solution);
      break;
    }
    case KNEIP:
    {
      solutions = opengv::absolute_pose::p3p_kneip(_adapter,indices);

      //transform solution into body frame (case of single shifted cam)
      translation_t t_bc = _adapter.getCamOffset(indices[0]);
      rotation_t R_bc = _adapter.getCamRotation(indices[0]);

      for(size_t i = 0; i < solutions.size(); i++)
      {
        translation_t translation = solutions[i].col(3);
        rotation_t rotation = solutions[i].block<3,3>(0,0);
        solutions[i].col(3) = translation - rotation * R_bc.transpose() * t_bc;
        solutions[i].block<3,3>(0,0) = rotation * R_bc.transpose();
      }

      break;
    }
    case GAO:
    {
      solutions = opengv::absolute_pose::p3p_gao(_adapter,indices);

      //transform solution into body frame (case of single shifted cam)
      translation_t t_bc = _adapter.getCamOffset(indices[0]);
      rotation_t R_bc = _adapter.getCamRotation(indices[0]);

      for(size_t i = 0; i < solutions.size(); i++)
      {
        translation_t translation = solutions[i].col(3);
        rotation_t rotation = solutions[i].block<3,3>(0,0);
        solutions[i].col(3) = translation - rotation * R_bc.transpose() * t_bc;
        solutions[i].block<3,3>(0,0) = rotation * R_bc.transpose();
      }

      break;
    }
    case EPNP:
    {
      transformation_t solution = opengv::absolute_pose::epnp(_adapter,indices);

      //transform solution into body frame (case of single shifted cam)
      translation_t t_bc = _adapter.getCamOffset(indices[0]);
      rotation_t R_bc = _adapter.getCamRotation(indices[0]);

      translation_t translation = solution.col(3);
      rotation_t rotation = solution.block<3,3>(0,0);
      solution.col(3) = translation - rotation * R_bc.transpose() * t_bc;
      solution.block<3,3>(0,0) = rotation * R_bc.transpose();

      solutions.push_back(solution);
      break;
    }
    case GP3P_KNEIP:
    {
      solutions = opengv::absolute_pose::gp3p(_adapter,indices);
      break;
    }
    case GP3P_LEE:
    {
      solutions = opengv::absolute_pose::gp3p_lee(_adapter,indices);
      break;
    }
    case GP3P_KUKELOVA:
    {
      solutions = opengv::absolute_pose::gp3p_kukelova(_adapter,indices);
      break;
    }
    case GPNP:
    {
      transformation_t solution = opengv::absolute_pose::gpnp(_adapter,indices);
      solutions.push_back(solution);
      break;
    }
    case UPNP:
    {
      solutions = opengv::absolute_pose::upnp(_adapter,indices);
      break;
    }
  }

  if( solutions.size() == 1 )
  {
    outModel = solutions[0];
    return true;
  }

  //now compute reprojection error of fourth point, in order to find the right one
  double minScore = 1000000.0;
  int minIndex = -1;
  for(size_t i = 0; i < solutions.size(); i++)
  {
    //compute inverse transformation
    model_t inverseSolution;
    inverseSolution.block<3,3>(0,0) = solutions[i].block<3,3>(0,0).transpose();
    inverseSolution.col(3) = -inverseSolution.block<3,3>(0,0)*solutions[i].col(3);

    //get the fourth point in homogeneous form
    Eigen::Matrix<double,4,1> p_hom;
    p_hom.block<3,1>(0,0) = _adapter.getPoint(indices[3]);
    p_hom[3] = 1.0;

    //compute the reprojection (this is working for both central and
    //non-central case)
    point_t bodyReprojection = inverseSolution * p_hom;
    point_t reprojection =
        _adapter.getCamRotation(indices[3]).transpose() *
        (bodyReprojection - _adapter.getCamOffset(indices[3]));
    reprojection = reprojection / reprojection.norm();

    //compute the score
    double score =
        1.0 - (reprojection.transpose() * _adapter.getBearingVector(indices[3]));

    //check for best solution
    if( score < minScore )
    {
      minScore = score;
      minIndex = i;
    }
  }

  if(minIndex == -1)
    return false;
  outModel = solutions[minIndex];
  return true;
}

bool
opengv::sac_problems::
    absolute_pose::AbsolutePoseSacProblem::computeModelCoefficients(
    const std::vector<int> &indices,
    model_t & outModel) const
{
  algorithm_t algorithm = _algorithm;
  return computeModelCoefficients(indices, algorithm, outModel);
}

bool
opengv::sac_problems::
absolute_pose::AbsolutePoseSacProblem::computeLoModelCoefficients(
    const std::vector<int> &indices,
    model_t & outModel) const
{
  algorithm_t algorithm = _lo_algorithm;
  return computeModelCoefficients(indices, algorithm, outModel);
}

void
opengv::sac_problems::
    absolute_pose::AbsolutePoseSacProblem::getSelectedDistancesToModel(
    const model_t & model,
    const std::vector<int> & indices,
    std::vector<double> & scores) const
{
  //compute the reprojection error of all points

  //compute inverse transformation
  model_t inverseSolution;
  inverseSolution.block<3,3>(0,0) = model.block<3,3>(0,0).transpose();
  inverseSolution.col(3) = -inverseSolution.block<3,3>(0,0)*model.col(3);

  Eigen::Matrix<double,4,1> p_hom;
  p_hom[3] = 1.0;

  for(size_t i = 0; i < indices.size(); i++)
  {
    //get point in homogeneous form
    p_hom.block<3,1>(0,0) = _adapter.getPoint(indices[i]);

    //compute the reprojection (this is working for both central and
    //non-central case)
    point_t bodyReprojection = inverseSolution * p_hom;
    point_t reprojection =
        _adapter.getCamRotation(indices[i]).transpose() *
        (bodyReprojection - _adapter.getCamOffset(indices[i]));
    reprojection = reprojection / reprojection.norm();

    //compute the score
    scores.push_back(
        1.0 - (reprojection.transpose() * _adapter.getBearingVector(indices[i])));
  }
}

void
opengv::sac_problems::
    absolute_pose::AbsolutePoseSacProblem::optimizeModelCoefficients(
    const std::vector<int> & inliers,
    const model_t & model,
    model_t & optimized_model)
{
  _adapter.sett(model.col(3));
  _adapter.setR(model.block<3,3>(0,0));
  optimized_model = opengv::absolute_pose::optimize_nonlinear(_adapter,inliers);
}

int
opengv::sac_problems::
    absolute_pose::AbsolutePoseSacProblem::getSampleSize() const
{
  int sampleSize = 4;
  if(_algorithm == TWOPT)
    sampleSize = 2;
  if(_algorithm == EPNP || _algorithm == GPNP)
    sampleSize = 6;
  return sampleSize;
}

int
opengv::sac_problems::
absolute_pose::AbsolutePoseSacProblem::getLoSampleSize() const
{
  int sampleSize = 4;
  if(_lo_algorithm == TWOPT)
    sampleSize = 2;
  if(_lo_algorithm == EPNP || _lo_algorithm == GPNP)
    sampleSize = 6;
  return sampleSize;
}
