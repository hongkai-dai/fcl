/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2018. Toyota Research Institute
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of CNRS-LAAS and AIST nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/** @author Hongkai Dai (hongkai.dai@tri.global) */

/** Tests the Expanded Polytope Algorithm (EPA) implementation inside FCL.
 * This test captures the test reported by Wolfgang Merkt in github issue
 * https://github.com/flexible-collision-library/fcl/issues/319
 */

#include "fcl/narrowphase/detail/convexity_based_algorithm/gjk_libccd-inl.h"

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include "fcl/narrowphase/detail/gjk_solver_libccd.h"

namespace fcl {
namespace detail {
template <typename S>
void TestBoxCylinderDistance() {
  fcl::Box<S> box(0.05, 0.39, 1.47);
  fcl::Cylinder<S> cylinder(0.18, 0.18);
  fcl::Transform3<S> tf1, tf2;
  tf1.setIdentity();
  tf2.setIdentity();
  tf1.translation() << 0.8, -0.385, 0.735;
  tf1.linear() = Eigen::Quaternion<S>(std::sqrt(2) / 2, 0, 0, std::sqrt(2) / 2)
                     .toRotationMatrix();
  tf2.translation() << 0.4, -0.4, 0.93;
  tf2.linear() =
      Eigen::Quaternion<S>(0.991445, 0, 0, -0.130526).toRotationMatrix();
  void* o1 = GJKInitializer<S, fcl::Box<S>>::createGJKObject(box, tf1);
  void* o2 =
      GJKInitializer<S, fcl::Cylinder<S>>::createGJKObject(cylinder, tf2);

  S dist;
  Vector3<S> p_FN1, p_FN2;
  GJKSolver_libccd<S> gjkSolver;
  gjkSolver.distance_tolerance = 1E-6;
  bool res = GJKSignedDistance(
      o1, detail::GJKInitializer<S, Box<S>>::getSupportFunction(), o2,
      detail::GJKInitializer<S, Cylinder<S>>::getSupportFunction(),
      gjkSolver.max_distance_iterations, gjkSolver.distance_tolerance, &dist,
      &p_FN1, &p_FN2);
  std::cout << "distance: " << dist << "\n";

  EXPECT_TRUE(res);
}

GTEST_TEST(FCL_GJK_EPA, box_cylinder) { TestBoxCylinderDistance<double>(); }
}  // namespace detail
}  // namespace fcl

//==============================================================================
int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
