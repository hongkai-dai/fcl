#include <gtest/gtest.h>

#include "eigen_matrix_compare.h"
#include "fcl/narrowphase/detail/gjk_solver_libccd.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"
#include "fcl/narrowphase/distance.h"
#include "fcl_resources/config.h"
#include "test_fcl_utility.h"

using namespace fcl;

bool verbose = false;

GTEST_TEST(FCL_BOX_BOX, EPA) {
  fcl::Box<double> box1{0.05, 0.02, 0.02};
  fcl::Box<double> box2{0.49, 0.05, 0.21};
  fcl::Transform3<double> tf1{};
  tf1.translation() << -0.4042785610953507, -0.5755893081793902,
      0.4672156674002889;
  tf1.linear() =
      fcl::Transform3<double>(
          Eigen::Quaternion<double>(0.7548086743258852, -0.1998958435501273,
                                    0.09822090472099809, 0.6169750163252735))
          .linear();
  fcl::Transform3<double> tf2{};
  tf2.translation() << -0.4399999999999999, -0.61499999999858, 0.355;
  tf2.linear() =
      fcl::Transform3<double>(Eigen::Quaternion<double>(0.7071067811882787, 0,
                                                        0, 0.7071067811848163))
          .linear();

  fcl::detail::GJKSolver_libccd<double> solver;
  double distance;
  Eigen::Vector3d p1, p2;
  solver.shapeSignedDistance(box1, tf1, box2, tf2, &distance, &p1, &p2);
}

//==============================================================================
int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
