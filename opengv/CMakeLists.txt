cmake_minimum_required( VERSION 3.23.1 )

# Versioning
set(OPENGV_VERSION_MAJOR 1)
set(OPENGV_VERSION_MINOR 0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)
include_directories(
        ${Eigen3_INCLUDE_DIRS}
        ${Eigen3_INCLUDE_DIRS}/eigen3/unsupported)

find_package(glog REQUIRED)
set(GOOGLE_LOG glog::glog)
include_directories(${glog_INCLUDE_DIRS})

set(WARNING_FLAGS "-Wall -Wextra -Werror -Wno-sign-compare -Wno-reorder -Wno-deprecated-declarations -Wno-unused-function -Wno-unused-variable -Wno-unused-parameter -Wno-format-security")
if(APPLE)
    set(WARNING_FLAGS "${WARNING_FLAGS} -Wno-int-in-bool-context -Wno-deprecated-copy -Wno-deprecated-copy-with-user-provided-copy -Wno-non-c-typedef-for-linkage -Wno-unused-but-set-variable -Wno-unused-value -Wno-unused-private-field -Wno-shorten-64-to-32")
elseif(WIN32)
    set(WARNING_FLAGS "${WARNING_FLAGS} -Wno-int-in-bool-context -Wno-deprecated-copy -Wno-unused-but-set-variable -Wno-unused-value -Wno-non-c-typedef-for-linkage -Wno-attributes -Wno-parentheses -Wno-cast-function-type -Wno-nonnull -Wno-array-bounds -Wno-implicit-const-int-float-conversion -Wno-unused-private-field")
elseif(ANDROID)
    set(WARNING_FLAGS "${WARNING_FLAGS} -Wno-unused-private-field -Wno-deprecated-copy -Wno-deprecated-copy-with-user-provided-copy -Wno-non-c-typedef-for-linkage -Wno-unused-but-set-variable -Wno-implicit-const-int-float-conversion -Wno-int-in-bool-context")
elseif(UNIX)
    set(WARNING_FLAGS "${WARNING_FLAGS} -Wno-int-in-bool-context -Wno-unused-but-set-variable -Wno-unused-value -Wno-maybe-uninitialized -Wno-implicit-fallthrough")
endif()

set(COMPILATION_FLAGS "${WARNING_FLAGS} -O3 -fPIC -fno-rtti")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILATION_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMPILATION_FLAGS}")

# Source and header files
set(OPENGV_SOURCE_FILES
  src/absolute_pose/modules/main.cpp
  src/absolute_pose/modules/gp3p/code.cpp
  src/absolute_pose/modules/gp3p/init.cpp
  src/absolute_pose/modules/gp3p/reductors.cpp
  src/absolute_pose/modules/gp3p/spolynomials.cpp
  src/absolute_pose/modules/gp3p_lee/code.cpp
  src/absolute_pose/modules/gp3p_lee/generalized_absolute_pose_coeffs.cc
  src/absolute_pose/modules/gp3p_kukelova/re3q3.cc
  src/absolute_pose/modules/Epnp.cpp
  src/absolute_pose/modules/gpnp1/code.cpp
  src/absolute_pose/modules/gpnp1/init.cpp
  src/absolute_pose/modules/gpnp1/reductors.cpp
  src/absolute_pose/modules/gpnp1/spolynomials.cpp
  src/absolute_pose/modules/gpnp2/code.cpp
  src/absolute_pose/modules/gpnp2/init.cpp
  src/absolute_pose/modules/gpnp2/reductors.cpp
  src/absolute_pose/modules/gpnp2/spolynomials.cpp
  src/absolute_pose/modules/gpnp3/code.cpp
  src/absolute_pose/modules/gpnp3/init.cpp
  src/absolute_pose/modules/gpnp3/reductors.cpp
  src/absolute_pose/modules/gpnp3/spolynomials.cpp
  src/absolute_pose/modules/gpnp4/code.cpp
  src/absolute_pose/modules/gpnp4/init.cpp
  src/absolute_pose/modules/gpnp4/reductors.cpp
  src/absolute_pose/modules/gpnp4/spolynomials.cpp
  src/absolute_pose/modules/gpnp5/code.cpp
  src/absolute_pose/modules/gpnp5/init.cpp
  src/absolute_pose/modules/gpnp5/reductors.cpp
  src/absolute_pose/modules/gpnp5/spolynomials.cpp
  src/absolute_pose/modules/upnp2.cpp
  src/absolute_pose/modules/upnp4.cpp
  src/relative_pose/modules/main.cpp
  src/relative_pose/modules/fivept_nister/modules.cpp
  src/relative_pose/modules/fivept_stewenius/modules.cpp
  src/relative_pose/modules/fivept_kneip/code.cpp
  src/relative_pose/modules/fivept_kneip/init.cpp
  src/relative_pose/modules/fivept_kneip/reductors.cpp
  src/relative_pose/modules/fivept_kneip/spolynomials.cpp
  src/relative_pose/modules/sixpt/modules2.cpp
  src/relative_pose/modules/eigensolver/modules.cpp
  src/relative_pose/modules/ge/modules.cpp
  src/math/cayley.cpp
  src/math/quaternion.cpp
  src/math/arun.cpp
  src/math/Sturm.cpp
  src/math/roots.cpp
  src/math/gauss_jordan.cpp
  src/math/polynomial.cc
  src/absolute_pose/methods.cpp
  src/absolute_pose/CentralAbsoluteAdapter.cpp
  src/absolute_pose/NoncentralAbsoluteAdapter.cpp
  src/absolute_pose/NoncentralAbsoluteMultiAdapter.cpp
  src/relative_pose/methods.cpp
  src/relative_pose/CentralRelativeAdapter.cpp
  src/relative_pose/CentralRelativeWeightingAdapter.cpp
  src/relative_pose/NoncentralRelativeAdapter.cpp
  src/relative_pose/CentralRelativeMultiAdapter.cpp
  src/relative_pose/NoncentralRelativeMultiAdapter.cpp
  src/sac/sprt.cpp
  src/triangulation/methods.cpp
  src/point_cloud/methods.cpp
  src/point_cloud/PointCloudAdapter.cpp
  src/sac_problems/absolute_pose/AbsolutePoseSacProblem.cpp
  src/sac_problems/absolute_pose/MultiNoncentralAbsolutePoseSacProblem.cpp
  src/sac_problems/relative_pose/CentralRelativePoseSacProblem.cpp
  src/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.cpp
  src/sac_problems/relative_pose/RotationOnlySacProblem.cpp
  src/sac_problems/relative_pose/TranslationOnlySacProblem.cpp
  src/sac_problems/relative_pose/EigensolverSacProblem.cpp
  src/sac_problems/relative_pose/MultiCentralRelativePoseSacProblem.cpp
  src/sac_problems/relative_pose/MultiNoncentralRelativePoseSacProblem.cpp
  src/sac_problems/point_cloud/PointCloudSacProblem.cpp
  src/absolute_pose/MACentralAbsolute.cpp
  src/absolute_pose/MANoncentralAbsolute.cpp
  src/relative_pose/MACentralRelative.cpp
  src/relative_pose/MANoncentralRelative.cpp
  src/relative_pose/MANoncentralRelativeMulti.cpp
  src/point_cloud/MAPointCloud.cpp )

set(OPENGV_HEADER_FILES
  include/opengv/types.hpp
  include/opengv/OptimizationFunctor.hpp
  include/opengv/absolute_pose/methods.hpp
  include/opengv/relative_pose/methods.hpp
  include/opengv/triangulation/methods.hpp
  include/opengv/point_cloud/methods.hpp
  include/opengv/math/cayley.hpp
  include/opengv/math/quaternion.hpp
  include/opengv/math/arun.hpp
  include/opengv/math/Sturm.hpp
  include/opengv/math/roots.hpp
  include/opengv/math/gauss_jordan.hpp
  include/opengv/absolute_pose/AbsoluteAdapterBase.hpp
  include/opengv/absolute_pose/CentralAbsoluteAdapter.hpp
  include/opengv/absolute_pose/NoncentralAbsoluteAdapter.hpp
  include/opengv/absolute_pose/NoncentralAbsoluteMultiAdapter.hpp
  include/opengv/absolute_pose/AbsoluteMultiAdapterBase.hpp
  include/opengv/relative_pose/RelativeAdapterBase.hpp
  include/opengv/relative_pose/RelativeMultiAdapterBase.hpp
  include/opengv/relative_pose/CentralRelativeAdapter.hpp
  include/opengv/relative_pose/CentralRelativeWeightingAdapter.hpp
  include/opengv/relative_pose/NoncentralRelativeAdapter.hpp
  include/opengv/relative_pose/CentralRelativeMultiAdapter.hpp
  include/opengv/relative_pose/NoncentralRelativeMultiAdapter.hpp
  include/opengv/point_cloud/PointCloudAdapterBase.hpp
  include/opengv/point_cloud/PointCloudAdapter.hpp
  include/opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp
  include/opengv/sac_problems/absolute_pose/MultiNoncentralAbsolutePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/CentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/MultiCentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/MultiNoncentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/EigensolverSacProblem.hpp
  include/opengv/sac_problems/relative_pose/RotationOnlySacProblem.hpp
  include/opengv/sac_problems/relative_pose/TranslationOnlySacProblem.hpp
  include/opengv/sac_problems/point_cloud/PointCloudSacProblem.hpp
  include/opengv/absolute_pose/MACentralAbsolute.hpp
  include/opengv/absolute_pose/MANoncentralAbsolute.hpp
  include/opengv/relative_pose/MACentralRelative.hpp
  include/opengv/relative_pose/MANoncentralRelative.hpp
  include/opengv/relative_pose/MANoncentralRelativeMulti.hpp
  include/opengv/point_cloud/MAPointCloud.hpp )

# Define the library and executables
add_library(opengv STATIC ${OPENGV_SOURCE_FILES} ${OPENGV_HEADER_FILES})
target_link_libraries(opengv glog::glog)
target_include_directories(opengv PUBLIC include)

install(DIRECTORY include/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(TARGETS opengv DESTINATION "${CMAKE_INSTALL_LIBDIR}")

#add_library( random_generators SHARED test/random_generators.cpp test/random_generators.hpp test/experiment_helpers.cpp test/experiment_helpers.hpp test/time_measurement.cpp test/time_measurement.hpp )
#target_link_libraries(random_generators ${PROJECT_NAME})
#
#add_executable( test_absolute_pose test/test_absolute_pose.cpp )
#target_link_libraries( test_absolute_pose ${PROJECT_NAME} random_generators )
#
#add_executable( test_absolute_pose_sac test/test_absolute_pose_sac.cpp )
#target_link_libraries( test_absolute_pose_sac ${PROJECT_NAME} random_generators )
#
#add_executable( test_noncentral_absolute_pose test/test_noncentral_absolute_pose.cpp )
#target_link_libraries( test_noncentral_absolute_pose ${PROJECT_NAME} random_generators )
#
#add_executable( test_noncentral_absolute_pose_sac test/test_noncentral_absolute_pose_sac.cpp )
#target_link_libraries( test_noncentral_absolute_pose_sac ${PROJECT_NAME} random_generators )
#
#add_executable( test_relative_pose test/test_relative_pose.cpp )
#target_link_libraries( test_relative_pose ${PROJECT_NAME} random_generators )
#
#add_executable( test_relative_pose_rotationOnly test/test_relative_pose_rotationOnly.cpp )
#target_link_libraries( test_relative_pose_rotationOnly ${PROJECT_NAME} random_generators )
#
#add_executable( test_relative_pose_rotationOnly_sac test/test_relative_pose_rotationOnly_sac.cpp )
#target_link_libraries( test_relative_pose_rotationOnly_sac ${PROJECT_NAME} random_generators )
#
#add_executable( test_relative_pose_sac test/test_relative_pose_sac.cpp )
#target_link_libraries( test_relative_pose_sac ${PROJECT_NAME} random_generators )
#
#add_executable( test_noncentral_relative_pose test/test_noncentral_relative_pose.cpp )
#target_link_libraries( test_noncentral_relative_pose ${PROJECT_NAME} random_generators )
#
#add_executable( test_noncentral_relative_pose_sac test/test_noncentral_relative_pose_sac.cpp )
#target_link_libraries( test_noncentral_relative_pose_sac ${PROJECT_NAME} random_generators )
#
#add_executable( test_multi_noncentral_relative_pose_sac test/test_multi_noncentral_relative_pose_sac.cpp )
#target_link_libraries( test_multi_noncentral_relative_pose_sac ${PROJECT_NAME} random_generators )
#
#add_executable( test_eigensolver_sac test/test_eigensolver_sac.cpp )
#target_link_libraries( test_eigensolver_sac ${PROJECT_NAME} random_generators )
#
#add_executable( test_triangulation test/test_triangulation.cpp )
#target_link_libraries( test_triangulation ${PROJECT_NAME} random_generators )
#
#add_executable( test_eigensolver test/test_eigensolver.cpp )
#target_link_libraries( test_eigensolver ${PROJECT_NAME} random_generators )
#
#add_executable( test_point_cloud test/test_point_cloud.cpp )
#target_link_libraries( test_point_cloud ${PROJECT_NAME} random_generators )
#
#add_executable( test_point_cloud_sac test/test_point_cloud_sac.cpp )
#target_link_libraries( test_point_cloud_sac ${PROJECT_NAME} random_generators )

