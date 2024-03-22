cmake_minimum_required( VERSION 3.23.1 )

project(opengv)

# Define relative dir for ${CMAKE_SOURCE_DIR}
set(PROJECT_SRC_DIR "opengv")

set(EXTERNAL_PROJECT_DIR ${PROJECT_SOURCE_DIR}/external/project)
list( APPEND CMAKE_MODULE_PATH "${EXTERNAL_PROJECT_DIR}/cmake" )

option(STATIC_LINK_PLUGINS "Static link plugins. Just for release" OFF)

set(ENABLE_REL_WITH_DEB_INFO OFF)
# After set environment variable, we includ project cmake.
include(Project)

# Set git hook
execute_process(COMMAND git config core.hooksPath external/project/githook)
