cmake_minimum_required(VERSION 2.8)
project(planer_utils)

find_package(catkin REQUIRED COMPONENTS
  std_msgs tf eigen_conversions cmake_modules collision_convex_model kin_dyn_model
  )

find_package(Eigen3 REQUIRED)
find_package(OMPL REQUIRED)

# Export package information (replaces catkin_package() macro) 
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES planer_utils
  DEPENDS roscpp geometry_msgs visualization_msgs message_runtime Eigen3
  collision_convex_model
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

link_directories(${OMPL_LIBRARY_DIRS})

include_directories(${catkin_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS} ${OMPL_INCLUDE_DIRS})

include_directories(include)

add_library(planer_utils
    src/marker_publisher.cpp
    src/random_uniform.cpp
    src/reachability_map.cpp
    src/task_col.cpp
    src/task_hand.cpp
    src/task_wcc.cpp
    src/task_jlc.cpp
    src/task_joint.cpp
    src/utilities.cpp
    src/rrt_star.cpp
    src/simulator.cpp
    src/activation_function.cpp
    src/double_joint_collision_checker.cpp)

target_link_libraries(planer_utils ${catkin_LIBRARIES} ${Boost_LIBRARIES})

### Orocos Package Exports and Install Targets ###

