cmake_minimum_required(VERSION 2.8)
project(planer_utils)

find_package(catkin REQUIRED COMPONENTS
  std_msgs tf eigen_conversions cmake_modules collision_convex_model kin_dyn_model rcprg_ros_utils
  )

find_package(Eigen3 REQUIRED)
find_package(ompl REQUIRED)

# Export package information (replaces catkin_package() macro) 
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES planer_utils
  CATKIN_DEPENDS roscpp geometry_msgs visualization_msgs message_runtime
  collision_convex_model rcprg_ros_utils kin_dyn_model
  DEPENDS EIGEN3
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

link_directories(${OMPL_LIBRARY_DIRS})

include_directories(${catkin_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS} ${OMPL_INCLUDE_DIRS})

include_directories(include)

add_library(planer_utils
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

install(TARGETS planer_utils
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)

