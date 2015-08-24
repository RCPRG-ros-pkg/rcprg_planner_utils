cmake_minimum_required(VERSION 2.8)
project(planer_utils)

find_package(catkin REQUIRED COMPONENTS
  rtt_ros rtt_rosclock std_msgs tf eigen_conversions cmake_modules collision_convex_model kin_dyn_model
  )

find_package(Eigen REQUIRED)
find_package(orocos_kdl REQUIRED)
find_package(OMPL REQUIRED)

link_directories(${orocos_kdl_LIBRARY_DIRS} ${OMPL_LIBRARY_DIRS})

include_directories(${catkin_INCLUDE_DIRS} ${Eigen_INCLUDE_DIRS} ${USE_OROCOS_INCLUDE_DIRS} ${OMPL_INCLUDE_DIRS})

include_directories(include)

add_library(planer_utils src/marker_publisher.cpp src/random_uniform.cpp src/reachability_map.cpp src/task_col.cpp src/task_hand.cpp src/task_wcc.cpp src/task_jlc.cpp src/utilities.cpp src/rrt_star.cpp src/simulator.cpp)

target_link_libraries(planer_utils ${catkin_LIBRARIES} ${USE_OROCOS_LIBRARIES} ${orocos_kdl_LIBRARIES})
add_dependencies(planer_utils ${PROJECT_NAME}_gencpp)

### Orocos Package Exports and Install Targets ###

# Export package information (replaces catkin_package() macro) 
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES planer_utils
  DEPENDS roscpp geometry_msgs visualization_msgs message_runtime Eigen orocos_kdl
)

