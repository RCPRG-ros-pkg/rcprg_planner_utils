cmake_minimum_required(VERSION 2.8)
project(planer_utils)

find_package(catkin REQUIRED COMPONENTS
  std_msgs tf eigen_conversions cmake_modules collision_convex_model kin_dyn_model
  )

find_package(Eigen REQUIRED)
find_package(OMPL REQUIRED)

# Export package information (replaces catkin_package() macro) 
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES planer_utils
  DEPENDS roscpp geometry_msgs visualization_msgs message_runtime Eigen 
  collision_convex_model
)

link_directories(${OMPL_LIBRARY_DIRS})

include_directories(${catkin_INCLUDE_DIRS} ${Eigen_INCLUDE_DIRS} ${OMPL_INCLUDE_DIRS})

include_directories(include)

add_library(planer_utils src/marker_publisher.cpp src/random_uniform.cpp src/reachability_map.cpp src/task_col.cpp src/task_hand.cpp src/task_wcc.cpp src/task_jlc.cpp src/task_joint.cpp src/utilities.cpp src/rrt_star.cpp src/simulator.cpp src/activation_function.cpp src/velma_q5q6_collision.cpp)

target_link_libraries(planer_utils ${catkin_LIBRARIES} ${Boost_LIBRARIES})
add_dependencies(planer_utils ${PROJECT_NAME}_gencpp)

### Orocos Package Exports and Install Targets ###

