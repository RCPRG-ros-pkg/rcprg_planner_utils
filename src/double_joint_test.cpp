#include <iostream>

#include <ros/ros.h>
#include <visualization_msgs/MarkerArray.h>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include "planer_utils/activation_function.h"
#include "planer_utils/double_joint_collision_checker.h"

using namespace cv;

int main (int argc, char** argv) {
    ros::init(argc, argv, "double_joint_test");
    ros::NodeHandle m_nh("~");

    ros::Publisher markerPub = m_nh.advertise<visualization_msgs::MarkerArray>("kernel", 10, true);
    ros::Duration(0.5).sleep();

    const double d0 = 0.25;

    const std::vector<double > constraint_polygon({
            0.397855401039 , 2.90307354927 ,
            -1.79 , 2.91307354927 ,
            -1.78 , 1.43123674393 ,
            -0.77621114254 , 1.39720571041 ,
            -0.36 , 1.00585031509 ,
            -0.35 , -0.414940297604 ,
            -0.8 , -0.942419290543 ,
            -1.8 , -1.01898884773 ,
            -1.81 , -2.88 ,
            0.4 , -2.89 ,
            0.81 , -2.27813267708 ,
            1.82 , -2.29514837265 ,
            1.83 , 1.66945314407 ,
            0.84 , 1.73751521111 ,
            0.423378348351 , 2.09483933449});

    DoubleJointCC cc(d0, constraint_polygon);

    const int sx = 1000;
    const int sy = 1000;
    Mat image = Mat::zeros( sx, sy, CV_8UC3 );

    double range = 3.2;

    // Draw the border (red)
    for (int i = 0; i < constraint_polygon.size()/2; ++i) {
        double x1 = constraint_polygon[i*2];
        double y1 = constraint_polygon[i*2+1];
        double x2 = constraint_polygon[(i*2+2)%constraint_polygon.size()];
        double y2 = constraint_polygon[(i*2+3)%constraint_polygon.size()];
        int ix1 = int( sx * (x1+range)/(range*2) );
        int iy1 = int( sy * (y1+range)/(range*2) );
        Point pt1(ix1, iy1);
        int ix2 = int( sx * (x2+range)/(range*2) );
        int iy2 = int( sy * (y2+range)/(range*2) );
        Point pt2(ix2, iy2);
        line( image, pt1, pt2, Scalar( 0, 0, 255 ), 2, LINE_8 );
    }

    // Draw the field
    for (double x = -range; x < range; x += 0.1) {
        for (double y = -range; y < range; y += 0.1) {
            DoubleJointCC::Joints q2(x, y);

            if (cc.inCollision(q2)) {
            }

            int min_idx;
            int min_type;
            double min_dist;
            DoubleJointCC::Joints min_v;

            bool found = cc.getMinDistanceIn(q2, min_v, min_dist, min_idx, min_type);

            if (found) {
                int ix = int( sx * (x+range)/(range*2) );
                int iy = int( sy * (y+range)/(range*2) );
                Point pt1(ix, iy);
                //min_v.normalize();
                double mult = 50.0;
                int ix2 = ix + int( min_v(0) * mult );
                int iy2 = iy + int( min_v(1) * mult );
                Point pt2(ix2, iy2);
                circle( image, pt1, 2, Scalar( 0, 255, 0 ), FILLED, LINE_8 );
                line( image, pt1, pt2, Scalar( 0, 128, 0 ), 1, LINE_8 );
            }

            found = cc.getMinDistanceOut(q2, min_v, min_dist, min_idx, min_type);

            if (found) {
                int ix = int( sx * (x+range)/(range*2) );
                int iy = int( sy * (y+range)/(range*2) );
                Point pt1(ix, iy);
                //min_v.normalize();
                double mult = 50.0;
                int ix2 = ix + int( min_v(0) * mult );
                int iy2 = iy + int( min_v(1) * mult );
                Point pt2(ix2, iy2);
                circle( image, pt1, 2, Scalar( 255, 0, 0 ), FILLED, LINE_8 );
                line( image, pt1, pt2, Scalar( 128, 0, 0 ), 1, LINE_8 );
            }
        }
    }
    char title_window[] = "double_joint_test";
    imshow( title_window, image );
    moveWindow( title_window, 0, 200 );
    waitKey( 0 );

    return 0;
}
