// Copyright (c) 2015, Robot Control and Pattern Recognition Group,
// Institute of Control and Computation Engineering
// Warsaw University of Technology
//
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Warsaw University of Technology nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYright HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Author: Dawid Seredynski
//

#include "planer_utils/utilities.h"
#include "planer_utils/random_uniform.h"
#include <boost/math/special_functions/bessel.hpp>
#include <sstream>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

static const double PI(3.141592653589793);

    void publishJointState(ros::Publisher &joint_state_pub, const Eigen::VectorXd &q, const std::vector<std::string > &joint_names) {
        sensor_msgs::JointState js;
        js.header.stamp = ros::Time::now();
        int q_idx = 0;
        for (std::vector<std::string >::const_iterator it = joint_names.begin(); it != joint_names.end(); it++, q_idx++) {
            js.name.push_back(*it);
            js.position.push_back(q[q_idx]);
        }
        joint_state_pub.publish(js);
    }

    void publishJointState(ros::Publisher &joint_state_pub, const Eigen::VectorXd &q, const std::vector<std::string > &joint_names, const Eigen::VectorXd &ign_q, const std::vector<std::string > &ign_joint_names) {
        sensor_msgs::JointState js;
        js.header.stamp = ros::Time::now();
        int q_idx = 0;
        for (std::vector<std::string >::const_iterator it = joint_names.begin(); it != joint_names.end(); it++, q_idx++) {
            js.name.push_back(*it);
            js.position.push_back(q[q_idx]);
        }
        q_idx = 0;
        for (std::vector<std::string >::const_iterator it = ign_joint_names.begin(); it != ign_joint_names.end(); it++, q_idx++) {
            js.name.push_back(*it);
            js.position.push_back(ign_q[q_idx]);
        }
        joint_state_pub.publish(js);
    }


    void publishTransform(tf::TransformBroadcaster &br, const KDL::Frame &T_B_F, const std::string &frame_id, const std::string &base_frame_id) {
        tf::Transform transform;
        transform.setOrigin( tf::Vector3(T_B_F.p.x(), T_B_F.p.y(), T_B_F.p.z()) );
        tf::Quaternion q;
        double qx, qy, qz, qw;
        T_B_F.M.GetQuaternion(q[0], q[1], q[2], q[3]);
        transform.setRotation(q);
        br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), base_frame_id, frame_id));
    }

    int addRobotModelVis(MarkerPublisher &markers_pub, int m_id, const boost::shared_ptr<self_collision::CollisionModel> &col_model, const std::vector<KDL::Frame > &T) {
        for (self_collision::CollisionModel::VecPtrLink::const_iterator l_it = col_model->getLinks().begin(); l_it != col_model->getLinks().end(); l_it++) {
            KDL::Frame T_B_L = T[(*l_it)->index_];
            for (self_collision::Link::VecPtrCollision::const_iterator it = (*l_it)->collision_array.begin(); it != (*l_it)->collision_array.end(); it++) {
                KDL::Frame T_B_O = T_B_L * (*it)->origin;
                double cr, cg, cb, ca;
                (*it)->geometry->getColor(cr, cg, cb, ca);
                if ((*it)->geometry->getType() == self_collision::Geometry::CONVEX) {
                    self_collision::Convex *convex = static_cast<self_collision::Convex* >((*it)->geometry.get());

                    if (convex->visualisation_hint_ == "lines") {
                        const std::vector<KDL::Vector >& points = convex->getPoints();
                        const std::vector<int >& polygons = convex->getPolygons();
                        std::vector<KDL::Vector > pts;

	                    int poly_idx = 0;
	                    while (poly_idx < polygons.size())
	                    {
		                    int points_in_poly = polygons[poly_idx];
		                    for (int j=0; j<points_in_poly; j++)
		                    {
                                int idx1 = poly_idx + 1 + j;
                                int idx2 = poly_idx + 1 + (j+1)%points_in_poly;
			                    pts.push_back( points[polygons[idx1]] );
                                pts.push_back( points[polygons[idx2]] );
		                    }
		                    poly_idx += points_in_poly+1;
	                    }
                        m_id = markers_pub.addLineListMarker(m_id, pts, T_B_O, cr, cg, cb, ca, 0.02, "world");
                    }
                    else if (convex->visualisation_hint_.find("box") == 0) {
                        double dx, dy, dz;
                        std::istringstream(convex->visualisation_hint_.substr(4)) >> dx >> dy >> dz;
                        m_id = markers_pub.addSinglePointMarkerCube(m_id, T_B_O.p, cr, cg, cb, ca, dx, dy, dz, "world");
                    }
                }
                else if ((*it)->geometry->getType() == self_collision::Geometry::SPHERE) {
                    self_collision::Sphere *sphere = static_cast<self_collision::Sphere* >((*it)->geometry.get());
                    m_id = markers_pub.addSinglePointMarker(m_id, T_B_O.p, cr, cg, cb, ca, sphere->getRadius()*2, "world");
                }
                else if ((*it)->geometry->getType() == self_collision::Geometry::CAPSULE) {
                    self_collision::Capsule *capsule = static_cast<self_collision::Capsule* >((*it)->geometry.get());
                    m_id = markers_pub.addCapsule(m_id, T_B_O, cr, cg, cb, ca, capsule->getLength(), capsule->getRadius(), "world");
                }
            }
        }
        return m_id;
    }

    void getPointOnPath(const std::list<Eigen::VectorXd > &path, double f, Eigen::VectorXd &x) {

        if (path.size() == 0) {
            std::cout << "ERROR: getPointOnPath: path size is 0" << std::endl;
            return;
        }
        else if (path.size() == 1 || f < 0.0) {
            x = (*path.begin());
            return;
        }

        if (f > 1.0) {
            x = (*(--path.end()));
            return;
        }

        double length = 0.0;
        for (std::list<Eigen::VectorXd >::const_iterator it1 = path.begin(), it2=++path.begin(); it2 != path.end(); it1++, it2++) {
            double dist = ((*it1) - (*it2)).norm();
            length += dist;
        }

        double pos = length * f;

        for (std::list<Eigen::VectorXd >::const_iterator it1 = path.begin(), it2=++path.begin(); it2 != path.end(); it1++, it2++) {
            Eigen::VectorXd v = ((*it2) - (*it1));
            double dist = v.norm();
            if (pos - dist > 0) {
                pos -= dist;
            }
            else {
                x = (*it1) + pos * v / dist;
                return;
            }
        }
        x = (*(--path.end()));
    }

    double getPathLength(const std::list<Eigen::VectorXd > &path) {

        if (path.size() < 2) {
            return 0.0;
        }

        double length = 0.0;
        for (std::list<Eigen::VectorXd >::const_iterator it1 = path.begin(), it2=++path.begin(); it2 != path.end(); it1++, it2++) {
            double dist = ((*it1) - (*it2)).norm();
            length += dist;
        }

        return length;
    }

void printFrameKDL(const KDL::Frame &f) {
    double qx, qy, qz, qw, px(f.p.x()), py(f.p.y()), pz(f.p.z());
    f.M.GetQuaternion(qx, qy, qz, qw);
    std::cout << "KDL::Frame(KDL::Rotation::Quaternion(" << qx << ", " << qy << ", " << qz << ", " << qw << "), KDL::Vector(" << px << ", " << py << ", " << pz << "))" << std::endl;
}

double getAngle(const KDL::Vector &v1, const KDL::Vector &v2) {
    return std::atan2((v1*v2).Norm(), KDL::dot(v1,v2));
}

void EigenTfToKDL(const Eigen::Isometry3d &tf, KDL::Frame &kdlT) {
    kdlT = KDL::Frame( KDL::Rotation(tf(0,0),tf(0,1),tf(0,2), tf(1,0), tf(1,1), tf(1,2), tf(2,0), tf(2,1), tf(2,2)), KDL::Vector(tf(0,3), tf(1,3), tf(2,3)) );
}

void KDLToEigenTf(const KDL::Frame &kdlT, Eigen::Isometry3d &tf) {
    double qx, qy, qz, qw;
    kdlT.M.GetQuaternion(qx, qy, qz, qw);
    tf.fromPositionOrientationScale(Eigen::Vector3d(kdlT.p.x(), kdlT.p.y(), kdlT.p.z()), Eigen::Quaterniond(qw, qx, qy, qz), Eigen::Vector3d(1.0, 1.0, 1.0));
}

std::ostream& operator<< (std::ostream& stream, const KDL::Frame& f) {
    double qx, qy, qz, qw;
    f.M.GetQuaternion(qx, qy, qz, qw);
    stream << f.p.x() << " " << f.p.y() << " " << f.p.z() << " " << qx << " " << qy << " " << qz << " " << qw;
    return stream;
}

std::istream& operator>> (std::istream& stream, KDL::Frame& f) {
    double x, y, z, qx, qy, qz, qw;
    Eigen::Vector4d q;
    stream >> x >> y >> z >> q(0) >> q(1) >> q(2) >> q(3);
    q.normalize();
    f = KDL::Frame(KDL::Rotation::Quaternion(q(0), q(1), q(2), q(3)), KDL::Vector(x, y, z));
    return stream;
}

std::string frameKdl2string(const KDL::Frame &f) {
    std::ostringstream strs;
    strs << f;
    return strs.str();
}

KDL::Frame string2frameKdl(const std::string &str) {
    std::istringstream strs(str);
    KDL::Frame ret;
    strs >> ret;
    return ret;
}

std::string double2string(double d) {
    std::ostringstream strs;
    strs << d;
    return strs.str();
}

std::string int2string(int i) {
    std::ostringstream strs;
    strs << i;
    return strs.str();
}

double string2double(const std::string &str) {
    std::istringstream strs(str);
    double ret;
    strs >> ret;
    return ret;
}

int string2int(const std::string &str) {
    std::istringstream strs(str);
    int ret;
    strs >> ret;
    return ret;
}

double deg2rad(double deg) {
    return deg / 180.0 * PI;
}

double rad2deg(double rad) {
    return rad / PI * 180.0;
}

/******************************************************************************************************************************************************************/

double triVariateIsotropicGaussianKernel(const Eigen::Vector3d &x, const Eigen::Vector3d &mean, double sigma) {
    static const double sqrt_2pi_pow_3 = std::sqrt(std::pow(2.0 * PI, 3.0));
    return std::exp( - ((x-mean).squaredNorm()) / (2.0 * sigma * sigma) ) / (sqrt_2pi_pow_3 * sigma);
}

double biVariateIsotropicGaussianKernel(const Eigen::Vector2d &x, const Eigen::Vector2d &mean, double sigma) {
    static const double sqrt_2pi_pow_2 = std::sqrt(std::pow(2.0 * PI, 2.0));
    return std::exp( - ((x-mean).squaredNorm()) / (2.0 * sigma * sigma) ) / (sqrt_2pi_pow_2 * sigma);
}

double uniVariateIsotropicGaussianKernel(double x, double mean, double sigma) {
    static const double sqrt_2pi = std::sqrt(2.0 * PI);
    return std::exp( - ((x-mean)*(x-mean)) / (2.0 * sigma * sigma) ) / (sqrt_2pi * sigma);
}

double uniVariateIsotropicGaussianKernelIntegral(double x1, double x2) {
    static const double x_min = -5.0;
    static const double x_max = 5.0;
    static const double cdf[201] = {0.0, 1.69666e-07, 2.91614e-07, 4.47223e-07, 6.45288e-07, 8.96763e-07, 1.21525e-06, 1.61761e-06, 2.12466e-06, 2.76202e-06, 3.56121e-06, 4.5608e-06, 5.80792e-06, 7.35999e-06, 9.28675e-06, 1.16727e-05, 1.46198e-05, 1.82511e-05, 2.27142e-05, 2.81859e-05, 3.48774e-05, 4.30403e-05, 5.2973e-05, 6.50294e-05, 7.96267e-05, 9.72565e-05, 0.000118495, 0.000144019, 0.000174614, 0.000211197, 0.000254831, 0.000306745, 0.000368356, 0.000441293, 0.000527421, 0.000628874, 0.000748078, 0.000887791, 0.00105113, 0.00124162, 0.00146321, 0.00172034, 0.00201797, 0.00236161, 0.00275738, 0.00321206, 0.00373311, 0.00432872, 0.00500787, 0.00578033, 0.00665675, 0.00764862, 0.00876834, 0.0100293, 0.0114456, 0.0130326, 0.0148063, 0.0167838, 0.018983, 0.0214226, 0.0241222, 0.0271019, 0.0303827, 0.0339859, 0.0379334, 0.0422473, 0.0469497, 0.052063, 0.057609, 0.0636095, 0.0700854, 0.0770569, 0.0845433, 0.0925624, 0.101131, 0.110263, 0.119973, 0.130269, 0.141162, 0.152656, 0.164755, 0.177458, 0.190762, 0.204661, 0.219146, 0.234203, 0.249815, 0.265964, 0.282625, 0.299772, 0.317376, 0.335402, 0.353815, 0.372577, 0.391647, 0.41098, 0.430532, 0.450256, 0.470104, 0.490026, 0.509973, 0.529896, 0.549743, 0.569467, 0.589019, 0.608353, 0.627422, 0.646184, 0.664598, 0.682624, 0.700227, 0.717374, 0.734036, 0.750184, 0.765797, 0.780854, 0.795338, 0.809238, 0.822542, 0.835245, 0.847343, 0.858837, 0.86973, 0.880027, 0.889736, 0.898869, 0.907437, 0.915456, 0.922943, 0.929914, 0.93639, 0.94239, 0.947937, 0.95305, 0.957752, 0.962066, 0.966014, 0.969617, 0.972898, 0.975877, 0.978577, 0.981016, 0.983216, 0.985193, 0.986967, 0.988554, 0.98997, 0.991231, 0.992351, 0.993343, 0.994219, 0.994992, 0.995671, 0.996266, 0.996787, 0.997242, 0.997638, 0.997982, 0.998279, 0.998536, 0.998758, 0.998948, 0.999112, 0.999251, 0.999371, 0.999472, 0.999558, 0.999631, 0.999693, 0.999745, 0.999788, 0.999825, 0.999855, 0.999881, 0.999902, 0.99992, 0.999934, 0.999947, 0.999956, 0.999965, 0.999971, 0.999977, 0.999981, 0.999985, 0.999988, 0.99999, 0.999992, 0.999994, 0.999995, 0.999996, 0.999997, 0.999997, 0.999998, 0.999998, 0.999999, 0.999999, 0.999999, 0.999999, 0.999999, 0.999999, 1.0};
    int x1_idx = int(201 * (x1-x_min) / (x_max-x_min));
    if (x1_idx < 0) {
        x1_idx = 0;
    }
    if (x1_idx >= 201) {
        x1_idx = 200;
    }

    int x2_idx = int(201 * (x2-x_min) / (x_max-x_min));
    if (x2_idx < 0) {
        x2_idx = 0;
    }
    if (x2_idx >= 201) {
        x2_idx = 200;
    }
    return cdf[x2_idx] - cdf[x1_idx];
}

double misesFisherKernelConstant(double sigma, int dimensions) {
    double kappa = sigma;
    double p_div_2 = static_cast<double >(dimensions) / 2.0;
    double Cp = std::pow(kappa, p_div_2 - 1.0) / (std::pow(2.0 * PI, p_div_2) * boost::math::cyl_bessel_i(p_div_2 - 1.0, kappa));
    return Cp;
}

double misesFisherKernel(const Eigen::Vector3d &q, const Eigen::Vector3d &mean, double sigma, double Cp) {
    double dot_prod = q.dot(mean);
    return Cp * (std::exp(sigma * dot_prod) + std::exp(-sigma * dot_prod)) / 2.0;
}

int vonMisesFisherSample(const Eigen::Vector3d &mean, double pdf_mean, double sigma, double Cp, Eigen::Vector3d &x) {
    for (int i = 0; i < 100000; i++) {
        randomUnitSphere(x);
        double pdf = misesFisherKernel(x, mean, sigma, Cp);
        double rand_pdf = randomUniform(0.0, pdf_mean);
        if (rand_pdf < pdf) {
            return i;
        }
    }
    return -1;
}

double misesFisherKernel(const Eigen::Vector4d &q, const Eigen::Vector4d &mean, double sigma, double Cp) {
    double dot_prod = q.dot(mean);
    return Cp * (std::exp(sigma * dot_prod) + std::exp(-sigma * dot_prod)) / 2.0;
}

int vonMisesFisherSample(const Eigen::Vector4d &mean, double pdf_mean, double sigma, double Cp, Eigen::Vector4d &x) {
    for (int i = 0; i < 100000; i++) {
        randomUnitQuaternion(x);
        double pdf = misesFisherKernel(x, mean, sigma, Cp);
        double rand_pdf = randomUniform(0.0, pdf_mean);
        if (rand_pdf < pdf) {
            return i;
        }
    }
    return -1;
}

double orientationNormalKernel(const Eigen::Vector4d &q, const Eigen::Vector4d &mean, double sigma) {
    KDL::Rotation r_q(KDL::Rotation::Quaternion(q(0), q(1), q(2), q(3)));
    KDL::Rotation r_mean(KDL::Rotation::Quaternion(mean(0), mean(1), mean(2), mean(3)));

    KDL::Vector diff = KDL::diff(r_q, r_mean);
    return uniVariateIsotropicGaussianKernel( diff.Norm(), 0.0, sigma );
}

int orientationNormalSample(const Eigen::Vector4d &mean, double sigma, Eigen::Vector4d &x) {

    static boost::mt19937 rng(rand());
    boost::normal_distribution<> nd(0.0, sigma);
    boost::variate_generator<boost::mt19937&, 
                           boost::normal_distribution<> > var_nor(rng, nd);

    double angle = var_nor();

    Eigen::Vector3d vec3;
    randomUnitSphere(vec3);

    KDL::Rotation r_mean(KDL::Rotation::Quaternion(mean(0), mean(1), mean(2), mean(3)));
    KDL::Rotation r_diff = KDL::Rotation::Rot(KDL::Vector(vec3(0), vec3(1), vec3(2)),angle);
    KDL::Rotation r_q = r_mean * r_diff;

    r_q.GetQuaternion(x(0), x(1), x(2), x(3));

    return 1;
}

bool checkSubtreeCollision( const boost::shared_ptr<self_collision::CollisionModel> &col_model, const boost::shared_ptr<KinematicModel > &kin_model,
                            const std::map<std::string, double > &q_map, const std::string &root_name, const KDL::Frame &T_W_T,
                            const boost::shared_ptr< self_collision::Link > &link2, const KDL::Frame &T_B_L2, double &min_dist) {

    Eigen::VectorXd q( kin_model->getJointCount() );
    Eigen::VectorXd ign_q( kin_model->getIgnoredJointCount() );
    q.setZero();
    ign_q.setZero();

    for (std::map<std::string, double >::const_iterator it = q_map.begin(); it != q_map.end(); it++) {
        int idx = kin_model->getJointIndex( it->first );
        if (idx >= 0) {
            q(idx) = it->second;
        }
        else {
            idx = kin_model->getIgnoredJointIndex( it->first );
            if (idx >= 0) {
                ign_q(idx) = it->second;
            }
            else {
                std::cout << "ERROR: checkSubtreeCollision: could not find joint " << it->first << std::endl;
                return false;
            }
        }
    }

//    std::cout << "checkSubtreeCollision: " << q.transpose() << " " << ign_q.transpose() << std::endl;
    KDL::Frame T_W_R;
    kin_model->calculateFk(T_W_R, root_name, q, ign_q);
    KDL::Frame T_R_W = T_W_R.Inverse();

    std::list<std::string > link_names;
    kin_model->getSubtreeLinks(root_name, link_names);

    min_dist = -1.0;

    for (std::list<std::string >::const_iterator it = link_names.begin(); it != link_names.end(); it++) {
        KDL::Frame T_W_Lf;
        kin_model->calculateFk(T_W_Lf, (*it), q, ign_q);
        KDL::Frame T_W_L = T_W_T * T_R_W * T_W_Lf;
        double dist;
        bool collision = self_collision::checkCollision(col_model->getLink( (*it) ), T_W_L, link2, T_B_L2, &dist);
        if (collision) {
            min_dist = 0.0;
            return true;
        }
        std::cout << (*it) << " " << dist << std::endl;
        if (min_dist < 0 || dist < min_dist) {
            min_dist = dist;
        }
    }
    return false;
}

