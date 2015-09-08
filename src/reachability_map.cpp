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

#include <kdl/frames.hpp>
#include "Eigen/Dense"

#include "planer_utils/reachability_map.h"
#include "planer_utils/random_uniform.h"

    ReachabilityMap::ReachabilityMap(double voxel_size, int dim) :
        voxel_size_(voxel_size),
        dim_(dim),
        ep_min_(dim),
        ep_max_(dim)
    {
        if (dim_ == 2) {
            for (int y = -1; y <= 1; y++ ) {
                for (int x = -1; x <= 1; x++ ) {
                    if (x == 0 && y == 0) {
                        continue;
                    }
                    std::vector<int > vec;
                    vec.push_back(x);
                    vec.push_back(y);
                    neighbours_.push_back(vec);
                }
            }
        }
        else if (dim_ == 3) {
            for (int z = -1; z <= 1; z++ ) {
                for (int y = -1; y <= 1; y++ ) {
                    for (int x = -1; x <= 1; x++ ) {
                        if (x == 0 && y == 0 && z == 0) {
                            continue;
                        }
                        std::vector<int > vec;
                        vec.push_back(x);
                        vec.push_back(y);
                        vec.push_back(z);
                        neighbours_.push_back(vec);
                    }
                }
            }
        }
        else {
            std::cout << "ERROR: ReachabilityMap: dim should be 2 or 3" << std::endl;
        }
    }

    ReachabilityMap::~ReachabilityMap() {
    }

    void ReachabilityMap::generate(const boost::shared_ptr<KinematicModel> &kin_model, const boost::shared_ptr<self_collision::CollisionModel> &col_model, const std::string &effector_name, int ndof, const Eigen::VectorXd &lower_limit, const Eigen::VectorXd &upper_limit) {
        std::list<Eigen::VectorXd > ep_B_list;
        for (int dim_idx = 0; dim_idx < dim_; dim_idx++) {
            ep_min_(dim_idx) = 1000000.0;
            ep_max_(dim_idx) = -1000000.0;
        }

        std::set<int> excluded_link_idx;
        excluded_link_idx.insert(col_model->getLinkIndex("env_link"));

        int effector_idx = col_model->getLinkIndex(effector_name);
        std::vector<KDL::Frame > links_fk(col_model->getLinksCount());

        for (int i = 0; i < 100000; i++) {
            Eigen::VectorXd tmp_q(ndof);
            for (int q_idx = 0; q_idx < ndof; q_idx++) {
                tmp_q(q_idx) = randomUniform(lower_limit(q_idx), upper_limit(q_idx));
            }

            // calculate forward kinematics for all links
            for (int l_idx = 0; l_idx < col_model->getLinksCount(); l_idx++) {
                kin_model->calculateFk(links_fk[l_idx], col_model->getLinkName(l_idx), tmp_q);
            }

            if (self_collision::checkCollision(col_model, links_fk, excluded_link_idx)) {
                continue;
            }

            const KDL::Frame &T_B_E = links_fk[effector_idx];
            Eigen::VectorXd x(dim_);
            for (int dim_idx = 0; dim_idx < dim_; dim_idx++) {
                x(dim_idx) = T_B_E.p[dim_idx];
                if (ep_min_(dim_idx) > x(dim_idx)) {
                    ep_min_(dim_idx) = x(dim_idx);
                }
                if (ep_max_(dim_idx) < x(dim_idx)) {
                    ep_max_(dim_idx) = x(dim_idx);
                }
            }
            ep_B_list.push_back(x);
        }

        steps_.clear();
        int map_size = 1;
        for (int dim_idx = 0; dim_idx < dim_; dim_idx++) {
            int steps = static_cast<int>( ceil( ( ep_max_(dim_idx) - ep_min_(dim_idx) ) / voxel_size_ ) );
            steps_.push_back( steps );
            map_size *= steps;
        }

        r_map_.resize(map_size, 0);
        p_map_.resize(map_size, 0);

        max_value_ = 0;
        for (std::list<Eigen::VectorXd >::const_iterator it = ep_B_list.begin(); it != ep_B_list.end(); it++) {
            int idx = getIndex( (*it) );
            if (idx < 0) {
                std::cout << "ERROR: ReachabilityMap::generate: idx < 0" << std::endl;
            }
            r_map_[idx]++;
            if (r_map_[idx] > max_value_) {
                max_value_ = r_map_[idx];
            }
        }
    }

    void ReachabilityMap::generateForArm(const boost::shared_ptr<KinematicModel> &kin_model, const std::string &base_name, const std::string &effector_name) {
        std::list<Eigen::VectorXd > ep_B_list;
        for (int dim_idx = 0; dim_idx < dim_; dim_idx++) {
            ep_min_(dim_idx) = 1000000.0;
            ep_max_(dim_idx) = -1000000.0;
        }

        Eigen::VectorXd tmp_q(kin_model->getDofCount());
        for (int i = 0; i < 1000000; i++) {
            for (int q_idx = 0; q_idx < kin_model->getDofCount(); q_idx++) {
                tmp_q(q_idx) = randomUniform(kin_model->getLowerLimit(q_idx), kin_model->getUpperLimit(q_idx));
            }

            KDL::Frame T_W_A, T_W_E;
            kin_model->calculateFk(T_W_A, base_name, tmp_q);
            kin_model->calculateFk(T_W_E, effector_name, tmp_q);

            KDL::Frame T_A_E = T_W_A.Inverse() * T_W_E;

            Eigen::VectorXd x(dim_);
            for (int dim_idx = 0; dim_idx < dim_; dim_idx++) {
                x(dim_idx) = T_A_E.p[dim_idx];
                if (ep_min_(dim_idx) > x(dim_idx)) {
                    ep_min_(dim_idx) = x(dim_idx);
                }
                if (ep_max_(dim_idx) < x(dim_idx)) {
                    ep_max_(dim_idx) = x(dim_idx);
                }
            }
            ep_B_list.push_back(x);
        }

        steps_.clear();
        int map_size = 1;
        for (int dim_idx = 0; dim_idx < dim_; dim_idx++) {
            int steps = static_cast<int>( ceil( ( ep_max_(dim_idx) - ep_min_(dim_idx) ) / voxel_size_ ) );
            steps_.push_back( steps );
            map_size *= steps;
        }

        r_map_.resize(map_size, 0);
        p_map_.resize(map_size, 0);
//        r_map_rot_.resize(map_size, 0);

        max_value_ = 0;
        for (std::list<Eigen::VectorXd >::const_iterator it = ep_B_list.begin(); it != ep_B_list.end(); it++) {
            int idx = getIndex( (*it) );
            if (idx < 0) {
                std::cout << "ERROR: ReachabilityMap::generate: idx < 0" << std::endl;
            }
            r_map_[idx]++;
            if (r_map_[idx] > max_value_) {
                max_value_ = r_map_[idx];
            }
        }
    }

    void ReachabilityMap::generate(const Eigen::VectorXd &lower_bound, const Eigen::VectorXd &upper_bound) {
        ep_min_ = lower_bound;
        ep_max_ = upper_bound;

        steps_.clear();
        int map_size = 1;
        for (int dim_idx = 0; dim_idx < dim_; dim_idx++) {
            int steps = static_cast<int>( ceil( ( ep_max_(dim_idx) - ep_min_(dim_idx) ) / voxel_size_ ) );
            steps_.push_back( steps );
            map_size *= steps;
        }

        r_map_.resize(map_size, 0);
        p_map_.resize(map_size, 0);
        d_map_.resize(map_size, 0);
        max_value_ = 0;
    }



    double ReachabilityMap::getValue(const Eigen::VectorXd &x) const {
        int idx = getIndex(x);
        if (idx < 0) {
            return 0;
        }
        // TODO: check what happens if the score is below 0
        return static_cast<double >(r_map_[idx] - p_map_[idx]) / static_cast<double >(max_value_);
    }

    void ReachabilityMap::setValue(const Eigen::VectorXd &x, int value) {
        int idx = getIndex(x);
        if (idx < 0) {
            return;
        }
        r_map_[idx] = value;
        if (r_map_[idx] > max_value_) {
            max_value_ = r_map_[idx];
        }
    }

    void ReachabilityMap::clear() {
        for (int idx = 0; idx < r_map_.size(); idx++) {
            r_map_[idx] = 0;
        }
        max_value_ = 0;
    }

    void ReachabilityMap::recurenceGrow(const std::list<Eigen::Vector3i > &states_to_expand, boost::function<bool(const KDL::Vector &x)> collision_func, const KDL::Vector &lower_bound, const KDL::Vector &upper_bound) {
        std::list<Eigen::Vector3i > added_states;

        for (std::list<Eigen::Vector3i >::const_iterator it = states_to_expand.begin(); it != states_to_expand.end(); it++) {
            int current_idx = composeIndex( (*it) );
            double current_val = d_map_[current_idx];
            Eigen::Vector3i indices[6] = {
            Eigen::Vector3i(-1, 0, 0) + (*it),
            Eigen::Vector3i(+1, 0, 0) + (*it),
            Eigen::Vector3i(0, -1, 0) + (*it),
            Eigen::Vector3i(0, +1, 0) + (*it),
            Eigen::Vector3i(0, 0, -1) + (*it),
            Eigen::Vector3i(0, 0, +1) + (*it),
            };

            for (int i = 0; i < 6; i++) {
                bool wrong_idx = false;
                for (int dim_idx = 0; dim_idx < 3; dim_idx++) {
                    if (indices[i][dim_idx] < 0 || indices[i][dim_idx] >= steps_[dim_idx]) {
                        wrong_idx = true;
                        break;
                    }
                }
                if (wrong_idx) {
                    continue;
                }

                int pt_idx = composeIndex(indices[i][0], indices[i][1], indices[i][2]);
                if (d_map_[pt_idx] == -1.0) {
                    KDL::Vector pt;
                    getIndexCenter(indices[i][0], indices[i][1], indices[i][2], pt);
                    if (getIndex(pt) != pt_idx) {
                        std::cout << "ERROR: ReachabilityMap::recurenceGrow: getIndex(pt) != pt_idx " << getIndex(pt) << " != " << pt_idx << std::endl;
                    }
                    if (collision_func(pt)) {
                        d_map_[pt_idx] = -2.0;
                    }
                    else {
                        d_map_[pt_idx] = current_val + voxel_size_;
                        added_states.push_back( Eigen::Vector3i(indices[i][0], indices[i][1], indices[i][2]) );
                    }
                }
                else if (d_map_[pt_idx] >= 0.0) {
                    // sanity check
                    if (d_map_[pt_idx] > current_val + voxel_size_) {
                        std::cout << "ERROR: ReachabilityMap::recurenceGrow" << std::endl;
                        return;
                    }
                }
            }
        }

        if (!added_states.empty()) {
            recurenceGrow(added_states, collision_func, lower_bound, upper_bound);
        }
    }

    bool ReachabilityMap::createDistanceMap(const KDL::Vector &origin, boost::function<bool(const KDL::Vector &x)> collision_func, const KDL::Vector &lower_bound, const KDL::Vector &upper_bound) {
        Eigen::VectorXd l_bound(3), u_bound(3);
        for (int i = 0; i < 3; i++) {
            l_bound(i) = lower_bound[i];
            u_bound(i) = upper_bound[i];
        }
        generate(l_bound, u_bound);

        std::cout << "ReachabilityMap::createDistanceMap: distance map size: " << d_map_.size() << std::endl;

        for (int idx = 0; idx < d_map_.size(); idx++) {
            d_map_[idx] = -1.0;
        }

        if (getIndex(origin) < 0) {
            return false;
        }

        origin_ = origin;

        int ix = getIndexDim(origin[0], 0);
        int iy = getIndexDim(origin[1], 1);
        int iz = getIndexDim(origin[2], 2);

        // start at the origin
        d_map_[composeIndex(ix, iy, iz)] = 0.0;

        {
            std::list<Eigen::Vector3i > states_to_expand;
            states_to_expand.push_back(Eigen::Vector3i(ix, iy, iz));
            recurenceGrow(states_to_expand, collision_func, lower_bound, upper_bound);
        }

        bounduary_set_.clear();

        for (int idx = 0; idx < d_map_.size(); idx++) {
            if (d_map_[idx] < 0.0) {
                continue;
            }
            int ix, iy, iz;
            decomposeIndex(idx, ix, iy, iz);
            bool found = false;
            // we are in the obstacle, get the smallest distance in the neighbouring cells
            for (int iix = std::max(0,ix-1); iix < std::min(steps_[0], ix+2) && !found; iix++) {
                for (int iiy = std::max(0,iy-1); iiy < std::min(steps_[1], iy+2) && !found; iiy++) {
                    for (int iiz = std::max(0,iz-1); iiz < std::min(steps_[2], iz+2) && !found; iiz++) {
                        int pt_idx = composeIndex(iix, iiy, iiz);
                        double pt_val = d_map_[pt_idx];
                        if (pt_val < 0.0) {
                            bounduary_set_.insert(idx);
                            found = true;
                        }
                    }
                }
            }
        }
        std::cout << "bounduary_set_.size() " << bounduary_set_.size() << std::endl;
/*
        // increase the value for the bounduary points
        // THIS DOES NOT WORK!
        for (std::set<int >::const_iterator it = bounduary_set_.begin(); it != bounduary_set_.end(); it++) {
            int idx = (*it);

            // get the highest value of non-bounduary neighbours
            double min_value = d_map_[idx];
            for (int iix = std::max(0,ix-1); iix < std::min(steps_[0], ix+2); iix++) {
                for (int iiy = std::max(0,iy-1); iiy < std::min(steps_[1], iy+2); iiy++) {
                    for (int iiz = std::max(0,iz-1); iiz < std::min(steps_[2], iz+2); iiz++) {
                        int pt_idx = composeIndex(iix, iiy, iiz);
                        double pt_val = d_map_[pt_idx];
                        if (pt_val >= 0.0 && pt_val < min_value && bounduary_set_.find(pt_idx) == bounduary_set_.end()) {
                            min_value = d_map_[pt_idx];
                        }
                    }
                }
            }
            min_value += voxel_size_;
            d_map_[idx] = min_value;
        }
//*/
        return true;
    }

    bool ReachabilityMap::getDistnace(const KDL::Vector &x, double &distance) const {
        int idx = getIndex(x);
        if (idx < 0) {
            return false;
        }
        distance = d_map_[idx];
        if (distance < 0.0) {
//            return true;
            int ix = getIndexDim(x.x(), 0);
            int iy = getIndexDim(x.y(), 1);
            int iz = getIndexDim(x.z(), 2);
            distance = 1000000.0;
            // we are in the obstacle, get the smallest distance int the neighbouring cells
            for (int iix = std::max(0,ix-1); iix < std::min(steps_[0], ix+2); iix++) {
                for (int iiy = std::max(0,iy-1); iiy < std::min(steps_[1], iy+2); iiy++) {
                    for (int iiz = std::max(0,iz-1); iiz < std::min(steps_[2], iz+2); iiz++) {
                        int pt_idx = composeIndex(iix, iiy, iiz);
                        double pt_val = d_map_[pt_idx];
                        if (pt_val >= 0.0 && distance > pt_val) {
                            distance = pt_val;
                        }
                    }
                }
            }
        }
        return true;
    }

    bool ReachabilityMap::collisionFreeLine(int ix1, int iy1, int iz1, int ix2, int iy2, int iz2) const {
        KDL::Vector pt1, pt2;
        getIndexCenter(ix1, iy1, iz1, pt1);
        getIndexCenter(ix2, iy2, iz2, pt2);
        KDL::Vector n = pt2-pt1;
        double dist = n.Norm();
        n.Normalize();
        for (double d = 0.0; d < dist; d += voxel_size_) {
            KDL::Vector pt = pt1 + n * d;
            int idx = getIndex(pt);
            if (d_map_[idx] < 0.0) {
                return false;
            }
        }
        return true;
    }

    bool ReachabilityMap::collisionFreeLine(KDL::Vector pt1, int ix2, int iy2, int iz2) const {
        KDL::Vector pt2;
        getIndexCenter(ix2, iy2, iz2, pt2);
        KDL::Vector n = pt2-pt1;
        double dist = n.Norm();
        n.Normalize();
        for (double d = 0.0; d < dist; d += voxel_size_) {
            KDL::Vector pt = pt1 + n * d;
            int idx = getIndex(pt);
            if (d_map_[idx] < 0.0) {
                return false;
            }
        }
        return true;
    }

    bool ReachabilityMap::getGradient(const KDL::Vector &x, KDL::Vector &gradient) const {
        int idx = getIndex(x);
        if (idx < 0) {
//            std::cout << "ReachabilityMap::getGradient: point is outside the map" << std::endl;
            return false;
        }

        double min_value = d_map_[idx];
        bool obstacle = false;
        int search_space = 1;
        if (min_value < 0.0) {
            // we are in the obstacle
            min_value = 100000.0;
            obstacle = true;
//            std::cout << "obstacle" << std::endl;
        }
        int min_ix=-1, min_iy=-1, min_iz=-1;
        int ix = getIndexDim(x.x(), 0);
        int iy = getIndexDim(x.y(), 1);
        int iz = getIndexDim(x.z(), 2);
        for (int iix = std::max(0,ix-search_space); iix < std::min(steps_[0], ix+search_space+1); iix++) {
            for (int iiy = std::max(0,iy-search_space); iiy < std::min(steps_[1], iy+search_space+1); iiy++) {
                for (int iiz = std::max(0,iz-search_space); iiz < std::min(steps_[2], iz+search_space+1); iiz++) {
                    if (ix == iix && iy == iiy && iz == iiz) {
                        continue;
                    }
                    int pt_idx = composeIndex(iix, iiy, iiz);
                    double pt_val = d_map_[pt_idx];
                    if (pt_val >= 0.0 && min_value > pt_val) {
//                    if (pt_val >= 0.0 && (obstacle || collisionFreeLine(x, iix, iiy, iiz)) && min_value > pt_val) {
//                    if (pt_val >= 0.0 && bounduary_set_.find(pt_idx) == bounduary_set_.end() && min_value > pt_val) {
                            min_value = pt_val;
                            min_ix = iix;
                            min_iy = iiy;
                            min_iz = iiz;
                    }
                }
            }
        }

        if (min_ix == -1) {
            std::cout << "ReachabilityMap::getGradient: could not find gradient, min_value: " << min_value << std::endl;
            return false;
        }
        KDL::Vector min_pt, pt;
        getIndexCenter(min_ix, min_iy, min_iz, min_pt);
        getIndexCenter(ix, iy, iz, pt);
        gradient = min_pt - pt;
        gradient.Normalize();
        return true;
    }

    bool ReachabilityMap::getAllGradients(const KDL::Vector &x, std::vector<ReachabilityMap::GradientInfo > &gradients) const {
        if (gradients.size() != 27) {
            std::cout << "ReachabilityMap::getAllGradients: wrong vector size: " << gradients.size() << std::endl;
            return false;
        }

        for (int gradient_idx = 0; gradient_idx < 27; gradient_idx++) {
            gradients[gradient_idx].valid_ = false;
        }

        int idx = getIndex(x);
        if (idx < 0) {
//            std::cout << "ReachabilityMap::getGradient: point is outside the map" << std::endl;
            return false;
        }

        bool isBounduary = bounduary_set_.find(idx) != bounduary_set_.end();

        double min_value = d_map_[idx];
        bool obstacle = false;
        int search_space = 1;
        if (min_value < 0.0) {
//            // we are in the obstacle
//            min_value = 100000.0;
            obstacle = true;
        }
//        int min_ix=-1, min_iy=-1, min_iz=-1;
        int ix = getIndexDim(x.x(), 0);
        int iy = getIndexDim(x.y(), 1);
        int iz = getIndexDim(x.z(), 2);
        int gradient_idx = 0;
//        for (int iix = std::max(0,ix-search_space); iix < std::min(steps_[0], ix+search_space+1); iix++) {
//            for (int iiy = std::max(0,iy-search_space); iiy < std::min(steps_[1], iy+search_space+1); iiy++) {
//                for (int iiz = std::max(0,iz-search_space); iiz < std::min(steps_[2], iz+search_space+1); iiz++) {
        for (int iix = ix-search_space; iix < ix+search_space+1; iix++) {
            for (int iiy = iy-search_space; iiy < iy+search_space+1; iiy++) {
                for (int iiz = iz-search_space; iiz < iz+search_space+1; iiz++) {

                    if ((iix == ix && iiy == iy && iiz == iz) || iix < 0 || iix >= steps_[0] || iiy < 0 || iiy >= steps_[1] || iiz < 0 || iiz >= steps_[2]) {
                        gradient_idx++;
                        continue;
                    }
                    int pt_idx = composeIndex(iix, iiy, iiz);
                    double pt_val = d_map_[pt_idx];

                    if (obstacle && bounduary_set_.find(pt_idx) != bounduary_set_.end()) {
                        gradients[gradient_idx].direction_ = KDL::Vector(iix - ix, iiy - iy, iiz - iz);
                        gradients[gradient_idx].direction_.Normalize();
                        gradients[gradient_idx].value_ = pt_val;// - d_map_[idx];
                        gradients[gradient_idx].valid_ = true;
                    }
                    else if (pt_val >= 0.0 && bounduary_set_.find(pt_idx) == bounduary_set_.end()) {
                        gradients[gradient_idx].direction_ = KDL::Vector(iix - ix, iiy - iy, iiz - iz);
                        gradients[gradient_idx].direction_.Normalize();
                        gradients[gradient_idx].value_ = pt_val;// - d_map_[idx];
                        gradients[gradient_idx].valid_ = true;
                    }
                    gradient_idx++;
                }
            }
        }

        return true;
    }

    void ReachabilityMap::getNeighbourIndices(const std::vector<int> &d, std::list<int> &n_indices) {
        n_indices.clear();
        for (int dim_idx = 0; dim_idx < dim_; dim_idx++) {
            if (d[dim_idx] > 0) {
                int total_idx = 0;
                for (int dim_idx2 = 0; dim_idx2 < dim_; dim_idx2++) {
                    if (dim_idx == dim_idx2) {
                        total_idx = total_idx * steps_[dim_idx2] + d[dim_idx2] - 1;
                    }
                    else {
                        total_idx = total_idx * steps_[dim_idx2] + d[dim_idx2];
                    }
                }
                n_indices.push_back(total_idx);
            }
            if (d[dim_idx] < steps_[dim_idx]-1) {
                int total_idx = 0;
                for (int dim_idx2 = 0; dim_idx2 < dim_; dim_idx2++) {
                    if (dim_idx == dim_idx2) {
                        total_idx = total_idx * steps_[dim_idx2] + d[dim_idx2] + 1;
                    }
                    else {
                        total_idx = total_idx * steps_[dim_idx2] + d[dim_idx2];
                    }
                }
                n_indices.push_back(total_idx);
            }
            
        }

/*        if (d0 > 0) {
            n_indices.push_back(((d0-1) * steps_[1] + d1) * steps_[2] + d2);
        }
        if (d0 < steps_[0]-1) {
            n_indices.push_back(((d0+1) * steps_[1] + d1) * steps_[2] + d2);
        }
*/
    }

    void ReachabilityMap::grow() {
        std::vector<int > map_copy(r_map_);
        for (int idx = 0; idx < map_copy.size(); idx++) {
            if (map_copy[idx] > 1) {
                map_copy[idx] = 1;
            }
        }

        std::vector<int > d(dim_);

        if (dim_ == 2) {
            for (int d0=0; d0<steps_[0]; d0++) {
                for (int d1=0; d1<steps_[1]; d1++) {
                    int idx = d0 * steps_[1] + d1;
                    if (map_copy[idx] == 1) {
                        d[0] = d0;
                        d[1] = d1;
                        std::list<int > n_indices;
                        getNeighbourIndices(d, n_indices);
                        for (std::list<int >::const_iterator it = n_indices.begin(); it != n_indices.end(); it++) {
                            if (map_copy[(*it)] == 0) {
                                map_copy[(*it)] = 2;
                            }
                        }
                    }
                }
            }
        }
        else if (dim_ == 3) {
            for (int d0=0; d0<steps_[0]; d0++) {
                for (int d1=0; d1<steps_[1]; d1++) {
                    for (int d2=0; d2<steps_[2]; d2++) {
                        int idx = (d0 * steps_[1] + d1) * steps_[2] + d2;
                        if (map_copy[idx] == 1) {
                            d[0] = d0;
                            d[1] = d1;
                            d[2] = d2;
                            std::list<int > n_indices;
                            getNeighbourIndices(d, n_indices);
                            for (std::list<int >::const_iterator it = n_indices.begin(); it != n_indices.end(); it++) {
                                if (map_copy[(*it)] == 0) {
                                    map_copy[(*it)] = 2;
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            std::cout << "ReachabilityMap::grow not implemented for " << dim_ << " dimensions" << std::endl;
            return;
        }

        for (int idx = 0; idx < map_copy.size(); idx++) {
            if (map_copy[idx] > 1) {
                map_copy[idx] = 1;
            }
            r_map_[idx] += map_copy[idx];
            if (r_map_[idx] > max_value_) {
                max_value_ = r_map_[idx];
            }
        }
    }

    void ReachabilityMap::addMap(const ReachabilityMap &map) {
        for (int idx = 0; idx < r_map_.size(); idx++) {
            r_map_[idx] += map.r_map_[idx];
            if (r_map_[idx] > max_value_) {
                max_value_ = r_map_[idx];
            }
        }
    }

    void ReachabilityMap::addMap(const boost::shared_ptr<ReachabilityMap > &pmap) {
        addMap( (*(pmap.get())) );
    }

    double ReachabilityMap::getMaxValue() const {
        return static_cast<double >(max_value_);
    }

    void ReachabilityMap::addPenalty(const Eigen::VectorXd &x) {
        int idx = getIndex(x);
        if (idx >= 0) {
            p_map_[idx] += max_value_;
        }
    }

    void ReachabilityMap::resetPenalty() {
        for (int idx = 0; idx < p_map_.size(); idx++) {
            p_map_[idx] = 0;
        }
    }

    int ReachabilityMap::getIndex(const Eigen::VectorXd &x) const {
        int total_idx = 0;
        for (int dim_idx = 0; dim_idx < dim_; dim_idx++) {
            int idx = static_cast<int >( floor( (x(dim_idx) - ep_min_(dim_idx)) / voxel_size_ ) );
            if (idx < 0 || idx >= steps_[dim_idx]) {
                return -1;
            }
            total_idx = total_idx * steps_[dim_idx] + idx;
        }
        return total_idx;
    }

    int ReachabilityMap::getIndex(const KDL::Vector &x) const {
        int total_idx = 0;
        for (int dim_idx = 0; dim_idx < dim_; dim_idx++) {
            int idx = static_cast<int >( (x[dim_idx] - ep_min_(dim_idx)) / voxel_size_ );
            if (idx < 0 || idx >= steps_[dim_idx]) {
                return -1;
            }
            total_idx = total_idx * steps_[dim_idx] + idx;
        }
        return total_idx;
    }

    int ReachabilityMap::getIndexDim(double x, int dim_idx) const {
        int idx = static_cast<int >( (x - ep_min_(dim_idx)) / voxel_size_ );
        if (idx < 0 || idx >= steps_[dim_idx]) {
            return -1;
        }
        return idx;
    }

    int ReachabilityMap::composeIndex(const Eigen::Vector3i &i) const {
        return (i(0) * steps_[1] + i(1)) * steps_[2] + i(2);
    }

    int ReachabilityMap::composeIndex(int ix, int iy, int iz) const {
        return (ix * steps_[1] + iy) * steps_[2] + iz;
    }

    void ReachabilityMap::decomposeIndex(int idx, int &ix, int &iy, int &iz) const {
        iz = idx % steps_[2];
        idx = (idx - iz) / steps_[2];
        iy = idx % steps_[1];
        idx = (idx - iy) / steps_[1];
        ix = idx;
    }

    void ReachabilityMap::getIndexCenter(int ix, int iy, int iz, KDL::Vector &pt) const {
        pt = KDL::Vector( (((double)ix)+0.5) * voxel_size_ + ep_min_(0), (((double)iy)+0.5) * voxel_size_ + ep_min_(1), (((double)iz)+0.5) * voxel_size_ + ep_min_(2) );
    }

    const KDL::Vector &ReachabilityMap::getOrigin() const {
        return origin_;
    }

