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

#include "planer_utils/random_uniform.h"

#include <stdlib.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

double randomUniform(double min, double max) {
    return min + (max - min) * (double)rand() / (double)RAND_MAX;
}

void randomUnitQuaternion(Eigen::VectorXd &quat) {
    static boost::random::mt19937 rng(time(NULL));
    static boost::random::normal_distribution<> normal;

    do {
        quat(0) = normal(rng);
        quat(1) = normal(rng);
        quat(2) = normal(rng);
        quat(3) = normal(rng);
    } while (quat.norm() < 0.00001);
    quat.normalize();
}

void randomUnitQuaternion(Eigen::Vector4d &quat) {
    static boost::random::mt19937 rng(time(NULL));
    static boost::random::normal_distribution<> normal;

    do {
        quat(0) = normal(rng);
        quat(1) = normal(rng);
        quat(2) = normal(rng);
        quat(3) = normal(rng);
    } while (quat.norm() < 0.00001);
    quat.normalize();
}

void randomUnitSphere(Eigen::VectorXd &sp) {
    static boost::random::mt19937 rng(time(NULL));
    static boost::random::normal_distribution<> normal;

    do {
        sp(0) = normal(rng);
        sp(1) = normal(rng);
        sp(2) = normal(rng);
    } while (sp.norm() < 0.00001);
    sp.normalize();
}

void randomUnitSphere(Eigen::Vector3d &sp) {
    static boost::random::mt19937 rng(time(NULL));
    static boost::random::normal_distribution<> normal;

    do {
        sp(0) = normal(rng);
        sp(1) = normal(rng);
        sp(2) = normal(rng);
    } while (sp.norm() < 0.00001);
    sp.normalize();
}

