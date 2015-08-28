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

#include "planer_utils/activation_function.h"

ActivationFunction::ActivationFunction(double z1, double Nmax) {
    // parameters k were calculated in matlab by commands:
    // syms k1 k2 k3 k4 z1 z2 u Nmax
    // [a1,a2,a3,a4,a5]=solve(k1*z1^3 + k2*z1^2 + k3*z1 + k4 == 0, k1*(u+z1)^3 + k2*(u+z1)^2 + k3*(u+z1) + k4 == 1, 3*k1*z1^2 + 2*k2*z1 + k3 == 0, 3*k1*(u+z1)^2 + 2*k2*(u+z1) + k3 == 0, 3*k1*((u+2*z1)/2)^2 + 2*k2*((u+2*z1)/2) + k3 == Nmax, k1, k2, k3, k4, u)

    z1_ = z1;

    k1_ = -(16.0*Nmax*Nmax*Nmax)/27.0;
    k2_ = (16.0*z1_*Nmax*Nmax*Nmax)/9.0 + (4.0*Nmax*Nmax)/3.0;
    k3_ = - (8.0*Nmax*Nmax*z1_)/3.0 - (16.0*Nmax*Nmax*Nmax*z1_*z1_)/9.0;
    k4_ = (4.0*Nmax*Nmax*z1_*z1_*(4.0*Nmax*z1_ + 9.0))/27.0;
    double u = 3.0/(2.0*Nmax);
    z2_ = u + z1_;

//    std::cout << "z2: " << z2_ << std::endl;
//    for (double z = 0.0; z <= 1.0; z += 0.01) {
//        std::cout << "Ndes(" << (z * activation_dist_) << ") = " << func_Ndes(z * activation_dist_) << std::endl;
//    }
}

double ActivationFunction::func_g(double z) const {
    return z * (z * (z * k1_ + k2_) + k3_) + k4_;
}

double ActivationFunction::func_Ndes(double z) const {
    if (z < z1_) {
        return 0.0;
    }
    else if (z <= z2_) {
        return func_g(z);
    }
    return 1.0;
}

