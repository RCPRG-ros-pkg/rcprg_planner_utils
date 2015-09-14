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

//
// implements:
// Tricubic interpolation in three dimensions
// by F. Lekien and J. Marsden
//
// sources:
// https://github.com/nbigaouette/libtricubic/tree/master/tricubic-1.0/src/libtricubic
// https://github.com/danielguterding/pytricubic
//

int A[64][64] = {
{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 9,-9,-9, 9, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-6, 6, 6,-6, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-6, 6, 6,-6, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 4,-4,-4, 4, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
{-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 9,-9, 0, 0,-9, 9, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-6, 6, 0, 0, 6,-6, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9, 0, 0,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0},
{ 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0},
{-27,27,27,-27,27,-27,-27,27,-18,-9,18, 9,18, 9,-18,-9,-18,18,-9, 9,18,-18, 9,-9,-18,18,18,-18,-9, 9, 9,-9,-12,-6,-6,-3,12, 6, 6, 3,-12,-6,12, 6,-6,-3, 6, 3,-12,12,-6, 6,-6, 6,-3, 3,-8,-4,-4,-2,-4,-2,-2,-1},
{18,-18,-18,18,-18,18,18,-18, 9, 9,-9,-9,-9,-9, 9, 9,12,-12, 6,-6,-12,12,-6, 6,12,-12,-12,12, 6,-6,-6, 6, 6, 6, 3, 3,-6,-6,-3,-3, 6, 6,-6,-6, 3, 3,-3,-3, 8,-8, 4,-4, 4,-4, 2,-2, 4, 4, 2, 2, 2, 2, 1, 1},
{-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0},
{18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6, 9,-9, 9,-9,-9, 9,-9, 9,12,-12,-12,12, 6,-6,-6, 6, 6, 3, 6, 3,-6,-3,-6,-3, 8, 4,-8,-4, 4, 2,-4,-2, 6,-6, 6,-6, 3,-3, 3,-3, 4, 2, 4, 2, 2, 1, 2, 1},
{-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-6, 6,-6, 6, 6,-6, 6,-6,-8, 8, 8,-8,-4, 4, 4,-4,-3,-3,-3,-3, 3, 3, 3, 3,-4,-4, 4, 4,-2,-2, 2, 2,-4, 4,-4, 4,-2, 2,-2, 2,-2,-2,-2,-2,-1,-1,-1,-1},
{ 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-6, 6, 0, 0, 6,-6, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 4,-4, 0, 0,-4, 4, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4, 0, 0,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0},
{-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0},
{18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6,12,-12, 6,-6,-12,12,-6, 6, 9,-9,-9, 9, 9,-9,-9, 9, 8, 4, 4, 2,-8,-4,-4,-2, 6, 3,-6,-3, 6, 3,-6,-3, 6,-6, 3,-3, 6,-6, 3,-3, 4, 2, 2, 1, 4, 2, 2, 1},
{-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-8, 8,-4, 4, 8,-8, 4,-4,-6, 6, 6,-6,-6, 6, 6,-6,-4,-4,-2,-2, 4, 4, 2, 2,-3,-3, 3, 3,-3,-3, 3, 3,-4, 4,-2, 2,-4, 4,-2, 2,-2,-2,-1,-1,-2,-2,-1,-1},
{ 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0},
{-12,12,12,-12,12,-12,-12,12,-8,-4, 8, 4, 8, 4,-8,-4,-6, 6,-6, 6, 6,-6, 6,-6,-6, 6, 6,-6,-6, 6, 6,-6,-4,-2,-4,-2, 4, 2, 4, 2,-4,-2, 4, 2,-4,-2, 4, 2,-3, 3,-3, 3,-3, 3,-3, 3,-2,-1,-2,-1,-2,-1,-2,-1},
{ 8,-8,-8, 8,-8, 8, 8,-8, 4, 4,-4,-4,-4,-4, 4, 4, 4,-4, 4,-4,-4, 4,-4, 4, 4,-4,-4, 4, 4,-4,-4, 4, 2, 2, 2, 2,-2,-2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2, 2,-2, 2,-2, 2,-2, 2,-2, 1, 1, 1, 1, 1, 1, 1, 1}};


void tricubic_get_coeff_stacked(double a[64], double x[64]) {
  int i,j;
  for (i=0;i<64;i++) {
    a[i]=(double)(0.0);
    for (j=0;j<64;j++) {
      a[i]+=A[i][j]*x[j];
    }
  }
}

void tricubic_get_coeff(double a[64], double f[8], double dfdx[8], double dfdy[8], double dfdz[8], double d2fdxdy[8], double d2fdxdz[8], double d2fdydz[8], double d3fdxdydz[8]) {
  int i;
  double x[64];
  for (i=0;i<8;i++) {
    x[0+i]=f[i];
    x[8+i]=dfdx[i];
    x[16+i]=dfdy[i];
    x[24+i]=dfdz[i];
    x[32+i]=d2fdxdy[i];
    x[40+i]=d2fdxdz[i];
    x[48+i]=d2fdydz[i];
    x[56+i]=d3fdxdydz[i];
  }
  tricubic_get_coeff_stacked(a,x);
}

int ijk2n(int i, int j, int k) {
  return(i+4*j+16*k);
}

double tricubic_eval(double a[64], double x, double y, double z, int derx, int dery, int derz) {
  int i,j,k;
  double ret=(double)(0.0);
  double cont;
  int w;
  /* TRICUBIC_EVAL 
     The full version takes 3 extra integers args that allows to evaluate
     any partial derivative of f at the point
     derx=dery=derz=0 => f
     derx=2 dery=derz=0 => d2f/dx2
     derx=dery=derz=1 =? d3f/dxdydz
     NOTICE that (derx>3)||(dery>3)||(derz>3) => returns 0.0
     this computes   \frac{\partial ^{derx+dery+derz} d}{\partial x ^{derx} \partial y ^{dery} \partial z ^{derz}}
  */
  for (i=derx;i<4;i++) {
    for (j=dery;j<4;j++) {
      for (k=derz;k<4;k++) {
	cont=a[ijk2n(i,j,k)]*pow(x,i-derx)*pow(y,j-dery)*pow(z,k-derz);
	for (w=0;w<derx;w++) {
	  cont*=(i-w);
	}
	for (w=0;w<dery;w++) {
	  cont*=(j-w);
	}
	for (w=0;w<derz;w++) {
	  cont*=(k-w);
	}
	ret+=cont;
      }
    }
  }
  return(ret);
}

void ReachabilityMap::tricubic_get_coeff(double a[64], int xi, int yi, int zi) const {
    int i;

    double x[64] = {
      // values of f(x,y,z) at each corner.
      d_map_[composeIndex(xi,yi,zi)],d_map_[composeIndex(xi+1,yi,zi)],d_map_[composeIndex(xi,yi+1,zi)],
      d_map_[composeIndex(xi+1,yi+1,zi)],d_map_[composeIndex(xi,yi,zi+1)],d_map_[composeIndex(xi+1,yi,zi+1)],
      d_map_[composeIndex(xi,yi+1,zi+1)],d_map_[composeIndex(xi+1,yi+1,zi+1)],
      // values of df/dx at each corner.
      0.5*(d_map_[composeIndex(xi+1,yi,zi)]-d_map_[composeIndex(xi-1,yi,zi)]),
      0.5*(d_map_[composeIndex(xi+2,yi,zi)]-d_map_[composeIndex(xi,yi,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi)]-d_map_[composeIndex(xi-1,yi+1,zi)]),
      0.5*(d_map_[composeIndex(xi+2,yi+1,zi)]-d_map_[composeIndex(xi,yi+1,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi,zi+1)]-d_map_[composeIndex(xi-1,yi,zi+1)]),
      0.5*(d_map_[composeIndex(xi+2,yi,zi+1)]-d_map_[composeIndex(xi,yi,zi+1)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi-1,yi+1,zi+1)]),
      0.5*(d_map_[composeIndex(xi+2,yi+1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi+1)]),
      // values of df/dy at each corner.
      0.5*(d_map_[composeIndex(xi,yi+1,zi)]-d_map_[composeIndex(xi,yi-1,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi)]-d_map_[composeIndex(xi+1,yi-1,zi)]),
      0.5*(d_map_[composeIndex(xi,yi+2,zi)]-d_map_[composeIndex(xi,yi,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi+2,zi)]-d_map_[composeIndex(xi+1,yi,zi)]),
      0.5*(d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi,yi-1,zi+1)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi-1,zi+1)]),
      0.5*(d_map_[composeIndex(xi,yi+2,zi+1)]-d_map_[composeIndex(xi,yi,zi+1)]),
      0.5*(d_map_[composeIndex(xi+1,yi+2,zi+1)]-d_map_[composeIndex(xi+1,yi,zi+1)]),
      // values of df/dz at each corner.
      0.5*(d_map_[composeIndex(xi,yi,zi+1)]-d_map_[composeIndex(xi,yi,zi-1)]),
      0.5*(d_map_[composeIndex(xi+1,yi,zi+1)]-d_map_[composeIndex(xi+1,yi,zi-1)]),
      0.5*(d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi-1)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi+1,zi-1)]),
      0.5*(d_map_[composeIndex(xi,yi,zi+2)]-d_map_[composeIndex(xi,yi,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi,zi+2)]-d_map_[composeIndex(xi+1,yi,zi)]),
      0.5*(d_map_[composeIndex(xi,yi+1,zi+2)]-d_map_[composeIndex(xi,yi+1,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi+2)]-d_map_[composeIndex(xi+1,yi+1,zi)]),
      // values of d2f/dxdy at each corner.
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi)]-d_map_[composeIndex(xi-1,yi+1,zi)]-d_map_[composeIndex(xi+1,yi-1,zi)]+d_map_[composeIndex(xi-1,yi-1,zi)]),
      0.25*(d_map_[composeIndex(xi+2,yi+1,zi)]-d_map_[composeIndex(xi,yi+1,zi)]-d_map_[composeIndex(xi+2,yi-1,zi)]+d_map_[composeIndex(xi,yi-1,zi)]),
      0.25*(d_map_[composeIndex(xi+1,yi+2,zi)]-d_map_[composeIndex(xi-1,yi+2,zi)]-d_map_[composeIndex(xi+1,yi,zi)]+d_map_[composeIndex(xi-1,yi,zi)]),
      0.25*(d_map_[composeIndex(xi+2,yi+2,zi)]-d_map_[composeIndex(xi,yi+2,zi)]-d_map_[composeIndex(xi+2,yi,zi)]+d_map_[composeIndex(xi,yi,zi)]),
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi-1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi-1,zi+1)]+d_map_[composeIndex(xi-1,yi-1,zi+1)]),
      0.25*(d_map_[composeIndex(xi+2,yi+1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi+2,yi-1,zi+1)]+d_map_[composeIndex(xi,yi-1,zi+1)]),
      0.25*(d_map_[composeIndex(xi+1,yi+2,zi+1)]-d_map_[composeIndex(xi-1,yi+2,zi+1)]-d_map_[composeIndex(xi+1,yi,zi+1)]+d_map_[composeIndex(xi-1,yi,zi+1)]),
      0.25*(d_map_[composeIndex(xi+2,yi+2,zi+1)]-d_map_[composeIndex(xi,yi+2,zi+1)]-d_map_[composeIndex(xi+2,yi,zi+1)]+d_map_[composeIndex(xi,yi,zi+1)]),
      // values of d2f/dxdz at each corner.
      0.25*(d_map_[composeIndex(xi+1,yi,zi+1)]-d_map_[composeIndex(xi-1,yi,zi+1)]-d_map_[composeIndex(xi+1,yi,zi-1)]+d_map_[composeIndex(xi-1,yi,zi-1)]),
      0.25*(d_map_[composeIndex(xi+2,yi,zi+1)]-d_map_[composeIndex(xi,yi,zi+1)]-d_map_[composeIndex(xi+2,yi,zi-1)]+d_map_[composeIndex(xi,yi,zi-1)]),
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi-1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi+1,zi-1)]+d_map_[composeIndex(xi-1,yi+1,zi-1)]),
      0.25*(d_map_[composeIndex(xi+2,yi+1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi+2,yi+1,zi-1)]+d_map_[composeIndex(xi,yi+1,zi-1)]),
      0.25*(d_map_[composeIndex(xi+1,yi,zi+2)]-d_map_[composeIndex(xi-1,yi,zi+2)]-d_map_[composeIndex(xi+1,yi,zi)]+d_map_[composeIndex(xi-1,yi,zi)]),
      0.25*(d_map_[composeIndex(xi+2,yi,zi+2)]-d_map_[composeIndex(xi,yi,zi+2)]-d_map_[composeIndex(xi+2,yi,zi)]+d_map_[composeIndex(xi,yi,zi)]),
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi+2)]-d_map_[composeIndex(xi-1,yi+1,zi+2)]-d_map_[composeIndex(xi+1,yi+1,zi)]+d_map_[composeIndex(xi-1,yi+1,zi)]),
      0.25*(d_map_[composeIndex(xi+2,yi+1,zi+2)]-d_map_[composeIndex(xi,yi+1,zi+2)]-d_map_[composeIndex(xi+2,yi+1,zi)]+d_map_[composeIndex(xi,yi+1,zi)]),
      // values of d2f/dydz at each corner.
      0.25*(d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi,yi-1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi-1)]+d_map_[composeIndex(xi,yi-1,zi-1)]),
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi-1,zi+1)]-d_map_[composeIndex(xi+1,yi+1,zi-1)]+d_map_[composeIndex(xi+1,yi-1,zi-1)]),
      0.25*(d_map_[composeIndex(xi,yi+2,zi+1)]-d_map_[composeIndex(xi,yi,zi+1)]-d_map_[composeIndex(xi,yi+2,zi-1)]+d_map_[composeIndex(xi,yi,zi-1)]),
      0.25*(d_map_[composeIndex(xi+1,yi+2,zi+1)]-d_map_[composeIndex(xi+1,yi,zi+1)]-d_map_[composeIndex(xi+1,yi+2,zi-1)]+d_map_[composeIndex(xi+1,yi,zi-1)]),
      0.25*(d_map_[composeIndex(xi,yi+1,zi+2)]-d_map_[composeIndex(xi,yi-1,zi+2)]-d_map_[composeIndex(xi,yi+1,zi)]+d_map_[composeIndex(xi,yi-1,zi)]),
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi+2)]-d_map_[composeIndex(xi+1,yi-1,zi+2)]-d_map_[composeIndex(xi+1,yi+1,zi)]+d_map_[composeIndex(xi+1,yi-1,zi)]),
      0.25*(d_map_[composeIndex(xi,yi+2,zi+2)]-d_map_[composeIndex(xi,yi,zi+2)]-d_map_[composeIndex(xi,yi+2,zi)]+d_map_[composeIndex(xi,yi,zi)]),
      0.25*(d_map_[composeIndex(xi+1,yi+2,zi+2)]-d_map_[composeIndex(xi+1,yi,zi+2)]-d_map_[composeIndex(xi+1,yi+2,zi)]+d_map_[composeIndex(xi+1,yi,zi)]),
      // values of d3f/dxdydz at each corner.
      0.125*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi-1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi-1,zi+1)]+d_map_[composeIndex(xi-1,yi-1,zi+1)]-d_map_[composeIndex(xi+1,yi+1,zi-1)]+d_map_[composeIndex(xi-1,yi+1,zi-1)]+d_map_[composeIndex(xi+1,yi-1,zi-1)]-d_map_[composeIndex(xi-1,yi-1,zi-1)]),
      0.125*(d_map_[composeIndex(xi+2,yi+1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi+2,yi-1,zi+1)]+d_map_[composeIndex(xi,yi-1,zi+1)]-d_map_[composeIndex(xi+2,yi+1,zi-1)]+d_map_[composeIndex(xi,yi+1,zi-1)]+d_map_[composeIndex(xi+2,yi-1,zi-1)]-d_map_[composeIndex(xi,yi-1,zi-1)]),
      0.125*(d_map_[composeIndex(xi+1,yi+2,zi+1)]-d_map_[composeIndex(xi-1,yi+2,zi+1)]-d_map_[composeIndex(xi+1,yi,zi+1)]+d_map_[composeIndex(xi-1,yi,zi+1)]-d_map_[composeIndex(xi+1,yi+2,zi-1)]+d_map_[composeIndex(xi-1,yi+2,zi-1)]+d_map_[composeIndex(xi+1,yi,zi-1)]-d_map_[composeIndex(xi-1,yi,zi-1)]),
      0.125*(d_map_[composeIndex(xi+2,yi+2,zi+1)]-d_map_[composeIndex(xi,yi+2,zi+1)]-d_map_[composeIndex(xi+2,yi,zi+1)]+d_map_[composeIndex(xi,yi,zi+1)]-d_map_[composeIndex(xi+2,yi+2,zi-1)]+d_map_[composeIndex(xi,yi+2,zi-1)]+d_map_[composeIndex(xi+2,yi,zi-1)]-d_map_[composeIndex(xi,yi,zi-1)]),
      0.125*(d_map_[composeIndex(xi+1,yi+1,zi+2)]-d_map_[composeIndex(xi-1,yi+1,zi+2)]-d_map_[composeIndex(xi+1,yi-1,zi+2)]+d_map_[composeIndex(xi-1,yi-1,zi+2)]-d_map_[composeIndex(xi+1,yi+1,zi)]+d_map_[composeIndex(xi-1,yi+1,zi)]+d_map_[composeIndex(xi+1,yi-1,zi)]-d_map_[composeIndex(xi-1,yi-1,zi)]),
      0.125*(d_map_[composeIndex(xi+2,yi+1,zi+2)]-d_map_[composeIndex(xi,yi+1,zi+2)]-d_map_[composeIndex(xi+2,yi-1,zi+2)]+d_map_[composeIndex(xi,yi-1,zi+2)]-d_map_[composeIndex(xi+2,yi+1,zi)]+d_map_[composeIndex(xi,yi+1,zi)]+d_map_[composeIndex(xi+2,yi-1,zi)]-d_map_[composeIndex(xi,yi-1,zi)]),
      0.125*(d_map_[composeIndex(xi+1,yi+2,zi+2)]-d_map_[composeIndex(xi-1,yi+2,zi+2)]-d_map_[composeIndex(xi+1,yi,zi+2)]+d_map_[composeIndex(xi-1,yi,zi+2)]-d_map_[composeIndex(xi+1,yi+2,zi)]+d_map_[composeIndex(xi-1,yi+2,zi)]+d_map_[composeIndex(xi+1,yi,zi)]-d_map_[composeIndex(xi-1,yi,zi)]),
      0.125*(d_map_[composeIndex(xi+2,yi+2,zi+2)]-d_map_[composeIndex(xi,yi+2,zi+2)]-d_map_[composeIndex(xi+2,yi,zi+2)]+d_map_[composeIndex(xi,yi,zi+2)]-d_map_[composeIndex(xi+2,yi+2,zi)]+d_map_[composeIndex(xi,yi+2,zi)]+d_map_[composeIndex(xi+2,yi,zi)]-d_map_[composeIndex(xi,yi,zi)])
    };
    tricubic_get_coeff_stacked(a,x);
}
/*
fptype ReachabilityMap::ip(list xyz){
  
  fptype x = extract<fptype>(xyz[0]);
  fptype y = extract<fptype>(xyz[1]);
  fptype z = extract<fptype>(xyz[2]);
  
  fptype dx = fmod(x/_spacing, _n1), dy = fmod(y/_spacing, _n2), dz = fmod(z/_spacing, _n3); //determine the relative position in the box enclosed by nearest data points
  
  if(dx < 0) dx += _n1; //periodicity is built in
  if(dy < 0) dy += _n2;
  if(dz < 0) dz += _n3;
  
  int xi = (int)floor(dx); //calculate lower-bound grid indices
  int yi = (int)floor(dy);
  int zi = (int)floor(dz);
  
  // Check if we can re-use coefficients from the last interpolation.
  if(!_initialized || xi != _i1 || yi != _i2 || zi != _i3) {
  // Extract the local vocal values and calculate partial derivatives.
  Eigen::Matrix<fptype,64,1> x;
  x << 
      // values of f(x,y,z) at each corner.
      d_map_[composeIndex(xi,yi,zi)],d_map_[composeIndex(xi+1,yi,zi)],d_map_[composeIndex(xi,yi+1,zi)],
      d_map_[composeIndex(xi+1,yi+1,zi)],d_map_[composeIndex(xi,yi,zi+1)],d_map_[composeIndex(xi+1,yi,zi+1)],
      d_map_[composeIndex(xi,yi+1,zi+1)],d_map_[composeIndex(xi+1,yi+1,zi+1)],
      // values of df/dx at each corner.
      0.5*(d_map_[composeIndex(xi+1,yi,zi)]-d_map_[composeIndex(xi-1,yi,zi)]),
      0.5*(d_map_[composeIndex(xi+2,yi,zi)]-d_map_[composeIndex(xi,yi,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi)]-d_map_[composeIndex(xi-1,yi+1,zi)]),
      0.5*(d_map_[composeIndex(xi+2,yi+1,zi)]-d_map_[composeIndex(xi,yi+1,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi,zi+1)]-d_map_[composeIndex(xi-1,yi,zi+1)]),
      0.5*(d_map_[composeIndex(xi+2,yi,zi+1)]-d_map_[composeIndex(xi,yi,zi+1)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi-1,yi+1,zi+1)]),
      0.5*(d_map_[composeIndex(xi+2,yi+1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi+1)]),
      // values of df/dy at each corner.
      0.5*(d_map_[composeIndex(xi,yi+1,zi)]-d_map_[composeIndex(xi,yi-1,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi)]-d_map_[composeIndex(xi+1,yi-1,zi)]),
      0.5*(d_map_[composeIndex(xi,yi+2,zi)]-d_map_[composeIndex(xi,yi,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi+2,zi)]-d_map_[composeIndex(xi+1,yi,zi)]),
      0.5*(d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi,yi-1,zi+1)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi-1,zi+1)]),
      0.5*(d_map_[composeIndex(xi,yi+2,zi+1)]-d_map_[composeIndex(xi,yi,zi+1)]),
      0.5*(d_map_[composeIndex(xi+1,yi+2,zi+1)]-d_map_[composeIndex(xi+1,yi,zi+1)]),
      // values of df/dz at each corner.
      0.5*(d_map_[composeIndex(xi,yi,zi+1)]-d_map_[composeIndex(xi,yi,zi-1)]),
      0.5*(d_map_[composeIndex(xi+1,yi,zi+1)]-d_map_[composeIndex(xi+1,yi,zi-1)]),
      0.5*(d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi-1)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi+1,zi-1)]),
      0.5*(d_map_[composeIndex(xi,yi,zi+2)]-d_map_[composeIndex(xi,yi,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi,zi+2)]-d_map_[composeIndex(xi+1,yi,zi)]),
      0.5*(d_map_[composeIndex(xi,yi+1,zi+2)]-d_map_[composeIndex(xi,yi+1,zi)]),
      0.5*(d_map_[composeIndex(xi+1,yi+1,zi+2)]-d_map_[composeIndex(xi+1,yi+1,zi)]),
      // values of d2f/dxdy at each corner.
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi)]-d_map_[composeIndex(xi-1,yi+1,zi)]-d_map_[composeIndex(xi+1,yi-1,zi)]+d_map_[composeIndex(xi-1,yi-1,zi)]),
      0.25*(d_map_[composeIndex(xi+2,yi+1,zi)]-d_map_[composeIndex(xi,yi+1,zi)]-d_map_[composeIndex(xi+2,yi-1,zi)]+d_map_[composeIndex(xi,yi-1,zi)]),
      0.25*(d_map_[composeIndex(xi+1,yi+2,zi)]-d_map_[composeIndex(xi-1,yi+2,zi)]-d_map_[composeIndex(xi+1,yi,zi)]+d_map_[composeIndex(xi-1,yi,zi)]),
      0.25*(d_map_[composeIndex(xi+2,yi+2,zi)]-d_map_[composeIndex(xi,yi+2,zi)]-d_map_[composeIndex(xi+2,yi,zi)]+d_map_[composeIndex(xi,yi,zi)]),
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi-1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi-1,zi+1)]+d_map_[composeIndex(xi-1,yi-1,zi+1)]),
      0.25*(d_map_[composeIndex(xi+2,yi+1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi+2,yi-1,zi+1)]+d_map_[composeIndex(xi,yi-1,zi+1)]),
      0.25*(d_map_[composeIndex(xi+1,yi+2,zi+1)]-d_map_[composeIndex(xi-1,yi+2,zi+1)]-d_map_[composeIndex(xi+1,yi,zi+1)]+d_map_[composeIndex(xi-1,yi,zi+1)]),
      0.25*(d_map_[composeIndex(xi+2,yi+2,zi+1)]-d_map_[composeIndex(xi,yi+2,zi+1)]-d_map_[composeIndex(xi+2,yi,zi+1)]+d_map_[composeIndex(xi,yi,zi+1)]),
      // values of d2f/dxdz at each corner.
      0.25*(d_map_[composeIndex(xi+1,yi,zi+1)]-d_map_[composeIndex(xi-1,yi,zi+1)]-d_map_[composeIndex(xi+1,yi,zi-1)]+d_map_[composeIndex(xi-1,yi,zi-1)]),
      0.25*(d_map_[composeIndex(xi+2,yi,zi+1)]-d_map_[composeIndex(xi,yi,zi+1)]-d_map_[composeIndex(xi+2,yi,zi-1)]+d_map_[composeIndex(xi,yi,zi-1)]),
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi-1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi+1,zi-1)]+d_map_[composeIndex(xi-1,yi+1,zi-1)]),
      0.25*(d_map_[composeIndex(xi+2,yi+1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi+2,yi+1,zi-1)]+d_map_[composeIndex(xi,yi+1,zi-1)]),
      0.25*(d_map_[composeIndex(xi+1,yi,zi+2)]-d_map_[composeIndex(xi-1,yi,zi+2)]-d_map_[composeIndex(xi+1,yi,zi)]+d_map_[composeIndex(xi-1,yi,zi)]),
      0.25*(d_map_[composeIndex(xi+2,yi,zi+2)]-d_map_[composeIndex(xi,yi,zi+2)]-d_map_[composeIndex(xi+2,yi,zi)]+d_map_[composeIndex(xi,yi,zi)]),
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi+2)]-d_map_[composeIndex(xi-1,yi+1,zi+2)]-d_map_[composeIndex(xi+1,yi+1,zi)]+d_map_[composeIndex(xi-1,yi+1,zi)]),
      0.25*(d_map_[composeIndex(xi+2,yi+1,zi+2)]-d_map_[composeIndex(xi,yi+1,zi+2)]-d_map_[composeIndex(xi+2,yi+1,zi)]+d_map_[composeIndex(xi,yi+1,zi)]),
      // values of d2f/dydz at each corner.
      0.25*(d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi,yi-1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi-1)]+d_map_[composeIndex(xi,yi-1,zi-1)]),
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi-1,zi+1)]-d_map_[composeIndex(xi+1,yi+1,zi-1)]+d_map_[composeIndex(xi+1,yi-1,zi-1)]),
      0.25*(d_map_[composeIndex(xi,yi+2,zi+1)]-d_map_[composeIndex(xi,yi,zi+1)]-d_map_[composeIndex(xi,yi+2,zi-1)]+d_map_[composeIndex(xi,yi,zi-1)]),
      0.25*(d_map_[composeIndex(xi+1,yi+2,zi+1)]-d_map_[composeIndex(xi+1,yi,zi+1)]-d_map_[composeIndex(xi+1,yi+2,zi-1)]+d_map_[composeIndex(xi+1,yi,zi-1)]),
      0.25*(d_map_[composeIndex(xi,yi+1,zi+2)]-d_map_[composeIndex(xi,yi-1,zi+2)]-d_map_[composeIndex(xi,yi+1,zi)]+d_map_[composeIndex(xi,yi-1,zi)]),
      0.25*(d_map_[composeIndex(xi+1,yi+1,zi+2)]-d_map_[composeIndex(xi+1,yi-1,zi+2)]-d_map_[composeIndex(xi+1,yi+1,zi)]+d_map_[composeIndex(xi+1,yi-1,zi)]),
      0.25*(d_map_[composeIndex(xi,yi+2,zi+2)]-d_map_[composeIndex(xi,yi,zi+2)]-d_map_[composeIndex(xi,yi+2,zi)]+d_map_[composeIndex(xi,yi,zi)]),
      0.25*(d_map_[composeIndex(xi+1,yi+2,zi+2)]-d_map_[composeIndex(xi+1,yi,zi+2)]-d_map_[composeIndex(xi+1,yi+2,zi)]+d_map_[composeIndex(xi+1,yi,zi)]),
      // values of d3f/dxdydz at each corner.
      0.125*(d_map_[composeIndex(xi+1,yi+1,zi+1)]-d_map_[composeIndex(xi-1,yi+1,zi+1)]-d_map_[composeIndex(xi+1,yi-1,zi+1)]+d_map_[composeIndex(xi-1,yi-1,zi+1)]-d_map_[composeIndex(xi+1,yi+1,zi-1)]+d_map_[composeIndex(xi-1,yi+1,zi-1)]+d_map_[composeIndex(xi+1,yi-1,zi-1)]-d_map_[composeIndex(xi-1,yi-1,zi-1)]),
      0.125*(d_map_[composeIndex(xi+2,yi+1,zi+1)]-d_map_[composeIndex(xi,yi+1,zi+1)]-d_map_[composeIndex(xi+2,yi-1,zi+1)]+d_map_[composeIndex(xi,yi-1,zi+1)]-d_map_[composeIndex(xi+2,yi+1,zi-1)]+d_map_[composeIndex(xi,yi+1,zi-1)]+d_map_[composeIndex(xi+2,yi-1,zi-1)]-d_map_[composeIndex(xi,yi-1,zi-1)]),
      0.125*(d_map_[composeIndex(xi+1,yi+2,zi+1)]-d_map_[composeIndex(xi-1,yi+2,zi+1)]-d_map_[composeIndex(xi+1,yi,zi+1)]+d_map_[composeIndex(xi-1,yi,zi+1)]-d_map_[composeIndex(xi+1,yi+2,zi-1)]+d_map_[composeIndex(xi-1,yi+2,zi-1)]+d_map_[composeIndex(xi+1,yi,zi-1)]-d_map_[composeIndex(xi-1,yi,zi-1)]),
      0.125*(d_map_[composeIndex(xi+2,yi+2,zi+1)]-d_map_[composeIndex(xi,yi+2,zi+1)]-d_map_[composeIndex(xi+2,yi,zi+1)]+d_map_[composeIndex(xi,yi,zi+1)]-d_map_[composeIndex(xi+2,yi+2,zi-1)]+d_map_[composeIndex(xi,yi+2,zi-1)]+d_map_[composeIndex(xi+2,yi,zi-1)]-d_map_[composeIndex(xi,yi,zi-1)]),
      0.125*(d_map_[composeIndex(xi+1,yi+1,zi+2)]-d_map_[composeIndex(xi-1,yi+1,zi+2)]-d_map_[composeIndex(xi+1,yi-1,zi+2)]+d_map_[composeIndex(xi-1,yi-1,zi+2)]-d_map_[composeIndex(xi+1,yi+1,zi)]+d_map_[composeIndex(xi-1,yi+1,zi)]+d_map_[composeIndex(xi+1,yi-1,zi)]-d_map_[composeIndex(xi-1,yi-1,zi)]),
      0.125*(d_map_[composeIndex(xi+2,yi+1,zi+2)]-d_map_[composeIndex(xi,yi+1,zi+2)]-d_map_[composeIndex(xi+2,yi-1,zi+2)]+d_map_[composeIndex(xi,yi-1,zi+2)]-d_map_[composeIndex(xi+2,yi+1,zi)]+d_map_[composeIndex(xi,yi+1,zi)]+d_map_[composeIndex(xi+2,yi-1,zi)]-d_map_[composeIndex(xi,yi-1,zi)]),
      0.125*(d_map_[composeIndex(xi+1,yi+2,zi+2)]-d_map_[composeIndex(xi-1,yi+2,zi+2)]-d_map_[composeIndex(xi+1,yi,zi+2)]+d_map_[composeIndex(xi-1,yi,zi+2)]-d_map_[composeIndex(xi+1,yi+2,zi)]+d_map_[composeIndex(xi-1,yi+2,zi)]+d_map_[composeIndex(xi+1,yi,zi)]-d_map_[composeIndex(xi-1,yi,zi)]),
      0.125*(d_map_[composeIndex(xi+2,yi+2,zi+2)]-d_map_[composeIndex(xi,yi+2,zi+2)]-d_map_[composeIndex(xi+2,yi,zi+2)]+d_map_[composeIndex(xi,yi,zi+2)]-d_map_[composeIndex(xi+2,yi+2,zi)]+d_map_[composeIndex(xi,yi+2,zi)]+d_map_[composeIndex(xi+2,yi,zi)]-d_map_[composeIndex(xi,yi,zi)])
    ;
    // Convert voxel values and partial derivatives to interpolation coefficients.
    _coefs = _C * x;
    // Remember this voxel for next time.
    _i1 = xi;
    _i2 = yi;
    _i3 = zi;
    _initialized = true;
  }
  // Evaluate the interpolation within this grid voxel.
  dx -= xi;
  dy -= yi;
  dz -= zi;
  int ijkn(0);
  fptype dzpow(1);
  fptype result(0);
  for(int k = 0; k < 4; ++k) {
    fptype dypow(1);
    for(int j = 0; j < 4; ++j) {
      result += dypow*dzpow*(_coefs[ijkn] + dx*(_coefs[ijkn+1] + dx*(_coefs[ijkn+2] + dx*_coefs[ijkn+3])));
      ijkn += 4;
      dypow *= dy;
    }
    dzpow *= dz;
  }
  return result;
}
*/
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
        dd_map_.resize(map_size);
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

//        std::cout << "ReachabilityMap::createDistanceMap: distance map size: " << d_map_.size() << std::endl;

        for (int idx = 0; idx < d_map_.size(); idx++) {
            d_map_[idx] = -1.0;
        }

        if (getIndex(origin) < 0) {
            std::cout << "ReachabilityMap::createDistanceMap: getIndex(origin) < 0" << std::endl;
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

        std::set<int > obstacle_ids;
        for (int idx = 0; idx < d_map_.size(); idx++) {
            if (d_map_[idx] < 0.0) {
                obstacle_ids.insert(idx);
            }
        }

        while (!obstacle_ids.empty()) {
            std::set<std::pair<int, double> > neighbouring_ids;
            for (std::set<int>::const_iterator it = obstacle_ids.begin(); it != obstacle_ids.end(); it++) {
                int ix, iy, iz;
                int idx = (*it);
                decomposeIndex(idx, ix, iy, iz);

                double max_value = -0.01;
                for (int iix = std::max(0,ix-1); iix < std::min(steps_[0], ix+2); iix++) {
                    for (int iiy = std::max(0,iy-1); iiy < std::min(steps_[1], iy+2); iiy++) {
                        for (int iiz = std::max(0,iz-1); iiz < std::min(steps_[2], iz+2); iiz++) {
                            if (ix == iix || iy == iiy || iy == iiy) {
                                continue;
                            }
                            int pt_idx = composeIndex(iix, iiy, iiz);
                            double pt_val = d_map_[pt_idx];
                            if (pt_val >= 0.0 && pt_val > max_value) {
                                max_value = pt_val;
                            }
                        }
                    }
                }
                if (max_value >= 0.0) {
                    neighbouring_ids.insert( std::make_pair(idx, max_value) );
                }
            }

            for (std::set<std::pair<int, double> >::const_iterator it = neighbouring_ids.begin(); it != neighbouring_ids.end(); it++) {
                int idx = it->first;
                double max_value = it->second;
                d_map_[idx] = max_value + voxel_size_;
                obstacle_ids.erase(idx);
            }
            neighbouring_ids.clear();
        }

        return true;

        for (int idx = 0; idx < d_map_.size(); idx++) {
            if (d_map_[idx] < 0.0) {
                double max_value = 0.0;
                for (int iix = std::max(0,ix-1); iix < std::min(steps_[0], ix+2); iix++) {
                    for (int iiy = std::max(0,iy-1); iiy < std::min(steps_[1], iy+2); iiy++) {
                        for (int iiz = std::max(0,iz-1); iiz < std::min(steps_[2], iz+2); iiz++) {
                            if (ix == iix || iy == iiy || iy == iiy) {
                                continue;
                            }
                            int pt_idx = composeIndex(iix, iiy, iiz);
                            double pt_val = d_map_[pt_idx];
                            if (pt_val >= 0.0 && pt_val > max_value) {
                                max_value = pt_val;
                            }
                        }
                    }
                }


                d_map_[idx] = max_value + voxel_size_;
                continue;
            }

            continue;

            int ix, iy, iz;
            decomposeIndex(idx, ix, iy, iz);

            if (ix == 0 || ix >= steps_[0] - 1) {
                continue;
            }
            if (iy == 0 || iy >= steps_[1] - 1) {
                continue;
            }
            if (iz == 0 || iz >= steps_[2] - 1) {
                continue;
            }

            double min_value = d_map_[idx];
            KDL::Vector min_gr;
            for (int iix = std::max(0,ix-1); iix < std::min(steps_[0], ix+2); iix++) {
                for (int iiy = std::max(0,iy-1); iiy < std::min(steps_[1], iy+2); iiy++) {
                    for (int iiz = std::max(0,iz-1); iiz < std::min(steps_[2], iz+2); iiz++) {
                        if (ix == iix || iy == iiy || iy == iiy) {
                            continue;
                        }
                        int pt_idx = composeIndex(iix, iiy, iiz);
                        double pt_val = d_map_[pt_idx];
                        if (pt_val >= 0.0 && pt_val < min_value) {
                            min_value = pt_val;
                            min_gr = KDL::Vector(iix - ix, iiy - iy, iiz - iz);
                        }
                    }
                }
            }
            min_gr.Normalize();
            dd_map_[idx].dx_ = min_gr.x();
            dd_map_[idx].dy_ = min_gr.y();
            dd_map_[idx].dz_ = min_gr.z();
            continue;



            double deriv[16] = {
                      // p n
                0,    // 0 0
                -1,   // 0 1
                -1,   // 0 2
                1,    // 0 3
                1,    // 1 0
                0,    // 1 1
                1,    // 1 2
                1,    // 1 3
                1,    // 2 0
                -1,   // 2 1
                0,    // 2 2
                1,    // 2 3
                -1,   // 3 0
                -1,   // 3 1
                -1,   // 3 2
                0,    // 3 3
            };

/*
            double deriv[16][2] = {
                          // np
                {0,0},    // 00
                {1,0},    // 01
                {0,0},    // 02
                {-1,0},   // 03
                {-1,0},   // 10
                {0,-2},   // 11
                {0,-1},   // 12
                {-1,0},   // 13
                {0,0},    // 20
                {0,-1},   // 21
                {0,0},    // 22
                {0,1},    // 23
                {1,0},    // 30
                {1,0},    // 31
                {0,1},    // 32
                {0,2},    // 33
            };

            double deriv[16][2] = {
                          // np
                {0,0},    // 00
                {0,1},    // 01
                {0,0},    // 02
                {0,-1},   // 03
                {0,-1},   // 10
                {-2,0},   // 11
                {-1,0},   // 12
                {0,-1},   // 13
                {0,0},    // 20
                {-1,0},   // 21
                {0,0},    // 22
                {1,0},    // 23
                {0,1},    // 30
                {0,1},    // 31
                {1,0},    // 32
                {2,0},    // 33
            };
//*/
            int idxp, idxn;
            idxp = composeIndex(ix-1, iy, iz);
            idxn = composeIndex(ix+1, iy, iz);
            // calculate partial derivatives for x
            int vidx = 0;
            if (ix == 0 || d_map_[idxp] < 0.0) {
                vidx = 0;
            }
            else if (d_map_[idxp] < d_map_[idx]) {
                vidx = 1 * 4;
            }
            else if (d_map_[idxp] == d_map_[idx]) {
                vidx = 2 * 4;
            }
            else {
                vidx = 3 * 4;
            }
            vidx = vidx;
            if (ix == steps_[0]-1 || d_map_[idxn] < 0.0) {
                vidx |= 0;
            }
            else if (d_map_[idxn] < d_map_[idx]) {
                vidx |= 1;
            }
            else if (d_map_[idxn] == d_map_[idx]) {
                vidx |= 2;
            }
            else {
                vidx |= 3;
            }

            dd_map_[idx].dx_ = deriv[vidx];

            idxp = composeIndex(ix, iy-1, iz);
            idxn = composeIndex(ix, iy+1, iz);

            // calculate partial derivatives for y
            vidx = 0;
            if (iy == 0 || d_map_[idxp] < 0.0) {
                vidx = 0;
            }
            else if (d_map_[idxp] < d_map_[idx]) {
                vidx = 1 * 4;
            }
            else if (d_map_[idxp] == d_map_[idx]) {
                vidx = 2 * 4;
            }
            else {
                vidx = 3 * 4;
            }
            vidx = vidx;
            if (iy == steps_[1]-1 || d_map_[idxn] < 0.0) {
                vidx |= 0;
            }
            else if (d_map_[idxn] < d_map_[idx]) {
                vidx |= 1;
            }
            else if (d_map_[idxn] == d_map_[idx]) {
                vidx |= 2;
            }
            else {
                vidx |= 3;
            }

            dd_map_[idx].dy_ = deriv[vidx];

            idxp = composeIndex(ix, iy, iz-1);
            idxn = composeIndex(ix, iy, iz+1);

            // calculate partial derivatives for z
            vidx = 0;
            if (iz == 0 || d_map_[idxp] < 0.0) {
                vidx = 0;
            }
            else if (d_map_[idxp] < d_map_[idx]) {
                vidx = 1 * 4;
            }
            else if (d_map_[idxp] == d_map_[idx]) {
                vidx = 2 * 4;
            }
            else {
                vidx = 3 * 4;
            }
            vidx = vidx;
            if (iz == steps_[2]-1 || d_map_[idxn] < 0.0) {
                vidx |= 0;
            }
            else if (d_map_[idxp] < d_map_[idx]) {
                vidx |= 1;
            }
            else if (d_map_[idxp] == d_map_[idx]) {
                vidx |= 2;
            }
            else {
                vidx |= 3;
            }

            dd_map_[idx].dz_ = deriv[vidx];

/*

            idxp = composeIndex(ix-1, iy, iz);
            idxn = composeIndex(ix+1, iy, iz);
            if ((d_map_[idxp] >= d_map_[idx] && d_map_[idxn] >= d_map_[idx]) || 
                (d_map_[idxp] <= d_map_[idx] && d_map_[idxn] <= d_map_[idx])) {
                dd_map_[idx].dx_ = 0.0;
            }
            else if (d_map_[idxp] >= d_map_[idx]) {
                dd_map_[idx].dx_ = 1.0;
            }
            else {
                dd_map_[idx].dx_ = -1.0;
            }

            idxp = composeIndex(ix, iy-1, iz);
            idxn = composeIndex(ix, iy+1, iz);
            if ((d_map_[idxp] >= d_map_[idx] && d_map_[idxn] >= d_map_[idx]) || 
                (d_map_[idxp] <= d_map_[idx] && d_map_[idxn] <= d_map_[idx])) {
                dd_map_[idx].dy_ = 0.0;
            }
            else if (d_map_[idxp] >= d_map_[idx]) {
                dd_map_[idx].dy_ = 1.0;
            }
            else {
                dd_map_[idx].dy_ = -1.0;
            }

            idxp = composeIndex(ix, iy, iz-1);
            idxn = composeIndex(ix, iy, iz+1);
            if ((d_map_[idxp] >= d_map_[idx] && d_map_[idxn] >= d_map_[idx]) || 
                (d_map_[idxp] <= d_map_[idx] && d_map_[idxn] <= d_map_[idx])) {
                dd_map_[idx].dz_ = 0.0;
            }
            else if (d_map_[idxp] >= d_map_[idx]) {
                dd_map_[idx].dz_ = 1.0;
            }
            else {
                dd_map_[idx].dz_ = -1.0;
            }
*/
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
//        std::cout << "bounduary_set_.size() " << bounduary_set_.size() << std::endl;
/*
        // increase the value for the bounduary points
        // THIS DOES NOT WORK!
        for (std::set<int >::const_iterator it = bounduary_set_.begin(); it != bounduary_set_.end(); it++) {
            int idx = (*it);
//            d_map_[idx] = -1;
//            continue;
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
        return true;
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

    bool ReachabilityMap::getGradient(int idx, KDL::Vector &gradient) const {
        double min_value = d_map_[idx];

//        int ix = getIndexDim(x.x(), 0);
//        int iy = getIndexDim(x.y(), 1);
//        int iz = getIndexDim(x.z(), 2);
        int ix, iy, iz;
        decomposeIndex(idx, ix, iy, iz);

        bool obstacle = false;
        int search_space = 1;
        if (min_value < 0.0) {
            // we are in the obstacle
            min_value = 100000.0;
            obstacle = true;
//            std::cout << "obst. " << pt[0] << " " << pt[1] << " " << pt[2] << " " << d_map_[idx] << std::endl;
            return false;
//            std::cout << "obstacle" << std::endl;
        }
//        int min_ix=-1, min_iy=-1, min_iz=-1;
        KDL::Vector min_vec;
        bool found = false;
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
                        min_vec = KDL::Vector(iix - ix, iiy - iy, iiz - iz);
                        found = true;
//                            min_ix = iix;
//                            min_iy = iiy;
//                            min_iz = iiz;
                    }
                }
            }
        }

        if (!found) {// min_ix == -1) {
//            std::cout << "ReachabilityMap::getGradient: could not find gradient, min_value: " << min_value << std::endl;
//            std::cout << "no min " << pt[0] << " " << pt[1] << " " << pt[2] << " " << d_map_[idx] << std::endl;
            return false;
        }
        min_vec.Normalize();
        gradient = min_vec;
        return true;
/*
        KDL::Vector min_pt;
        KDL::Vector pt;
        getIndexCenter(ix, iy, iz, pt);
        getIndexCenter(min_ix, min_iy, min_iz, min_pt);
//        std::cout << " ok " << pt[0] << " " << pt[1] << " " << pt[2] << " " << d_map_[idx] << " to " << min_pt[0] << " " << min_pt[1] << " " << min_pt[2] << " " << min_value << std::endl;
        gradient = min_pt - pt;
        gradient.Normalize();
*/
        return true;
    }

    bool ReachabilityMap::getGradient(const KDL::Vector &x, KDL::Vector &gradient) const {

        int ix0 = std::floor((x.x() - ep_min_(0)) / voxel_size_);
        int iy0 = std::floor((x.y() - ep_min_(1)) / voxel_size_);
        int iz0 = std::floor((x.z() - ep_min_(2)) / voxel_size_);
        int ix1 = ix0+1;
        int iy1 = iy0+1;
        int iz1 = iz0+1;

        if (ix0 < 1 || iy0 < 1 || iz0 < 1 || ix1 >= steps_[0]-2 || iy1 >= steps_[1]-2 || iz1 >= steps_[2]-2) {
            return false;
        }

        double x0 = ep_min_(0) + ix0 * voxel_size_;
        double y0 = ep_min_(1) + iy0 * voxel_size_;
        double z0 = ep_min_(2) + iz0 * voxel_size_;
        double x1 = ep_min_(0) + ix1 * voxel_size_;
        double y1 = ep_min_(1) + iy1 * voxel_size_;
        double z1 = ep_min_(2) + iz1 * voxel_size_;

/*        int indices[8] = {
            composeIndex(ix0, iy0, iz0),
            composeIndex(ix1, iy0, iz0),
            composeIndex(ix0, iy1, iz0),
            composeIndex(ix1, iy1, iz0),
            composeIndex(ix0, iy0, iz1),
            composeIndex(ix1, iy0, iz1),
            composeIndex(ix0, iy1, iz1),
            composeIndex(ix1, iy1, iz1),
        };
*/
        int indices[8] = {
            composeIndex(ix0, iy0, iz0),
            composeIndex(ix0, iy0, iz1),
            composeIndex(ix0, iy1, iz0),
            composeIndex(ix0, iy1, iz1),
            composeIndex(ix1, iy0, iz0),
            composeIndex(ix1, iy0, iz1),
            composeIndex(ix1, iy1, iz0),
            composeIndex(ix1, iy1, iz1),
        };

        double a[64];
        double f[8], dfdx[8], dfdy[8], dfdz[8], d2fdxdy[8], d2fdxdz[8], d2fdydz[8], d3fdxdydz[8];
        for (int i = 0; i < 8; i++) {
            int nidx = indices[i];
            f[i] = d_map_[nidx];
            dfdx[i] = -dd_map_[nidx].dx_ * voxel_size_;
            dfdy[i] = -dd_map_[nidx].dy_ * voxel_size_;
            dfdz[i] = -dd_map_[nidx].dz_ * voxel_size_;
            d2fdxdy[i] = 0.0;//dd_map_[nidx].dx_ * dd_map_[nidx].dy_;// * voxel_size_ * voxel_size_;
            d2fdxdz[i] = 0.0;//dd_map_[nidx].dx_ * dd_map_[nidx].dz_;// * voxel_size_ * voxel_size_;
            d2fdydz[i] = 0.0;//dd_map_[nidx].dy_ * dd_map_[nidx].dz_;// * voxel_size_ * voxel_size_;
            d3fdxdydz[i] = 0.0;//dd_map_[nidx].dx_ * dd_map_[nidx].dy_ * dd_map_[nidx].dz_;// * voxel_size_ * voxel_size_ * voxel_size_;
        }

        tricubic_get_coeff(a, ix0, iy0, iz0);

        //tricubic_get_coeff(a, f, dfdx, dfdy, dfdz, d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz);
//std::cout << (x.x() - x0) / voxel_size_ << " " << (x.y() - y0) / voxel_size_ << " " << (x.z() - z0) / voxel_size_ << std::endl;
        double dx, dy, dz;
        dx = tricubic_eval(a, (x.x() - x0) / voxel_size_, (x.y() - y0) / voxel_size_, (x.z() - z0) / voxel_size_, 1, 0, 0)/ voxel_size_;
        dy = tricubic_eval(a, (x.x() - x0) / voxel_size_, (x.y() - y0) / voxel_size_, (x.z() - z0) / voxel_size_, 0, 1, 0)/ voxel_size_;
        dz = tricubic_eval(a, (x.x() - x0) / voxel_size_, (x.y() - y0) / voxel_size_, (x.z() - z0) / voxel_size_, 0, 0, 1)/ voxel_size_;
        gradient = -KDL::Vector(dx, dy, dz);
        gradient.Normalize();

//        gradient = KDL::Vector(dd_map_[indices[0]].dx_, dd_map_[indices[0]].dy_, dd_map_[indices[0]].dz_);
//        gradient.Normalize();

//        gradient = KDL::Vector(0, 4.0 * tricubic_eval(a, (x.x() - x0) / voxel_size_, (x.y() - y0) / voxel_size_, (x.z() - z0) / voxel_size_, 0, 0, 0), 0);

        return true;


        bool valid_gradients[8];
        KDL::Vector gradients[8];

        double max_value = -0.0001;
//        KDL::Vector mean_vec = 0.0;
        for (int i = 0; i < 8; i++) {
            if (getGradient(indices[i], gradients[i])) {
//                mean_vec += 
            }

            if (d_map_[indices[i]] > max_value) {
                max_value = d_map_[indices[i]];
            }
        }

//        if (d_map_[n_idx] < 0.0) {
//            return false;
//        }

        KDL::Vector V[8];
        int obst_count = 0;
        for (int i = 0; i < 8; i++) {
            if (!getGradient(indices[i], V[i])) {
                V[i] = KDL::Vector();
                obst_count++;
            }
/*            if (d_map_[indices[i]] < 0.0) {
                V[i] = max_value + voxel_size_;
            }
            else {
                V[i] = d_map_[indices[i]];
            }
*/
        }
        if (obst_count == 8) {
            return false;
        }
/*
        double x0 = ep_min_(0) + ix0 * voxel_size_;
        double y0 = ep_min_(1) + iy0 * voxel_size_;
        double z0 = ep_min_(2) + iz0 * voxel_size_;
        double x1 = ep_min_(0) + ix1 * voxel_size_;
        double y1 = ep_min_(1) + iy1 * voxel_size_;
        double z1 = ep_min_(2) + iz1 * voxel_size_;
*/
        double xd = (x.x() - x0) / (x1 - x0);
        double yd = (x.y() - y0) / (y1 - y0);
        double zd = (x.z() - z0) / (z1 - z0);

        KDL::Vector c00 = V[0] * (1.0 - xd) + V[1] * xd;
        KDL::Vector c10 = V[2] * (1.0 - xd) + V[3] * xd;
        KDL::Vector c01 = V[4] * (1.0 - xd) + V[5] * xd;
        KDL::Vector c11 = V[6] * (1.0 - xd) + V[7] * xd;

        KDL::Vector c0 = c00 * (1.0 - yd) + c10 * yd;
        KDL::Vector c1 = c01 * (1.0 - yd) + c11 * yd;

        KDL::Vector c = c0 * (1.0 - zd) + c1 * zd;
/*
c = ((V[0] * (1.0 - xd) + V[1] * xd) * (1.0 - yd) + (V[2] * (1.0 - xd) + V[3] * xd) * yd) * (1.0 - zd) + ((V[4] * (1.0 - xd) + V[5] * xd) * (1.0 - yd) + (V[6] * (1.0 - xd) + V[7] * xd) * yd) * zd;


c = ((V[0] * (1.0 - xd) + V[1] * xd) * (1.0 - yd) + (V[2] * (1.0 - xd) + V[3] * xd) * yd) * (1.0 - zd) + ((V[4] * (1.0 - xd) + V[5] * xd) * (1.0 - yd) + (V[6] * (1.0 - xd) + V[7] * xd) * yd) * zd;



dc/dx
*/
//        c.Normalize();
        gradient = c;
        return true;

/*
        int idx = getIndex(x);
        if (idx < 0) {
//            std::cout << "ReachabilityMap::getGradient: point is outside the map" << std::endl;
            return false;
        }

        double min_value = d_map_[idx];

        int ix = getIndexDim(x.x(), 0);
        int iy = getIndexDim(x.y(), 1);
        int iz = getIndexDim(x.z(), 2);
        KDL::Vector pt;
        getIndexCenter(ix, iy, iz, pt);

        bool obstacle = false;
        int search_space = 1;
        if (min_value < 0.0) {
            // we are in the obstacle
            min_value = 100000.0;
            obstacle = true;
//            std::cout << "obst. " << pt[0] << " " << pt[1] << " " << pt[2] << " " << d_map_[idx] << std::endl;
//            return false;
//            std::cout << "obstacle" << std::endl;
        }
        int min_ix=-1, min_iy=-1, min_iz=-1;

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
//            std::cout << "ReachabilityMap::getGradient: could not find gradient, min_value: " << min_value << std::endl;
//            std::cout << "no min " << pt[0] << " " << pt[1] << " " << pt[2] << " " << d_map_[idx] << std::endl;
            return false;
        }
        KDL::Vector min_pt;
        getIndexCenter(min_ix, min_iy, min_iz, min_pt);
//        std::cout << " ok " << pt[0] << " " << pt[1] << " " << pt[2] << " " << d_map_[idx] << " to " << min_pt[0] << " " << min_pt[1] << " " << min_pt[2] << " " << min_value << std::endl;
        gradient = min_pt - pt;
        gradient.Normalize();
        return true;
*/
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

    void ReachabilityMap::printDistanceMap() const {
        std::cout << steps_[0] << " " << steps_[1] << " " << steps_[2] << std::endl;
        for (int i = 0; i < d_map_.size(); i++) {
            std::cout << d_map_[i] << " ";
        }
        std::cout << std::endl;
    }

