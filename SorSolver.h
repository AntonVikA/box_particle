//
// Created by anton on 11/16/22.
//

#ifndef PARTICLE_IN_SQUARE_BOX_SORSOLVER_H
#define PARTICLE_IN_SQUARE_BOX_SORSOLVER_H

#include <cmath>
#include "DataType.h"
#include <iostream>

class SorSolver {
    size_t N = 0;
    double r_optimal = 0;
    const double e = 0.000001;
    double _dx;
public:
    SorSolver() = default;
    void init(boost_matrix& A, const boost_matrix& W, double dx);
};


#endif //PARTICLE_IN_SQUARE_BOX_SORSOLVER_H
