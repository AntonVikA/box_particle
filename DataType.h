//
// Created by anton on 11/12/22.
//

#ifndef PARTICLE_IN_SQUARE_BOX_DATATYPE_H
#define PARTICLE_IN_SQUARE_BOX_DATATYPE_H
#include <vector>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/io.hpp>

const double PI = 4.0 * std::atan(1.0);

using boost_vector = boost::numeric::ublas::vector<double>;
using boost_matrix = boost::numeric::ublas::matrix<double>;
using zero_matrix = boost::numeric::ublas::zero_matrix<double>;

enum class Diff {
    d_x,
    d_xx,
    d_y,
    d_yy
};

#endif //PARTICLE_IN_SQUARE_BOX_DATATYPE_H
