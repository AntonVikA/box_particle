//
// Created by anton on 11/13/22.
//

#ifndef PARTICLE_IN_SQUARE_BOX_INTERPOLATION2D_H
#define PARTICLE_IN_SQUARE_BOX_INTERPOLATION2D_H

#include <vector>
#include <map>
#include <iostream>


class Interpolation2D {
    double _L = 0;
    std::size_t _n = 0;
    double h = 0;
    std::size_t i_ind = 0, j_ind = 0;
    std::map<double, std::size_t> axis;
public:
    Interpolation2D (std::size_t n, double L) : _n{n}, _L{L} {
        h = L / (_n - 1);
        for (std::size_t i = 0; i < _n; ++i) {
            axis[i * h] = i;
        }
    }
    double operator() (const std::vector<std::vector<double>>& values, double x, double y);
};


#endif //PARTICLE_IN_SQUARE_BOX_INTERPOLATION2D_H
