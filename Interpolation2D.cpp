//
// Created by anton on 11/13/22.
//

#include "Interpolation2D.h"


double Interpolation2D::operator()(const std::vector<std::vector<double>>& values, double x, double y) {
    i_ind = axis.lower_bound(x)->second;
    j_ind = axis.lower_bound(y)->second;

    //std::cout << i_ind << " " << j_ind << std::endl;

    double x1 = i_ind * h;
    double y1 = j_ind * h;

    double x2 = (i_ind + 1) * h;
    double y2 = (j_ind + 1) * h;

    double C = 1.0 / ((x2 - x1) * (y2 - y1));
    double C_11 = values[i_ind][j_ind] * (x2 - x) * (y2 - y);
    double C_21 = values[i_ind + 1][j_ind] * (x - x1) * (y2 - y);
    double C_12 = values[i_ind][j_ind + 1] * (x2 - x) * (y - y1);
    double C_22 = values[i_ind + 1][j_ind + 1] * (x - x1) * (y - y1);

    return C * (C_11 + C_21 + C_12 + C_22);
}