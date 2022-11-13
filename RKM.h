//
// Created by anton on 10/15/22.
//
#ifndef COLD_SPRAY_PROJECT_RKM_H
#define COLD_SPRAY_PROJECT_RKM_H
#include <cassert>
#include <functional>
#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include <fstream>
#include <ios>


using function_type = double(double, const std::vector<double>&);

class RKM {
    size_t number_of_equations = 0;
    double t_ = 0;
    std::vector<double> K1, K2, K3, K4;
    std::vector<double> Y, dY;
    std::vector<std::function<function_type>> functions_;
    std::vector<double> operator() (double t, const std::vector<double>& y);
    std::ofstream file;
    void check_boundary_collision();
    const double L = 25;

public:
    RKM() = default;
    RKM(const std::vector<std::function<function_type>>& functions);
    void init(const std::vector<std::function<function_type>>& functions);
    void solve(double t_begin, double t_end, const std::vector<double>& init_conditions, double h);
};


#endif //COLD_SPRAY_PROJECT_RKM_H
