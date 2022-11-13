#include <iostream>
#include <vector>
#include <map>
#include "RKM.h"
#include "DataType.h"
#include <cmath>
#include "Interpolation2D.h"

constexpr double PI = 4.0 * std::atan(1.0);

const double C = 15.1, Re = 10.0, L = 50;
size_t n = 25;
Interpolation2D interpolator(n, L);
std::vector<std::vector<double>> U_flow = std::vector<std::vector<double>>(n, std::vector<double>(n, 0)),
                                 V_flow = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));

const std::map<std::string, size_t> var_map = {{"x", 0},
                                                {"u", 1},
                                                {"y", 2},
                                                {"v", 3}};
double var(const std::string& var_name, const std::vector<double>& data) {
    return data[var_map.at(var_name)];
}

const double xw = 20.0,
             yw = 20.0;

double U(double x, double y) {
    return 1.0 * (y - yw) / std::sqrt(std::pow(x - xw, 2.0) + std::pow(y - yw, 2.0)); //1.0;
}

double V(double x, double y) {
    return -1.0 * (x - xw) / std::sqrt(std::pow(x - xw, 2.0) + std::pow(y - yw, 2.0));//y > 0 ? -1.0 : 1.0;
}

double U_desc(double x, double y) {
    return interpolator(U_flow, x, y);
}

double V_desc(double x, double y) {
    return interpolator(V_flow, x, y);
}

double f_dx(double t, const std::vector<double>& y) {
    double x_p = y[0];
    double u = y[1];
    double y_p = y[2];
    double v = y[3];

    return u;
}

double f_du(double t, const std::vector<double>& y) {
    double x_p = y[0];
    double u = y[1];
    double y_p = y[2];
    double v = y[3];

    return C * std::sqrt(std::pow(u - U_desc(x_p, y_p), 2.0) +
                         std::pow(v - V_desc(x_p, y_p), 2.0)) *
                            (U_desc(x_p, y_p) - u);
}

double f_dy(double t, const std::vector<double>& y) {
    double x_p = y[0];
    double u = y[1];
    double y_p = y[2];
    double v = y[3];

    return v;
}

double f_dv(double t, const  std::vector<double>& y) {
    double x_p = y[0];
    double u = y[1];
    double y_p = y[2];
    double v = y[3];

    return C * std::sqrt(std::pow(u - U_desc(x_p, y_p), 2.0) +
                         std::pow(v - V_desc(x_p, y_p), 2.0)) *
                            (V_desc(x_p, y_p) - v);
}


int main() {

    const std::vector<std::function<double(double, const std::vector<double>&)>> functions = {f_dx, f_du, f_dy, f_dv};
    std::vector<double> initial_condition{25, 0.0, 25, 0.0};

    double t_begin = 0.0;
    double t_end = 135.00;
    double th = 0.001;

    RKM rkm;

    rkm.init(functions);

    double h = L / (n - 1);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double x = i * h;
            double y = j * h;
            U_flow[i][j] = U(x, y);
            V_flow[i][j] = V(x, y);
        }
    }

    size_t iter = 0;
    rkm.solve(t_begin, t_end, initial_condition, th);


    return 0;
}