#include <iostream>
#include <vector>
#include <map>
#include "RKM.h"
#include "DataType.h"
#include <cmath>
#include <fstream>
#include "Interpolation2D.h"
#include "SorSolver.h"


constexpr double C = 15.1, Re = 10.0, L = 1;
const size_t n = 21;
constexpr double dx = 1.0 / (n - 1);
constexpr double dt = 0.5 * (1.0 / 4.0 * Re * dx * dx);

Interpolation2D interpolator(n, L);
std::vector<std::vector<double>> U_flow = std::vector<std::vector<double>>(n, std::vector<double>(n, 0)),
                                 V_flow = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));

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

double f_dx(double t, const boost_vector& y) {
    double x_p = y[0];
    double u = y[1];
    double y_p = y[2];
    double v = y[3];

    return u;
}

double f_du(double t, const boost_vector& y) {
    double x_p = y[0];
    double u = y[1];
    double y_p = y[2];
    double v = y[3];

    return C * std::sqrt(std::pow(u - U_desc(x_p, y_p), 2.0) +
                         std::pow(v - V_desc(x_p, y_p), 2.0)) *
                            (U_desc(x_p, y_p) - u);
}

double f_dy(double t, const boost_vector& y) {
    double x_p = y[0];
    double u = y[1];
    double y_p = y[2];
    double v = y[3];

    return v;
}

double f_dv(double t, const  boost_vector& y) {
    double x_p = y[0];
    double u = y[1];
    double y_p = y[2];
    double v = y[3];

    return C * std::sqrt(std::pow(u - U_desc(x_p, y_p), 2.0) +
                         std::pow(v - V_desc(x_p, y_p), 2.0)) *
                            (V_desc(x_p, y_p) - v);
}


double diff(const boost_matrix& A, size_t i, size_t j, const Diff& type) {
    switch(type) {
        case Diff::d_x:
        {
            return (A(i + 1, j) - A(i - 1, j)) / (2.0 * dx);
        } break;
        case Diff::d_y:
        {
            return (A(i, j + 1) - A(i, j - 1)) / (2.0 * dx);
        } break;
        case Diff::d_xx:
        {
            return (A(i + 1, j) -2.0 * A(i, j) + A(i - 1, j)) / (dx * dx);
        } break;
        case Diff::d_yy:
        {
            return (A(i, j + 1) -2.0 * A(i, j) + A(i, j - 1)) / (dx * dx);
        } break;
    }
}

void get_vorticity(const boost_matrix& u, const boost_matrix& v, boost_matrix& w) {
    for (size_t i = 1; i < n - 1; ++i) {
        for (size_t j = 1; j < n - 1; ++j) {
            double dv_dx = diff(v, i, j, Diff::d_x);
            double du_dy = diff(u, i, j, Diff::d_y);

            w(i, j) = dv_dx - du_dy;
        }
    }
}

void fill_boundary_vorticity(boost_matrix& w, const boost_matrix& psi, const boost_matrix& u, const boost_matrix& v) {
    // boundary vorticity
    //upper
    for (size_t i = 0; i < n; ++i) {
        size_t j = n - 1;
        double u_half = (psi(i, j) - psi(i, j - 1)) / dx;
        double w_quarter = -(u(i, j) - u_half) / (dx * 1.0 / 2.0);
        w(i, j) = (4.0 * w_quarter - w(i, j - 1) / 3);
    }
    //bottom
    for (size_t i = 0; i < n; ++i) {
        size_t j = 0;
        double u_half = (psi(i, j + 1) - psi(i, j)) / dx;
        double w_quarter = -(u_half - u(i, j)) / (dx * 1.0 / 2.0);
        w(i, j) = (4.0 * w_quarter - w(i, j + 1)) / 3;
    }

    //left
    for (size_t j = 0; j < n; ++j) {
        size_t i = 0;
        double v_half = -(psi(i + 1, j) - psi(i, j)) / dx;
        double w_quarter = (v_half - v(i, j)) / (dx * 1.0 / 2.0);
        w(i, j) = (4.0 * w_quarter - w(i + 1, j)) / 3;
    }

    //right
    for (size_t j = 0; j < n; ++j) {
        size_t i = n - 1;
        double v_half = -(psi(i, j) - psi(i - 1, j)) / dx;
        double w_quarter = (v(i, j) - v_half) / (dx * 1.0 / 2.0);
        w(i, j) = (4.0 * w_quarter - w(i - 1, j)) / 3;
    }
}

void get_v_from_psi(const boost_matrix& psi, boost_matrix& u, boost_matrix& v) {
    for (size_t i = 1; i < n - 1; ++i) {
        for (size_t j = 1; j < n - 1; ++j) {
            double dpsi_dy = diff(psi, i, j, Diff::d_y);
            double dpsi_dx = diff(psi, i, j, Diff::d_x);

            u(i, j) = dpsi_dy;
            v(i, j) = -dpsi_dx;
        }
    }
}

void right (boost_matrix& A, std::function<double(double, double)> boundary_function) {

    size_t i = n - 1;
    for (size_t j = 0; j < n; ++j) {
        double x = i * dx;
        double y = j * dx;
        A(i, j) = boundary_function(x, y);
    }
}

void upper (boost_matrix& A, std::function<double(double, double)> boundary_function) {

    size_t j = n - 1;
    for (size_t i = 0; i < n; ++i) {
        double x = i * dx;
        double y = j * dx;
        A(i, j) = boundary_function(x, y);
    }
}

void left (boost_matrix& A, std::function<double(double, double)> boundary_function) {

    size_t i = 0;
    for (size_t j = 0; j < n; ++j) {
        double x = i * dx;
        double y = j * dx;
        A(i, j) = boundary_function(x, y);
    }
}

void bottom (boost_matrix& A, std::function<double(double, double)> boundary_function) {

    size_t j = 0;
    for (size_t i = 0; i < n; ++i) {
        double x = i * dx;
        double y = j * dx;
        A(i, j) = boundary_function(x, y);
    }
}

boost_matrix vorticity_evolution(const boost_matrix& w, const boost_matrix psi) {
    boost_matrix w_new = zero_matrix(n, n);

    for (size_t i = 1; i < n - 1; ++i) {
        for (size_t j = 1; j < n - 1; ++j) {

            double w_y = diff(w, i, j, Diff::d_y);
            double psi_x = diff(psi, i, j, Diff::d_x);
            double w_x = diff(w, i, j, Diff::d_x);
            double psi_y = diff(psi, i, j, Diff::d_y);

            double w_xx = diff(w, i, j, Diff::d_xx);
            double w_yy = diff(w, i, j, Diff::d_yy);

            double jacob = (w_x * psi_y - w_y * psi_x);
            double laplace_w = (w_xx + w_yy);


            w_new(i, j) = w(i, j) + dt * (-jacob + 1.0 / Re * laplace_w);
        }
    }

    return w_new;
}

int main() {


    const std::vector<std::function<double(double, const boost_vector&)>> functions = {f_dx, f_du, f_dy, f_dv};
    boost_vector initial_condition(4);
    initial_condition <<= 25, 0., 25., 0;

    double t_begin = 0.0;
    double t_end = 135.00;
    double th = 0.001;
    std::ofstream file;

    RKM rkm;
    SorSolver sor_solver;

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

    boost_matrix    u = zero_matrix(n, n),
                    v = zero_matrix(n, n),
                    psi = zero_matrix(n, n),
                    w = zero_matrix(n, n),
                    w_new = zero_matrix(n, n);

    //init_w(w);

    upper(u, [](double x, double y) {return std::pow(std::sin(PI * x), 2.0);});

    for (size_t i = 0; i < n; ++i)
        diff(u, i, n - 1, Diff::dx);

    double t_end_vorticity = 3.0;
    double t = 0;
    get_vorticity(u, v, w);
    int point = (n - 1) / 2;
    file.open("Psi.txt");

    while(t < t_end_vorticity) {
        std::cout << "t: " << t << std::endl;
        file << t << " " << w(point, point) << std::endl;
        sor_solver.init(psi, w, dx);
        fill_boundary_vorticity(w, psi, u, v);
        w = vorticity_evolution(w, psi);
        t += dt;
    }
    fill_boundary_vorticity(w, psi, u, v);
    get_v_from_psi(psi, u, v);

    std::cout << w((n - 1) / 2, (n - 1) / 2) << std::endl;

    //for (size_t i = 0; i < n; ++i) {
        size_t i = (n - 1) / 2;
        for (size_t j = 0; j < n; ++j) {
            double x = i * dx;
            double y = j * dx;
            //file << y << " " << u(i, j) << std::endl;
            //file << x << " " << y << " " << psi(i, j) << " " << w(i, j) << " " << u(i, j) << " " << v(i, j) << std::endl;
        }
    //}
    file.close();

    /*file.open("Psi.txt");

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double x = i * dx;
            double y = j * dx;

            file << x << " " << y << " " << psi(i, j) << std::endl;
        }
    }*/


    return 0;
}