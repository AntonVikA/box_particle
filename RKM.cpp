//
// Created by anton on 10/15/22.
//

#include "RKM.h"

boost_vector RKM::operator() (double t, const boost_vector& y) {
    boost_vector res(number_of_equations);
    for (size_t i = 0; i < number_of_equations; ++i) {
        res[i] = functions_[i](t, y);
    }

    return res;
}

RKM::RKM(const std::vector<std::function<function_type>>& functions) : functions_{functions}, number_of_equations{functions_.size()} {

}

void RKM::init(const std::vector<std::function<function_type>>& functions) {
    functions_ = functions;
    number_of_equations = functions_.size();
    Y = boost_vector (number_of_equations);
    dY = boost_vector (number_of_equations);
}

void RKM::check_boundary_collision() {
    if (std::abs(Y[2] - 0.) < 0.2 && Y[3] < 0) {
        Y[3] *= -1.0;
    }

    if (std::abs(Y[2] - L) < 0.2 && Y[3] > 0) {
        Y[3] *= -1.0;
    }

    if (std::abs(Y[0] - .0) < 0.2 && Y[1] < 0) {
        Y[1] *= -1.0;
    }
    if (std::abs(Y[0] - L) < 0.2 && Y[1] > 0) {
        Y[1] *= -1.0;
    }

}


void RKM::solve(double t_begin, double t_end, const boost_vector& init_conditions, double h) {
    assert(init_conditions.size() == functions_.size());

    t_ = t_begin;

    Y = init_conditions;

    file.open("box_move.txt");
    size_t iter = 0;
    size_t n = t_end / h;
    for (size_t i = 0; i < n; ++i) {
        K1 = h * (*this)(t_, Y);
        K2 = h * (*this)(t_ + h / 2.0, Y + K1 / 2.0);
        K3 = h * (*this)(t_ + h / 2.0, Y + K2 / 2.0);
        K4 = h * (*this)(t_ + h, Y + K3);
        t_ += h;
        dY = 1.0 / 6.0 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
        Y = Y + dY;
        file << t_ << " ";
        for (size_t i = 0; i < Y.size(); ++i)
            file << Y[i] << " ";
        file << std::endl;

        //check_boundary_collision();
    }

    file.close();
}