//
// Created by anton on 11/16/22.
//

#include "SorSolver.h"

void SorSolver::init(boost_matrix& A, const boost_matrix& W, double dx) {
    N = A.size1();
    r_optimal = 2.0 / (1.0 + PI * N);
    _dx = dx;
    //boost_matrix sol = zero_matrix(N, N);

    size_t iter = 0;
    while(true) {
        double max_error = 0;
        for (size_t i = 1; i < N - 1; ++i) {
            for (size_t j = 1; j < N - 1; ++j) {
                double A_ij_old = A(i, j);
                double A_ij_new =
                        1.0 / 4.0 * (A(i + 1,j) + A(i - 1, j) + A(i, j + 1) + A(i, j - 1) + dx * dx * W(i, j));

                A(i, j) = (1.0 - r_optimal) * A(i, j) + r_optimal * A_ij_new;
                //std::cout << A[i][j] << " " << A_ij_old << " " << std::abs(A[i][j] - A_ij_old) << std::endl;
                if (std::abs(A(i, j) - A_ij_old) > max_error)
                    max_error = std::abs(A(i, j) - A_ij_old);
            }
        }

        ++iter;
        //std::cout << max_error << std::endl;
        if (max_error < e || iter > 10000)
            break;
    }
}