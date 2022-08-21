#include <iostream>

#include "sdqp/sdqp.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    int m = 5;
    Eigen::Matrix<double, 3, 3> Q;
    Eigen::Matrix<double, 3, 1> c;
    Eigen::Matrix<double, 3, 1> x; // decision variables
    Eigen::Matrix<double, 3, 1> x_est;
    Eigen::Matrix<double, 3, 1> dx;
    Eigen::Matrix<double, -1, 3> A(m, 3); // constraint matrix
    Eigen::VectorXd b(m);                 // constraint bound

    Q << 8.0, -6.0, 2.0, -6.0, 6.0, -3.0, 2.0, -3.0, 2.0;
    c << 1.0, 3.0, -2.0;

    A << 0.0, -1.0, -2.0,
        -1.0, 1.0, -3.0,
        1.0, -2.0, 0.0,
        -1.0, -2.0, -1.0,
        3.0, 5.0, 1.0;
    b << -1.0, 2.0, 7.0, 2.0, -1.0;

    size_t iter_num = 0;
    double minobj = 0.0;
    double pho = 100000;

    Eigen::Matrix<double, 3, 3> Q_reg = Q + 1.0 / pho * Eigen::MatrixXd::Identity(3, 3);
    x.setZero();
    x_est.setZero();
    dx.setZero();
    while (iter_num == 0 || dx.norm() > 1e-10)
    {
        x_est = x;
        Eigen::Matrix<double, 3, 1> c_reg = c + 1.0 / pho * x_est;
        x.setZero();
        minobj = sdqp::sdqp<3>(Q_reg, c_reg, A, b, x);
        dx = x - x_est;
        ++iter_num;
    }

    std::cout << "iteration number: " << iter_num << std::endl;
    std::cout << "optimal sol: " << x.transpose() << std::endl;
    std::cout << "optimal obj: " << 0.5 * x.transpose() * Q_reg * x + c.transpose() * x << std::endl;
    std::cout << "reg optimal obj: " << minobj << std::endl;
    std::cout << "cons precision: " << (A * x - b).maxCoeff() << std::endl;

    return 0;
}
