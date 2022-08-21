#include <iostream>

#include "sdqp/sdqp.hpp"

using namespace std;
using namespace Eigen;

Eigen::VectorXd VecProject(const Eigen::VectorXd &vec)
{
    int size = vec.size();
    Eigen::VectorXd proj_vec = Eigen::VectorXd::Zero(size);
    Eigen::VectorXd v1 = vec.segment(1, size - 1);
    double v0 = vec[0];
    double v1_norm = v1.norm();

    if (v0 <= -v1_norm)
    {
        return proj_vec;
    }
    else if (v0 >= v1_norm)
    {
        return vec;
    }
    else
    {
        proj_vec = vec;
        proj_vec[0] = v1_norm;
        proj_vec = (v0 + v1_norm) / (2 * v1_norm) * proj_vec;
    }
    return proj_vec;
}

Eigen::MatrixXd BSubDiff(const Eigen::VectorXd &vec)
{
    int size = vec.size();
    Eigen::MatrixXd bsub_diff = Eigen::MatrixXd::Zero(size, size);
    Eigen::VectorXd v1 = vec.segment(1, size - 1);
    double v0 = vec[0];
    double v1_norm = v1.norm();

    if (v0 <= -v1_norm)
    {
        return bsub_diff;
    }
    else if (v0 >= v1_norm)
    {
        return Eigen::MatrixXd::Identity(size, size);
    }
    else
    {
        bsub_diff(0, 0) = 0.5;
        bsub_diff.block(0, 1, 1, size - 1) = 0.5 / v1_norm * v1.transpose();
        bsub_diff.block(1, 0, size - 1, 1) = 0.5 / v1_norm * v1;
        bsub_diff.block(1, 1, size - 1, size - 1) = (v0 + v1_norm) / (2 * v1_norm) * Eigen::MatrixXd::Identity(size - 1, size - 1) - v0 / (2 * std::pow(v1_norm, 3)) * v1 * v1.transpose();
    }
    return bsub_diff;
}

int main(int argc, char **argv)
{
    constexpr int m = 7;
    Eigen::Matrix<double, m + 1, m> A;
    Eigen::Matrix<double, m + 1, 1> b;
    Eigen::Matrix<double, m, 1> x; // decision variables
    Eigen::Matrix<double, m, 1> c;
    Eigen::Matrix<double, m + 1, 1> mu;
    double pho = 1;
    double belta = 1e3;
    double gamma = 1;

    A.setZero();
    b.setZero();
    x.setZero();
    c.setZero();
    mu.setZero();

    A(0, 0) = 1.0;
    b(0) = 1.0;
    for (int index = 0; index < m; ++index)
    {
        A(index + 1, index) = m - index;
        b(index + 1) = 2 * index + 1;
        c(index) = index + 1;
    }

    int outer_loop_iter_num = 0;
    double grad_alm_infinity_norm = 0.0;
    double conic_constraint_infinity_norm = 0.0;

    while (outer_loop_iter_num == 0 || !(conic_constraint_infinity_norm < 1e-4 && grad_alm_infinity_norm < 1e-4))
    {
        int in_loop_iter_num = 0;
        Eigen::Matrix<double, m, 1> grad_alm;
        grad_alm.setZero();
        while (in_loop_iter_num == 0 || grad_alm_infinity_norm > 1e-4)
        {
            // gradient of alm conic problem
            Eigen::Matrix<double, m + 1, 1> conic_vec = mu - pho * (A * x + b);
            grad_alm = c - A.transpose() * VecProject(conic_vec);
            grad_alm_infinity_norm = std::numeric_limits<double>::lowest();
            for (int i = 0; i < grad_alm.size(); ++i)
            {
                grad_alm_infinity_norm = std::max<double>(grad_alm_infinity_norm, std::fabs(grad_alm[i]));
            }

            // semi-smooth hessian of alm conic problem
            double epsilon = std::min<double>(1, grad_alm_infinity_norm) / 10.0;
            Eigen::Matrix<double, m, m> ss_hession_alm = pho * A.transpose() * BSubDiff(conic_vec) * A + epsilon * Eigen::MatrixXd::Identity(m, m);
            x = x - ss_hession_alm.inverse() * grad_alm;
            in_loop_iter_num++;
        }

        Eigen::Matrix<double, m + 1, 1> conic_constraint_vec = mu / pho - VecProject(mu / pho - A * x - b);
        conic_constraint_infinity_norm = std::numeric_limits<double>::lowest();
        for (int i = 0; i < conic_constraint_vec.size(); ++i)
        {
            conic_constraint_infinity_norm = std::max<double>(conic_constraint_infinity_norm, std::fabs(conic_constraint_vec[i]));
        }

        mu = VecProject(mu - pho * (A * x + b));
        pho = std::min<double>((1 + gamma) * pho, belta);
        ++outer_loop_iter_num;
    }

    std::cout << "outer loop iteration number: " << outer_loop_iter_num << std::endl;
    std::cout << "optimal sol: " << x.transpose() << std::endl;
    std::cout << "optimal obj: " << c.transpose() * x << std::endl;
    return 0;
}
