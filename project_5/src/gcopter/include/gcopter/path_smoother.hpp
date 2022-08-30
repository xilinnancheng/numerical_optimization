#ifndef PATH_SMOOTHER_HPP
#define PATH_SMOOTHER_HPP

#include "cubic_spline.hpp"
#include "lbfgs.hpp"
#include "sdqp.hpp"
#include <Eigen/Eigen>

#include <cmath>
#include <cfloat>
#include <iostream>
#include <vector>

namespace path_smoother
{

    class PathSmoother
    {
    private:
        cubic_spline::CubicSpline cubSpline;

        int pieceN;
        Eigen::Matrix3Xd diskObstacles;
        std::vector<Eigen::MatrixX3d> polygon_obstacles;
        double penaltyWeight;
        Eigen::Vector2d headP;
        Eigen::Vector2d tailP;
        Eigen::Matrix2Xd points;
        Eigen::Matrix2Xd gradByPoints;

        lbfgs::lbfgs_parameter_t lbfgs_params;

    private:
        inline static double costFunction(void *ptr,
                                          const Eigen::VectorXd &x,
                                          Eigen::VectorXd &g)
        {
            PathSmoother *path_smoother = (PathSmoother *)ptr;
            // Cubic spline.
            double cost = 0.0;
            int max_index = path_smoother->pieceN - 1;
            path_smoother->points.row(0) = x.segment(0, max_index);
            path_smoother->points.row(1) = x.segment(max_index, max_index);
            path_smoother->cubSpline.setInnerPoints(path_smoother->points);
            path_smoother->cubSpline.getGrad(path_smoother->gradByPoints);

            // Energy gradient.
            Eigen::VectorXd grad_energy;
            grad_energy.resize(2 * (max_index));
            grad_energy.setZero();
            grad_energy.segment(0, max_index) = path_smoother->gradByPoints.row(0).transpose();
            grad_energy.segment(max_index, max_index) = path_smoother->gradByPoints.row(1).transpose();

            // Potential gradient.
            Eigen::VectorXd grad_potential;
            grad_potential.resize(2 * (max_index));
            grad_potential.setZero();
            double potential_value = 0.0;
            potential_value = path_smoother->ValueAndGradientOfPotentialFunc(x, grad_potential);

            // Total gradient.
            g.setZero();
            g = grad_energy + grad_potential;

            // Cost value.
            double energy_value = 0.0;
            path_smoother->cubSpline.getStretchEnergy(energy_value);
            cost = energy_value + potential_value;
            return cost;
        }

        inline double ValueAndGradientOfPotentialFunc(const Eigen::VectorXd &x, Eigen::VectorXd &g)
        {
            double potential_value = 0.0;

            for (int point_index = 0; point_index < pieceN - 1; ++point_index)
            {
                const double x_0 = x[point_index];
                const double y_0 = x[point_index + pieceN - 1];
                const Eigen::Vector2d pos(x_0, y_0);
                // const Eigen::Vector2d c = -2 * pos;
                for (size_t obs_index = 0; obs_index < polygon_obstacles.size(); ++obs_index)
                {
                    const Eigen::MatrixX3d polygon_obstalce = polygon_obstacles[obs_index];
                    Eigen::Vector2d min_x(0.0, 0.0);
                    Eigen::MatrixXd A = polygon_obstalce.leftCols(2);
                    Eigen::VectorXd b = -polygon_obstalce.rightCols(1);
                    Eigen::VectorXd diff = A * pos - b;
                    bool is_out_of_obstacle = false;
                    for(int index = 0;index < diff.size();++index){
                        if(diff[index] > 0){
                            is_out_of_obstacle = true;
                            break;
                        }
                    }
                    // double min_dist_to_obs = std::max<double>(sdqp::sdqp<2>(2 * Eigen::Matrix2d::Identity(), c, polygon_obstalce.leftCols(2), -polygon_obstalce.rightCols(1), min_x) + std::pow(pos.norm(), 2), 0.0);
                    if (!is_out_of_obstacle)
                    {
                        double min_dist = std::numeric_limits<double>::max();
                        Eigen::Vector2d min_dist_vec;
                        for (int row = 0; row < polygon_obstalce.rows(); ++row)
                        {
                            double A = polygon_obstalce.row(row)(0);
                            double B = polygon_obstalce.row(row)(1);
                            double C = polygon_obstalce.row(row)(2);
                            double A_square = A * A;
                            double B_square = B * B;
                            
                            double pen_x = (B_square * x_0 - A * B * y_0 - A * C) / (A_square + B_square);
                            double pen_y = (A_square * y_0 - A * B * x_0 - B * C) / (A_square + B_square);
                            double pen_dist = std::hypot(pen_x - x_0, pen_y - y_0);
                            if (pen_dist < min_dist)
                            {
                                min_dist = pen_dist;
                                min_dist_vec = Eigen::Vector2d{x_0 - pen_x, y_0 - pen_y} / min_dist;
                            }
                        }
                        g[point_index] += penaltyWeight * min_dist_vec[0];
                        g[point_index + pieceN - 1] += penaltyWeight * min_dist_vec[1];
                        potential_value += penaltyWeight * min_dist;
                    }
                }
                // std::cout << x_0 << " " << y_0 << std::endl;
                // std::cout << 1111111111111111111 << std::endl;
                // std::cout <<  g << std::endl;
            }
            return potential_value;
        }

    public:
        inline bool setup(const Eigen::Vector2d &initialP,
                          const Eigen::Vector2d &terminalP,
                          const int &pieceNum,
                          const Eigen::Matrix3Xd &diskObs,
                          const std::vector<Eigen::MatrixX4d> &polgon_obs,
                          const double penaWeight)
        {
            pieceN = pieceNum;
            diskObstacles = diskObs;
            penaltyWeight = penaWeight;
            headP = initialP;
            tailP = terminalP;

            cubSpline.setConditions(headP, tailP, pieceN);

            points.resize(2, pieceN - 1);
            gradByPoints.resize(2, pieceN - 1);
            for (size_t index = 0; index < polgon_obs.size(); ++index)
            {
                Eigen::MatrixX3d polygon(polgon_obs[index].rows() - 2, 3);
                for (int col = 0; col < polgon_obs[index].cols(); ++col)
                {
                    if(col == 2) continue;
                    int col_tmp = col >= 2? col-1:col;
                    polygon.col(col_tmp) = polgon_obs[index].col(col).head(polgon_obs[index].rows() - 2);
                }
                polygon_obstacles.emplace_back(polygon);
            }
            return true;
        }

        inline double optimize(CubicCurve &curve,
                               const Eigen::Matrix2Xd &iniInPs,
                               const double &relCostTol)
        {
            Eigen::VectorXd x(pieceN * 2 - 2);
            x.segment(0, pieceN - 1) = iniInPs.row(0).transpose();
            x.segment(pieceN - 1, pieceN - 1) = iniInPs.row(1).transpose();

            double minCost;
            lbfgs_params.mem_size = 64;
            lbfgs_params.past = 3;
            lbfgs_params.min_step = 1.0e-32;
            lbfgs_params.g_epsilon = 1e-4;
            lbfgs_params.delta = relCostTol;
            // lbfgs_params.max_iterations = 0;

            int ret = lbfgs::lbfgs_optimize(x,
                                            minCost,
                                            &PathSmoother::costFunction,
                                            nullptr,
                                            nullptr,
                                            this,
                                            lbfgs_params);

            if (ret >= 0)
            {
                std::cout << "Congratulation!! Optimize successfully!!!" << std::endl;
                points.row(0) = x.segment(0, pieceN - 1);
                points.row(1) = x.segment(pieceN - 1, pieceN - 1);
                cubSpline.setInnerPoints(points);
                cubSpline.getCurve(curve);
            }
            else
            {
                curve.clear();
                minCost = INFINITY;
                std::cout << "Oh, sorry, optimization Failed: "
                          << lbfgs::lbfgs_strerror(ret)
                          << std::endl;
            }
            return minCost;
        }
    };

}

#endif
