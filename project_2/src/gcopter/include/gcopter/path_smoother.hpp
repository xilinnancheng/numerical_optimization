#ifndef PATH_SMOOTHER_HPP
#define PATH_SMOOTHER_HPP

#include "cubic_spline.hpp"
#include "lbfgs.hpp"

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
            PathSmoother* path_smoother = (PathSmoother*)ptr;
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
            path_smoother->GradOfPotentialFunc(x, grad_potential);

            // Total gradient.
            g.setZero();
            g = grad_energy + grad_potential;

            // Cost value.
            double energy_value = 0.0;
            path_smoother->cubSpline.getStretchEnergy(energy_value);
            double potential_value = 0.0;
            potential_value = path_smoother->ValueOfPotentialFunc(x);
            cost = energy_value + potential_value;
            return cost;
        }

        inline double ValueOfPotentialFunc(const Eigen::VectorXd &x)
        {
            double potential_value = 0.0;

            for (int point_index = 0; point_index < pieceN - 1; ++point_index)
            {
                for (int obs_index = 0; obs_index < diskObstacles.cols(); ++obs_index)
                {
                    Eigen::Vector3d obs_vec = diskObstacles.col(obs_index);
                    potential_value += penaltyWeight * std::max<double>(0.0, obs_vec[2] - std::hypot(obs_vec[0] - x[point_index], obs_vec[1] - x[point_index + pieceN - 1]));
                }
            }
            return potential_value;
        }

        inline void GradOfPotentialFunc(const Eigen::VectorXd &x, Eigen::VectorXd &g)
        {
            for (int point_index = 0; point_index < pieceN - 1; ++point_index)
            {
                for (int obs_index = 0; obs_index < diskObstacles.cols(); ++obs_index)
                {
                    Eigen::Vector3d obs_vec = diskObstacles.col(obs_index);
                    double denomitor_value = std::hypot(x[point_index] - obs_vec[0], x[point_index + pieceN - 1] - obs_vec[1]);
                    if(denomitor_value < obs_vec[2]){
                        g[point_index] += -penaltyWeight * (x[point_index] - obs_vec[0]) / denomitor_value;
                        g[point_index + pieceN - 1] += -penaltyWeight * (x[point_index + pieceN - 1] - obs_vec[1]) / denomitor_value;
                    }else{
                        g[point_index] += 0.0;
                        g[point_index + pieceN - 1] += 0.0;
                    }
                }
            }
        }

    public:
        inline bool setup(const Eigen::Vector2d &initialP,
                          const Eigen::Vector2d &terminalP,
                          const int &pieceNum,
                          const Eigen::Matrix3Xd &diskObs,
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

            int ret = lbfgs::lbfgs_optimize(x,
                                            minCost,
                                            &PathSmoother::costFunction,
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
