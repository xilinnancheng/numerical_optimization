# Cubic Spline Optimization Using LBFGS
## Problem Description
<img src=https://github.com/xilinnancheng/numerical_optimization/blob/main/project_2/problem_description.png width = "600" height="400"/><br/>

## Code Review
```
// To calculate cost and gradient based on cubic spline and obstacle info.
inline static double costFunction(void *ptr, const Eigen::VectorXd &x, Eigen::VectorXd &g)
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
```

```
    inline int line_search_lewisoverton(Eigen::VectorXd &x,
                                        double &f,
                                        Eigen::VectorXd &g,
                                        double &stp,
                                        const Eigen::VectorXd &s,
                                        const Eigen::VectorXd &xp,
                                        const Eigen::VectorXd &gp,
                                        const double stpmin,
                                        const double stpmax,
                                        const callback_data_t &cd,
                                        const lbfgs_parameter_t &param)
    {
        double upper_bound = std::numeric_limits<double>::max();
        double lower_bound = 0.0;

        double f_next = 0.0;
        Eigen::VectorXd g_next = g;
        int ls = 0;
        f = cd.proc_evaluate(cd.instance, xp, g);
        double slope = s.dot(gp);

        while (1)
        {
            f_next = cd.proc_evaluate(cd.instance, xp + stp * s, g_next);
            // armijo condition
            // f(x_k) - f(x_k + alpha * d) >= -c1 * alpha * df(x_k)
            if (f - f_next < -param.f_dec_coeff * stp * slope)
            {
                upper_bound = stp;
            }
            else if (s.dot(g_next) < param.s_curv_coeff * slope)
            {
                // weak wolfe condition
                // d_t * d(x_k + alpha * d) >= c2 * d_t * df(x_k)
                lower_bound = stp;
            }
            else
            {
                break;
            }

            if (upper_bound < std::numeric_limits<double>::max())
            {
                stp = 0.5 * (upper_bound + lower_bound);
            }
            else
            {
                stp = 2 * lower_bound;
            }
            ls++;
        }

        x = xp + stp * s;
        g = g_next;
        f = f_next;
        return ls;
    }
```
## Experimental Result
<img src=https://github.com/xilinnancheng/numerical_optimization/blob/main/project_2/lbfgs_1.png width = "600" height="400"/><br/>

<img src=https://github.com/xilinnancheng/numerical_optimization/blob/main/project_2/lbfgs_2.png width = "600" height="400"/><br/>
