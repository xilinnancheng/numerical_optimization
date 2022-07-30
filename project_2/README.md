# Cubic Spline Optimization Using LBFGS
## Problem Description
<img src=https://github.com/xilinnancheng/numerical_optimization/blob/main/project_2/problem_description.png width = "600" height="400"/><br/>

## Code Review
```
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
```
## Experimental Result
<img src=https://github.com/xilinnancheng/numerical_optimization/blob/main/project_2/lbfgs_1.png width = "600" height="400"/><br/>

<img src=https://github.com/xilinnancheng/numerical_optimization/blob/main/project_2/lbfgs_2.png width = "600" height="400"/><br/>
