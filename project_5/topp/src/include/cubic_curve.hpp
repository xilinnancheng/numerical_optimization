#ifndef CUBIC_CURVE_HPP
#define CUBIC_CURVE_HPP

#include <Eigen/Eigen>

#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>

class CubicPolynomial
{
private:
    double duration;
    Eigen::Matrix<double, 1, 4> coeffMat;

public:
    CubicPolynomial() = default;

    CubicPolynomial(double dur, const Eigen::Matrix<double, 1, 4> &cMat)
        : duration(dur), coeffMat(cMat) {}

    inline int getDim() const
    {
        return 2;
    }

    inline int getDegree() const
    {
        return 3;
    }

    inline double getDuration() const
    {
        return duration;
    }

    inline const Eigen::Matrix<double, 1, 4> &getCoeffMat() const
    {
        return coeffMat;
    }

    inline double getPos(const double &t) const
    {
        return coeffMat(3) + t * (coeffMat(2) + t * (coeffMat(1) + t * coeffMat(0)));
    }

    inline double getVel(const double &t) const
    {
        return coeffMat(2) + t * (2.0 * coeffMat(1) + 3.0 * t * coeffMat(0));
    }

    inline double getAcc(const double &t) const
    {
        return 2.0 * coeffMat(1) + 6.0 * t * coeffMat(0);
    }
};
#endif
