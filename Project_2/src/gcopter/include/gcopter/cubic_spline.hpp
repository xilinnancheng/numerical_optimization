#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include "cubic_curve.hpp"

#include <Eigen/Eigen>

#include <cmath>
#include <vector>

namespace cubic_spline
{

    // The banded system class is used for solving
    // banded linear system Ax=b efficiently.
    // A is an N*N band matrix with lower band width lowerBw
    // and upper band width upperBw.
    // Banded LU factorization has O(N) time complexity.
    class BandedSystem
    {
    public:
        // The size of A, as well as the lower/upper
        // banded width p/q are needed
        inline void create(const int &n, const int &p, const int &q)
        {
            // In case of re-creating before destroying
            destroy();
            N = n;
            lowerBw = p;
            upperBw = q;
            int actualSize = N * (lowerBw + upperBw + 1);
            ptrData = new double[actualSize];
            std::fill_n(ptrData, actualSize, 0.0);
            return;
        }

        inline void destroy()
        {
            if (ptrData != nullptr)
            {
                delete[] ptrData;
                ptrData = nullptr;
            }
            return;
        }

    private:
        int N;
        int lowerBw;
        int upperBw;
        // Compulsory nullptr initialization here
        double *ptrData = nullptr;

    public:
        // Reset the matrix to zero
        inline void reset(void)
        {
            std::fill_n(ptrData, N * (lowerBw + upperBw + 1), 0.0);
            return;
        }

        // The band matrix is stored as suggested in "Matrix Computation"
        inline const double &operator()(const int &i, const int &j) const
        {
            return ptrData[(i - j + upperBw) * N + j];
        }

        inline double &operator()(const int &i, const int &j)
        {
            return ptrData[(i - j + upperBw) * N + j];
        }

        // This function conducts banded LU factorization in place
        // Note that NO PIVOT is applied on the matrix "A" for efficiency!!!
        inline void factorizeLU()
        {
            int iM, jM;
            double cVl;
            for (int k = 0; k <= N - 2; ++k)
            {
                iM = std::min(k + lowerBw, N - 1);
                cVl = operator()(k, k);
                for (int i = k + 1; i <= iM; ++i)
                {
                    if (operator()(i, k) != 0.0)
                    {
                        operator()(i, k) /= cVl;
                    }
                }
                jM = std::min(k + upperBw, N - 1);
                for (int j = k + 1; j <= jM; ++j)
                {
                    cVl = operator()(k, j);
                    if (cVl != 0.0)
                    {
                        for (int i = k + 1; i <= iM; ++i)
                        {
                            if (operator()(i, k) != 0.0)
                            {
                                operator()(i, j) -= operator()(i, k) * cVl;
                            }
                        }
                    }
                }
            }
            return;
        }

        // This function solves Ax=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template <typename EIGENMAT>
        inline void solve(EIGENMAT &b) const
        {
            int iM;
            for (int j = 0; j <= N - 1; ++j)
            {
                iM = std::min(j + lowerBw, N - 1);
                for (int i = j + 1; i <= iM; ++i)
                {
                    if (operator()(i, j) != 0.0)
                    {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; --j)
            {
                b.row(j) /= operator()(j, j);
                iM = std::max(0, j - upperBw);
                for (int i = iM; i <= j - 1; ++i)
                {
                    if (operator()(i, j) != 0.0)
                    {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            return;
        }

        // This function solves ATx=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template <typename EIGENMAT>
        inline void solveAdj(EIGENMAT &b) const
        {
            int iM;
            for (int j = 0; j <= N - 1; ++j)
            {
                b.row(j) /= operator()(j, j);
                iM = std::min(j + upperBw, N - 1);
                for (int i = j + 1; i <= iM; ++i)
                {
                    if (operator()(j, i) != 0.0)
                    {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; --j)
            {
                iM = std::max(0, j - lowerBw);
                for (int i = iM; i <= j - 1; ++i)
                {
                    if (operator()(j, i) != 0.0)
                    {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
        }
    };

    class CubicSpline
    {
    public:
        CubicSpline() = default;
        ~CubicSpline() { A.destroy(); }

    private:
        int N;
        Eigen::Vector2d headP;
        Eigen::Vector2d tailP;
        BandedSystem A;
        Eigen::MatrixX2d b;

    public:
        inline void setConditions(const Eigen::Vector2d &headPos,
                                  const Eigen::Vector2d &tailPos,
                                  const int &pieceNum)
        {
            headP = headPos;
            tailP = tailPos;
            N = pieceNum;
            A.create(4 * N, 4, 4);
            b.resize(4 * N, 2);

            A.reset();
            A(0, 0) = 1.0;
            A(1, 1) = 1.0;
            for (int i = 0; i < N - 1; ++i)
            {
                A(4 * i + 2, 4 * i + 2) = 2.0;
                A(4 * i + 2, 4 * i + 3) = 6.0;
                A(4 * i + 2, 4 * i + 6) = -2.0;
                A(4 * i + 3, 4 * i) = 1.0;
                A(4 * i + 3, 4 * i + 1) = 1.0;
                A(4 * i + 3, 4 * i + 2) = 1.0;
                A(4 * i + 3, 4 * i + 3) = 1.0;
                A(4 * i + 4, 4 * i) = 1.0;
                A(4 * i + 4, 4 * i + 1) = 1.0;
                A(4 * i + 4, 4 * i + 2) = 1.0;
                A(4 * i + 4, 4 * i + 3) = 1.0;
                A(4 * i + 4, 4 * i + 4) = -1.0;
                A(4 * i + 5, 4 * i + 1) = 1.0;
                A(4 * i + 5, 4 * i + 2) = 2.0;
                A(4 * i + 5, 4 * i + 3) = 3.0;
                A(4 * i + 5, 4 * i + 5) = -1.0;
            }
            A(4 * N - 2, 4 * N - 4) = 1.0;
            A(4 * N - 2, 4 * N - 3) = 1.0;
            A(4 * N - 2, 4 * N - 2) = 1.0;
            A(4 * N - 2, 4 * N - 1) = 1.0;
            A(4 * N - 1, 4 * N - 3) = 1.0;
            A(4 * N - 1, 4 * N - 2) = 2.0;
            A(4 * N - 1, 4 * N - 1) = 3.0;
            A.factorizeLU();
            return;
        }

        inline void setInnerPoints(const Eigen::Ref<const Eigen::Matrix2Xd> &inPs)
        {
            b.setZero();
            b.row(0) = headP.transpose();
            for (int i = 0; i < N - 1; ++i)
            {
                b.row(4 * i + 3) = inPs.col(i).transpose();
            }
            b.row(4 * N - 2) = tailP.transpose();
            A.solve(b);
            return;
        }

        inline void getCurve(CubicCurve &curve) const
        {
            curve.clear();
            curve.reserve(N);
            for (int index = 0; index < N; ++index)
            {
                curve.emplace_back(1, b.block<4, 2>(4*index, 0).transpose().rowwise().reverse());
            }
            return;
        }

        inline void getStretchEnergy(double &energy) const
        {
            energy = 0.0;
            for (int index = 0; index < N; ++index)
            {
                energy += 4 * b.row(4 * index + 2).squaredNorm() + 12.0 * b.row(4 * index + 3).squaredNorm() + 12.0 * b.row(4 * index + 2).dot(b.row(4 * index + 3));
            }
            return;
        }

        inline const Eigen::MatrixX2d &getCoeffs(void) const
        {
            return b;
        }

        inline void getGrad(Eigen::Ref<Eigen::Matrix2Xd> gradByPoints) const
        {
            // TODO
            Eigen::MatrixX2d gdC(4 * N, 2);
            gdC.setZero();
            for (int i = 0; i < N; ++i)
            {
                gdC.row(4 * i + 3) = 12.0 * b.row(4 * i + 2) +
                                     24.0 * b.row(4 * i + 3);
                gdC.row(4 * i + 2) = 8.0 * b.row(4 * i + 2) +
                                     12.0 * b.row(4 * i + 3);
                gdC.block<2, 2>(4 * i, 0).setZero();
            }

            A.solveAdj(gdC);

            gradByPoints.resize(2, N - 1);
            gradByPoints.setZero();
            for (int i = 0; i < N - 1; ++i)
            {
                gradByPoints.col(i) = gdC.row(4 * i + 3).transpose();
            }
        }
    };
}

#endif
