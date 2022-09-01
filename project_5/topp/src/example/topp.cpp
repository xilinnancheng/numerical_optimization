#include <iostream>

#include "cubic_curve.hpp"

using namespace std;
using namespace Eigen;

class TOPP {
 public:
  TOPP(const CubicPolynomial& cubic, double delta_s, double max_vel_,
       double max_acc_, int step, int total_variable_num)
      : cubic_(cubic),
        delta_s_(delta_s),
        max_vel_(max_vel_),
        max_acc_(max_acc_),
        total_variable_num_(total_variable_num),
        step_(step) {
    lambda_ = Eigen::VectorXd(step_ + 1);
    lambda_.setZero();

    mu_conic_1_.resize(step_ - 1, Eigen::Vector3d{0.0, 0.0, 0.0});
    mu_conic_2_.resize(step_, Eigen::Vector3d{0.0, 0.0, 0.0});
    mu_conic_3_.resize(step_, Eigen::Vector2d{0.0, 0.0});
    mu_conic_4_.resize(step_, Eigen::Vector2d{0.0, 0.0});
    mu_conic_5_.resize(step_, Eigen::Vector2d{0.0, 0.0});
    mu_conic_6_.resize(step_, Eigen::Vector2d{0.0, 0.0});
    mu_conic_7_.resize(step_, Eigen::Vector2d{0.0, 0.0});

    SetUp();
  }

  Eigen::VectorXd VecProject(const Eigen::VectorXd& vec) {
    int size = vec.size();
    Eigen::VectorXd proj_vec = Eigen::VectorXd::Zero(size);
    Eigen::VectorXd v1 = vec.segment(1, size - 1);
    double v0 = vec[0];
    double v1_norm = v1.norm();

    if (v0 <= -v1_norm) {
      return proj_vec;
    } else if (v0 >= v1_norm) {
      return vec;
    } else {
      proj_vec = vec;
      proj_vec[0] = v1_norm;
      proj_vec = (v0 + v1_norm) / (2 * v1_norm) * proj_vec;
    }
    return proj_vec;
  }

  Eigen::MatrixXd BSubDiff(const Eigen::VectorXd& vec) {
    int size = vec.size();
    Eigen::MatrixXd bsub_diff = Eigen::MatrixXd::Zero(size, size);
    Eigen::VectorXd v1 = vec.segment(1, size - 1);
    double v0 = vec[0];
    double v1_norm = v1.norm();

    if (v0 <= -v1_norm) {
      return bsub_diff;
    } else if (v0 >= v1_norm) {
      return Eigen::MatrixXd::Identity(size, size);
    } else {
      bsub_diff(0, 0) = 0.5;
      bsub_diff.block(0, 1, 1, size - 1) = 0.5 / v1_norm * v1.transpose();
      bsub_diff.block(1, 0, size - 1, 1) = 0.5 / v1_norm * v1;
      bsub_diff.block(1, 1, size - 1, size - 1) =
          (v0 + v1_norm) / (2 * v1_norm) *
              Eigen::MatrixXd::Identity(size - 1, size - 1) -
          v0 / (2 * std::pow(v1_norm, 3)) * v1 * v1.transpose();
    }
    return bsub_diff;
  }

  double MaxConicConstraitInfinityNorm(const Eigen::VectorXd& x, double pho) {
    double conic_constraint_infinity_norm =
        std::numeric_limits<double>::lowest();

    auto MaxConicConstraitInfinityNormImpl =
        [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
            const Eigen::VectorXd& mu_conic) {
          Eigen::VectorXd vec =
              mu_conic / pho - VecProject(mu_conic / pho - A * x - b);
          for (int index = 0; index < vec.size(); ++index) {
            conic_constraint_infinity_norm =
                std::max(conic_constraint_infinity_norm, std::fabs(vec[index]));
          }
        };

    for (int index = 0; index < (step_ - 1); ++index) {
      const Eigen::MatrixXd& A_conic_1 = A_conic_1_[index];
      const Eigen::VectorXd& b_conic_1 = b_conic_1_[index];
      MaxConicConstraitInfinityNormImpl(A_conic_1, b_conic_1,
                                        mu_conic_1_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_2 = A_conic_2_[index];
      const Eigen::VectorXd& b_conic_2 = b_conic_2_[index];
      MaxConicConstraitInfinityNormImpl(A_conic_2, b_conic_2,
                                        mu_conic_2_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_3 = A_conic_3_[index];
      const Eigen::VectorXd& b_conic_3 = b_conic_3_[index];
      MaxConicConstraitInfinityNormImpl(A_conic_3, b_conic_3,
                                        mu_conic_3_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_4 = A_conic_4_[index];
      const Eigen::VectorXd& b_conic_4 = b_conic_4_[index];
      MaxConicConstraitInfinityNormImpl(A_conic_4, b_conic_4,
                                        mu_conic_4_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_5 = A_conic_5_[index];
      const Eigen::VectorXd& b_conic_5 = b_conic_5_[index];
      MaxConicConstraitInfinityNormImpl(A_conic_5, b_conic_5,
                                        mu_conic_5_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_6 = A_conic_6_[index];
      const Eigen::VectorXd& b_conic_6 = b_conic_6_[index];
      MaxConicConstraitInfinityNormImpl(A_conic_6, b_conic_6,
                                        mu_conic_6_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_7 = A_conic_7_[index];
      const Eigen::VectorXd& b_conic_7 = b_conic_7_[index];
      MaxConicConstraitInfinityNormImpl(A_conic_7, b_conic_7,
                                        mu_conic_7_[index]);
    }
    return conic_constraint_infinity_norm;
  }

  void UpdateMuConic(const Eigen::VectorXd& x, double pho) {
    auto UpdateMuConicImpl = [&](const Eigen::MatrixXd& A,
                                 const Eigen::VectorXd& b,
                                 Eigen::VectorXd& mu_conic) {
      mu_conic = VecProject(mu_conic - pho * (A * x + b));
    };

    for (int index = 0; index < (step_ - 1); ++index) {
      const Eigen::MatrixXd& A_conic_1 = A_conic_1_[index];
      const Eigen::VectorXd& b_conic_1 = b_conic_1_[index];
      UpdateMuConicImpl(A_conic_1, b_conic_1, mu_conic_1_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_2 = A_conic_2_[index];
      const Eigen::VectorXd& b_conic_2 = b_conic_2_[index];
      UpdateMuConicImpl(A_conic_2, b_conic_2, mu_conic_2_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_3 = A_conic_3_[index];
      const Eigen::VectorXd& b_conic_3 = b_conic_3_[index];
      UpdateMuConicImpl(A_conic_3, b_conic_3, mu_conic_3_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_4 = A_conic_4_[index];
      const Eigen::VectorXd& b_conic_4 = b_conic_4_[index];
      UpdateMuConicImpl(A_conic_4, b_conic_4, mu_conic_4_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_5 = A_conic_5_[index];
      const Eigen::VectorXd& b_conic_5 = b_conic_5_[index];
      UpdateMuConicImpl(A_conic_5, b_conic_5, mu_conic_5_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_6 = A_conic_6_[index];
      const Eigen::VectorXd& b_conic_6 = b_conic_6_[index];
      UpdateMuConicImpl(A_conic_6, b_conic_6, mu_conic_6_[index]);
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_7 = A_conic_7_[index];
      const Eigen::VectorXd& b_conic_7 = b_conic_7_[index];
      UpdateMuConicImpl(A_conic_7, b_conic_7, mu_conic_7_[index]);
    }
  }

  void Solve() {
    Eigen::VectorXd x(total_variable_num_);
    x.setZero();
    int outer_loop_iter_num = 0;
    double grad_alm_infinity_norm = 0.0;
    double conic_constraint_infinity_norm = 0.0;
    double eq_constraint_infinity_norm = 0.0;

    double pho = 1.0;
    double belta = 1000;
    double gamma = 1;

    while (outer_loop_iter_num == 0 ||
           !(std::max(conic_constraint_infinity_norm,
                      eq_constraint_infinity_norm) < 1e-4 &&
             grad_alm_infinity_norm < 1e-4)) {
      int in_loop_iter_num = 0;
      Eigen::VectorXd grad_alm(total_variable_num_);
      grad_alm.setZero();
      while (in_loop_iter_num == 0 || grad_alm_infinity_norm > 1e-4) {
        // gradient of alm conic problem
        Eigen::VectorXd eq_vec = lambda_ + pho * (G_ * x - h_);
        Eigen::VectorXd grad_soc(total_variable_num_);
        Eigen::MatrixXd hessian_soc(total_variable_num_, total_variable_num_);
        GradientAndHessianOfSecondOrderCone(x, pho, grad_soc, hessian_soc);
        grad_alm = c_ + G_.transpose() * eq_vec - grad_soc;
        grad_alm_infinity_norm = std::numeric_limits<double>::lowest();
        for (int i = 0; i < grad_alm.size(); ++i) {
          grad_alm_infinity_norm =
              std::max<double>(grad_alm_infinity_norm, std::fabs(grad_alm[i]));
        }

        // semi-smooth hessian of alm conic problem
        double epsilon = std::min<double>(1, grad_alm_infinity_norm) / 10.0;
        Eigen::MatrixXd ss_hession_alm =
            pho * (G_.transpose() * G_ + hessian_soc) +
            epsilon * Eigen::MatrixXd::Identity(total_variable_num_,
                                                total_variable_num_);
        x = x - ss_hession_alm.inverse() * grad_alm;
        in_loop_iter_num++;
        // std::cout << x << std::endl;
        // std::cout << 1111111111111111 << std::endl;
        // std::cout << G_.transpose() * eq_vec << std::endl;
        // std::cout << 2222222222222222 << std::endl;
        // std::cout << grad_soc << std::endl;
        // std::cout << ss_hession_alm.inverse() * grad_alm << std::endl;
        // if (in_loop_iter_num == 5) break;
      }

      conic_constraint_infinity_norm = MaxConicConstraitInfinityNorm(x, pho);

      Eigen::VectorXd eq_constraint_vec = G_ * x - h_;
      eq_constraint_infinity_norm = std::numeric_limits<double>::lowest();
      for (int i = 0; i < eq_constraint_vec.size(); ++i) {
        eq_constraint_infinity_norm = std::max<double>(
            eq_constraint_infinity_norm, std::fabs(eq_constraint_vec[i]));
      }

      UpdateMuConic(x, pho);
      lambda_ = lambda_ + pho * eq_constraint_vec;
      pho = std::min<double>((1 + gamma) * pho, belta);
      ++outer_loop_iter_num;
      if(outer_loop_iter_num == 3)
          break;
    }

    std::cout << "outer loop iteration number: " << outer_loop_iter_num
              << std::endl;
    std::cout << "optimal sol: " << x.transpose() << std::endl;
  }

 private:
  void GradientAndHessianOfSecondOrderCone(const Eigen::VectorXd& x, double pho,
                                           Eigen::VectorXd& grad,
                                           Eigen::MatrixXd& hessian) {
    grad.setZero();
    hessian.setZero();

    auto GradImpl = [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& conic_vec){
        return A.transpose() * VecProject(conic_vec);
    };

    for (int index = 0; index < (step_ - 1); ++index) {
      const Eigen::MatrixXd& A_conic_1 = A_conic_1_[index];
      const Eigen::VectorXd& b_conic_1 = b_conic_1_[index];

      Eigen::VectorXd conic_vec =
          mu_conic_1_[index] - pho * (A_conic_1 * x + b_conic_1);
      grad += A_conic_1.transpose() * VecProject(conic_vec);
      hessian += A_conic_1.transpose() * BSubDiff(conic_vec) * A_conic_1;
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_2 = A_conic_2_[index];
      const Eigen::VectorXd& b_conic_2 = b_conic_2_[index];
      Eigen::VectorXd conic_vec =
          mu_conic_2_[index] - pho * (A_conic_2 * x + b_conic_2);
      grad += A_conic_2.transpose() * VecProject(conic_vec);
      hessian += A_conic_2.transpose() * BSubDiff(conic_vec) * A_conic_2;
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_3 = A_conic_3_[index];
      const Eigen::VectorXd& b_conic_3 = b_conic_3_[index];
      Eigen::VectorXd conic_vec =
          mu_conic_3_[index] - pho * (A_conic_3 * x + b_conic_3);
      grad += A_conic_3.transpose() * VecProject(conic_vec);
      hessian += A_conic_3.transpose() * BSubDiff(conic_vec) * A_conic_3;
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_4 = A_conic_4_[index];
      const Eigen::VectorXd& b_conic_4 = b_conic_4_[index];
      Eigen::VectorXd conic_vec =
          mu_conic_4_[index] - pho * (A_conic_4 * x + b_conic_4);
      grad += A_conic_4.transpose() * VecProject(conic_vec);
      hessian += A_conic_4.transpose() * BSubDiff(conic_vec) * A_conic_4;
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_5 = A_conic_5_[index];
      const Eigen::VectorXd& b_conic_5 = b_conic_5_[index];
      Eigen::VectorXd conic_vec =
          mu_conic_5_[index] - pho * (A_conic_5 * x + b_conic_5);
      grad += A_conic_5.transpose() * VecProject(conic_vec);
      hessian += A_conic_5.transpose() * BSubDiff(conic_vec) * A_conic_5;
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_6 = A_conic_6_[index];
      const Eigen::VectorXd& b_conic_6 = b_conic_6_[index];
      Eigen::VectorXd conic_vec =
          mu_conic_6_[index] - pho * (A_conic_6 * x + b_conic_6);
      grad += A_conic_6.transpose() * VecProject(conic_vec);
      hessian += A_conic_6.transpose() * BSubDiff(conic_vec) * A_conic_6;
    }

    for (int index = 0; index < step_; ++index) {
      const Eigen::MatrixXd& A_conic_7 = A_conic_7_[index];
      const Eigen::VectorXd& b_conic_7 = b_conic_7_[index];
      Eigen::VectorXd conic_vec =
          mu_conic_7_[index] - pho * (A_conic_7 * x + b_conic_7);
      grad += A_conic_7.transpose() * VecProject(conic_vec);
      hessian += A_conic_7.transpose() * BSubDiff(conic_vec) * A_conic_7;
    }
  }

  void SetUp() {
    c_ = Eigen::VectorXd(total_variable_num_);
    c_.setZero();

    for (int index = 0; index < step_ - 1; ++index) {
      c_(3 * step_ + index) = 2.0 * delta_s_;
    }

    Eigen::MatrixXd G_1(step_ - 1, total_variable_num_);
    Eigen::VectorXd h_1(step_ - 1);
    G_1.setZero();
    h_1.setZero();
    for (int index = 0; index < step_ - 1; ++index) {
      G_1(index, step_ + index + 1) = 1.0;
      G_1(index, step_ + index) = -1.0;
      G_1(index, index) = -2.0 * delta_s_;
    }

    Eigen::MatrixXd G_2(2, total_variable_num_);
    Eigen::VectorXd h_2(2);
    G_2.setZero();
    h_2.setZero();
    G_2(0, step_) = 1;
    G_2(1, 2 * step_ - 1) = 1;
    h_2(0) = 0.5;
    h_2(1) = 0.5;

    G_ = Eigen::MatrixXd(step_ + 1, total_variable_num_);
    h_ = Eigen::VectorXd(step_ + 1);
    G_.block(0, 0, step_ - 1, total_variable_num_) = G_1;
    G_.block(step_ - 1, 0, 2, total_variable_num_) = G_2;
    h_.segment(0, step_ - 1) = h_1;
    h_.segment(step_ - 1, 2) = h_2;

    for (int index = 0; index < (step_ - 1); ++index) {
      Eigen::MatrixXd A_conic_1(3, total_variable_num_);
      Eigen::VectorXd b_conic_1(3);
      A_conic_1.setZero();
      b_conic_1.setZero();
      A_conic_1(0, 2 * step_ + index) = 1.0;
      A_conic_1(0, 2 * step_ + index + 1) = 1.0;
      A_conic_1(0, 3 * step_ + index) = 1.0;

      A_conic_1(2, 2 * step_ + index) = 1.0;
      A_conic_1(2, 2 * step_ + index + 1) = 1.0;
      A_conic_1(2, 3 * step_ + index) = -1.0;
      b_conic_1(1) = 2.0;
      A_conic_1_.emplace_back(std::move(A_conic_1));
      b_conic_1_.emplace_back(std::move(b_conic_1));
    }

    for (int index = 0; index < step_; ++index) {
      Eigen::MatrixXd A_conic_2(3, total_variable_num_);
      Eigen::VectorXd b_conic_2(3);
      A_conic_2.setZero();
      b_conic_2.setZero();
      A_conic_2(0, step_ + index) = 1.0;
      A_conic_2(1, 2 * step_ + index) = 2.0;
      A_conic_2(2, step_ + index) = 1.0;
      b_conic_2(0) = 1.0;
      b_conic_2(2) = -1.0;
      A_conic_2_.emplace_back(std::move(A_conic_2));
      b_conic_2_.emplace_back(std::move(b_conic_2));
    }

    for (int index = 0; index < step_; ++index) {
      Eigen::MatrixXd A_conic_3(2, total_variable_num_);
      Eigen::VectorXd b_conic_3(2);
      A_conic_3.setZero();
      b_conic_3.setZero();
      A_conic_3(0, step_ + index) = 1.0;
      b_conic_3(0) = 0.0;
      b_conic_3(1) = 0.0;
      A_conic_3_.emplace_back(std::move(A_conic_3));
      b_conic_3_.emplace_back(std::move(b_conic_3));
    }

    for (int index = 0; index < step_; ++index) {
      Eigen::MatrixXd A_conic_4(2, total_variable_num_);
      Eigen::VectorXd b_conic_4(2);
      A_conic_4.setZero();
      b_conic_4.setZero();
      A_conic_4(1, 2 * step_ + index) = cubic_.getVel(index * delta_s_);
      b_conic_4(0) = max_vel_;
      b_conic_4(1) = 0.0;
      A_conic_4_.emplace_back(std::move(A_conic_4));
      b_conic_4_.emplace_back(std::move(b_conic_4));
    }

    for (int index = 0; index < step_; ++index) {
      Eigen::MatrixXd A_conic_5(2, total_variable_num_);
      Eigen::VectorXd b_conic_5(2);
      A_conic_5.setZero();
      b_conic_5.setZero();
      A_conic_5(0, 2 * step_ + index) = cubic_.getVel(index * delta_s_);
      b_conic_5(0) = 0.0;
      b_conic_5(1) = -max_vel_;
      A_conic_5_.emplace_back(std::move(A_conic_5));
      b_conic_5_.emplace_back(std::move(b_conic_5));
    }

    for (int index = 0; index < step_; ++index) {
      Eigen::MatrixXd A_conic_6(2, total_variable_num_);
      Eigen::VectorXd b_conic_6(2);
      A_conic_6.setZero();
      b_conic_6.setZero();
      A_conic_6(1, index) = cubic_.getVel(index * delta_s_);
      A_conic_6(1, step_ + index) = cubic_.getAcc(index * delta_s_);
      b_conic_6(0) = max_acc_;
      b_conic_6(1) = 0.0;
      A_conic_6_.emplace_back(std::move(A_conic_6));
      b_conic_6_.emplace_back(std::move(b_conic_6));
    }

    for (int index = 0; index < step_; ++index) {
      Eigen::MatrixXd A_conic_7(2, total_variable_num_);
      Eigen::VectorXd b_conic_7(2);
      A_conic_7.setZero();
      b_conic_7.setZero();
      A_conic_7(0, index) = cubic_.getVel(index * delta_s_);
      A_conic_7(0, step_ + index) = cubic_.getAcc(index * delta_s_);
      b_conic_7(0) = 0.0;
      b_conic_7(1) = -max_acc_;
      A_conic_7_.emplace_back(std::move(A_conic_7));
      b_conic_7_.emplace_back(std::move(b_conic_7));
    }
  }

  CubicPolynomial cubic_;
  double delta_s_;
  double max_vel_;
  double max_acc_;
  int total_variable_num_;
  int step_;
  Eigen::VectorXd lambda_;
  Eigen::VectorXd c_;
  Eigen::MatrixXd G_;
  Eigen::VectorXd h_;
  std::vector<Eigen::VectorXd> mu_conic_1_;
  std::vector<Eigen::VectorXd> mu_conic_2_;
  std::vector<Eigen::VectorXd> mu_conic_3_;
  std::vector<Eigen::VectorXd> mu_conic_4_;
  std::vector<Eigen::VectorXd> mu_conic_5_;
  std::vector<Eigen::VectorXd> mu_conic_6_;
  std::vector<Eigen::VectorXd> mu_conic_7_;

  std::vector<Eigen::MatrixXd> A_conic_1_;
  std::vector<Eigen::MatrixXd> A_conic_2_;
  std::vector<Eigen::MatrixXd> A_conic_3_;
  std::vector<Eigen::MatrixXd> A_conic_4_;
  std::vector<Eigen::MatrixXd> A_conic_5_;
  std::vector<Eigen::MatrixXd> A_conic_6_;
  std::vector<Eigen::MatrixXd> A_conic_7_;

  std::vector<Eigen::VectorXd> b_conic_1_;
  std::vector<Eigen::VectorXd> b_conic_2_;
  std::vector<Eigen::VectorXd> b_conic_3_;
  std::vector<Eigen::VectorXd> b_conic_4_;
  std::vector<Eigen::VectorXd> b_conic_5_;
  std::vector<Eigen::VectorXd> b_conic_6_;
  std::vector<Eigen::VectorXd> b_conic_7_;
};

int main(int argc, char** argv) {
  double max_vel = 5.0;
  double max_acc = 1.0;
  double delta_s = 1.0;
  double cubic_s = 5;
  int step = cubic_s / delta_s + 1;
  int total_variable_num = step * 4;

  Eigen::Matrix<double, 1, 4> x_coeffMat;
  Eigen::Matrix<double, 1, 4> y_coeffMat;
  x_coeffMat << 1.0, 2.0, 3.0, 4.0;
  y_coeffMat << 1.0, 2.0, 3.0, 4.0;
  CubicPolynomial x_cubic(cubic_s, x_coeffMat);
//   CubicPolynomial y_cubic(cubic_s, y_coeffMat);

  TOPP x_topp(x_cubic, delta_s, max_vel, max_acc, step, total_variable_num);
  x_topp.Solve();
  return 0;
}
