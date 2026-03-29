/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_FINITE_DIFFERENCES_HH
#define OPTIMIST_FINITE_DIFFERENCES_HH

#include "Optimist.hh"

namespace Optimist {

  namespace FiniteDifferences {

    /**
     * Helper class for computing finite differences epsilons.
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
      requires TypeTrait<Scalar>::IsScalar
    class Epsilon {
      OPTIMIST_BASIC_CONSTANTS(Scalar)

      Scalar m_epsilon_1{SQRT_EPSILON}; /**< Square root of machine epsilon \f$
                                           \sqrt{\epsilon} \f$. */
      Scalar m_epsilon_2{
        std::pow(EPSILON, 0.75)};       /**< Epsilon to the power of 3/4. */
      Scalar m_epsilon_3{
        std::pow(EPSILON, 0.25)};       /**< Epsilon to the power of 1/4. */
     public:
      /**
       * Returns \f$ \epsilon_1 \f$ scaled by |v|+1.
       * \param[in] v Value to scale epsilon.
       */
      Scalar epsilon_1(const Scalar v) const {
        return (std::abs(v) + 1.0) * this->m_epsilon_1;
      }

      /**
       * Returns \f$ \epsilon_2 \f$ scaled by |v|+1.
       * \param[in] v Value to scale epsilon.
       */
      Scalar epsilon_2(const Scalar v) const {
        return (std::abs(v) + 1.0) * this->m_epsilon_2;
      }

      /**
       * Returns \f$ \epsilon_3 \f$ scaled by |v|+1.
       * \param[in] v Value to scale epsilon.
       */
      Scalar epsilon_3(const Scalar v) const {
        return (std::abs(v) + 1.0) * this->m_epsilon_3;
      }
    };

    /**
     * Compute one-sided finite differences.
     * \tparam Scalar Floating-point number type.
     * \param[in] f_0 Function value at \f$ x \f$.
     * \param[in] f_1 Function value at \f$ x+h_1 \f$.
     * \param[in] f_2 Function value at \f$ x+h_2 \f$.
     * \param[in] h_1 Step size 1.
     * \param[in] h_2 Step size 2.
     * \param[in] i Index of the variable.
     * \param[out] out The approximated derivative.
     */
    template <typename Vector, typename Scalar>
      requires TypeTrait<Vector>::IsEigen && TypeTrait<Scalar>::IsScalar
    inline void SideFiniteDifferences(const Scalar f_0,
                                      const Scalar f_1,
                                      const Scalar f_2,
                                      const Scalar h_1,
                                      const Scalar h_2,
                                      const Integer i,
                                      Vector &out) {
      const Scalar d_f_1{f_1 - f_0}, d_f_2{f_2 - f_0}, d_h{h_1 - h_2};
      if constexpr (TypeTrait<Vector>::IsFixed ||
                    TypeTrait<Vector>::IsDynamic) {
        out(i) = ((h_1 / h_2) * d_f_2 - (h_2 / h_1) * d_f_1) / d_h;
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        out.coeffRef(i) = ((h_1 / h_2) * d_f_2 - (h_2 / h_1) * d_f_1) / d_h;
      }
    }

    /**
     * Compute centered finite differences.
     * \tparam Scalar Floating-point number type.
     * \param[in] f_l Function value at \f$ x-h \f$.
     * \param[in] f_c Function value at \f$ x \f$.
     * \param[in] f_r Function value at \f$ x+h \f$.
     * \param[in] h Step size.
     * \param[in] i Index of the variable.
     * \param[out] out The approximated derivative.
     */
    template <typename Vector, typename Scalar>
      requires TypeTrait<Vector>::IsEigen && TypeTrait<Scalar>::IsScalar
    inline void CenteredFiniteDifferences(const Scalar f_l,
                                          const Scalar f_c,
                                          const Scalar f_r,
                                          const Scalar h,
                                          const Integer i,
                                          Vector &out) {
      constexpr Scalar EPSILON{std::numeric_limits<Scalar>::epsilon()};
      const Scalar diff_r{(f_r - f_c) / h}, diff_l{(f_c - f_l) / h};
      Scalar weight_r{std::abs(diff_r) + EPSILON},
          weight_l{std::abs(diff_l) + EPSILON};
      const Scalar weight_max{std::max(weight_r, weight_l)};
      weight_r = std::sqrt(weight_r / weight_max);
      weight_l = std::sqrt(weight_l / weight_max);
      if constexpr (TypeTrait<Vector>::IsFixed ||
                    TypeTrait<Vector>::IsDynamic) {
        out(i) =
            (diff_r * weight_l + diff_l * weight_r) / (weight_r + weight_l);
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        out.coeffRef(i) =
            (diff_r * weight_l + diff_l * weight_r) / (weight_r + weight_l);
      }
    }

    /**
     * Compute one-sided finite differences (for dense Eigen vectors).
     * \tparam Vector Dense or sparse Eigen vector type.
     * \tparam Matrix Dense or sparse Eigen matrix type.
     * \tparam Scalar Floating-point number type.
     * \param[in] f_0 Vector at \f$ x \f$.
     * \param[in] f_1 Vector at \f$ x+h_1 \f$.
     * \param[in] f_2 Vector at \f$ x+h_2 \f$.
     * \param[in] h_1 Step size 1.
     * \param[in] h_2 Step size 2.
     * \param[in] i Index of the variable.
     * \param[out] out The approximated derivative.
     */
    template <typename Vector, typename Matrix, typename Scalar>
      requires TypeTrait<Vector>::IsEigen && TypeTrait<Scalar>::IsScalar
    inline void SideFiniteDifferences(const Vector &f_0,
                                      const Vector &f_1,
                                      const Vector &f_2,
                                      const Scalar h_1,
                                      const Scalar h_2,
                                      const Integer i,
                                      Matrix &out) {
      const Scalar d_h{h_1 - h_2};
      const Scalar t_1{-(h_2 / h_1) / d_h}, t_2{(h_1 / h_2) / d_h};
      if constexpr (TypeTrait<Vector>::IsFixed ||
                    TypeTrait<Vector>::IsDynamic) {
        out.col(i) = t_1 * (f_1 - f_0) + t_2 * (f_2 - f_0);
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        const Vector tmp(t_1 * (f_1 - f_0) + t_2 * (f_2 - f_0));
        for (Eigen::Index j{0}; j < tmp.size(); ++j) {
          out.coeffRef(j, i) = tmp.coeff(j);
        }
      }
    }

    /**
     * Compute centered finite differences (for dense and Eigen vectors).
     * \tparam Vector Dense or sparse Eigen vector type.
     * \tparam Matrix Dense or sparse Eigen matrix type.
     * \tparam Scalar Floating-point number type.
     * \param[in] f_l Vector at \f$ x-h \f$.
     * \param[in] f_c Vector at \f$ x \f$.
     * \param[in] f_r Vector at \f$ x+h \f$.
     * \param[in] h Step size.
     * \param[in] i Index of the variable.
     * \param[out] out The approximated derivative vector.
     */
    template <typename Vector, typename Matrix, typename Scalar>
      requires TypeTrait<Vector>::IsEigen && TypeTrait<Scalar>::IsScalar
    inline void CenteredFiniteDifferences(const Vector &f_l,
                                          const Vector &f_c,
                                          const Vector &f_r,
                                          const Scalar h,
                                          const Integer i,
                                          Matrix &&out) {
      constexpr Scalar EPSILON{std::numeric_limits<Scalar>::epsilon()};
      const Vector diff_r((f_r - f_c) / h), diff_l((f_c - f_l) / h);
      if constexpr (TypeTrait<Vector>::IsFixed ||
                    TypeTrait<Vector>::IsDynamic) {
        Vector weight_r(diff_r.array().abs() + EPSILON),
            weight_l(diff_l.array().abs() + EPSILON);
        const Vector weight_max(weight_r.array().max(weight_l.array()));
        weight_r   = (weight_r.array() / weight_max.array()).sqrt();
        weight_l   = (weight_l.array() / weight_max.array()).sqrt();
        out.col(i) = (diff_r.array() * weight_l.array() +
                      diff_l.array() * weight_r.array()) /
                     (weight_r.array() + weight_l.array());
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        Scalar weight_r, weight_l, weight_max;
        for (Eigen::Index j{0}; j < diff_r.size(); ++j) {
          weight_max =
              std::max(std::abs(diff_r.coeff(j)), std::abs(diff_l.coeff(j)));
          weight_r = std::abs(diff_r.coeff(j)) + EPSILON;
          weight_l = std::abs(diff_l.coeff(j)) + EPSILON;
          weight_r = std::sqrt(weight_r / weight_max);
          weight_l = std::sqrt(weight_l / weight_max);
          out.coeffRef(j, i) =
              (diff_r.coeff(j) * weight_l + diff_l.coeff(j) * weight_r) /
              (weight_r + weight_l);
        }
      }
    }

    /**
     * Compute finite differences gradient for a scalar function (for dense
     * Eigen vectors).
     * \tparam Function Callable type <tt>bool(Vector const &, Scalar &)</tt>.
     * \tparam Vector Dense Eigen vector type.
     * \tparam Scalar Floating-point number type.
     * \param[in] function Function to differentiate.
     * \param[in] x Point at which to compute gradient.
     * \param[out] out Gradient vector.
     * \return True if successful, false otherwise.
     */
    template <typename Function,
              typename Vector,
              typename Scalar = typename Vector::Scalar>
      requires std::
                   is_invocable_r_v<bool, Function, const Vector &, Scalar &> &&
               TypeTrait<Scalar>::IsScalar && TypeTrait<Vector>::IsEigen
    inline bool Gradient(Function &&function, const Vector &x, Vector &out) {
      Epsilon<Scalar> eps;
      Scalar v_c{0.0};
      if (!function(x, v_c) || !std::isfinite(v_c)) {
        return false;
      }
      Eigen::Index dim_x{x.size()};
      if constexpr (TypeTrait<Vector>::IsDynamic ||
                    TypeTrait<Vector>::IsSparse) {
        out.resize(dim_x);
      }
      if constexpr (!TypeTrait<Vector>::IsSparse) {
        out.setZero();
      } else {
        out.reserve(dim_x);
      }
      for (Eigen::Index i{0}; i < dim_x; ++i) {
        Vector v_x(x);
        Scalar tmp;
        if constexpr (TypeTrait<Vector>::IsFixed ||
                      TypeTrait<Vector>::IsDynamic) {
          tmp = x(i);
        } else if constexpr (TypeTrait<Vector>::IsSparse) {
          tmp = x.coeff(i);
        }
        Scalar h_1{eps.epsilon_1(tmp)}, h_2{eps.epsilon_2(tmp)};
        if constexpr (TypeTrait<Vector>::IsFixed ||
                      TypeTrait<Vector>::IsDynamic) {
          v_x(i) = tmp + h_1;
        } else if constexpr (TypeTrait<Vector>::IsSparse) {
          v_x.coeffRef(i) = tmp + h_1;
        }
        Scalar v_r{0.0};
        bool is_finite_r{function(v_x, v_r) && std::isfinite(v_r)};
        if constexpr (TypeTrait<Vector>::IsFixed ||
                      TypeTrait<Vector>::IsDynamic) {
          v_x(i) = tmp - h_1;
        } else if constexpr (TypeTrait<Vector>::IsSparse) {
          v_x.coeffRef(i) = tmp - h_1;
        }
        Scalar v_l{0.0};
        bool is_finite_l{function(v_x, v_l) && std::isfinite(v_l)};
        Eigen::Index ic{(is_finite_r && is_finite_l)
                            ? 0
                            : (is_finite_r ? 1 : (is_finite_l ? -1 : -2))};
        switch (ic) {
          case 0: {
            CenteredFiniteDifferences(v_l, v_c, v_r, h_1, i, out);
            break;
          }
          case 1: {
            if constexpr (TypeTrait<Vector>::IsFixed ||
                          TypeTrait<Vector>::IsDynamic) {
              v_x(i) = tmp + h_2;
            } else if constexpr (TypeTrait<Vector>::IsSparse) {
              v_x.coeffRef(i) = tmp + h_2;
            }
            Scalar v_rr{0.0};
            bool is_finite_rr{function(v_x, v_rr) && std::isfinite(v_rr)};
            if (is_finite_rr) {
              SideFiniteDifferences(v_c, v_r, v_rr, h_1, h_2, i, out);
            } else {
              if constexpr (TypeTrait<Vector>::IsFixed ||
                            TypeTrait<Vector>::IsDynamic) {
                out(i) = (v_r - v_c) / h_1;
              } else if constexpr (TypeTrait<Vector>::IsSparse) {
                out.coeffRef(i) = (v_r - v_c) / h_1;
              }
            }
            break;
          }
          case -1: {
            if constexpr (TypeTrait<Vector>::IsFixed ||
                          TypeTrait<Vector>::IsDynamic) {
              v_x(i) = tmp - h_2;
            } else if constexpr (TypeTrait<Vector>::IsSparse) {
              v_x.coeffRef(i) = tmp - h_2;
            }
            Scalar v_ll{0.0};
            bool is_finite_ll{function(v_x, v_ll) && std::isfinite(v_ll)};
            if (is_finite_ll) {
              SideFiniteDifferences(v_c, v_l, v_ll, -h_1, -h_2, i, out);
            } else {
              if constexpr (TypeTrait<Vector>::IsFixed ||
                            TypeTrait<Vector>::IsDynamic) {
                out(i) = (v_c - v_l) / h_1;
              } else if constexpr (TypeTrait<Vector>::IsSparse) {
                out.coeffRef(i) = (v_c - v_l) / h_1;
              }
            }
            break;
          }
          case -2: {
            if constexpr (TypeTrait<Vector>::IsFixed ||
                          TypeTrait<Vector>::IsDynamic) {
              out(i) = 0.0;
            }
            return false;
          }
        }
      }
      if constexpr (TypeTrait<Vector>::IsFixed ||
                    TypeTrait<Vector>::IsDynamic) {
        return out.allFinite();
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        for (Eigen::Index j{0}; j < out.size(); ++j) {
          if (!std::isfinite(out.coeff(j))) {
            return false;
          }
        }
        return true;
      }
    }

    /**
     * Compute finite differences Jacobian for a vector function (for dense
     * Eigen vectors).
     * \tparam Function Callable type <tt>bool(Vector const &, Vector &)</tt>.
     * \tparam Vector Dense Eigen vector type.
     * \tparam Matrix Dense Eigen matrix type.
     * \tparam Scalar Floating-point number type.
     * \param[in] fun Function to differentiate.
     * \param[in] x Point at which to compute Jacobian.
     * \param[out] out Jacobian matrix.
     * \return True if successful, false otherwise.
     */
    template <typename Function,
              typename Vector,
              typename Matrix,
              typename Scalar = typename Vector::Scalar>
      requires std::
                   is_invocable_r_v<bool, Function, const Vector &, Vector &> &&
               TypeTrait<Scalar>::IsScalar && TypeTrait<Vector>::IsEigen &&
               TypeTrait<Matrix>::IsEigen
    inline bool Jacobian(Function &&function, const Vector &x, Matrix &out) {
      Epsilon<Scalar> eps;
      Vector v_c;
      if (!function(x, v_c) || !v_c.allFinite()) {
        return false;
      }
      Eigen::Index dim_x{x.size()};
      if constexpr (TypeTrait<Matrix>::IsDynamic ||
                    TypeTrait<Matrix>::IsSparse) {
        out.resize(dim_x, dim_x);
      }
      if constexpr (!TypeTrait<Matrix>::IsSparse) {
        out.setZero();
      } else {
        out.reserve(dim_x * dim_x);
      }
      for (Eigen::Index j{0}; j < dim_x; ++j) {
        Vector v_x(x);
        Scalar tmp{x(j)}, h_1{eps.epsilon_1(tmp)}, h_2{eps.epsilon_2(tmp)};
        v_x(j) = tmp + h_1;
        Vector v_r;
        bool is_finite_r{function(v_x, v_r) && v_r.allFinite()};
        v_x(j) = tmp - h_1;
        Vector v_l;
        bool is_finite_l{function(v_x, v_l) && v_l.allFinite()};
        Eigen::Index ic{(is_finite_r && is_finite_l)
                            ? 0
                            : (is_finite_r ? 1 : (is_finite_l ? -1 : -2))};
        switch (ic) {
          case 0:
            CenteredFiniteDifferences(v_l, v_c, v_r, h_1, j, out);
            break;
          case 1: {
            v_x(j) = tmp + h_2;
            Vector v_rr;
            bool is_finite_rr{function(v_x, v_rr) && v_rr.allFinite()};
            if (is_finite_rr) {
              SideFiniteDifferences(v_c, v_r, v_rr, h_1, h_2, j, out);
            } else {
              out.col(j) = (v_r - v_c) / h_1;
            }
            break;
          }
          case -1: {
            v_x(j) = tmp - h_2;
            Vector v_ll;
            bool is_finite_ll{function(v_x, v_ll) && v_ll.allFinite()};
            if (is_finite_ll) {
              SideFiniteDifferences(v_c, v_l, v_ll, -h_1, -h_2, j, out);
            } else {
              out.col(j) = (v_c - v_l) / h_1;
            }
            break;
          }
          case -2: {
            out.col(j).setZero();
            return false;
          }
        }
      }
      if constexpr (TypeTrait<Vector>::IsFixed ||
                    TypeTrait<Vector>::IsDynamic) {
        return out.allFinite();
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        for (Eigen::Index r{0}; r < out.rows(); ++r) {
          for (Eigen::Index c{0}; c < out.cols(); ++c) {
            if (!std::isfinite(out.coeff(r, c))) {
              return false;
            }
          }
        }
        return true;
      }
    }

    /**
     * Compute finite differences Hessian for a scalar function (for dense
     * Eigen vectors).
     * \tparam Function Callable type <tt>bool(Vector const &, Real &)</tt>.
     * \tparam Vector Dense Eigen vector type.
     * \tparam Matrix Dense Eigen matrix type.
     * \tparam Scalar Floating-point number type.
     * \param[in] function Function to differentiate.
     * \param[in] x Point at which to compute Hessian.
     * \param[out] out Hessian matrix.
     * \return True if successful, false otherwise.
     */
    template <typename Function,
              typename Vector,
              typename Matrix,
              typename Scalar = typename Vector::Scalar>
      requires std::
                   is_invocable_r_v<bool, Function, const Vector &, Scalar &> &&
               TypeTrait<Scalar>::IsScalar && TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsSparse) && TypeTrait<Matrix>::IsEigen &&
               (!TypeTrait<Matrix>::IsSparse)
    inline bool Hessian(Function &&function, const Vector &x, Matrix &out) {
      Epsilon<Scalar> eps;
      Eigen::Index dim_x{x.size()};
      out.resize(dim_x, dim_x);
      out.setZero();
      Scalar fc{0.0};
      if (!function(x, fc) || !std::isfinite(fc)) {
        return false;
      }
      for (Eigen::Index j{0}; j < dim_x; ++j) {
        Scalar tmp_j{x(j)}, h_j{eps.epsilon_3(tmp_j)};
        Vector v_x(x);
        v_x(j) = tmp_j + h_j;
        Scalar fp;
        if (!function(v_x, fp) || !std::isfinite(fp)) {
          return false;
        }
        v_x(j) = tmp_j - h_j;
        Scalar fm;
        if (!function(v_x, fm) || !std::isfinite(fm)) {
          return false;
        }
        out(j, j) = ((fp + fm) - 2.0 * fc) / (h_j * h_j);
        for (Eigen::Index i{j + 1}; i < dim_x; ++i) {
          Scalar tmp_i{x(i)}, h_i{eps.epsilon_3(tmp_i)};
          v_x(i) = tmp_i + h_i;
          v_x(j) = tmp_j + h_j;
          Scalar fpp;
          if (!function(v_x, fpp) || !std::isfinite(fpp)) {
            return false;
          }
          v_x(i) = tmp_i - h_i;
          Scalar fmp;
          if (!function(v_x, fmp) || !std::isfinite(fmp)) {
            return false;
          }
          v_x(j) = tmp_j - h_j;
          Scalar fmm;
          if (!function(v_x, fmm) || !std::isfinite(fmm)) {
            return false;
          }
          v_x(i) = tmp_i + h_i;
          Scalar fpm;
          if (!function(v_x, fpm) || !std::isfinite(fpm)) {
            return false;
          }
          Scalar h_ij{4.0 * h_i * h_j},
              value{((fpp + fmm) - (fpm + fmp)) / h_ij};
          out(j, i) = out(i, j) = value;
          v_x(i)                = tmp_i;
        }
        v_x(j) = tmp_j;
      }
      if constexpr (TypeTrait<Vector>::IsFixed ||
                    TypeTrait<Vector>::IsDynamic) {
        return out.allFinite();
      } else if constexpr (TypeTrait<Vector>::IsSparse) {
        for (Eigen::Index r{0}; r < out.rows(); ++r) {
          for (Eigen::Index c{0}; c < out.cols(); ++c) {
            if (!std::isfinite(out.coeff(r, c))) {
              return false;
            }
          }
        }
        return true;
      }
    }

  }  // namespace FiniteDifferences

}  // namespace Optimist

#endif  // OPTIMIST_FINITE_DIFFERENCES_HH
