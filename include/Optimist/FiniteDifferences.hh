/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_FINITE_DIFFERENCES_HH
#define OPTIMIST_FINITE_DIFFERENCES_HH

#include "Optimist.hh"

namespace Optimist {

  namespace FiniteDifferences {

    /**
    * Helper class for computing finite differences epsilons.
    * \tparam Real Scalar number type
    */
    template<typename Real>
    class Epsilon {
      OPTIMIST_BASIC_CONSTANTS(Real)
      Real m_epsilon_1{std::sqrt(EPSILON)}; /**< Square root of machine epsilon \f$ \sqrt{\epsilon} \f$. */
      Real m_epsilon_2{std::pow(EPSILON, 0.75)}; /**< Epsilon to the power of 3/4. */
      Real m_epsilon_3{std::pow(EPSILON, 0.25)}; /**< Epsilon to the power of 1/4. */
    public:

      /** Returns \f$ \epsilon_1 \f$ scaled by |v|+1. */
      Real epsilon_1(Real const v) const {return (std::abs(v) + 1.0)*this->m_epsilon_1;}

      /** Returns \f$ \epsilon_2 \f$ scaled by |v|+1. */
      Real epsilon_2(Real const v) const {return (std::abs(v) + 1.0)*this->m_epsilon_2;}

      /** Returns \f$ \epsilon_3 \f$ scaled by |v|+1. */
      Real epsilon_3(Real const v) const {return (std::abs(v) + 1.0)*this->m_epsilon_3;}
    };


    /**
    * Compute one-sided finite differences (for scalars).
    * \param[in] f_0 Function value at \f$ x \f$.
    * \param[in] f_1 Function value at \f$ x+h_1 \f$.
    * \param[in] f_2 Function value at \f$ x+h_2 \f$.
    * \param[in] h_1 Step size 1.
    * \param[in] h_2 Step size 2.
    * \tparam Real Scalar number type.
    * \return Approximated derivative.
    */
    template<typename Real>
    inline Real SideFiniteDifferences(Real const f_0, Real const f_1, Real const f_2, Real h_1,
      Real h_2)
    {
      Real d_f_1{f_1 - f_0};
      Real d_f_2{f_2 - f_0};
      Real d_h{h_1 - h_2};
      return ((h_1/h_2)*d_f_2 - (h_2/h_1)*d_f_1) / d_h;
    }

    /**
    * Compute centered finite differences (for scalars).
    * \param[in] f_l Function value at \f$ x-h \f$.
    * \param[in] f_c Function value at \f$ x \f$.
    * \param[in] f_r Function value at \f$ x+h \f$.
    * \param[in] h Step size.
    * \tparam Real Scalar number type.
    * \return Approximated derivative.
    */
    template<typename Real>
    inline Real CenteredFiniteDifferences(Real const f_l, Real const f_c, Real const f_r, Real h)
    {
      constexpr Real EPSILON{std::numeric_limits<Real>::epsilon()};
      Real diff_r{(f_r-f_c)/h};
      Real diff_l{(f_c-f_l)/h};
      Real weight_r{std::abs(diff_r) + EPSILON};
      Real weight_l{std::abs(diff_l) + EPSILON};
      Real weight_max{std::max(weight_r, weight_l)};
      weight_r = std::sqrt(weight_r/weight_max);
      weight_l = std::sqrt(weight_l/weight_max);
      return (diff_r*weight_l + diff_l*weight_r) / (weight_r+weight_l);
    }

    /**
    * Compute one-sided finite differences (for fixed-size Eigen vectors).
    * \param[in] f_0 Vector at \f$ x \f$.
    * \param[in] f_1 Vector at \f$ x+h_1 \f$.
    * \param[in] f_2 Vector at \f$ x+h_2 \f$.
    * \param[in] h_1 Step size 1.
    * \param[in] h_2 Step size 2.
    * \tparam Vector Fixed-size Eigen vector type.
    * \tparam Real Scalar number type.
    * \return Approximated derivative vector.
    */
    template<typename Vector, typename Real = typename Vector::Scalar>
    inline Vector SideFiniteDifferences(Vector const & f_0, Vector const & f_1,
      Vector const & f_2, Real h_1, Real h_2)
    {
      Real d_h{h_1-h_2}, t_1{-(h_2/h_1)/d_h}, t_2{(h_1/h_2)/d_h};
      return t_1*(f_1-f_0) + t_2*(f_2-f_0);
    }

    /**
    * Compute centered finite differences (for fixed-size Eigen vectors).
    * \param[in] f_l Vector at \f$ x-h \f$.
    * \param[in] f_c Vector at \f$ x \f$.
    * \param[in] f_r Vector at \f$ x+h \f$.
    * \param[in] h Step size.
    * \tparam Vector Fixed-size Eigen vector type.
    * \tparam Real Scalar number type.
    * \return Approximated derivative vector.
    */
    template<typename Vector, typename Real = typename Vector::Scalar>
    inline Vector CenteredFiniteDifferences(Vector const & f_l, Vector const & f_c, Vector const & f_r,
      Real h)
    {
      constexpr Real EPSILON{std::numeric_limits<Real>::epsilon()};
      Vector diff_r((f_r - f_c) / h);
      Vector diff_l((f_c - f_l) / h);
      Vector weight_r(diff_r.array().abs() + EPSILON);
      Vector weight_l(diff_l.array().abs() + EPSILON);
      Vector weight_max(weight_r.array().max(weight_l.array()));
      weight_r = (weight_r.array()/weight_max.array()).sqrt();
      weight_l = (weight_l.array()/weight_max.array()).sqrt();
      return (diff_r.array()*weight_l.array() + diff_l.array()*weight_r.array()) / (weight_r.array()+weight_l.array());
    }

    /**
    * Compute finite differences gradient for a scalar function (for fixed-size Eigen vectors).
    * \tparam Function Callable type <tt>bool(Vector const &, Real &)</tt>.
    * \tparam Vector Fixed-size Eigen vector type.
    * \tparam Real Scalar number type.
    * \param[in] x Point at which to compute gradient.
    * \param[in] function Function to differentiate.
    * \param[out] out Gradient vector.
    * \return True if successful, false otherwise.
    */
    template <typename Function, typename Vector, typename Real = typename Vector::Scalar>
    bool Gradient(Vector const & x, Function && function, Vector & out)
    {
      Epsilon<Real> eps;
      Real v_c{0.0};
      if (!function(x, v_c) || !std::isfinite(v_c)) {return false;}
      Eigen::Index dim_x{x.size()};
      out.resize(dim_x);
      out.setZero();
      for (Eigen::Index i{0}; i < dim_x; ++i) {
        Vector v_x(x);
        Real tmp{x[i]}, h_1{eps.epsilon_1(tmp)}, h_2{eps.epsilon_2(tmp)};
        v_x[i] = tmp + h_1; Real v_r{0.0}; bool is_finite_r{function(v_x, v_r) && std::isfinite(v_r)};
        v_x[i] = tmp - h_1; Real v_l{0.0}; bool is_finite_l{function(v_x, v_l) && std::isfinite(v_l)};
        Integer ic{(is_finite_r && is_finite_l) ? 0 : (is_finite_r ? 1 : (is_finite_l ? -1 : -2))};
        switch (ic) {
          case 0:
            out[i] = CenteredFiniteDifferences(v_l, v_c, v_r, h_1);
            break;
          case 1: {
            v_x[i] = tmp + h_2; Real v_rr{0.0}; bool is_finite_rr{function(v_x, v_rr) && std::isfinite(v_rr)};
            out[i] = is_finite_rr ? SideFiniteDifferences(v_c, v_r, v_rr, h_1, h_2) : (v_r-v_c)/h_1;
            break;
          }
          case -1: {
            v_x[i] = tmp - h_2; Real v_ll{0.0}; bool is_finite_ll{function(v_x, v_ll) && std::isfinite(v_ll)};
            out[i] = is_finite_ll ? SideFiniteDifferences(v_c, v_l, v_ll, -h_1, -h_2) : (v_c-v_l)/h_1;
            break;
          }
          case -2: {
            out[i] = 0.0;
            return false;
          }
        }
      }
      return out.allFinite();
    }

    /**
    * Compute finite differences Jacobian for a vector function (for fixed-size Eigen vectors).
    * \tparam Function Callable type <tt>bool(Vector const &, Vector &)</tt>.
    * \tparam Vector Fixed-size Eigen vector type.
    * \tparam Matrix Fixed-size Eigen matrix type.
    * \tparam Real Scalar number type.
    * \param[in] x Point at which to compute Jacobian.
    * \param[in] fun Function to differentiate.
    * \param[out] out Jacobian matrix.
    * \return True if successful, false otherwise.
    */
    template <typename Function, typename Vector, typename Matrix, typename Real = typename Vector::Scalar>
    bool Jacobian(Vector const & x, Function && function, Matrix & out)
    {
      Epsilon<Real> eps;
      Vector v_c;
      if (!function(x, v_c) || !v_c.allFinite()) {return false;}
      Eigen::Index dim_x{x.size()};
      out.resize(dim_x, dim_x);
      out.setZero();
      for (Eigen::Index j{0}; j < dim_x; ++j) {
        Vector v_x(x);
        Real tmp{x[j]}, h_1{eps.epsilon_1(tmp)}, h_2{eps.epsilon_2(tmp)};
        v_x[j] = tmp + h_1; Vector v_r; bool is_finite_r{function(v_x, v_r) && v_r.allFinite()};
        v_x[j] = tmp - h_1; Vector v_l; bool is_finite_l{function(v_x, v_l) && v_l.allFinite()};
        Integer ic{(is_finite_r && is_finite_l) ? 0 : (is_finite_r ? 1 : (is_finite_l ? -1 : -2))};
        switch (ic) {
          case 0:
            out.col(j) = CenteredFiniteDifferences(v_l, v_c, v_r, h_1);
            break;
          case 1: {
            v_x[j] = tmp + h_2; Vector v_rr; bool is_finite_rr{function(v_x, v_rr) && v_rr.allFinite()};
            out.col(j) = is_finite_rr ? SideFiniteDifferences(v_c, v_r, v_rr, h_1, h_2) : (v_r-v_c)/h_1;
            break;
          }
          case -1: {
            v_x[j] = tmp - h_2; Vector v_ll; bool is_finite_ll{function(v_x, v_ll) && v_ll.allFinite()};
            out.col(j) = is_finite_ll ? SideFiniteDifferences(v_c, v_l, v_ll, -h_1, -h_2) : (v_c-v_l)/h_1;
            break;
          }
          case -2: {
            out.col(j).setZero();
            return false;
          }
        }
      }
      return out.allFinite();
    }

    /**
    * Compute finite differences Hessian for a scalar function (for fixed-size Eigen vectors).
    * \tparam Function Callable type <tt>bool(Vector const &, Real &)</tt>.
    * \tparam Vector Fixed-size Eigen vector type.
    * \tparam Matrix Fixed-size Eigen matrix type.
    * \tparam Real Scalar number type.
    * \param[in] x Point at which to compute Hessian.
    * \param[in] function Function to differentiate.
    * \param[out] out Hessian matrix.
    * \return True if successful, false otherwise.
    */
    template <typename Function, typename Vector, typename Matrix, typename Real = typename Vector::Scalar>
    bool Hessian(Vector const & x, Function && function, Matrix & out)
    {
      Epsilon<typename Vector::Scalar> eps;
      Eigen::Index dim_x{x.size()};
      out.resize(dim_x, dim_x);
      out.setZero();
      Real fc{0.0};
      if (!function(x, fc) || !std::isfinite(fc)) {return false;}
      for (Eigen::Index j{0}; j < dim_x; ++j) {
        Real tmp_j{x[j]}, h_j{eps.epsilon_3(tmp_j)};
        Vector v_x(x);
        v_x[j] = tmp_j + h_j; Real fp; if (!function(v_x, fp) || !std::isfinite(fp)) {return false;}
        v_x[j] = tmp_j - h_j; Real fm; if (!function(v_x, fm) || !std::isfinite(fm)) {return false;}
        out(j, j) = ((fp + fm) - 2.0 * fc) / (h_j * h_j);
        for (Eigen::Index i{j + 1}; i < dim_x; ++i) {
          Real tmp_i{x[i]}, h_i{eps.epsilon_3(tmp_i)};
          v_x[i] = tmp_i + h_i;
          v_x[j] = tmp_j + h_j; Real fpp; if (!function(v_x, fpp) || !std::isfinite(fpp)) {return false;}
          v_x[i] = tmp_i - h_i; Real fmp; if (!function(v_x, fmp) || !std::isfinite(fmp)) {return false;}
          v_x[j] = tmp_j - h_j; Real fmm; if (!function(v_x, fmm) || !std::isfinite(fmm)) {return false;}
          v_x[i] = tmp_i + h_i; Real fpm; if (!function(v_x, fpm) || !std::isfinite(fpm)) {return false;}
          Real h_ij{4.0 * h_i * h_j}, value{((fpp + fmm) - (fpm + fmp)) / h_ij};
          out(j, i) = out(i, j) = value;
          v_x[i] = tmp_i;
        }
        v_x[j] = tmp_j;
      }
      return out.allFinite();
    }

  } // namespace FiniteDifferences

} // namespace Optimist

#endif // OPTIMIST_FINITE_DIFFERENCES_HH
