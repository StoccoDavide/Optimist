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

#ifndef OPTIMIST_FINITE_DIFFERENCE_HH
#define OPTIMIST_FINITE_DIFFERENCE_HH

#include "Optimist.hh"

namespace Optimist {

  namespace FiniteDifference {

    /**
    * Helper class for computing finite difference epsilons.
    * \tparam Real Scalar number type
    */
    template<typename Real>
    class Epsilon {
      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */
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
    * Compute one-sided finite difference (for scalars).
    * \param[in] fun_0 Function value at \f$ x \f$.
    * \param[in] fun_1 Function value at \f$ x+h_1 \f$.
    * \param[in] fun_2 Function value at \f$ x+h_2 \f$.
    * \param[in] h_1 Step size 1.
    * \param[in] h_2 Step size 2.
    * \tparam Real Scalar number type.
    * \return Approximated derivative.
    */
    template<typename Real>
    inline Real SideFiniteDifference(Real const fun_0, Real const fun_1, Real const fun_2, Real h_1,
      Real h_2)
    {
      Real d_fun_1{fun_1 - fun_0};
      Real d_fun_2{fun_2 - fun_0};
      Real d_h{h_1 - h_2};
      return ((h_1/h_2)*d_fun_2 - (h_2/h_1)*d_fun_1) / d_h;
    }

    /**
    * Compute centered finite difference (for scalars).
    * \param[in] fun_l Function value at \f$ x-h \f$.
    * \param[in] fun_c Function value at \f$ x \f$.
    * \param[in] fun_r Function value at \f$ x+h \f$.
    * \param[in] h Step size.
    * \tparam Real Scalar number type.
    * \return Approximated derivative.
    */
    template<typename Real>
    inline Real CenteredFiniteDifference(Real const fun_l, Real const fun_c, Real const fun_r, Real h)
    {
      constexpr Real EPSILON{std::numeric_limits<Real>::epsilon()};
      Real diff_r{(fun_r-fun_c)/h};
      Real diff_l{(fun_c-fun_l)/h};
      Real weight_r{std::abs(diff_r) + EPSILON};
      Real weight_l{std::abs(diff_l) + EPSILON};
      Real weight_max{std::max(weight_r, weight_l)};
      weight_r = std::sqrt(weight_r/weight_max);
      weight_l = std::sqrt(weight_l/weight_max);
      return (diff_r*weight_l + diff_l*weight_r) / (weight_r+weight_l);
    }

    /*\
     |   ____                              _
     |  |  _ \ _   _ _ __   __ _ _ __ ___ (_) ___
     |  | | | | | | | '_ \ / _` | '_ ` _ \| |/ __|
     |  | |_| | |_| | | | | (_| | | | | | | | (__
     |  |____/ \__, |_| |_|\__,_|_| |_| |_|_|\___|
     |         |___/
    \*/

    /**< Template alias for <tt>VectorXd<Real></tt>. */
    template<typename Real>
    using VectorXd = Eigen::Vector<Real, Eigen::Dynamic>;

    /**< Template alias for <tt>MatrixXd<Real></tt>. */
    template<typename Real>
    using MatrixXd = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

    /**
    * Compute one-sided finite difference (for dynamic-size Eigen vectors).
    * \param[in] fun_0 Vector at \f$ x \f$.
    * \param[in] fun_1 Vector at \f$ x+h_1 \f$.
    * \param[in] fun_2 Vector at \f$ x+h_2 \f$.
    * \param[in] h_1 Step size 1.
    * \param[in] h_2 Step size 2.
    * \tparam Real Scalar number type.
    * \return Approximated derivative vector.
    */
    template<typename Real>
    inline VectorXd<Real> SideFiniteDifference(VectorXd<Real> const & fun_0, VectorXd<Real> const & fun_1,
      VectorXd<Real> const & fun_2, Real h_1, Real h_2)
    {
      Real d_h{h_1 - h_2};
      Real t_1{-(h_2/h_1) / d_h};
      Real t_2{(h_1/h_2) / d_h};
      return t_1*(fun_1 - fun_0) + t_2*(fun_2 - fun_0);
    }

    /**
    * Compute centered finite difference (for dynamic-size Eigen vectors).
    * \param[in] fun_l Vector at \f$ x-h \f$.
    * \param[in] fun_c Vector at \f$ x \f$.
    * \param[in] fun_r Vector at \f$ x+h \f$.
    * \param[in] h Step size.
    * \tparam Real Scalar number type.
    * \return Approximated derivative vector.
    */
    template<typename Real>
    inline VectorXd<Real> CenteredFiniteDifference(VectorXd<Real> const & fun_l, VectorXd<Real> const & fun_c,
      VectorXd<Real> const & fun_r, Real h)
    {
      constexpr Real EPSILON{std::numeric_limits<Real>::epsilon()};
      VectorXd<Real> diff_r((fun_r - fun_c) / h);
      VectorXd<Real> diff_l((fun_c - fun_l) / h);
      VectorXd<Real> weight_r(diff_r.array().abs() + EPSILON);
      VectorXd<Real> weight_l(diff_l.array().abs() + EPSILON);
      VectorXd<Real> weight_max(weight_r.array().max(weight_l.array()));
      weight_r = (weight_r.array()/weight_max.array()).sqrt();
      weight_l = (weight_l.array()/weight_max.array()).sqrt();
      return (diff_r.array()*weight_l.array() + diff_l.array()*weight_r.array()) / (weight_r.array()+weight_l.array());
    }

    /**
    * Compute finite difference gradient of a scalar function (for dynamic-size Eigen vectors).
    * \param[in] x Point at which to compute gradient.
    * \param[in] fun Function to differentiate.
    * \param[out] grad Gradient vector.
    * \tparam Function Callable type <tt>bool(VectorXd const &, Real &)</tt>.
    * \tparam Real Scalar number type.
    * \return True if successful, false otherwise.
    */
    template <typename Function, typename Real>
    bool Gradient(VectorXd<Real> const & x, Function fun, VectorXd<Real> & grad)
    {
      Epsilon<Real> eps;
      Real v_c{0.0};
      if (!fun(x, v_c) || !std::isfinite(v_c)) {return false;}
      grad.resize(x.size());
      for (Integer i{0}; i < x.size(); ++i) {
        VectorXd<Real> v_x(x);
        Real tmp{x[i]};
        Real h_1{eps.epsilon_1(tmp)};
        Real h_2{eps.epsilon_2(tmp)};
        v_x[i] = tmp + h_1; Real v_r{0.0}; bool is_finite_r{fun(v_x, v_r) && std::isfinite(v_r)};
        v_x[i] = tmp - h_1; Real v_l{0.0}; bool is_finite_l{fun(v_x, v_l) && std::isfinite(v_l)};
        Integer ic{(is_finite_r && is_finite_l) ? 0 : (is_finite_r ? 1 : (is_finite_l ? -1 : -2))};
        switch (ic) {
          case 0:
            grad[i] = CenteredFiniteDifference<Real>(v_l, v_c, v_r, h_1);
            break;
          case 1: {
            v_x[i] = tmp + h_2; Real v_rr{0.0}; bool is_finite_rr{fun(v_x, v_rr) && std::isfinite(v_rr)};
            grad[i] = is_finite_rr ? SideFiniteDifference<Real>(v_c, v_r, v_rr, h_1, h_2) : (v_r-v_c)/h_1;
            break;
          }
          case -1: {
            v_x[i] = tmp - h_2; Real v_ll{0.0}; bool is_finite_ll{fun(v_x, v_ll) && std::isfinite(v_ll)};
            grad[i] = is_finite_ll ? SideFiniteDifference<Real>(v_c, v_l, v_ll, -h_1, -h_2) : (v_c-v_l)/h_1;
            break;
          }
          case -2: {
            grad[i] = 0.0;
            return false;
          }
        }
      }
      return grad.allFinite();
    }

    /**
    * Compute finite difference Jacobian of a vector function (for dynamic-size Eigen vectors).
    * \param[in] x Point at which to compute Jacobian.
    * \param[in] fun Function to differentiate.
    * \param[out] jac Output Jacobian matrix.
    * \tparam Function Callable type <tt>bool(VectorXd<Real> const &, VectorXd<Real> &)</tt>.
    * \tparam Real Scalar number type.
    * \return True if successful, false otherwise.
    */
    template <typename Function, typename Real>
    bool Jacobian(VectorXd<Real> const & x, Function fun, MatrixXd<Real> & jac)
    {
      Epsilon<Real> eps;
      VectorXd<Real> v_c;
      if (!fun(x, v_c) || !v_c.allFinite()) {return false;}
      Integer dim_x{static_cast<Integer>(x.size())};
      Integer dim_f{static_cast<Integer>(v_c.size())};
      jac.resize(dim_f, dim_x);
      for (Integer j{0}; j < dim_x; ++j) {
        VectorXd<Real> v_x(x);
        Real tmp{x[j]};
        Real h_1{eps.epsilon_1(tmp)};
        Real h_2{eps.epsilon_2(tmp)};
        v_x[j] = tmp + h_1; VectorXd<Real> v_r; bool is_finite_r{fun(v_x, v_r) && v_r.allFinite()};
        v_x[j] = tmp - h_1; VectorXd<Real> v_l; bool is_finite_l{fun(v_x, v_l) && v_l.allFinite()};
        Integer ic{(is_finite_r && is_finite_l) ? 0 : (is_finite_r ? 1 : (is_finite_l ? -1 : -2))};
        switch (ic) {
          case 0:
            jac.col(j) = CenteredFiniteDifference<Real>(v_l, v_c, v_r, h_1);
            break;
          case 1: {
            v_x[j] = tmp + h_2; VectorXd<Real> v_rr; bool is_finite_rr{fun(v_x, v_rr) && v_rr.allFinite()};
            jac.col(j) = is_finite_rr ? SideFiniteDifference<Real>(v_c, v_r, v_rr, h_1, h_2) : (v_r-v_c)/h_1;
            break;
          }
          case -1: {
            v_x[j] = tmp - h_2; VectorXd<Real> v_ll; bool is_finite_ll{fun(v_x, v_ll) && v_ll.allFinite()};
            jac.col(j) = is_finite_ll ? SideFiniteDifference<Real>(v_c, v_l, v_ll, -h_1, -h_2) : (v_c-v_l)/h_1;
            break;
          }
          case -2: {
            jac.col(j).setZero();
            return false;
          }
        }
      }
      return jac.allFinite();
    }

    /**
    * Compute finite difference Hessian of a scalar function (for dynamic-size Eigen vectors).
    * \param[in] x Point at which to compute Hessian.
    * \param[in] fun Function to differentiate.
    * \param[out] hes Hessian matrix.
    * \tparam Function Callable type <tt>bool(VectorXd<Real> const &, Real&)</tt>.
    * \tparam Real Scalar number type.
    * \return True if successful, false otherwise.
    */
    template <typename Function, typename Real>
    bool Hessian(VectorXd<Real> const & x, Function fun, MatrixXd<Real> & hes)
    {
      Epsilon<Real> eps;
      Integer dim_x{x.size()};
      hes.resize(dim_x, dim_x);
      Real fc{0.0};
      if (!fun(x, fc) || !std::isfinite(fc)) {return false;}
      for (Integer j{0}; j < dim_x; ++j) {
        Real tmp_j{x[j]};
        Real h_j{eps.epsilon_3(tmp_j)};
        VectorXd<Real> v_x(x);
        v_x[j] = tmp_j + h_j; Real fp{0.0}; if (!fun(v_x, fp) || !std::isfinite(fp)) {return false;}
        v_x[j] = tmp_j - h_j; Real fm{0.0}; if (!fun(v_x, fm) || !std::isfinite(fm)) {return false;}
        hes(j, j) = ((fp + fm) - 2.0 * fc) / (h_j * h_j);
        for (Integer i{j + 1}; i < dim_x; ++i) {
          Real tmp_i{x[i]};
          Real h_i{eps.epsilon_3(tmp_i)};
          v_x[i] = tmp_i + h_i; v_x[j] = tmp_j + h_j; Real fpp{0.0}; if (!fun(v_x, fpp) || !std::isfinite(fpp)) {return false;}
          v_x[i] = tmp_i - h_i; Real fmp{0.0}; if (!fun(v_x, fmp) || !std::isfinite(fmp)) {return false;}
          v_x[j] = tmp_j - h_j; Real fmm{0.0}; if (!fun(v_x, fmm) || !std::isfinite(fmm)) {return false;}
          v_x[i] = tmp_i + h_i; Real fpm{0.0}; if (!fun(v_x, fpm) || !std::isfinite(fpm)) {return false;}
          Real h_ij{4.0 * h_i * h_j};
          Real value{((fpp + fmm) - (fpm + fmp)) / h_ij};
          hes(j, i) = hes(i, j) = value;
          v_x[i] = tmp_i;
        }
        v_x[j] = tmp_j;
      }
      return hes.allFinite();
    }

    /*\
     |   _____ _              _
     |  |  ___(_)_  _____  __| |
     |  | |_  | \ \/ / _ \/ _` |
     |  |  _| | |>  <  __/ (_| |
     |  |_|   |_/_/\_\___|\__,_|
     |
    \*/

    /**
    * Compute one-sided finite difference (for fixed-size Eigen vectors).
    * \param[in] fun_0 Vector at \f$ x \f$.
    * \param[in] fun_1 Vector at \f$ x+h_1 \f$.
    * \param[in] fun_2 Vector at \f$ x+h_2 \f$.
    * \param[in] h_1 Step size 1.
    * \param[in] h_2 Step size 2.
    * \tparam Vector Fixed-size Eigen vector type.
    * \return Approximated derivative vector.
    */
    template<typename Vector>
    inline Vector SideFiniteDifference(Vector const & fun_0, Vector const & fun_1,
      Vector const & fun_2, typename Vector::Scalar h_1, typename Vector::Scalar h_2)
    {
      using Real = typename Vector::Scalar;
      Real d_h{h_1 - h_2};
      Real t_1{-(h_2/h_1) / d_h};
      Real t_2{(h_1/h_2) / d_h};
      return t_1*(fun_1 - fun_0) + t_2*(fun_2 - fun_0);
    }

    /**
    * Compute centered finite difference (for fixed-size Eigen vectors).
    * \param[in] fun_l Vector at \f$ x-h \f$.
    * \param[in] fun_c Vector at \f$ x \f$.
    * \param[in] fun_r Vector at \f$ x+h \f$.
    * \param[in] h Step size.
    * \tparam Vector Fixed-size Eigen vector type.
    * \return Approximated derivative vector.
    */
    template<typename Vector>
    inline Vector CenteredFiniteDifference(Vector const & fun_l, Vector const & fun_c, Vector const & fun_r,
      typename Vector::Scalar h)
    {
      using Real = typename Vector::Scalar;
      constexpr Real EPSILON{std::numeric_limits<Real>::epsilon()};
      Vector diff_r((fun_r - fun_c) / h);
      Vector diff_l((fun_c - fun_l) / h);
      Vector weight_r(diff_r.array().abs() + EPSILON);
      Vector weight_l(diff_l.array().abs() + EPSILON);
      Vector weight_max(weight_r.array().max(weight_l.array()));
      weight_r = (weight_r.array()/weight_max.array()).sqrt();
      weight_l = (weight_l.array()/weight_max.array()).sqrt();
      return (diff_r.array()*weight_l.array() + diff_l.array()*weight_r.array()) / (weight_r.array()+weight_l.array());
    }

    /**
    * Compute finite difference gradient for a scalar function (for fixed-size Eigen vectors).
    * \param[in] x Point at which to compute gradient.
    * \param[in] fun Function to differentiate.
    * \param[out] grad Gradient vector.
    * \tparam Function Callable type <tt>bool(Vector const &, Real &)</tt>.
    * \tparam Vector Fixed-size Eigen vector type.
    * \return True if successful, false otherwise.
    */
    template <typename Function, typename Vector>
    bool Gradient(Vector const & x, Function fun, Vector & grad)
    {
      using Real = typename Vector::Scalar;
      Epsilon<Real> eps;
      Real v_c{0.0};
      if (!fun(x, v_c) || !std::isfinite(v_c)) {return false;}
      constexpr Integer dim_x{Vector::RowsAtCompileTime};
      grad.setZero();
      for (Integer i = 0; i < dim_x; ++i) {
        Vector v_x(x);
        Real tmp{x(i)};
        Real h_1{eps.epsilon_1(tmp)};
        Real h_2{eps.epsilon_2(tmp)};
        v_x(i) = tmp + h_1; Real v_r{0.0}; bool is_finite_r = fun(v_x, v_r) && std::isfinite(v_r);
        v_x(i) = tmp - h_1; Real v_l{0.0}; bool is_finite_l = fun(v_x, v_l) && std::isfinite(v_l);
        Integer ic{(is_finite_r && is_finite_l) ? 0 : (is_finite_r ? 1 : (is_finite_l ? -1 : -2))};
        switch (ic) {
          case 0:
            grad(i) = CenteredFiniteDifference<Real>(v_l, v_c, v_r, h_1);
            break;
          case 1: {
            v_x(i) = tmp + h_2; Real v_rr{0.0}; bool is_finite_rr = fun(v_x, v_rr) && std::isfinite(v_rr);
            grad(i) = is_finite_rr ? SideFiniteDifference<Real>(v_c, v_r, v_rr, h_1, h_2) : (v_r-v_c)/h_1;
            break;
          }
          case -1: {
            v_x(i) = tmp - h_2; Real v_ll{0.0}; bool is_finite_ll = fun(v_x, v_ll) && std::isfinite(v_ll);
            grad(i) = is_finite_ll ? SideFiniteDifference<Real>(v_c, v_l, v_ll, -h_1, -h_2) : (v_c-v_l)/h_1;
            break;
          }
          case -2: {
            grad(i) = 0.0;
            return false;
          }
        }
      }
      return grad.allFinite();
    }

    /**
    * Compute finite difference Jacobian for a vector function (for fixed-size Eigen vectors).
    * \param[in] x Point at which to compute Jacobian.
    * \param[in] fun Function to differentiate.
    * \param[out] jac Jacobian matrix.
    * \tparam Function Callable type <tt>bool(Vector const &, Vector &)</tt>.
    * \tparam Vector Fixed-size Eigen vector type.
    * \tparam Matrix Fixed-size Eigen matrix type.
    * \return True if successful, false otherwise.
    */
    template <typename Function, typename Vector, typename Matrix>
    bool Jacobian(Vector const & x, Function fun, Matrix & jac)
    {
      using Real = typename Vector::Scalar;
      Epsilon<Real> eps;
      Vector v_c;
      if (!fun(x, v_c) || !v_c.allFinite()) {return false;}
      constexpr Integer dim_x{Vector::RowsAtCompileTime};
      jac.setZero();
      for (Integer j{0}; j < dim_x; ++j) {
        Vector v_x(x);
        Real tmp{x(j)};
        Real h_1{eps.epsilon_1(tmp)};
        Real h_2{eps.epsilon_2(tmp)};
        v_x(j) = tmp + h_1; Vector v_r; bool is_finite_r = fun(v_x, v_r) && v_r.allFinite();
        v_x(j) = tmp - h_1; Vector v_l; bool is_finite_l = fun(v_x, v_l) && v_l.allFinite();
        Integer ic{(is_finite_r && is_finite_l) ? 0 : (is_finite_r ? 1 : (is_finite_l ? -1 : -2))};
        switch (ic) {
          case 0:
            jac.col(j) = CenteredFiniteDifference<Real>(v_l, v_c, v_r, h_1);
            break;
          case 1: {
            v_x(j) = tmp + h_2; Vector v_rr; bool is_finite_rr = fun(v_x, v_rr) && v_rr.allFinite();
            jac.col(j) = is_finite_rr ? SideFiniteDifference<Real>(v_c, v_r, v_rr, h_1, h_2) : (v_r-v_c)/h_1;
            break;
          }
          case -1: {
            v_x(j) = tmp - h_2; Vector v_ll; bool is_finite_ll = fun(v_x, v_ll) && v_ll.allFinite();
            jac.col(j) = is_finite_ll ? SideFiniteDifference<Real>(v_c, v_l, v_ll, -h_1, -h_2) : (v_c-v_l)/h_1;
            break;
          }
          case -2: {
            jac.col(j).setZero();
            return false;
          }
        }
      }
      return jac.allFinite();
    }

    /**
    * Compute finite difference Hessian for a scalar function (for fixed-size Eigen vectors).
    * \param[in] x Point at which to compute Hessian.
    * \param[in] fun Function to differentiate.
    * \param[out] hes Hessian matrix.
    * \tparam Function Callable type <tt>bool(Vector const &, Real &)</tt>.
    * \tparam Vector Fixed-size Eigen vector type.
    * \tparam Matrix Fixed-size Eigen matrix type.
    * \return True if successful, false otherwise.
    */
    template <typename Function, typename Vector, typename Matrix>
    bool Hessian(Vector const & x, Function fun, Matrix & hes)
    {
      using Real = typename Vector::Scalar;
      Epsilon<Real> eps;
      constexpr Integer dim_x = Vector::RowsAtCompileTime;
      hes.setZero();
      Real fc{0.0};
      if (!fun(x, fc) || !std::isfinite(fc)) {return false;}
      for (Integer j{0}; j < dim_x; ++j) {
        Real tmp_j{x(j)};
        Real h_j{eps.epsilon_3(tmp_j)};
        Vector v_x = x;
        v_x(j) = tmp_j + h_j; Real fp{0.0}; if (!fun(v_x, fp) || !std::isfinite(fp)) {return false;}
        v_x(j) = tmp_j - h_j; Real fm{0.0}; if (!fun(v_x, fm) || !std::isfinite(fm)) {return false;}
        hes(j, j) = ((fp + fm) - 2.0 * fc) / (h_j * h_j);
        for (Integer i{j + 1}; i < dim_x; ++i) {
          Real tmp_i{x(i)};
          Real h_i{eps.epsilon_3(tmp_i)};
          v_x(i) = tmp_i + h_i; v_x(j) = tmp_j + h_j; Real fpp{0.0}; if (!fun(v_x, fpp) || !std::isfinite(fpp)) {return false;}
          v_x(i) = tmp_i - h_i; Real fmp{0.0}; if (!fun(v_x, fmp) || !std::isfinite(fmp)) {return false;}
          v_x(j) = tmp_j - h_j; Real fmm{0.0}; if (!fun(v_x, fmm) || !std::isfinite(fmm)) {return false;}
          v_x(i) = tmp_i + h_i; Real fpm{0.0}; if (!fun(v_x, fpm) || !std::isfinite(fpm)) {return false;}
          Real h_ij{4.0 * h_i * h_j};
          Real value{((fpp + fmm) - (fpm + fmp)) / h_ij};
          hes(j, i) = hes(i, j) = value;
          v_x(i) = tmp_i;
        }
        v_x(j) = tmp_j;
      }
      return hes.allFinite();
    }

  } // namespace FiniteDifference

} // namespace Optimist

#endif // OPTIMIST_FINITE_DIFFERENCE_HH
