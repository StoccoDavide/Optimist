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

#ifndef OPTIMIST_FINITE_DIFFERENCE_HXX
#define OPTIMIST_FINITE_DIFFERENCE_HXX

namespace Optimist
{

  namespace FiniteDifference
  {

    class Epsilon
    {
      Real const m_epsilon{std::numeric_limits<Real>::epsilon()};
      Real const m_epsilon1{std::sqrt(std::numeric_limits<Real>::epsilon())};
      Real const m_epsilon2{std::pow(std::numeric_limits<Real>::epsilon(), 0.75)};
      Real const m_epsilon3{std::pow(std::numeric_limits<Real>::epsilon(), 0.25)};
    public:
      Real epsilon1( Real v ) const { return (abs(v)+1)*m_epsilon1; }
      Real epsilon2( Real v ) const { return (abs(v)+1)*m_epsilon2; }
      Real epsilon3( Real v ) const { return (abs(v)+1)*m_epsilon3; }
    };

    /*\
     |   ____  _     _
     |  / ___|(_) __| | ___
     |  \___ \| |/ _` |/ _ \
     |   ___) | | (_| |  __/
     |  |____/|_|\__,_|\___|
     |
    \*/

    static
    Real
    finite_difference_side(
      Real f0,
      Real f1,
      Real f2,
      Real h1,
      Real h2
    ) {
      Real df1{f1 - f0};
      Real df2{f2 - f0};
      Real dH{h1 - h2};
      return ((h1/h2)*df2 - (h2/h1)*df1) / dH;
    }

    template <Integer N>
    static
    bool
    finite_difference_side(
      const Vector<N> & f0,
      const Vector<N> & f1,
      const Vector<N> & f2,
      const Real        h1,
      const Real        h2,
      Vector<N> &       diff
    ) {
      Real dH{h1 - h2};
      Real t1{-(h2/h1) / dH};
      Real t2{(h1/h2) / dH};
      diff = t1 * (f1 - f0) + t2 * (f2 - f0);
      return diff.allFinite();
    }

    /*\
     |    ____           _                    _
     |   / ___|___ _ __ | |_ ___ _ __ ___  __| |
     |  | |   / _ \ '_ \| __/ _ \ '__/ _ \/ _` |
     |  | |__|  __/ | | | ||  __/ | |  __/ (_| |
     |   \____\___|_| |_|\__\___|_|  \___|\__,_|
     |
    \*/

    static
    Real
    finite_difference_centered(
      const Real fun_l,
      const Real fun_c,
      const Real fun_r,
      const Real h
    ) {
      Real diff_r{(fun_r-fun_c)/h};
      Real diff_l{(fun_c-fun_l)/h};
      Real weight_r{std::abs(diff_r) + EPSILON};
      Real weight_l{std::abs(diff_l) + EPSILON};
      Real weight_max{std::max(weight_r, weight_l)};
      weight_r = std::sqrt(weight_r/weight_max);
      weight_l = std::sqrt(weight_l/weight_max);
      return (diff_r*weight_l + diff_l*weight_r) / (weight_r+weight_l);
    }

    static
    Integer
    finite_difference_centered(
      const Real fun_l, const bool is_finite_l,
      const Real fun_c, const bool is_finite_c,
      const Real fun_r, const bool is_finite_r,
      const Real h,
      Real & diff
    ) {

      if (is_finite_r && is_finite_l && is_finite_c) {
        diff = finite_difference_centered(fun_l, fun_c, fun_r, h);
        return 0;
      } else if (is_finite_c && is_finite_r) {
        diff = (fun_r-fun_c)/h;
        return 1;
      } else if (is_finite_c && is_finite_l) {
        diff = (fun_c-fun_l)/h;
        return -1;
      }
      diff = 0;
      return -2;
    }

    template <Integer N>
    static
    Integer
    finite_difference_centered(
      const Vector<N> & fun_l, bool is_finite_l,
      const Vector<N> & fun_c, bool is_finite_c,
      const Vector<N> & fun_r, bool is_finite_r,
      Real h,
      Real diff
    ) {
      if ( is_finite_r && is_finite_l && is_finite_c  ) {
        Vector<N> diff_r((fun_r - fun_c) / h);
        Vector<N> diff_l((fun_c - fun_l) / h);
        Vector<N> weight_r{diff_r.array().abs() + EPSILON};
        Vector<N> weight_l{diff_l.array().abs() + EPSILON};
        Vector<N> weight_max = weight_r.array().max(weight_l.array());
        weight_r = (weight_r/weight_max).array().sqrt();
        weight_l = (weight_l/weight_max).array().sqrt();
        diff = (diff_r*weight_l + diff_l*weight_r) / (weight_r+weight_l);
        return 0;
      } else if ( is_finite_c && is_finite_r ) {
        diff = (fun_r-fun_c)/h;
        return 1;
      } else if ( is_finite_c && is_finite_l ) {
        diff = (fun_c-fun_l)/h;
        return -1;
      }
      diff.setZero();
      return -2;
    }

    /*\
     |    ____               _ _            _
     |   / ___|_ __ __ _  __| (_) ___ _ __ | |_
     |  | |  _| '__/ _` |/ _` | |/ _ \ '_ \| __|
     |  | |_| | | | (_| | (_| | |  __/ | | | |_
     |   \____|_|  \__,_|\__,_|_|\___|_| |_|\__|
     |
    \*/

    template <typename FUNCTION, typename Real>
    static
    bool
    finite_difference_gradient(
      Real const x[],
      Integer      dim_x,
      FUNCTION     fun,
      Real       grad[]
    ) {

      static FiniteDifferenceEpsilon<Real> EPS;

      auto FUN = [&fun]( Real const x[], Real & f ) -> bool {
        bool ok{fun( x, f )};
        if ( ok ) ok = std::is_finite(f);
        return ok;
      };

      Real vC{0}, vR{0}, vL{0}; // only to stop warning
      bool is_finite_c = FUN( x, vC );
      if ( !is_finite_c ) return false;

      Real * X = const_cast<Real*>(x);

      for ( Integer i{0}; i < dim_x; ++i ) {
        Real temp{x[i]};
        Real h1{EPS.epsilon1(temp)};
        Real h2{EPS.epsilon2(temp)};
        X[i] = temp+h1; bool is_finite_r{FUN( X, vR )};
        X[i] = temp-h1; bool is_finite_l{FUN( X, vL )};
        Integer ic = finite_difference_centered( vL, is_finite_l, vC, is_finite_c, vR, is_finite_r, h1, grad[i] );
        switch ( ic ) {
        case  0:
          break;
        case  1:
          {
            Real & vRR = vL;
            X[i] = temp+h2; // modify the vector only at i position
            bool is_finite_rr{FUN( X, vRR )};
            if ( is_finite_rr ) grad[i] = finite_difference_side( vC, vR, vRR, h1, h2 );
            if ( ! (is_finite_rr&&std::is_finite(grad[i])) ) grad[i] = (vR-vC)/h1; // low precision FD
          }
          break;
        case -1:
          {
            Real & vLL = vR;
            X[i] = temp-h2; // modify the vector only at i position
            bool is_finite_ll{FUN( X, vLL )};
            if ( is_finite_ll ) grad[i] = finite_difference_side( vC, vL, vLL, -h1, -h2 );
            if ( ! (is_finite_ll&&std::is_finite(grad[i])) ) grad[i] = (vC-vL)/h1; // low precision FD
          }
          break;
        case -2:
          return false;
          break;
        }
        X[i] = temp; // restore i position
        if ( !std::is_finite(grad[i]) ) return false;
      }
      return true;
    }

    /*\
    |       _                 _     _
    |      | | __ _  ___ ___ | |__ (_) __ _ _ __
    |   _  | |/ _` |/ __/ _ \| '_ \| |/ _` | '_ \
    |  | |_| | (_| | (_| (_) | |_) | | (_| | | | |
    |   \___/ \__,_|\___\___/|_.__/|_|\__,_|_| |_|
    |
    \*/

    template <typename FUNCTION, typename Real>
    static
    bool
    finite_difference_jacobian(
      Real const x[],
      Integer      dim_x,
      FUNCTION     fun,
      Integer      dim_f,
      Real       Jac[],
      Integer      ldJ,
      Real *     work,
      Integer      lwork
    ) {

      static FiniteDifferenceEpsilon<Real> EPS;

      UTILS_ASSERT(
        lwork >= 3*dim_f,
        "finite_difference_jacobian(...,dim_f={},...,lwork={}), lwork must be >= 3*dim_f\n",
        dim_f, lwork
      );

      Real * vC{work};
      Real * vR{work+dim_f};
      Real * vL{work+2*dim_f};

      auto FUN = [&fun,dim_f]( Real const x[], Real f[] ) -> bool {
        bool ok{fun( x, f )};
        for ( Integer i{0}; ok && i < dim_f; ++i ) ok = std::is_finite(f[i]);
        return ok;
      };

      bool is_finite_c{FUN( x, vC )};
      if ( !is_finite_c ) return false;

      Real * X{const_cast<Real*>(x)};
      Real * pjac{Jac};

      for ( Integer j{0}; j < dim_x; ++j ) {
        Real temp{x[j]};
        Real h1{EPS.epsilon1(temp)};
        Real h2{EPS.epsilon2(temp)};

        X[j] = temp+h1; bool is_finite_r{FUN( X, vR )};
        X[j] = temp-h1; bool is_finite_l{FUN( X, vL )};

        Integer ic = finite_difference_centered( vL, is_finite_l, vC, is_finite_c, vR, is_finite_r, h1, pjac, dim_f );

        switch ( ic ) {
        case  0:
          break;
        case  1:
          {
            Real * vRR{vL};
            X[j] = temp+h2; // modify the vector only at i position
            bool is_finite_rr{FUN( X, vRR )};
            if ( is_finite_rr ) is_finite_rr = finite_difference_side( vC, vR, vRR, h1, h2, pjac, dim_f );
            if ( !is_finite_rr ) {
              for ( Integer i{0}; i < dim_f; ++i ) {
                pjac[i] = (vR[i]-vC[i])/h1; // low precision FD
                if ( !std::is_finite( pjac[i] ) ) return false;
              }
            }
          }
          break;
        case -1:
          {
            Real * vLL{vR};
            X[j] = temp-h2; // modify the vector only at i position
            bool is_finite_ll{FUN( X, vLL )};
            if ( is_finite_ll ) is_finite_ll = finite_difference_side( vC, vL, vLL, -h1, -h2, pjac, dim_f );
            if ( !is_finite_ll ) {
              for ( Integer i{0}; i < dim_f; ++i ) {
                pjac[i] = (vC[i]-vL[i])/h1; // low precision FD
                if ( !std::is_finite( pjac[i] ) ) return false;
              }
            }
          }
          break;
        case -2:
          return false;
          break;
        }
        X[j] = temp; // restore i position
        pjac += ldJ;
      }
      return true;
    }

    /*\
    |   _   _               _
    |  | | | | ___  ___ ___(_) __ _ _ __
    |  | |_| |/ _ \/ __/ __| |/ _` | '_ \
    |  |  _  |  __/\__ \__ \ | (_| | | | |
    |  |_| |_|\___||___/___/_|\__,_|_| |_|
    |
    \*/

    template <typename FUNCTION, typename Real>
    static
    bool
    finite_difference_hessian(
      Real const x[],
      Integer      dim_x,
      FUNCTION     fun,
      Real       Hess[],
      Integer      ldH
    ) {
      static FiniteDifferenceEpsilon<Real> EPS;

      auto FUN = [&fun]( Real const x[], Real & f ) -> bool {
        bool ok{fun( x, f )};
        if ( ok ) ok = std::is_finite(f);
        return ok;
      };

      Real fp{0}, fm{0}, fc{0}, hij{0}, tempi{0}, tempj{0}, hi{0}, hj{0},
            fpp{0}, fpm{0}, fmp{0}, fmm{0};
      bool ok{true};

      Real * X{const_cast<Real*>(x)};
      for ( Integer j{0}; j < dim_x && ok; ++j ) {
        tempj = x[j];
        hj    = EPS.epsilon3(tempj);
        ok    = FUN( X, fc );
        if ( !ok ) goto skip;
        X[j] = tempj+hj;
        ok   = FUN( X, fp );
        if ( !ok ) goto skip;
        X[j] = tempj-hj;
        ok   = FUN( X, fm );
        if ( !ok ) goto skip;
        Hess[j*(ldH+1)] = ((fp+fm)-2*fc)/(hj*hj);
        for ( Integer i{j+1}; i < dim_x && ok; ++i ) {
          tempi = X[i];
          hi    = EPS.epsilon3(tempi);
          X[i]  = tempi+hi;
          X[j]  = tempj+hj;
          ok    = FUN( X, fpp );
          if ( !ok ) goto skip2;
          X[i] = tempi-hi;
          ok   = FUN( X, fmp );
          if ( !ok ) goto skip2;
          X[j] = tempj-hj;
          ok   = FUN( X, fmm );
          if ( !ok ) goto skip2;
          X[i] = tempi+hi;
          ok   = FUN( X, fpm );
          if ( !ok ) goto skip2;
          hij = 4*hi*hj;
          Hess[j+i*ldH] = Hess[i+j*ldH] = ( (fpp+fmm) - (fpm+fmp) )/hij;
        skip2:
          X[i] = tempi;
        }
      skip:
        X[j] = tempj;
      }
      return ok;
    }

  } // namespace FiniteDifference

} // namespace Optimist

#endif // OPTIMIST_FINITE_DIFFERENCE_HXX
