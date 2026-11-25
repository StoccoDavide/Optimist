/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_TESTSET_DENSE_FUNCTION_HH
#define OPTIMIST_TESTSET_DENSE_FUNCTION_HH

#include "Optimist/Function.hh"
#include "Optimist/Utils.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Helper base class for scalar cost functions \f$ f : \mathbb{R}^N
     * \to \mathbb{R} \f$ written once in terms of dense vectors.
     *
     * The derived class implements
     * - <tt>Scalar value(const Dense &x) const</tt>,
     * - <tt>void gradient(const Dense &x, Dense &g) const</tt>,
     * - <tt>void hessian(const Dense &x, DenseMatrix &H) const</tt>,
     *
     * and this class takes care of the conversion from/to fixed-size,
     * dynamic-size and sparse Eigen vectors.
     *
     * \tparam Vector Eigen vector type.
     * \tparam N Problem dimension.
     * \tparam Derived Derived function class.
     */
    template <typename Vector, Integer N, typename Derived>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension == N)
    class DenseCostFunction
        : public Function<Vector, typename Vector::Scalar, Derived> {
     public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using Dense       = Utils::DenseVector<Vector>;
      using DenseMatrix = Utils::DenseMatrix<Vector>;
      using typename Function<Vector, Scalar, Derived>::FirstDerivative;
      using typename Function<Vector, Scalar, Derived>::SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Problem dimension.
       */
      static constexpr Integer dimension() {
        return N;
      }

      /**
       * Evaluate the cost function.
       */
      bool evaluate_impl(const Vector &x, Scalar &out) const {
        out = this->derived().value(Utils::to_dense(x));
        return std::isfinite(out);
      }

      /**
       * Evaluate the gradient of the cost function.
       */
      bool first_derivative_impl(const Vector &x, FirstDerivative &out) const {
        Dense g(N);
        this->derived().gradient(Utils::to_dense(x), g);
        out = Utils::from_dense<FirstDerivative>(g);
        return g.allFinite();
      }

      /**
       * Evaluate the Hessian of the cost function.
       */
      bool second_derivative_impl(const Vector &x,
                                  SecondDerivative &out) const {
        DenseMatrix H(N, N);
        this->derived().hessian(Utils::to_dense(x), H);
        out = Utils::from_dense_matrix<SecondDerivative>(H);
        return H.allFinite();
      }

     protected:
      /**
       * Add a known solution given its dense coefficients.
       */
      void add_solution(const Dense &x) {
        this->m_solutions.push_back(Utils::from_dense<Vector>(x));
      }

      /**
       * Add a suggested initial guess given its dense coefficients.
       */
      void add_guess(const Dense &x) {
        this->m_guesses.push_back(Utils::from_dense<Vector>(x));
      }

      /**
       * Build a dense vector from a list of coefficients.
       */
      static Dense make(std::initializer_list<Scalar> values) {
        Dense x(static_cast<Integer>(values.size()));
        Integer i{0};
        for (const auto v : values) {
          x(i++) = v;
        }
        return x;
      }

     private:
      const Derived &derived() const {
        return *static_cast<const Derived *>(this);
      }

    };  // class DenseCostFunction

    /**
     * \brief Helper base class for nonlinear systems \f$ \mathbf{F} :
     * \mathbb{R}^N \to \mathbb{R}^M \f$ written once in terms of dense vectors
     * (square by default).
     *
     * The derived class implements
     * - <tt>void residual(const Dense &x, Dense &F) const</tt>,
     * - <tt>void jacobian(const Dense &x, DenseMatrix &J) const</tt>,
     * - <tt>void hessian(const Dense &x, Integer i, DenseMatrix &H) const</tt>,
     *   the Hessian of the \f$ i \f$-th component.
     *
     * \tparam Vector Eigen input vector type.
     * \tparam N Number of unknowns.
     * \tparam Derived Derived function class.
     * \tparam OutVector Eigen output vector type.
     * \tparam M Number of equations.
     */
    template <typename Vector,
              Integer N,
              typename Derived,
              typename OutVector = Vector,
              Integer M          = N>
      requires TypeTrait<Vector>::IsEigen && TypeTrait<OutVector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed ||
                TypeTrait<Vector>::Dimension == N) &&
               (!TypeTrait<OutVector>::IsFixed ||
                TypeTrait<OutVector>::Dimension == M)
    class DenseSystemFunction : public Function<Vector, OutVector, Derived> {
     public:
      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using Dense       = Utils::DenseVector<Vector>;
      using DenseMatrix = Utils::DenseMatrix<Vector>;
      using typename Function<Vector, OutVector, Derived>::FirstDerivative;
      using typename Function<Vector, OutVector, Derived>::SecondDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Number of unknowns.
       */
      static constexpr Integer dimension() {
        return N;
      }

      /**
       * Number of equations.
       */
      static constexpr Integer equations() {
        return M;
      }

      /**
       * Evaluate the system residual.
       */
      bool evaluate_impl(const Vector &x, OutVector &out) const {
        Dense F(M);
        this->derived().residual(Utils::to_dense(x), F);
        out = Utils::from_dense<OutVector>(F);
        return F.allFinite();
      }

      /**
       * Evaluate the system Jacobian.
       */
      bool first_derivative_impl(const Vector &x, FirstDerivative &out) const {
        DenseMatrix J(M, N);
        this->derived().jacobian(Utils::to_dense(x), J);
        out = Utils::from_dense_matrix<FirstDerivative>(J);
        return J.allFinite();
      }

      /**
       * Evaluate the Hessians of the system components.
       */
      bool second_derivative_impl(const Vector &x,
                                  SecondDerivative &out) const {
        using Matrix = typename SecondDerivative::value_type;
        const Dense xd(Utils::to_dense(x));
        out.resize(M);
        bool finite{true};
        for (Integer i{0}; i < M; ++i) {
          DenseMatrix H(N, N);
          this->derived().hessian(xd, i, H);
          finite = finite && H.allFinite();
          out[i] = Utils::from_dense_matrix<Matrix>(H);
        }
        return finite;
      }

     protected:
      /**
       * Add a known solution given its dense coefficients.
       */
      void add_solution(const Dense &x) {
        this->m_solutions.push_back(Utils::from_dense<Vector>(x));
      }

      /**
       * Add a suggested initial guess given its dense coefficients.
       */
      void add_guess(const Dense &x) {
        this->m_guesses.push_back(Utils::from_dense<Vector>(x));
      }

      /**
       * Build a dense vector from a list of coefficients.
       */
      static Dense make(std::initializer_list<Scalar> values) {
        Dense x(static_cast<Integer>(values.size()));
        Integer i{0};
        for (const auto v : values) {
          x(i++) = v;
        }
        return x;
      }

     private:
      const Derived &derived() const {
        return *static_cast<const Derived *>(this);
      }

    };  // class DenseSystemFunction

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_DENSE_FUNCTION_HH
