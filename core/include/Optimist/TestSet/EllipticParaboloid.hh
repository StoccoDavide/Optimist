

#pragma once

#ifndef OPTIMIST_TESTSET_ELLIPTIC_PARABOLOID_HH
#define OPTIMIST_TESTSET_ELLIPTIC_PARABOLOID_HH

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Class container for the elliptic paraboloid.
     *
     * Convex quadratic cost function
     * \f[
     *   f(\mathbf{x}) = a x_1^2 + b x_2^2 \text{,}
     * \f]
     * with \f$ a = b = 1 \f$ and minimum \f$ f(\mathbf{0}) = 0 \f$. The
     * guesses lie on a regular grid of the square \f$ [-10, 10]^2 \f$.
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    class EllipticParaboloid
        : public DenseCostFunction<Vector, 2, EllipticParaboloid<Vector>> {
     public:
      using Base = DenseCostFunction<Vector, 2, EllipticParaboloid<Vector>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

     private:
      Scalar m_a{1.0}; /**< Coefficient \f$ a \f$. */
      Scalar m_b{1.0}; /**< Coefficient \f$ b \f$. */

     public:
      /**
       * Class constructor for the elliptic paraboloid.
       */
      EllipticParaboloid() {
        this->add_solution(Dense::Zero(2));
        const Scalar step{4};
        for (Scalar x{-10}; x < 10 + Base::EPSILON; x += step) {
          for (Scalar y{-10}; y < 10 + Base::EPSILON; y += step) {
            this->add_guess(this->make({x, y}));
          }
        }
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "EllipticParaboloid";
      }

      Scalar value(const Dense &x) const {
        return this->m_a * x(0) * x(0) + this->m_b * x(1) * x(1);
      }

      void gradient(const Dense &x, Dense &g) const {
        g << 2 * this->m_a * x(0), 2 * this->m_b * x(1);
      }

      void hessian([[maybe_unused]] const Dense &x, DenseMatrix &H) const {
        H << 2 * this->m_a, 0, 0, 2 * this->m_b;
      }

    };  // class EllipticParaboloid

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_ELLIPTIC_PARABOLOID_HH
