

#pragma once

#ifndef OPTIMIST_TESTSET_BOOTH_HH
#define OPTIMIST_TESTSET_BOOTH_HH

#include "Optimist/TestSet/DenseFunction.hh"

namespace Optimist {

  namespace TestSet {

    /**
     * \brief Class container for the Booth function.
     *
     * Linear system of two equations,
     * \f[
     *   \mathbf{f}(\mathbf{x}) = \begin{bmatrix}
     *     x_1 + 2x_2 - 7 \\
     *     2x_1 + x_2 - 5
     *   \end{bmatrix} \text{,}
     * \f]
     * with solution \f$ \mathbf{x} = (1, 3) \f$.
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    class Booth : public DenseSystemFunction<Vector, 2, Booth<Vector>> {
     public:
      using Base = DenseSystemFunction<Vector, 2, Booth<Vector>>;
      using typename Base::Dense;
      using typename Base::DenseMatrix;
      using typename Base::Scalar;

      /**
       * Class constructor for the Booth function.
       */
      Booth() {
        this->add_solution(this->make({1.0, 3.0}));
        this->add_guess(this->make({-10.0, -10.0}));
        this->add_guess(this->make({10.0, 10.0}));
      }

      /**
       * Get the function name.
       * \return The function name.
       */
      constexpr std::string name_impl() const {
        return "Booth";
      }

      void residual(const Dense &x, Dense &F) const {
        F(0) = x(0) + 2 * x(1) - 7;
        F(1) = 2 * x(0) + x(1) - 5;
      }

      void jacobian([[maybe_unused]] const Dense &x, DenseMatrix &J) const {
        J << 1, 2, 2, 1;
      }

      void hessian([[maybe_unused]] const Dense &x,
                   [[maybe_unused]] const Integer i,
                   DenseMatrix &H) const {
        H.setZero();
      }

    };  // class Booth

  }  // namespace TestSet

}  // namespace Optimist

#endif  // OPTIMIST_TESTSET_BOOTH_HH
