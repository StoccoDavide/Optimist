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

#ifndef OPTIMIST_TIC_TOC_HXX
#define OPTIMIST_TIC_TOC_HXX

namespace Optimist {

  /*\
   |   _____ _     _____
   |  |_   _(_) __|_   _|__   ___
   |    | | | |/ __|| |/ _ \ / __|
   |    | | | | (__ | | (_) | (__
   |    |_| |_|\___||_|\___/ \___|
   |
  \*/

  /**
  * \brief Class container for a timer.
  *
  * Class container for a timer that uses the C++ standard library chrono to measure the elapsed time
  * between two points in the code.
  */
  class TicToc {

    using Clock = std::chrono::high_resolution_clock; /**< Clock type. */
    using Time  = std::chrono::microseconds;          /**< Time type. */

    Clock::time_point m_start_time;   /**< Start time. */
    Clock::time_point m_stop_time;    /**< Stop time. */
    Time              m_elapsed_time; /**< Elapsed time. */

  public:
    /**
    * Copy constructor for the timer.
    */
    TicToc(const TicToc &) = delete;

    /**
    * Assignment operator for the timer.
    */
    TicToc & operator=(TicToc const &) = delete;

    /**
    * Class constructor for the timer.
    */
    TicToc() : m_elapsed_time(0) {}

    /**
    * Start the timer.
    */
    void tic() {this->m_start_time = Clock::now();}

    /**
    * Stop the timer.
    */
    void toc()
    {
      this->m_stop_time    = Clock::now();
      this->m_elapsed_time = std::chrono::duration_cast<Time>(
        this->m_stop_time - this->m_start_time
      );
    }

    /**
    * Get the elapsed time in seconds.
    * \return The elapsed time in seconds.
    */
    Real elapsed_s() const {return Real(1.0e-6)*this->m_elapsed_time.count();}

    /**
    * Get the elapsed time in milliseconds.
    * \return The elapsed time in milliseconds.
    */
    Real elapsed_ms() const {return Real(1.0e-3)*this->m_elapsed_time.count();}

    /**
    * Get the elapsed time in microseconds.
    * \return The elapsed time in microseconds.
    */
    Real elapsed_us() const {return this->m_elapsed_time.count();}

  }; // class TicToc

} // namespace Optimist

#endif // OPTIMIST_TIC_TOC_HXX
