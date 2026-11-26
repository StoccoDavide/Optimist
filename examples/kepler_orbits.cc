/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Satellites around Earth: two-body orbits propagated by solving Kepler's
// equation
//
//   M = E - e sin(E)
//
// for the eccentric anomaly E at every time step, with Optimist's
// Newton-Raphson and Halley root finders.
//
// Usage: kepler_orbits [output directory]
//
// Writes to the output directory (default: the current one):
//   satellites.csv      orbital elements of each satellite,
//   orbits.csv          positions and iterations every three minutes,
//   kepler_orbits.svg   the orbits drawn to scale around Earth,
//   kepler_altitude.svg the altitude of each satellite over time.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Optimist/RootFinder/Halley.hh"
#include "Optimist/RootFinder/NewtonRaphson.hh"

namespace {

  constexpr double MU           = 398600.4418;  // Earth's gravity, km^3/s^2
  constexpr double EARTH_RADIUS = 6378.137;     // Equatorial radius, km
  constexpr double HOURS        = 72.0;         // Propagation span, h
  constexpr double STEP_MINUTES = 3.0;          // Output step, min
  constexpr double PI           = 3.14159265358979323846;

  // Orbital elements in the orbital plane
  struct Satellite {
    std::string name;
    double a;          // Semi-major axis, km
    double e;          // Eccentricity
    double omega_deg;  // Argument of perigee, deg
    double M0;         // Mean anomaly at t = 0, rad
    std::string color;
  };

  // One propagated position
  struct Sample {
    double t_h, x, y, altitude, E;
    int newton_iterations, halley_iterations;
  };

  double period_hours(const Satellite &s) {
    return 2 * PI * std::sqrt(s.a * s.a * s.a / MU) / 3600;
  }

  // Propagate a satellite by solving Kepler's equation at every step
  std::vector<Sample> propagate(const Satellite &s) {
    Optimist::RootFinder::NewtonRaphson<double> newton;
    Optimist::RootFinder::Halley<double> halley;
    newton.tolerance(1.0e-12);
    halley.tolerance(1.0e-12);
    newton.disable_damped_mode();
    halley.disable_damped_mode();

    const double n{std::sqrt(MU / (s.a * s.a * s.a))};  // Mean motion, rad/s
    const double b{s.a * std::sqrt(1 - s.e * s.e)};      // Semi-minor axis, km
    const double w{s.omega_deg * PI / 180};

    std::vector<Sample> samples;
    for (double t{0}; t <= HOURS * 3600 + 1e-9; t += STEP_MINUTES * 60) {
      const double M{std::remainder(s.M0 + n * t, 2 * PI)};
      auto f   = [&](double E, double &v) { v = E - s.e * std::sin(E) - M; return true; };
      auto df  = [&](double E, double &v) { v = 1 - s.e * std::cos(E); return true; };
      auto ddf = [&](double E, double &v) { v = s.e * std::sin(E); return true; };

      // Standard starting guess, accurate for small eccentricities
      const double E0{M + s.e * std::sin(M)};
      double E_newton, E_halley;
      if (!newton.solve(f, df, E0, E_newton) ||
          !halley.solve(f, df, ddf, E0, E_halley)) {
        throw std::runtime_error("Kepler's equation did not converge for " +
                                 s.name);
      }

      // Position in the orbital plane, rotated by the argument of perigee
      const double xp{s.a * (std::cos(E_newton) - s.e)};
      const double yp{b * std::sin(E_newton)};
      const double x{xp * std::cos(w) - yp * std::sin(w)};
      const double y{xp * std::sin(w) + yp * std::cos(w)};
      samples.push_back({t / 3600, x, y, std::hypot(x, y) - EARTH_RADIUS,
                         E_newton, newton.iterations(), halley.iterations()});
    }
    return samples;
  }

  std::string fmt(const double v, const int digits = 1) {
    std::ostringstream os;
    os.setf(std::ios::fixed);
    os.precision(digits);
    os << v;
    return os.str();
  }

  // Shared SVG frame: white background, so the image reads in any viewer
  std::string svg_open(const int w, const int h) {
    return "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" + std::to_string(w) +
           "\" height=\"" + std::to_string(h) + "\" viewBox=\"0 0 " + std::to_string(w) +
           " " + std::to_string(h) +
           "\" font-family=\"Helvetica, Arial, sans-serif\" font-size=\"13\">\n"
           "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n";
  }

  // The orbits around Earth, to scale, with the satellites at t = 0
  void write_orbits_svg(const std::string &path,
                        const std::vector<Satellite> &sats,
                        const std::vector<std::vector<Sample>> &tracks) {
    // Fit the frame to the orbits and Earth, with equal scales on both axes
    double x_min{-EARTH_RADIUS}, x_max{EARTH_RADIUS}, y_min{-EARTH_RADIUS}, y_max{EARTH_RADIUS};
    for (const auto &tr : tracks) {
      for (const auto &p : tr) {
        x_min = std::min(x_min, p.x);
        x_max = std::max(x_max, p.x);
        y_min = std::min(y_min, p.y);
        y_max = std::max(y_max, p.y);
      }
    }
    const double pad{0.04 * std::max(x_max - x_min, y_max - y_min)};
    x_min -= pad, x_max += pad, y_min -= pad, y_max += pad;
    const int W{640}, head_h{40};
    const double scale{W / (x_max - x_min)};
    const int plot_h{static_cast<int>(std::ceil((y_max - y_min) * scale))};
    const int legend_h{16 + 22 * static_cast<int>(sats.size())};
    const int H{head_h + plot_h + legend_h};
    auto X = [&](double v) { return (v - x_min) * scale; };
    auto Y = [&](double v) { return head_h + (y_max - v) * scale; };

    std::ofstream out(path);
    out << svg_open(W, H);
    out << "<text x=\"16\" y=\"26\" font-size=\"16\" fill=\"#1d2330\">Satellites around Earth, "
        << fmt(HOURS, 0) << " h of two-body motion</text>\n";
    // Earth
    const double r{X(EARTH_RADIUS) - X(0)};
    out << "<circle cx=\"" << fmt(X(0)) << "\" cy=\"" << fmt(Y(0)) << "\" r=\"" << fmt(r)
        << "\" fill=\"#3a6ea5\" fill-opacity=\"0.85\"/>\n"
        << "<text x=\"" << fmt(X(0)) << "\" y=\"" << fmt(Y(0) + 4)
        << "\" text-anchor=\"middle\" fill=\"#ffffff\">Earth</text>\n";
    // Orbits and satellites
    for (std::size_t k{0}; k < sats.size(); ++k) {
      const auto &tr{tracks[k]};
      // Outline of the ellipse, sampled in the eccentric anomaly
      const Satellite &s{sats[k]};
      const double b{s.a * std::sqrt(1 - s.e * s.e)}, w{s.omega_deg * PI / 180};
      out << "<polygon fill=\"none\" stroke=\"" << s.color << "\" stroke-width=\"1.6\" points=\"";
      for (int i{0}; i < 360; ++i) {
        const double E{2 * PI * i / 360};
        const double xp{s.a * (std::cos(E) - s.e)}, yp{b * std::sin(E)};
        out << fmt(X(xp * std::cos(w) - yp * std::sin(w))) << ","
            << fmt(Y(xp * std::sin(w) + yp * std::cos(w))) << " ";
      }
      out << "\"/>\n";
      out << "<circle cx=\"" << fmt(X(tr[0].x)) << "\" cy=\"" << fmt(Y(tr[0].y))
          << "\" r=\"5\" fill=\"" << sats[k].color << "\" stroke=\"#ffffff\" stroke-width=\"1.5\"/>\n";
    }
    // Legend
    for (std::size_t k{0}; k < sats.size(); ++k) {
      const int y{head_h + plot_h + 12 + 22 * static_cast<int>(k)};
      out << "<circle cx=\"24\" cy=\"" << y << "\" r=\"5\" fill=\"" << sats[k].color << "\"/>\n"
          << "<text x=\"36\" y=\"" << y + 4 << "\" fill=\"#1d2330\">" << sats[k].name << "</text>\n"
          << "<text x=\"250\" y=\"" << y + 4 << "\" fill=\"#5b6475\">e = " << fmt(sats[k].e, 3)
          << ", period " << fmt(period_hours(sats[k]), 2) << " h, altitude "
          << fmt(sats[k].a * (1 - sats[k].e) - EARTH_RADIUS, 0) << "–"
          << fmt(sats[k].a * (1 + sats[k].e) - EARTH_RADIUS, 0) << " km</text>\n";
    }
    out << "</svg>\n";
  }

  // Altitude of each satellite over the whole propagation span
  void write_altitude_svg(const std::string &path,
                          const std::vector<Satellite> &sats,
                          const std::vector<std::vector<Sample>> &tracks) {
    const int W{760}, H{380};
    const int l{64}, r{16}, t{40}, b{44};
    double top{0};
    for (const auto &tr : tracks) {
      for (const auto &p : tr) top = std::max(top, p.altitude);
    }
    top = std::ceil(top / 10000) * 10000;
    auto X = [&](double v) { return l + v / HOURS * (W - l - r); };
    auto Y = [&](double v) { return H - b - v / top * (H - t - b); };

    std::ofstream out(path);
    out << svg_open(W, H);
    out << "<text x=\"16\" y=\"24\" font-size=\"16\" fill=\"#1d2330\">Altitude over "
        << fmt(HOURS, 0) << " hours</text>\n";
    // Grid and axes
    for (double v{0}; v <= top + 1; v += 10000) {
      out << "<line x1=\"" << l << "\" x2=\"" << W - r << "\" y1=\"" << fmt(Y(v)) << "\" y2=\""
          << fmt(Y(v)) << "\" stroke=\"#e3e6ec\"/>\n"
          << "<text x=\"" << l - 8 << "\" y=\"" << fmt(Y(v) + 4)
          << "\" text-anchor=\"end\" fill=\"#5b6475\">"
          << (v == 0 ? std::string("0") : fmt(v / 1000, 0) + "k") << " km</text>\n";
    }
    for (double h{0}; h <= HOURS + 1e-9; h += 12) {
      out << "<text x=\"" << fmt(X(h)) << "\" y=\"" << H - b + 18
          << "\" text-anchor=\"middle\" fill=\"#5b6475\">" << fmt(h, 0) << " h</text>\n";
    }
    // Curves
    for (std::size_t k{0}; k < sats.size(); ++k) {
      out << "<polyline fill=\"none\" stroke=\"" << sats[k].color
          << "\" stroke-width=\"1.5\" points=\"";
      for (const auto &p : tracks[k]) out << fmt(X(p.t_h)) << "," << fmt(Y(p.altitude)) << " ";
      out << "\"/>\n";
      const int lx{l + 10 + 170 * static_cast<int>(k)};
      out << "<circle cx=\"" << lx << "\" cy=\"" << H - 10 << "\" r=\"5\" fill=\"" << sats[k].color
          << "\"/>\n<text x=\"" << lx + 10 << "\" y=\"" << H - 6 << "\" fill=\"#1d2330\">"
          << sats[k].name << "</text>\n";
    }
    out << "</svg>\n";
  }

}  // namespace

int main(int argc, char **argv) {
  const std::string dir{argc > 1 ? std::string(argv[1]) + "/" : std::string()};
  const std::vector<Satellite> sats{
    {"Low Earth orbit", 6778.0, 0.001, 0.0, 0.0, "#1f77b4"},
    {"Eccentric orbit", 15000.0, 0.30, 120.0, 1.0, "#ff7f0e"},
    {"Molniya orbit", 26600.0, 0.74, 270.0, 2.0, "#2ca02c"},
    {"Transfer orbit", 24400.0, 0.73, 30.0, 3.0, "#d62728"},
  };

  try {
    std::vector<std::vector<Sample>> tracks;
    for (const auto &s : sats) tracks.push_back(propagate(s));

    std::ofstream elements(dir + "satellites.csv");
    elements << "satellite,a_km,e,omega_deg,period_h,perigee_alt_km,apogee_alt_km\n";
    std::ofstream orbits(dir + "orbits.csv");
    orbits << "satellite,t_h,x_km,y_km,alt_km,E_rad,newton_iterations,halley_iterations\n";
    std::printf("%-18s %6s %9s %12s %12s %8s %8s\n", "Satellite", "e", "Period h",
                "Perigee km", "Apogee km", "Newton", "Halley");
    for (std::size_t k{0}; k < sats.size(); ++k) {
      const auto &s{sats[k]};
      elements << s.name << "," << fmt(s.a) << "," << fmt(s.e, 3) << "," << fmt(s.omega_deg, 0)
               << "," << fmt(period_hours(s), 4) << "," << fmt(s.a * (1 - s.e) - EARTH_RADIUS)
               << "," << fmt(s.a * (1 + s.e) - EARTH_RADIUS) << "\n";
      double newton{0}, halley{0};
      for (const auto &p : tracks[k]) {
        orbits << s.name << "," << fmt(p.t_h, 4) << "," << fmt(p.x, 2) << "," << fmt(p.y, 2)
               << "," << fmt(p.altitude, 2) << "," << fmt(p.E, 9) << ","
               << p.newton_iterations << "," << p.halley_iterations << "\n";
        newton += p.newton_iterations;
        halley += p.halley_iterations;
      }
      std::printf("%-18s %6.3f %9.2f %12.0f %12.0f %8.2f %8.2f\n", s.name.c_str(), s.e,
                  period_hours(s), s.a * (1 - s.e) - EARTH_RADIUS,
                  s.a * (1 + s.e) - EARTH_RADIUS, newton / tracks[k].size(),
                  halley / tracks[k].size());
    }

    write_orbits_svg(dir + "kepler_orbits.svg", sats, tracks);
    write_altitude_svg(dir + "kepler_altitude.svg", sats, tracks);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  return 0;
}
