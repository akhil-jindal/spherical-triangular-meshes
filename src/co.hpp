// AMS 562 final project.

// spherical triangle mesh manupilation

// This is the interface for coordinates

#ifndef _CO_HPP
#define _CO_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

namespace ams562_final {

// Class that represents spherical coordinates
// the underlying data structure is a vector of length3 double vectors,
// where each of them is a cartisian coordinates.
// NOTE you can assume that the original is (0,0,0)
class SphCo {
 public:
  // constructor with number of points
  // default is 0
  explicit SphCo(unsigned np = 0u) : _pts(np) {}

  // get n points
  unsigned npoints() const { return _pts.size(); }

  // resize
  void resize(unsigned size) { _pts.resize(size); }

  // data accessing
  std::array<double, 3> &operator[](unsigned i) { return _pts.at(i); }
  const std::array<double, 3> &operator[](unsigned i) const {
    return _pts.at(i);
  }

  // normalize the spherical coordinates with unit Euclidian length 1
  void normalize() {
    // helper lambda that compute the distance, assuming (0,0,0) is the origin
    const auto dist = [](const std::array<double, 3> &pt) {
      return std::sqrt(pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]);
    };
    // reference pt and the dist lambda, define the callback that will be
    // triggered inside for_each
    auto do_nrm = [&](std::array<double, 3> &pt) {
      const double d_inv = 1.0 / dist(pt);
      pt[0] *= d_inv;
      pt[1] *= d_inv;
      pt[2] *= d_inv;
    };
    std::for_each(_pts.begin(), _pts.end(), do_nrm);
  }

  // get reference of underlying std::vector
  std::vector<std::array<double, 3>> &to_vector() { return _pts; }
  const std::vector<std::array<double, 3>> &to_vector() const { return _pts; }

  // given an output streamer, write data to it
  void _write(std::ostream &o_str) const {
    auto write_point = [&o_str](const std::array<double, 3> &pt) {
      o_str << pt[0] << ' ' << pt[1] << ' ' << pt[2] << '\n';
    };
    std::for_each(_pts.begin(), _pts.end(), write_point);
  }

  // given an input streamer, extract data from  it
  void _read(std::istream &i_str) {
    auto read_point = [&i_str](std::array<double, 3> &pt) {
      i_str >> pt[0] >> pt[1] >> pt[2];
    };
    std::for_each(_pts.begin(), _pts.end(), read_point);
  }

 private:
  // underlying datastructure
  std::vector<std::array<double, 3>> _pts;
};

inline std::ostream &operator<<(std::ostream &o_str, const SphCo &co) {
  co._write(o_str);
  return o_str;
}

inline std::istream &operator>>(std::istream &i_str, SphCo &co) {
  co._read(i_str);
  return i_str;
}

}  // namespace ams562_final

#endif
