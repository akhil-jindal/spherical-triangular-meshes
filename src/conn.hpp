// AMS 562 final project.

// spherical triangle mesh manupilation

// This is the interface for connectivity

#ifndef _CONN_HPP
#define _CONN_HPP

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

namespace ams562_final {

// triangle representation
// the connectivity table is defines as a vector of triangles
// where a triangle is defined as a length3 array of integers that
// are the indices of the coordinates
class Triangles {
 public:
  // constructor with number of triangles
  explicit Triangles(unsigned ne = 0u) : _conn(ne) {}

  // get number of triangles
  unsigned ntris() const { return _conn.size(); }

  // resize the connectivity table
  void resize(unsigned size) { _conn.resize(size); }

  // access to i-th triangle
  std::array<int, 3> &operator[](unsigned i) { return _conn.at(i); }

  // access to i-th triangle with constant specifier
  const std::array<int, 3> &operator[](unsigned i) const { return _conn.at(i); }

  // get the underlying std::vector
  std::vector<std::array<int, 3>> &to_vector() { return _conn; }
  const std::vector<std::array<int, 3>> &to_vector() const { return _conn; }

  // write the contents to streamer with white space delimiter
  void _write(std::ostream &o_str) const {
    // define the lambda that references o_str, for each of the triangle
    // write them with a newline
    auto write_tri = [&o_str](const std::array<int, 3> &tri) {
      o_str << tri[0] << ' ' << tri[1] << ' ' << tri[2] << '\n';
    };
    // invoke for_each
    std::for_each(_conn.begin(), _conn.end(), write_tri);
  }

  void _read(std::istream &i_str) {
    // read data from input streamer to the triangle. Note that the tri
    // is passed in with normal reference
    auto read_tri = [&i_str](std::array<int, 3> &tri) {
      i_str >> tri[0] >> tri[1] >> tri[2];
    };
    std::for_each(_conn.begin(), _conn.end(), read_tri);
  }

 private:
  // connectivity data structure
  std::vector<std::array<int, 3>> _conn;
};

inline std::ostream &operator<<(std::ostream &o_str, const Triangles &conn) {
  conn._write(o_str);
  return o_str;
}

inline std::istream &operator>>(std::istream &i_str, Triangles &conn) {
  conn._read(i_str);
  return i_str;
}

}  // namespace ams562_final

#endif