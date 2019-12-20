// AMS 562 final project.

// TASK 1:
// spherical triangle mesh manupilation

// This is the mesh option interface

#ifndef _MANIP_HPP
#define _MANIP_HPP

#include "co.hpp"
#include "conn.hpp"

namespace ams562_final {

/// \brief compute node two element/triangle adjacency
/// \param[in] n number of points
/// \param[in] conn triangle connectivity
/// \param[out] adj adjacency mapping
///
/// On output, the \a adj stores a mapping that maps a node to its adjacent
/// elements/triangles. \a adj has length of \a n.
void compute_n2e_adj(const unsigned n, const Triangles &conn,
                     std::vector<std::vector<int>> &adj);

/// \brief compute averaged outward normals
/// \param[in] points coordinates
/// \param[in] conn triangle connectivity table
/// \param[in] n2e_adj node to triangles adjacency mapping, given by compute_n2e_adj
/// \param[out] nrms output averaged normals
///
/// \a points are passed in as normalized coordinates, i.e. their Euclidian
/// distances with respect to the origin (0.0.0) is 1. You should do the same
/// thing for \a nrms before finishing this function.
void compute_avg_normals(const SphCo &points, const Triangles &conn,
                         const std::vector<std::vector<int>> &n2e_adj,
                         SphCo &nrms);

/// \brief compute the errors in terms of the angles (rad)
/// \param[in] exact analytic normal vectors, same as the coordinates
/// \param[in] num numerical normal vectors, computed by compute_avg_normals
/// \param[out] acos_theta the angle errors between them
///
/// It's worth noting that under our setting, the coordinates are nothing but
/// the analytic outward normal vectors. You need to compute the angles (in rad)
/// between the numerical and analytic vectors. Note that you can assume that
/// \a exact is normalized and \a num should be normalized inside
/// compute_avg_normals.
void compute_errors(const SphCo &exact, const SphCo &num,
                    std::vector<double> &acos_theta);

}  // namespace ams562_final

#endif
