#include <iostream>

// you will need cmath and algorithm headers.
#include <cmath>
#include <algorithm>

#include "manip.hpp"


namespace ams562_final {

void compute_n2e_adj(const unsigned n, const Triangles &conn,
                     std::vector<std::vector<int>> &adj) {
  // resize adj to n
  adj.resize(n);
  // reserve for each of them with a reasonable upper bound

  unsigned count = 0;

  auto add_to_adj = [&count,&adj](const std::array<int, 3> &tri) {
    adj[tri[0]].push_back(count);
    adj[tri[1]].push_back(count);
    adj[tri[2]].push_back(count);
    count++;
  };
  // for_each #1
  std::for_each(conn.to_vector().begin(), conn.to_vector().end(), add_to_adj);
}

void compute_avg_normals(const SphCo &points, const Triangles &conn,
                         const std::vector<std::vector<int>> &n2e_adj,
                         SphCo &nrms) {
  // resize the nrms
  nrms.resize(points.npoints());

  // number of triangles
  unsigned m = conn.ntris();
  // initialize workspace for facial norms (m x 3)
  ams562_final::SphCo nrmsf(m);

  // for each triangle in conn
  for (unsigned i = 0; i < m; i++)
  {
    //The cross product of two sides of the triangle equals the surface normal. 
    //So, if V = B - A and W = C - A, and N is the surface normal, then:
    //Nx=(Vy∗Wz)−(Vz∗Wy) 
    //Ny=(Vz∗Wx)−(Vx∗Wz)
    //Nz=(Vx∗Wy)−(Vy∗Wx)
    std::array<int, 3> arr = conn[i];
    std::array<double, 3> ptA = points[arr[0]];
    std::array<double, 3> ptB = points[arr[1]];
    std::array<double, 3> ptC = points[arr[2]];

    double Vx = ptB[0] - ptA[0];
    double Vy = ptB[1] - ptA[1];
    double Vz = ptB[2] - ptA[2];

    double Wx = ptC[0] - ptA[0];
    double Wy = ptC[1] - ptA[1];
    double Wz = ptC[2] - ptA[2];

    std::array<double, 3> &normfArr = nrmsf[i];
    normfArr[0] = (Vy*Wz) - (Vz*Wy);
    normfArr[1] = (Vz*Wx) - (Vx*Wz);
    normfArr[2] = (Vx*Wy) - (Vy*Wx);
  }
  //normalize the surface normals
  nrmsf.normalize();

  // perform averaging
  for (unsigned i = 0; i < points.npoints(); i++)
  {
    std::array<double, 3> &nrmsArr = nrms[i];
    unsigned k = n2e_adj[i].size();
    for (unsigned j = 0; j < k; j++)
    {
        nrmsArr[0] += nrmsf[n2e_adj[i][j]][0];
        nrmsArr[1] += nrmsf[n2e_adj[i][j]][1];
        nrmsArr[2] += nrmsf[n2e_adj[i][j]][2];
    }
    nrmsArr[0] /= k;
    nrmsArr[1] /= k;
    nrmsArr[2] /= k;
  }
  // hint don't forget normalizing the normals
  nrms.normalize();
}

void compute_errors(const SphCo &exact, const SphCo &num,
                    std::vector<double> &acos_theta) {
  // resize the error array
  acos_theta.resize(num.npoints());

  unsigned count = 0;
  auto compute_theta = [&count,&exact,&acos_theta](const std::array<double, 3> &pt) {
    acos_theta.push_back(acos(pt[0]*exact[count][0] +
                              pt[1]*exact[count][1] +
                              pt[2]*exact[count][2]));
    count++;
  };
  // for_each #2
  std::for_each(num.to_vector().begin(), num.to_vector().end(), compute_theta);
}

}  // namespace ams562_final
