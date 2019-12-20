#include "src/manip.hpp"

#include <fstream>

// rad to degree conversion
#define RAD2DEG 57.2958

// thresholds for three cases
const static double thresholds[3] = {4.0, 3.0, 2.0};

#define FAIL "\033[0;31mFAIL\033[0m"
#define PASS "\033[0;32mPASS\033[0m"

int main() {
  std::cout << "Akhil Jindal - AMS 562 Final Project" << std::endl;
  for (int kase = 1; kase <= 3; ++kase) {
    std::cout << "Enterring test case " << kase << "\n\n";
    const std::string filename = "test/mesh" + std::to_string(kase) + ".txt";
    std::ifstream file(filename.c_str());
    if (!file.is_open()) {
      std::cerr << "ERROR! Cannot open file " << filename << '\n';
      return 1;
    }
    // read in our number of points and triangles
    int np, ne;
    file >> np >> ne;
    std::cout << "Number of points: " << np << '\n'
              << "Number of triangles: " << ne << '\n';

    // step 1, load the data
    ams562_final::SphCo points(np);
    ams562_final::Triangles conn(ne);
    // normalize it
    points.normalize();
    std::cout << "Loading the geometric and topological data sets...\n";
    file >> points >> conn;

    // step 2, compute the node to cell adjacency
    std::vector<std::vector<int>> adj;
    std::cout << "Computing node to triangles adjacency...\n";
    ams562_final::compute_n2e_adj(np, conn, adj);

    // step 3, compute the numerical normals
    ams562_final::SphCo nrms;
    std::cout << "Computing the averaged outward normal vectors...\n";
    ams562_final::compute_avg_normals(points, conn, adj, nrms);

    // final step, compute the errors
    std::vector<double> errors;
    std::cout << "Computing the angle errors...\n";
    ams562_final::compute_errors(points, nrms, errors);

    std::cout << "All computations are done, doing error analysis...\n";
    // note that max_element returns iterators thus requiring dereferencing
    const double inf_err =
        *std::max_element(errors.begin(), errors.end()) * RAD2DEG;
    std::cout << "Inf error = " << inf_err << " degrees\n";
    if (inf_err < thresholds[kase - 1] && inf_err != 0.0)
      std::cout << PASS << " test case " << kase << "\n\n";
    else
      std::cerr << FAIL << " the error is larger than threshold "
                << thresholds[kase - 1] << ", or you haven\'t finished yet" << "\n\n";
  }
}
