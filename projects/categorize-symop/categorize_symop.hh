#ifndef CAT_SYMOP
#define CAT_SYMOP

#include "../../submodules/eigen-git-mirror/Eigen/Core"
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include "../avdv-factor-group/include.hpp"
#include "./symop_4d.hpp"

bool compare_diff_to_prec(Eigen::Vector3d difference, double tol);
bool has_translation(const Eigen::Vector3d translation, const Eigen::Matrix3d lattice, double tol);
Eigen::Vector3d project_translation_onto_vectors(const std::vector<Eigen::Vector3d>& eigen_vectors, Eigen::Vector3d translation);
bool lmost_equal(double LHS, double RHS, double tol);
std::vector<Eigen::Vector3d> eigenvectors_with_positive_unit_eigenvalues(const Eigen::Matrix3d& cart_matrix, double tol);
std::string check_op_type(const SymOp sym_op, const Lattice lattice, double tol);
Subspace find_invariant_subspace(Symop_4d symop4d);
Subspace find_invariant_subspace(SymOp symop);


#endif