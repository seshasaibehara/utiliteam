#ifndef SYMOP_4D
#define SYMOP_4D

#include "../../submodules/eigen-git-mirror/Eigen/Core"
#include "../avdv-factor-group/include.hpp"
#include "wyckoff.hpp"
#include <vector>

class Symop_4d
{
public:
    Symop_4d(const Eigen::Matrix3d& cart_matrix, const Eigen::Vector3d& translation);
    Symop_4d(const Eigen::Matrix4d& symop_matrix);
    Symop_4d(const SymOp& symop_3d);
    const Eigen::Matrix4d get_matrix() { return this->symop_matrix; };
    Subspace find_invariant_subspace();
    int subspace_dimension();
    Subspace invariant_subspace;

private:
    Eigen::Matrix4d symop_matrix;
    int dimension = -1;
};

#endif
