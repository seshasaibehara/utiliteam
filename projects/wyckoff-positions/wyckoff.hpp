#ifndef WYCKOFF_HH
#define WYCKOFF_HH

#include "../../submodules/eigen-git-mirror/Eigen/Core"
#include "../avdv-factor-group/include.hpp"
#include <vector>

class Subspace
{
    public:
        /* inherent from symop??*/
        Subspace(const Eigen::Matrix3d& basis_col_vectors, const Eigen::Vector3d& offset);
        std::string formula() const;//Sesha suggested we overload the << operator
        Eigen::Vector3d offset;
        Eigen::Matrix3d basis_col_vectors;/* [100    [x   
                                              010 dot y  + [offset] 
                                              001]    z]*/

};

typedef SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f> PeriodicGroup;

std::vector<SymOp> find_coset(PeriodicGroup factor_group,PeriodicGroup subgroup);

Eigen::Matrix4d make_symop_4dmatrix(SymOp input_op);

//Reynolds?
Eigen::Matrix4d make_reynolds_operator(PeriodicGroup subgroup);
/*template this-->(need begin and end iterator)? vector of symops?*/
Subspace find_invariant_subspace(SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f> subgroup);

Subspace operator*(const SymOp& lhs, const Subspace& rhs);

std::vector<Subspace> find_symmetrically_equivalent_wyckoff_positions(std::vector<SymOp> coset, Subspace wyckoff_position);


//coordinate.cpp has symop*vector, on the utiliteam/find_interstitial/ on  github






#endif
