#ifndef WYCOFF_HH
#define WYCOFF_HH

#include "../../submodules/eigen-git-mirror/Eigen/Core"
#include "../avdv-factor-group/include.hpp"
#include <vector>

class Subspace
{
    public:
        Subspace(const Eigen::Matrix3d& basis_col_vectors, const Eigen::Vector3d& offset);
        std::string formula() const;
        Eigen::Vector3d offset;
        Eigen::Matrix3d basis_col_vectors;/* [100    [x   
                                      010 dot y  + [offset] 
                                      001]    z]*/

    //Maybe a 4x4 matrix?
};

/* std::vector<Subspace> make_orbit(const std::vector<SymOp>& coset, const Subspace& prototype); */
/* Subspace operator*(const SymOp& op, const Subspace& prototype); */
/*everytime we see x in the wycoff position, we can replace it with a unit vector

/yyyyyyy?? Wyckoff?? */
/*     vector of vectors */
/*
 * Some notes:
 *  Wyckoff positions are vectors with some variable components.
 *  they are points lines or planes? how do we define each?
 *  ex: axis(x,0,0), pt(1/4, 14/, 1/4) plane(-z,1/2y, 0)
 *
 *  wyckoff positions are vectors of functions of new class called wyckoff function
 *      ex:
 *      wyxcoff_func(x,y,z){
 *          double a=0;
 *          double b=x-y;
 *          double c=0;
 *          }
 *       wyckoff2(x,y,z){
 *        double a=1;
 *        double b=0;
 *        double c=.5;
 *        }
 *find_all_equivalent_site({x,y,z}, wyckoff_1i, coset)
    {
        is_valid_wyckoff(wyckoff_1, {x,y,z})
    for( symop : coset){
            coset*{x,y,z}
            }
       i}
i*/

//TODO: Consider using a typedef for SymGroup
typedef SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f> PeriodicGroup;

std::vector<SymOp> find_coset(PeriodicGroup factor_group,PeriodicGroup subgroup);

Eigen::Matrix4d make_symop_4dmatrix(SymOp input_op);

//Reynolds?
Eigen::Matrix4d make_reynolds_operator(PeriodicGroup subgroup);

Subspace find_invariant_subspace(SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f> subgroup);

Subspace operator*(const SymOp& lhs, const Subspace& rhs);

std::vector<Subspace> find_symmetrically_equivalent_wyckoff_positions(std::vector<SymOp> coset, Subspace wyckoff_position);









#endif
