#include "wyckoff.hpp"
// is there a way to include all the headers in the avdv-factorgroup folder??
#include "../../submodules/eigen-git-mirror/Eigen/Core"
#include "../../submodules/eigen-git-mirror/Eigen/Dense"
#include "../avdv-factor-group/tests.hpp"

bool test_wyckoff_construction() { return true; }

PeriodicGroup make_identity_group(double tol)
{
    Eigen::Matrix3d identity_matrix=Eigen::Matrix3d::Identity();
    Lattice identity_lattice(identity_matrix.col(0), identity_matrix.col(1), identity_matrix.col(2));
    BinarySymOpPeriodicCompare_f symop_comparator(identity_lattice, tol);
    BinarySymOpPeriodicMultiplier_f symop_multiplier(identity_lattice, tol);
    PeriodicGroup identity_group({identity_matrix}, symop_comparator, symop_multiplier);
    return identity_group;
}

PeriodicGroup make_rotation_60_group(double tol)
{
    Eigen::Matrix3d identity_matrix=Eigen::Matrix3d::Identity();
    Lattice identity_lattice(identity_matrix.col(0), identity_matrix.col(1), identity_matrix.col(2));
    BinarySymOpPeriodicCompare_f symop_comparator(identity_lattice, tol);
    BinarySymOpPeriodicMultiplier_f symop_multiplier(identity_lattice, tol);
    PeriodicGroup rotation_60_group({make_z_rotation_matrix(60)}, symop_comparator, symop_multiplier);
    return rotation_60_group;
}

bool test_find_coset(double tol)
{
    /*    make a rotation group
     *    make identity group
     *    find coset of identity group in rotations group
     *    should be length rotation group -1, and not contain identity, and contain all other elements of rotation group.*/

    PeriodicGroup rotation_60_group = make_rotation_60_group(tol);
    PeriodicGroup identity_group = make_identity_group(tol);
    std::vector<SymOp> coset = find_coset(rotation_60_group, identity_group);

    // check if subgroup is not in coset:

    BinarySymOpPeriodicCompare_f symop_comparator = rotation_60_group.get_comparator();
    auto identity_compare = [identity_group, symop_comparator](SymOp input_symop) {
        for (SymOp group_op : identity_group.operations())
        {
            if (symop_comparator(group_op, input_symop))
            {
                return true;
            }
        }
        return false;
    };
    if (std::find_if(coset.begin(), coset.end(), identity_compare) != coset.end())
    {
        return false;
    }

    // Check if every coset element is in larger group(factor group)
    for (SymOp coset_symop : coset)
    {
        auto coset_compare = [coset_symop, symop_comparator](SymOp input_symop) { return symop_comparator(coset_symop, input_symop); };
        if (std::find_if(rotation_60_group.operations().begin(), rotation_60_group.operations().end(), coset_compare) ==
            rotation_60_group.operations().end())
        {
            return false;
        }
    }

    return true;
}

bool test_make_SymOp4dMatrix()
{
    /* make symop
     * make symop 4d, with translation
     * check all elements
     */
    Eigen::Matrix3d test_matrix;
    test_matrix<<1, 0, 0, 0, -1, 0, 0, 0, 1;
    Eigen::Vector3d test_translation;
    test_translation<<.25, .25, .25;
    SymOp test_input_symop(test_matrix, test_translation);
    Eigen::Matrix4d symop_4d = make_symop_4dmatrix(test_input_symop);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (test_matrix(i, j) != symop_4d(i, j))
            {
                return false;
            }
        }
    }
    for (int i = 0; i < 3; i++)
    {
        if (test_translation(i) != symop_4d(i,3))
        {
            return false;
        }
    }
    for (int i = 0; i < 3; i++)
    {
        if (i < 3 && symop_4d(3,i) != 0)
        {
            return false;
        }
    }
    if (symop_4d(3,3) != 1){ return false;}

    return true;
}

bool test_make_reynolds_operator(double tol)
{
    /*make roation group
     * find average
     * average should be
     */
    Eigen::Matrix4d expected_average;
    expected_average << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
    /*if the translation was not 000, would this move the axis from the 00z?
     */
    PeriodicGroup rotation_60_group = make_rotation_60_group(tol);
    Eigen::Matrix4d average_matrix = make_reynolds_operator(rotation_60_group);
    return average_matrix.isApprox(expected_average, tol);
}

bool test_find_invariant_subspace() { return true; }

bool test_find_symmetrically_equivalent_wyckoff_positions()
{
    /*make wyckoff position
     * make factor group?
     * find subgroup
     * find coset
     * find_symetrically_equvalent_poitns
     * check list vs expected
     */
    return true;
}

int main()
{
    double tol = 1e-6;
    /* all the tests above*/
   // EXPECT_TRUE(//test_wyckoff_construction, "Construction of Wyckoff class");
    EXPECT_TRUE(test_make_SymOp4dMatrix(), "Convert SymOp to 4d Matrix");
    EXPECT_TRUE(test_find_coset(tol), "Generate Expected Coset");
    EXPECT_TRUE(test_make_reynolds_operator(tol), "Generate Reynolds operator");
    //EXPECT_TRUE(test_find_invariant_subspace, "Find the expected invariant subspace of subgroup");

    return 0;
}

