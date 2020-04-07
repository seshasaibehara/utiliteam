#include "symmetry_operations.hpp"
#include <iostream>

std::vector<Eigen::Vector3d> calculate_gridpoints(const Lattice& my_lattice,
                                                  int radius) // Might take in a structure instead and get a lattice from it
{
    Eigen::Matrix3d lattice = my_lattice.row_vector_matrix();
    std::vector<Eigen::Vector3d> gridpoints; // left in int since I put radius in int
    Eigen::Vector3d pn;

    for (int m = -radius; m <= radius; m++)
    {
        for (int n = -radius; n <= radius; n++)
        {
            for (int j = -radius; j <= radius; j++)
            {
                pn = n * lattice.col(0) + m * lattice.col(1) + j * lattice.col(2);
                gridpoints.push_back(pn);
            }
        }
    }

    return gridpoints;
}
// change complete

// calculate a list of potentials Lprimes
std::vector<Eigen::Matrix3d> Calculate_Lprimes(const Lattice& my_lattice)
{
    Eigen::Matrix3d lattice = my_lattice.row_vector_matrix();
    int radius = 1;
    std::vector<Eigen::Matrix3d> Lprimes;
    auto PS = calculate_gridpoints(lattice, radius);
    Eigen::Matrix3d MakeMatrix;
    for (const auto& p1 : PS)
    {
        for (const auto& p2 : PS)
        {
            for (const auto& p3 : PS)
            {
                MakeMatrix << p1(0), p2(0), p3(0), p1(1), p2(1), p3(1), p1(2), p2(2),
                    p3(2); // I may have these flipped. If Im supposed to have [p1][p2]
                           // then this woulr be correct
                Lprimes.push_back(MakeMatrix.inverse());
            }
        }
    }
    return Lprimes;
}

// Determine whether or not the calculated Symmetry operation is valid through
// comparison of S^T*S=I
bool is_symop_valid(const Eigen::Matrix3d& SymMatrix)
{
    auto Matrixcheck = SymMatrix.transpose() * SymMatrix;
    if (!Matrixcheck.isIdentity(.005))
    {
        return false;
    }
    return true;
}
// This function calculates the symmetry operations that are valid for a given
// lattice
std::vector<SymOp> Calculate_point_group(const Lattice& my_lattice) // Is the type symops?
{
    Eigen::Matrix3d lattice = my_lattice.row_vector_matrix();
    int radius = 1;
    std::vector<SymOp> validsymops;
    Eigen::Vector3d trans;
    trans << 0, 0, 0;
    auto Lprimes = Calculate_Lprimes(lattice);
    Eigen::Matrix3d SymmetryOp;
    for (const auto& Lp : Lprimes)
    {
        SymmetryOp = lattice * Lp;
        if (is_symop_valid(SymmetryOp))
        {
            SymOp test_symop(SymmetryOp, trans);
            validsymops.push_back(test_symop);
        }
    }
    return validsymops;
}

//    bool is_SymOp_Unitary(SymOp SymMatrix); //Look up C++ general determinant
//    calcs

//
bool group_is_closed(const std::vector<Eigen::Matrix3d>& SymMatrix) //
{
    Eigen::Matrix3d GroupMultiplication;
    for (const auto& S1 : SymMatrix)
    {
        for (const auto& S2 : SymMatrix)
        {
            GroupMultiplication = S1 * S2;
            MatrixCompare_f compare_matrix(GroupMultiplication);

            if (std::find_if(SymMatrix.begin(), SymMatrix.end(), compare_matrix) != SymMatrix.end())
                return true;
        }
    }
    return false;
}

std::vector<Site> transform_basis(const SymOp& symop, const std::vector<Site>& basis)
{
    std::vector<Site> test_basis;
    for (const auto& basis_index : basis)
    {
        Eigen::Vector3d transformed_coordinate = symop.get_cart_matrix() * basis_index.get_coordinate() + symop.get_translation();
        test_basis.push_back(Site(basis_index.get_atom(), transformed_coordinate));
    }
    return test_basis;
}

bool basis_maps_onto_itself(const std::vector<Site>& original_basis, const std::vector<Site>& transformed_basis)
{
    for (const auto& basis : transformed_basis)
    {
        SiteCompare_f test_basis(basis, 1E-3);
        if (std::find_if(original_basis.begin(), original_basis.end(), test_basis) == original_basis.end())
        {
            return false;
        }
    }
    return true;
    // std::vector<Eigen::Vector3d> vector_of_coordinates;
    // for (const auto& site : original_basis)
    // {
    //     vector_of_coordinates.push_back(site.get_coordinate());
    // }

    // for (int l = 0; l < transformed_basis.size(); l++)
    // {
    //    VectorCompare_f compare(transformed_basis[l].get_coordinate());
    //   if (std::find_if(vector_of_coordinates.begin(), vector_of_coordinates.end(), compare) ==
    //       vector_of_coordinates.end())
    //   {
    //       return false;
    //   }
    // }
    // return true;
}

// find factor group using constructor
// std::vector<factor_group> find_factor_group(std::vector<Eigen::Matrix3d>
// ValidSymOps, Eigen::Matrix3d Lattice)
std::vector<SymOp> find_factor_group(Structure my_struc)
{

    std::vector<Site> Basis;
    std::vector<Eigen::Matrix3d> ValidCartMatricies;
    std::vector<SymOp> ValidSymOps = Calculate_point_group(my_struc.get_lattice());

    Eigen::Matrix3d Lattice = my_struc.get_lattice().col_vector_matrix();

    for (int j = 0; j < my_struc.get_sites().size(); j++)
    {
        Basis.push_back(my_struc.get_sites()[j]);
    }

    for (int i = 0; i < ValidSymOps.size(); i++)
    {
        ValidCartMatricies.push_back(ValidSymOps.at(i).get_cart_matrix());
    }

    std::vector<SymOp> tally;
    Eigen::Matrix3d test_cart_matrix;
    Eigen::Vector3d test_tau;
    test_tau << 0, 0, 0;
    SymOp symop(test_cart_matrix, test_tau);
    for (const auto& cart_matrix : ValidCartMatricies)
    {

        SymOp symop(cart_matrix, test_tau);
        //	symop.m_cart_matrix=cart_matrix;
        //      symop.m_translation << 0, 0, 0;
        // Eigen::Vector3d new_test_tau;
        // new_test_tau <<0,0,0;
        // symop(cart_matrix, new_test_tau);
        auto test_basis = transform_basis(symop, Basis);
        if (basis_maps_onto_itself(Basis, test_basis))
        {
            tally.push_back(symop); // necessary?
        }
        else
        {
            std::vector<Site> total_basis;
            Eigen::Vector3d trans;
            for (int j = 0; j < Basis.size(); j++)
            {
                for (int k = 0; k < Basis.size(); k++)
                {
                    trans = Basis[j].get_coordinate() - test_basis[k].get_coordinate();
                    for (int m = 0; m < test_basis.size(); m++)
                    {
                        Eigen::Vector3d changed_basis = test_basis[m].get_coordinate() + trans;
                        total_basis.push_back(Site(test_basis[m].get_atom(), changed_basis));
                    }
                }
            }

            if (basis_maps_onto_itself(Basis, total_basis))
            {
                // symop.get_translation() = trans;

                // tally.push_back(symop);
                tally.push_back(SymOp(cart_matrix, trans));
            }
        }
    }

    for (const auto& x : tally)
    {
        std::cout << x.get_cart_matrix() << std::endl;
        std::cout << x.get_translation() << std::endl;
        std::cout << "-----" << std::endl;
    }

    std::cout << tally.size() << std::endl;
    return tally;
}
