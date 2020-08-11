#include "./wyckoff.hpp"
#include "../eigen-git-mirror/Eigen/Eigenvalues"
#include <cmath>

Subspace::Subspace(const Eigen::Vector3d& input_vec_0, const Eigen::Vector3d& input_vec_1,
                   const Eigen::Vector3d& input_vec_2, const Eigen::Vector3d& offset)
    : m_offset(offset)
{
//    double tol=1e-6;
    Eigen::Vector3d basis_vec_0=input_vec_0.normalized();
    Eigen::Vector3d basis_vec_1=input_vec_1.normalized();
    Eigen::Vector3d basis_vec_2=input_vec_2.normalized();
//    int v0=0, v1=0, v2=0;
    for (int i = 0; i < 3; ++i)
    {
        m_basis_col_matrix(i, 0) = basis_vec_0(i);
//        if (abs(basis_vec_0(i))<tol && v0==0){
//            v0=1;}
        m_basis_col_matrix(i, 1) = basis_vec_1(i);
        m_basis_col_matrix(i, 2) = basis_vec_2(i);
    }

    //DONE: basis vectors are normalized
    //minimize the offset (orthormal to the basis)
    //This may still be a different subspace though order if given axis, or plane could be define in different orders
    //TODO: make sure planes are defined by orthoganol vectors?
}

Subspace::Subspace(const Eigen::Matrix3d& input_basis_vectors, const Eigen::Vector3d& input_offset)
    : m_basis_col_matrix(input_basis_vectors), m_offset(input_offset)
{
}

Subspace::Subspace():m_basis_col_matrix(Eigen::Matrix3d::Zero()), m_offset(Eigen::Vector3d::Zero()){}

Eigen::Matrix3d Subspace::basis_col_matrix() const { return this->m_basis_col_matrix; }

Eigen::Vector3d Subspace::offset() const { return this->m_offset; }

Subspace Subspace::operator*(const SymOp& lhs)
{
    Eigen::Matrix3d basis_vectors = lhs.get_cart_matrix() * this->m_basis_col_matrix;
    Eigen::Vector3d offset_vector = lhs.get_cart_matrix() * this->m_offset + lhs.get_translation();
    Subspace product_wycoff_position(basis_vectors, offset_vector);
    return product_wycoff_position;
}

std::string Subspace::formula() const
{
    double tol = 1e-5;
    std::string formula = "{";
    std::vector<std::string> xyz = {"x", "y", "z"};

    for (int i = 0; i < 3; ++i)
    {
        std::string temp_string;
        for (int j = 0; j < 3; ++j)
        {
            if (std::abs(this->m_basis_col_matrix(i,j)) > tol)
            {
                std::string sign = this->m_basis_col_matrix(i,j) > 0 ? "+" : "";
                temp_string += sign + std::to_string(this->m_basis_col_matrix(i,j)) + xyz[j];
            }
        }
        if (temp_string.size() == 0)
        {
            temp_string = "0";
        }
        formula += temp_string + ", ";
    }

    formula.pop_back();
    formula.pop_back();
    formula += "}";
    return formula;
}

std::vector<SymOp> find_coset(PeriodicGroup factor_group, PeriodicGroup subgroup)
{
    /*finds the coset to a subgroup in a factor group and returns it as
     * a std::vector of SymOps*/
    std::vector<SymOp> coset;
    auto binary_compare = factor_group.get_comparator();
    for (SymOp factor_group_symop : factor_group.operations())
    {
        auto compare = [factor_group_symop, binary_compare](SymOp subgroup_symop) {
            return binary_compare(factor_group_symop, subgroup_symop);
        };
        if (std::find_if(subgroup.operations().begin(), subgroup.operations().end(), compare) ==
            subgroup.operations().end())
        {
            coset.push_back(factor_group_symop);
        }
    }
    return coset;
}

Eigen::Matrix4d make_symop_4dmatrix(SymOp input_symop)
{
    Eigen::Matrix4d symop_4d;
    symop_4d.row(3) << 0, 0, 0, 1;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            symop_4d(i, j) = input_symop.get_cart_matrix()(i, j);
        }
    }
    for (int i = 0; i < 3; i++)
    {
        symop_4d(i, 3) = input_symop.get_translation()(i);
    }
    return symop_4d;
}

Eigen::Matrix4d make_4d_reynolds_operator(PeriodicGroup subgroup)
{
    Eigen::Matrix4d average_op = Eigen::Matrix4d::Zero();
    for (SymOp sym_op : subgroup.operations())
    {
        average_op += make_symop_4dmatrix(sym_op);
    }
    average_op = average_op / subgroup.operations().size();
    return average_op;
}

/*Subspace find_invariant_subspace(PeriodicGroup subgroup){
 return;
}*/

// Subspace operator*(const SymOp& lhs, const Subspace& rhs)
//{
//    Eigen::Matrix3d basis_vectors = lhs.get_cart_matrix() * rhs.basis_col_matrix();
//    Eigen::Vector3d offset_vector = lhs.get_cart_matrix() * rhs.offset() + lhs.get_translation();
//    Subspace product_wycoff_position(basis_vectors, offset_vector);
//    return product_wycoff_position;
//}

Eigen::Matrix3d make_3d_reynolds_operator(PeriodicGroup subgroup)
{
    Eigen::Matrix3d average_op = Eigen::Matrix3d::Zero();
    for (const auto& sym_op : subgroup.operations())
    {
        average_op += sym_op.get_cart_matrix();
    }

    average_op = average_op / subgroup.operations().size();
    return average_op;
}

Subspace find_invariant_subspace(PeriodicGroup subgroup)
{
    // TODO: How is tolerance is being dealt? Any global variable to use or should it be an input arg?
    double tol = 1e-5;
    Eigen::Matrix3d reynolds_operator = make_3d_reynolds_operator(subgroup);
    Eigen::EigenSolver<Eigen::Matrix3d> eigen_solver(reynolds_operator);
    Eigen::VectorXcd eigen_values = eigen_solver.eigenvalues();
    Eigen::MatrixXcd eigen_vectors = eigen_solver.eigenvectors();
    Eigen::Matrix3d basis_col_matrix = Eigen::Matrix3d::Zero();
    // TODO: Should compute the average translation in subgroup to get the offset
    Eigen::Vector3d offset(0, 0, 0);

    for (int i = 0; i < 3; ++i)
    {
        if (std::abs(eigen_values(i).real() - 1) < tol && std::abs(eigen_values(i).imag()) < tol)
        {
            for (int j = 0; j < 3; ++j)
            {
                basis_col_matrix(j, i) = eigen_vectors(j, i).real();
            }
        }
    }

    Subspace invariant_subspace(basis_col_matrix, offset);
    return invariant_subspace;
}

std::vector<Subspace> find_symmetrically_equivalent_wyckoff_positions(std::vector<SymOp> coset,
                                                                      Subspace wyckoff_position)
{
    std::vector<Subspace> equivalent_wyckoff_positions;
    equivalent_wyckoff_positions.push_back(wyckoff_position);
    for (SymOp symop : coset)
    {
        equivalent_wyckoff_positions.push_back(wyckoff_position * symop);
        //TODO: check if in list before push_back
    }

    return equivalent_wyckoff_positions;
}

bool subspaces_are_equal(Subspace lhs, Subspace rhs, double tol)
{
    return (lhs.basis_col_matrix().isApprox(rhs.basis_col_matrix(), tol) && lhs.offset().isApprox(rhs.offset(), tol));

    //TODO: insufficient comparison, multiple vectors can span the same subspace
    // is A is contained in B and B is contained in A then they are the same subspace (projections?)
}

bool wyckoff_positions_are_equal(std::vector<Subspace> lhs, std::vector<Subspace> rhs, double tol)
{ //TODO: 
    return false;
}
