#include "./wyckoff.hpp"

Subspace::Subspace(const Eigen::Matrix3d& input_basis_vectors, const Eigen::Vector3d& input_offset):basis_col_vectors(input_basis_vectors), offset(input_offset){};

std::vector<SymOp> find_coset(PeriodicGroup factor_group, PeriodicGroup subgroup){
    /*finds the coset to a subgroup in a factor group and returns it as 
     * a std::vector of SymOps*/
    std::vector<SymOp> coset;
    auto binary_compare=factor_group.get_comparator();
    for( SymOp factor_group_symop:factor_group.operations()){
        auto compare=[factor_group_symop, binary_compare](SymOp subgroup_symop){
            return binary_compare(factor_group_symop, subgroup_symop);
        };
        if(std::find_if(subgroup.operations().begin(), subgroup.operations().end(), compare)==subgroup.operations().end()){
            coset.push_back(factor_group_symop);
        }
    }
    return coset;
}

Eigen::Matrix4d make_symop_4dmatrix(SymOp input_symop){
    Eigen::Matrix4d symop_4d;
    symop_4d.row(3)<<0,0,0,1;
    for(int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            symop_4d(i,j)=input_symop.m_cart_matrix(i,j);
        }
    }
    for(int i=0; i<3;i++){
        symop_4d(i,3)=input_symop.get_translation()(i);
    }
    return symop_4d;
}

Eigen::Matrix4d make_reynolds_operator(PeriodicGroup subgroup){
    Eigen::Matrix4d average_op=Eigen::Matrix4d::Zero();
    for(SymOp sym_op:subgroup.operations()){
        average_op+=make_symop_4dmatrix(sym_op);
    }
    average_op=average_op/subgroup.operations().size();
    return average_op;
}

/*Subspace find_invariant_subspace(PeriodicGroup subgroup){
 return;
}*/

Subspace operator*(const SymOp& lhs, const Subspace& rhs)
{
    Eigen::Matrix3d basis_vectors= lhs.get_cart_matrix() * rhs.basis_col_vectors;
    Eigen::Vector3d offset_vector= lhs.get_cart_matrix()* rhs.offset+lhs.get_translation();
    Subspace product_wycoff_position(basis_vectors, offset_vector);
    return product_wycoff_position;
}

std::vector<Subspace> find_symmetrically_equivalent_wyckoff_positions(std::vector<SymOp> coset, Subspace wyckoff_position)
{
    std::vector<Subspace> equivalent_wyckoff_positions;
    equivalent_wyckoff_positions.push_back(wyckoff_position);
    for( SymOp symop : coset){
        equivalent_wyckoff_positions.push_back(symop*wyckoff_position);
    }
    return equivalent_wyckoff_positions;
}
