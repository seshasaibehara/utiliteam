#include "symopd_4d.hpp"
#include "wyckoff.hpp"

Symop_4d::Symop_4d(const Eigen::Matrix4d& input_matrix):symop_matrix(input_matrix)
{ 

}

Symop_4d::Symop_4d(const Eigen::Matrix3d& cart_matrix, const Eigen::Vector3d& translation){
    Eigen::Matrix4d symop_matrix=Eigen::Matrix4d::Zero();
    for(int i=0; i<3; i++){
        for (int j =0; j<3; j++){
            symop_matrix(i,j)=cart_matrix(i,j);
        }
        symop_matrix(i,3)=translation(i);
    }
    symop_matrix(3,3)=1;
    Eigen::Matrix3d empty_subspace=Eigen::Matrix3d::Zero();
    Eigen::Vector3d empty_offset=Eigen::Vector3d::Zero();
}

Symop_4d::Symop_4d(const SymOp& symop_3d){
    Symop_4d(symop_3d.get_cart_matrix(), symop_3d.get_translation());
}

Subspace Symop_4d::find_invariant_subspace()
{   double tol=1e-6;
    if (dimension==-1)
    {
        Eigen::Matrix3d basis_matrix=Eigen::Matrix3d::Zero();
        Eigen::Vector3d offset=Eigen::Vector3d::Zero();
        Eigen::EigenSolver<Eigen::Matrix4d> solver(symop_matrix, true);
        Eigen::Vector4<std::complex<double>>eigenvals=solver.eigenvalues();
        std::vector<Eigen::Vector4d> output_eigen_vectors;
        Eigen::Matrix4<std::complex<double>> eigenvectors = solver.eigenvectors();
        dimension++;
        for (int i=0; i<4;i++)
        {
            if ((1-abs(eigenvals(i).real()))<tol)
            { if (1-abs(eigenvectors.col(i)(3))<tol)
                {   for (int k=0; k<3; k++)
                    {offset(k)=eigenvectors.col(i)(k).real();}
                }
              else{ for (int k=0; k<3; k++)
                    {basis_matrix(k,dimension)=eigenvectors(k,i).real();}
                    dimension++;}
            }
        }
        invariant_subspace= Subspace(basis_matrix, offset);
    
    }
    return this->invariant_subspace;
}

int Symop_4d::subspace_dimension()
{
    if (dimension==1)
    { this->find_invariant_subspace();}
    return this->dimension;
}

