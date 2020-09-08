#include "categorize_symop.hh" 
#include "./symop_4d.hpp"
#include <stdexcept>

//#define tol 1e-6

bool compare_diff_to_prec(Eigen::Vector3d difference, double tol)
{
    int sig_diff = 0;
    for (int i = 0; i < 3; i++)
    {
        if (abs(difference(i)) > tol)
        {
            sig_diff++;
        }
    }
    if (sig_diff != 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool has_translation(const Eigen::Vector3d translation, const Lattice lattice_input, double tol)
{ // check if translation is 0 or integer multiple of lattice vectors-> false,
    // else truei
    Eigen::Matrix3d lattice = lattice_input.col_vector_matrix();
    Eigen::Vector3d lattice_displacement = lattice.inverse() * translation;
    Eigen::Vector3d modification(0.5, 0.5, 0.5);
    Eigen::Vector3d modif_displ = lattice_displacement + modification;
    Eigen::Vector3d trunc_displ;

    for (int i = 0; i < 3; i++)
    {
        trunc_displ(i) = (int)modif_displ(i);
    }
    Eigen::Vector3d difference = lattice_displacement - trunc_displ;
    bool check_trans = compare_diff_to_prec(difference, tol);

    if (check_trans == true)
    {
        return true;
    }
    else
    {
        return false;
    }
}

Eigen::Vector3d project_translation_onto_vectors(const std::vector<Eigen::Vector3d>& eigen_vectors, Eigen::Vector3d translation)
{
    // for each eigen vector this call the above project onto vector function and then added these vectors together
    Eigen::Vector3d projected_translation_vector = Eigen::Vector3d::Zero();
    for (const auto& vector : eigen_vectors)
    {
        double projected_component = vector.dot(translation);
        projected_translation_vector += (projected_component * vector);
    }
    return projected_translation_vector;
}

bool almost_equal(double LHS, double RHS, double tol) { return std::abs(LHS - RHS) < tol; }

/**
 * Returns a container of eigenvectors of the given matrix that have
 * positive unit eigenvalues
 */

std::vector<Eigen::Vector3d> eigenvectors_with_positive_unit_eigenvalues(const Eigen::Matrix3d& cart_matrix, double tol)
{
    // some function that evaluate eigne vlaue and counts how many ones
    // Now returns a matrix of the Eigen vectors with eigen values of 1, in column vectors
    // need to edit so that it returns only the real components of the eigenvector
    Eigen::EigenSolver<Eigen::Matrix3d> solver(cart_matrix, true);
    Eigen::Vector3<std::complex<double>> eigenvals = solver.eigenvalues();
    std::vector<Eigen::Vector3d> output_eigen_vectors;
    Eigen::Matrix3<std::complex<double>> eigenvectors = solver.eigenvectors();

    int ct=0;
    for (int i = 0; i < 3; i++)
    {
//      std::cout <<eigenvectors.col(i).real()<<std::endl;
        auto eigenval = eigenvals(i);
//        std::cout<<eigenval.real()<<std::endl;
        if (almost_equal(eigenval.real(), 1, tol))
        {
            ct++;
//            std::cout<<"EigenVal=1 Found!"<<std::endl;
            Eigen::Vector3d real_eigenvector = eigenvectors.col(i).real();
            output_eigen_vectors.push_back(real_eigenvector);
        }
    }
    if (ct==0){ 
        for (int i = 0; i < 3; i++)
        {
//        std::cout <<eigenvectors.col(i).real()<<std::endl;
          auto eigenval = eigenvals(i);
//          std::cout<<eigenval.real()<<std::endl;
          if (almost_equal(eigenval.real(), -1, tol))
          {
              ct++;
//              std::cout<<"EigenVal=1 Found!"<<std::endl;
              Eigen::Vector3d real_eigenvector = eigenvectors.col(i).real();
              output_eigen_vectors.push_back(real_eigenvector);
           }
         }
    }

    return output_eigen_vectors;
}

std::string check_op_type(const SymOp sym_op, const Lattice lattice, double tol)
{ // take in sym_op returns string of op type
    int trace = sym_op.get_cart_matrix().trace();
    std::string type;
    double det = sym_op.get_cart_matrix().determinant();

    if ((3 - trace < tol))
    {
        type = "Identity";
        return type;
    }
    if (trace + 3 < tol)
    {
        type = "Inversion";
        return type;
    }
    std::vector<Eigen::Vector3d> eigen_vectors = eigenvectors_with_positive_unit_eigenvalues(sym_op.get_cart_matrix(), tol);
//    std::cout<<eigen_vectors<<std::endl;
//    std::cout << eigen_vectors.cols() << std::endl;
    if (eigen_vectors.size() == 2)
    {
        if (has_translation(project_translation_onto_vectors(eigen_vectors, sym_op.get_translation()), lattice, tol))
        {
            type = "Glide";
            return type;
        }
        else
        {
            type = "Mirror";
            return type;
        }
    }
    else if (eigen_vectors.size() == 1)
    {
        if (abs(det - 1) < tol)
        {
            if (has_translation(project_translation_onto_vectors(eigen_vectors, sym_op.get_translation()), lattice, tol))
            {
                type = "Screw";
                return type;
            }
            else
            {
                type = "Rotation";
                return type;
            }
        }
        else if (abs(det + 1) < tol)
        {
            type = "Improper rotation";
            return type;
        }
        else
        {
            type = "Error type 1: Type not idenitified!!!";
            return type;
        }
    }
    else
    {
        type = "Error type 2: Type not idenitified!!!";
        return type;
    }
}

Subspace find_invariant_subspace(Symop_4d symop, Lattice lattice, double tol)
{
    SymOp symop_3d=symop_4d_to_3d(symop);
    if(has_translation(symop_3d.get_translation(), lattice, tol)){// throw error because it has not invariant subspace;
        throw std::runtime_error("This symop does not have an invariant subspace.");
    }
        
    return symop.invariant_subspace;
}

Subspace find_invariant_subspace(SymOp symop, Lattice lattice, double tol)
{
    if(has_translation(symop.get_translation(), lattice, tol)){// throw error because it has not invariant subspace;
        throw std::runtime_error("This symop does not have an invariant subspace.");
       
    }
    Symop_4d symop_4d(symop);    

    return symop_4d.invariant_subspace;
}
