#include "symmetry_operations.hpp"
// Initial Class set up
// Calculate the grid from which we will calculate potential Lprime values. Need 3 by 3 integers




/*bool MatrixComparison(Eigen::Matrix3d& Matrix1, Eigen::Matrix3d& Matrix2)
{
	return Matrix1.isApprox(Matrix2);
}


bool VectorComparison(Eigen::Vector3d& Vector1, Eigen::Vector3d& Vector2)
{
	return Vector1.isApprox(Vector2);
}
*/


std::vector<Eigen::Vector3d> calculate_gridpoints(Lattice my_lattice,
                                                  int radius) // Might take in a structure instead and get a lattice from it
{
    Eigen::Matrix3d lattice=my_lattice.col_vector_matrix();
    std::vector<Eigen::Vector3d> gridpoints; // left in int since I put radius in int
    Eigen::Vector3d pn;
    
    for (int m = -radius; m < (radius + 1); m++)
    {
        for (int n = -radius; n < (radius + 1); n++)
        {
            for (int j = -radius; j < (radius + 1); j++)
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
std::vector<Eigen::Matrix3d> Calculate_Lprimes(Lattice my_lattice)  
{
    Eigen::Matrix3d lattice=my_lattice.col_vector_matrix();
    int radius=1;
    std::vector<Eigen::Matrix3d> Lprimes;
    auto PS = calculate_gridpoints(lattice, radius);
    Eigen::Matrix3d MakeMatrix;
    for (auto p1 : PS)
    {
        for (auto p2 : PS)
        {
            for (auto p3 : PS)
            {
                MakeMatrix << p1(0), p2(0), p3(0), p1(1), p2(1), p3(1), p1(2), p2(2),
                    p3(2); // I may have these flipped. If Im supposed to have [p1][p2] then this woulr be correct
                Lprimes.push_back(MakeMatrix.inverse());
            }
        }
    }
    return Lprimes;
}

// Determine whether or not the calculated Symmetry operation is valid through comparison of S^T*S=I
bool is_symop_valid(Eigen::Matrix3d SymMatrix)
{
    auto Matrixcheck = SymMatrix.transpose() * SymMatrix;
    if (!Matrixcheck.isIdentity(.0005))
        return false;
    else
        return true;
}
// This function calculates the symmetry operations that are valid for a given lattice
std::vector<SymOp> Calculate_point_group(Lattice my_lattice) // Is the type symops?
{
    Eigen::Matrix3d lattice=my_lattice.col_vector_matrix();
    int radius=1;
    Eigen::Matrix3d dummy_matrix;
    Eigen::Vector3d dummy_tau;
    SymOp test_symop(dummy_matrix, dummy_tau);
    std::vector<SymOp> validsymops;
    auto Lprimes = Calculate_Lprimes(lattice);
    Eigen::Matrix3d SymmetryOp;
    for (auto Lp : Lprimes)
    {
        SymmetryOp = lattice * Lp;
        if (is_symop_valid(SymmetryOp))
        {
            test_symop.get_cart_matrix()=SymmetryOp;
	    test_symop.get_translation()<<0,0,0; //would I need to push this back as well?
	    validsymops.push_back(test_symop);
        }
    }
    return validsymops;
}

//    bool is_SymOp_Unitary(SymOp SymMatrix); //Look up C++ general determinant calcs

//
bool group_is_closed(std::vector<Eigen::Matrix3d> SymMatrix) //
{
    Eigen::Matrix3d GroupMultiplication;
    for (auto S1 : SymMatrix)
    {
        for (auto S2 : SymMatrix)
        {
            GroupMultiplication = S1 * S2;
            MatrixCompare_f compare_matrix(GroupMultiplication);

            if (std::find_if(SymMatrix.begin(), SymMatrix.end(), compare_matrix) != SymMatrix.end())
                return true;
            else
                return false;
        }
    }
}

std::vector<Site> transform_basis(SymOp symop, std::vector<Site> basis)
{
	std::vector<Site> test_basis;
	
	for (auto basis_index:basis)
	{
		test_basis.push_back(symop.get_cart_matrix()*basis_index(basis.get_atom(), basis.get_coordinate()+symop.get_translation());
	}
	return test_basis;
}


bool basis_maps_onto_itself(std::vector<Eigen::Vector3d> original_basis, std::vector<Eigen::Vector3d> transformed_basis)
{
//	std::cout<<"in basis maps onto itself";
	for (int l=0; l<transformed_basis.size(); l++)
	{
		VectorCompare_f compare(transformed_basis[l]);
		if (std::find_if(original_basis.begin(), original_basis.end(), compare)==original_basis.end())
		{
			return false;
		}
	}
	return true;
}


//find factor group using constructor
//std::vector<factor_group> find_factor_group(std::vector<Eigen::Matrix3d> ValidSymOps, Eigen::Matrix3d Lattice)
std::vector<SymOp> find_factor_group(std::vector<SymOp> ValidSymOps, Structure my_struc)
{ 
	std::vector<Eigen::Vector3d> Basis;
	std::vector<Eigen::Matrix3d> ValidCartMatricies;
	Eigen::Matrix3d Lattice= my_struc.get_lattice().col_vector_matrix();
	for (int j=0; j<my_struc.get_sites().size(); j++)
	{
		Basis.push_back(my_struc.get_sites().at(j).get_coordinate());
	}
	for (int i=0; i<ValidSymOps.size(); i++)
	{
		ValidCartMatricies.push_back(ValidSymOps.at(i).get_cart_matrix());
	}
	std::vector<SymOp> tally;
	Eigen::Matrix3d test_cart_matrix;
	Eigen::Vector3d test_tau;
	SymOp symop(test_cart_matrix, test_tau); 
	for (auto cart_matrix: ValidCartMatricies)
	{
		symop.get_cart_matrix()=cart_matrix;
		symop.get_translation()<<0,0,0; 
		auto test_basis=transform_basis(symop, Basis);
		if (basis_maps_onto_itself(Basis, test_basis))	
		{
			tally.push_back(symop);  //necessary?
		}
		else
		{
			std::vector<Eigen::Vector3d> total_basis;
			Eigen::Vector3d trans;
			for (int j=0; j<Basis.size(); j++)
			{
				for (int k=0; k<Basis.size(); k++)
				{
				trans=Basis[j]-test_basis[k];
					for (int m=0; m<test_basis.size(); m++)
					{
						Eigen::Vector3d changed_basis=test_basis[k]+trans;
						total_basis.push_back(changed_basis);
					}
				}	

			}

					if (basis_maps_onto_itself(Basis, total_basis))
					{
						symop.get_translation()=trans;
						tally.push_back(symop);

					
					}
		}
	}
	return tally;
}
