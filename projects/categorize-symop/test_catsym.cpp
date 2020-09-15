#include "./categorize_symop.hh"
#include "../avdv-factor-group/tests.hpp"

bool test_categorize_identity(double tol)
{
    Lattice lattice({1,0,0}, {0,1,0}, {0,0,1});
    Eigen::Matrix3d identity_matrix=Eigen::Matrix3d::Identity();
    SymOp identity_op(identity_matrix);
    if (check_op_type(identity_op, lattice, tol)==SYMOP_TYPE::IDENTITY)
    {return true;}
    return false;
}

int main()
{ 
    double tol=1e-6;
    EXPECT_TRUE(test_categorize_identity(tol), "test categorize identity");

    return 0;

}
