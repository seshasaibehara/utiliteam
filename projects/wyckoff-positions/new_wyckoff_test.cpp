#include "./wyckoff.hpp"
#include "../avdv-factor-group/include.hpp"

bool test_convert_symop_to_fractional()
{


    return false;
}

bool test_combine_symops()
{
  //  take a know group and combine comparing M-I matric and b matrix
  //
  return false;
}

bool test_calc_smith_normal_form()
{
    return false;
}

bool test_calc_wyckoff_prototypes()
{ 
    //test what the initial wyckoff prototypes are for given small group
    return false;
}

bool test_calc_wyckoff_prototypes_screw()
{ 
    //test what the initial wyckoff prototypes are for given small group with a screw operation
    //this should not generate a wyckoff position
    return false;
}

bool test_calc_wyckoff_prototypes_translation()
{
    //test with group which includes translations, such as an off center rotation axis
    //
    return false;
}

bool test_check_if_wyckoff_orbit_is_kept()
{   
    //compare size of coset with size of wyckoff orbit
    //if they are the same keep, if not discard
    return false;
}

bool test_generate_all_wyckoff_positions_pointgroup()
{   
    // use a simple point group and lattice to genreate wyckoff positions. 
    // compare to those known first checking size and then checking actual wyckoff positions
    
    return false;
}

bool test_generate_all_wyckoff_positions_from_POSCAR()
{
    //start with simple poscar, read in, genreate factor group, generate wyckoff positions
    return false;
}

bool test_genreate_all_wyckoff_positions_from_hexagonal_POSCAR()
{
    return false;
}

int main()
{
    //list all above tests with EXEPECT_TRUE()//
    return 0;
}


