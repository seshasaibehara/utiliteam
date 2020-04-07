#include "xtal_classes.hpp"
#include "symmetry_operations.hpp"
#include "interstitials.hpp"
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

Structure read_poscar(const std::string &poscar_path) {
  std::ifstream file;
  file.open(poscar_path, std::ios::in);
  std::string line;

  // Comment line
  std::getline(file, line);

  // Scaling line
  std::getline(file, line);
  std::istringstream ss0(line);
  double scaling;
  ss0 >> scaling;

  // Eigen Matrix line
  Eigen::Matrix3d lat_row_matrix;
  for (int i = 0; i < 3; ++i) {
    std::getline(file, line);
    std::istringstream ss1(line);
    ss1 >> lat_row_matrix(i, 0) >> lat_row_matrix(i, 1) >> lat_row_matrix(i, 2);
  }

  // species line
  std::getline(file, line);
  std::vector<std::string> species;
  std::istringstream ss2(line);
  std::string specie;
  while (ss2 >> specie) {
    species.push_back(specie);
  }

  // Number of species line
  std::getline(file, line);
  std::vector<int> NumSpecies;
  std::istringstream ss3(line);
  int spec;
  while (ss3 >> spec) {
    NumSpecies.push_back(spec);
  }

  // Coord-type line
  std::getline(file, line);
  std::string coord_type;
  std::istringstream ss4(line);
  ss4 >> coord_type;

  // Raw coordinates
  Eigen::Vector3d coord;
  std::vector<Eigen::Vector3d> raw_coordinate_values;
  while (std::getline(file, line)) {
    std::istringstream ss5(line);
    ss5 >> coord(0) >> coord(1) >> coord(2);
    raw_coordinate_values.push_back(coord);
  }

  // make Lattice
  Lattice latt(lat_row_matrix);
  // making the sites, creating one site per atom
  std::vector<Site> Sites;
  int t = 0;
  for (int i = 0; i < species.size(); i++)
    for (int j = 0; j < NumSpecies[i]; j++) {
      Coordinate temp(raw_coordinate_values[t]);
      Site Single_site(species[i], temp);
      Sites.push_back(Single_site);
      t = t + 1;
    }

  // making the structure
  // for instance printing first row
  return Structure(latt, Sites);
};

void EXPECT_T(bool true_statement, std::string test_type) {
  if (true_statement) {
    std::cout << "PASS" << std::endl;
  } else {
    std::cout << "FAILED: " << test_type << std::endl;
  }
}

int main() {

  // Test Coordinate class
  //
  // Test Coordinate comparator
  //
  // Test Site class
  //
  // Test Site Comparator

  // First Read POSCAR and get Structure
  Structure my_structure = read_poscar("POSCAR");

  // get_sites from structure
  std::vector<Site> my_sites = my_structure.get_sites();
  EXPECT_T(my_sites.size() == 5, "Wrong number of sites in structure");
  EXPECT_T(my_sites[0].get_atom() == "Fe", "Expected Fe for first atom");
  EXPECT_T(my_sites.back().get_atom() == "Se", "Expected Se for last atom");
  SiteCompare_f test_site(my_sites[0], 1e-5);
  EXPECT_T(test_site(my_sites[0]), "Site doesn't recognize itself");
  EXPECT_T(!test_site(my_sites[1]),
           "Site incorrectly believes it is equal to another site");
  // get_lattice from structure
  Lattice my_lattice = my_structure.get_lattice();
  Eigen::Matrix3d raw_poscar_lattice;
  raw_poscar_lattice << 3.8511262576863907, 0.0000000000000000,
      0.0000000000000000, -1.9255631288431954, 3.3351731755587712,
      0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
      6.1326753166606398;

  // Stop using isApprox, build Coordinates/Sites manually,then use comparators
  EXPECT_T(my_lattice.col_vector_matrix().isApprox(raw_poscar_lattice, 1e-5),
           "Lattice matrix does not match POSCAR");
  EXPECT_T(my_lattice.lattice_vector(0).isApprox(
               Eigen::Vector3d(raw_poscar_lattice.col(0)), 1e-5),
           "Lattice vector does not match");

  // Coordinate test
  // Is site coordinate the same as Coordinate coordinate?
  Eigen::MatrixXd raw_coordinate_rows(5, 3);
  raw_coordinate_rows << 0.3333330000000032, 0.6666669999999968,
      0.6257012643865139, 0.0000000000000000, 0.0000000000000000,
      0.3495550691578657, 0.6666669999999968, 0.3333330000000032,
      0.0268432865670718, 0.0000000000000000, 0.0000000000000000,
      0.7508303712375408, 0.3333330000000032, 0.6666669999999968,
      0.2470700086510149;

  for (int i = 0; i < my_sites.size(); ++i) {
    EXPECT_T(Eigen::Vector3d(raw_coordinate_rows.row(i))
                 .isApprox(my_sites.at(i).get_coordinate()),
             "Coordinate mismatch for " + std::to_string(i) + "th site");
  }

  // Make test cluster
  std::vector<Site> cluster_sites;
  for (int i = 0; i < 3; i++) {
    cluster_sites.push_back(my_sites[i]);
  }
  Cluster test_cluster(cluster_sites);

  cluster_sites.push_back(my_sites[3]);
  Cluster other_test_cluster(cluster_sites);
  // test test cluster
  EXPECT_T(test_cluster.cluster_size() == 3, "Wrong Cluster Size");
  ClusterCompare_f first_cluster(test_cluster, 1e-5);
  EXPECT_T(first_cluster(test_cluster),
           "Cluster Doesn't even recognize itself");
  EXPECT_T(!first_cluster(other_test_cluster),
           "Cluster incorrectly believes it is another site");
  ////////////////////////////////////
  // EXPECT_T(find_factor_group(Valid_symops, my_lattice), "Doesn't have factor
  // group over 0");
  //
  // Test SitePeriodicCompare_f
  //
  // Test point group
  
  // Test transform_basis
  //
  // Test basis_maps_onto_itself

  // Test SymOp
 
  // Test make_factor_group
  auto factorgroup=find_factor_group(my_structure);
  EXPECT_T(factorgroup.size()==6, "Wrong number of factor Group Operations");

  //Test asymmetric unit without factor group routine
  //Start with an asymmetric unit (make it up)
  //Start with a small symmetry group (make it up)
  //Make a basis by applying symmetry group to asymmetric unit
  //call make_asymmetric_unit on basis
  //make sure you get the original asymmetric unit back

  auto asym_units=make_asymmetric_unit(cluster_sites, factorgroup);
  EXPECT_T(asym_units.size()>=3, "Wrong number of asymmetric units");
  auto geometric_center	=find_geometric_center(test_cluster);

  return 0;

}
