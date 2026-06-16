#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <GSOCoefficientMatrix.h>
#include <ModelBuilder.h>
#include <SystematicGaugeBuilder.h>

using namespace std;

string Make_Output_File_Name(const string& Base_Output_File_Name,
                             const vector<int>& Extension_Orders) {
  stringstream ss;
  ss << Base_Output_File_Name;
  for (size_t a = 0; a < Extension_Orders.size(); ++a) {
    ss << "_O" << Extension_Orders.at(a);
    ss << "_L" << (a + 1);
  }
  ss << ".txt";
  return ss.str();
}

int main(int argc, char* argv[]) {
  if (argc == 1) {
    cout << "Add an order or orders to specify search parameters." << endl;
    return 0;
  }

  const int Large_ST_Dimensions = 4;
  const string Base_Output_File_Name =
      "AI_GENERATED/benchmarks/modern_d4_gauge_wrapper/D4_Gauge";

  vector<int> Extension_Orders;
  for (int a = 1; a < argc; ++a) {
    Extension_Orders.push_back(atoi(argv[a]));
  }

  const string Output_File_Name =
      Make_Output_File_Name(Base_Output_File_Name, Extension_Orders);

  cout << "Output_File_Name: " << Output_File_Name << endl;

  ModelBuilder Initial_Conditions(Large_ST_Dimensions);
  Initial_Conditions.Load_S_Vector(Large_ST_Dimensions);
  Initial_Conditions.Load_Default_k_ij();

  GSOCoefficientMatrix Initial_k_ij =
      Initial_Conditions.FFHS_Model().k_ij();

  if (!Initial_Conditions.Check_Model_Consistency()) {
    cout << "Bad initial list." << endl;
    return 1;
  }

  Initial_Conditions.rFFHS_Model().Set_k_ij(Initial_k_ij);

  SystematicGaugeBuilder Bob_The_Builder(
      Initial_Conditions,
      Large_ST_Dimensions,
      Extension_Orders,
      Output_File_Name);

  Bob_The_Builder.Perform_Search();

  cout << "Search finished." << endl;
  return 0;
}

