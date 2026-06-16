#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <GSOCoefficientMatrix.h>
#include <ModelBuilder.h>
#include <SystematicGaugeBuilder.h>

namespace {

int ParseInt(const char* value) {
  return std::atoi(value);
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 6 || std::string(argv[1]) != "--dimension" ||
      std::string(argv[3]) != "--output") {
    std::cerr << "Usage: HumanGaugeBaseline --dimension 4 --output <path> <order>\n";
    return 1;
  }

  const int dimension = ParseInt(argv[2]);
  if (dimension != 4) {
    std::cerr << "HumanGaugeBaseline supports only dimension 4\n";
    return 1;
  }

  std::vector<int> orders;
  orders.push_back(ParseInt(argv[5]));
  const std::string output = argv[4];

  ModelBuilder initial(dimension);
  initial.Load_S_Vector(dimension);
  initial.Load_Default_k_ij();
  GSOCoefficientMatrix initial_kij = initial.FFHS_Model().k_ij();
  if (!initial.Check_Model_Consistency()) {
    std::cerr << "Bad initial list.\n";
    return 1;
  }
  initial.rFFHS_Model().Set_k_ij(initial_kij);

  SystematicGaugeBuilder builder(initial, dimension, orders, output);
  builder.Perform_Search();
  return 0;
}
