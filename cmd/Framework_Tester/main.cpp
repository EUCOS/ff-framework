#include <fstream>
#include <iostream>
#include <vector>

#include <BasisVector.h>
#include <LEEFT.h>
#include <LatexWriter.h>
#include <ModelBuilder.h>

using namespace std;

/* This program builds individual models and outputs the basis vectors,
 * k_ij matrix, and particle content to a LaTeX file.
 */
int main() {
  cout << "Begin program ... forgive me." << endl;
  const int Large_ST_Dimensions = 4;
  const int Order = 2;
  const int BV_Size = 64;
  int BV_Array[BV_Size] = {
    1,1,                              //ST1, ST2
    1,0,0,1,0,0,                      //x12,y12,w12
    1,0,0,1,0,0,                      //x34,y34,w34
    1,1,1,1,1,1,                      //x56,y56,w56
    0,0,1,1,1,1,1,1,1,1,              //bar-psi1-10
    1,1,1,1,0,0,                      //bar-eta1-6
    1,1,1,1,1,1,                      //bar-y1-6
    0,0,0,0,1,1,                      //bar-w1-6
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; //bar-phi1-16

  int BV_Array_2[BV_Size] = {
    1,1,                              //ST1, ST2
    1,0,0,1,0,0,                      //x12,y12,w12
    1,1,1,1,1,1,                      //x34,y34,w34
    1,0,0,1,0,0,                      //x56,y56,w56
    0,0,1,1,1,1,1,1,1,1,              //bar-psi1-10
    1,1,0,0,1,1,                      //bar-eta1-6
    0,0,1,1,1,1,                      //bar-y1-6
    1,1,1,1,0,0,                      //bar-w1-6
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; //bar-phi1-16

  ModelBuilder NAHE_Builder(Large_ST_Dimensions);
  NAHE_Builder.Load_S_Vector(Large_ST_Dimensions);
  NAHE_Builder.Load_NAHE_Set();

  vector<int> k_ij_Row(5,0);
  //k_ij_Row.at(1) = 1;
  //k_ij_Row.at(2) = 1;
  //k_ij_Row.at(3) = 1;
  //k_ij_Row.at(4) = 1;

  NAHE_Builder.Load_Basis_Vector(BasisVector(std::vector<int>(BV_Array, BV_Array + BV_Size), Order));
  NAHE_Builder.Load_k_ij_Row(k_ij_Row);

  if (NAHE_Builder.Check_Linear_Independence()) {
    cout << "Linearly independent." << endl;
  } else {
    cout << "Linearly dependent - not consistent." << endl;
  }

  if (NAHE_Builder.Check_k_ij_Consistency()) {
    cout << "k_ij good." << endl;
  } else {
    cout << "k_ij inconsistent - inconsistent model." << endl;
  }

  if (!NAHE_Builder.Check_Modular_Invariance()) {
    cout << "Not modular invariant - inconsistent model." << endl;
    return 0;
  }

  NAHE_Builder.Build_Model();
  NAHE_Builder.FFHS_Model().Display_Particle_Content();
  cout << endl;
  NAHE_Builder.FFHS_Model().Display_k_ij();
  cout << endl;

  ModelBuilder NAHE_Builder_2(Large_ST_Dimensions);
  NAHE_Builder_2.Load_S_Vector(Large_ST_Dimensions);
  NAHE_Builder_2.Load_NAHE_Set();
  vector<int> k_ij_Row_2(5,0);
  //k_ij_Row.at(1) = 1;
  //k_ij_Row.at(2) = 1;
  //k_ij_Row_2.at(3) = 1;
  //k_ij_Row_2.at(4) = 1;

  NAHE_Builder_2.Load_Basis_Vector(BasisVector(std::vector<int>(BV_Array_2, BV_Array_2 + BV_Size), Order));
  NAHE_Builder_2.Load_k_ij_Row(k_ij_Row_2);

  if (NAHE_Builder_2.Check_Linear_Independence()) {
    cout << "Linearly independent." << endl;
  } else {
    cout << "Linearly dependent - not consistent." << endl;
  }

  if (NAHE_Builder_2.Check_k_ij_Consistency()) {
    cout << "k_ij good." << endl;
  } else {
    cout << "k_ij inconsistent - inconsistent model." << endl;
  }

  if (!NAHE_Builder_2.Check_Modular_Invariance()) {
    cout << "Not modular invariant - inconsistent model." << endl;
    return 0;
  }

  NAHE_Builder_2.Build_Model();
  NAHE_Builder_2.FFHS_Model().Display_Particle_Content();
  cout << endl;
  NAHE_Builder_2.FFHS_Model().Display_k_ij();
  cout << endl;

  LEEFT New_LEEFT(NAHE_Builder.FFHS_Model());
  cout << "New_LEEFT matter representation classes: " << endl;
  New_LEEFT.Display_Matter_Rep_Classes();

  LEEFT New_LEEFT_2(NAHE_Builder_2.FFHS_Model());
  cout << "New_LEEFT_2 matter representation classes: " << endl;
  New_LEEFT_2.Display_Matter_Rep_Classes();

  if (New_LEEFT<New_LEEFT_2) {
    cout << "New_LEEFT less than New_LEEFT_2" << endl;
  } else if (New_LEEFT_2<New_LEEFT) {
    cout << "New_LEEFT greater than New_LEEFT_2" << endl;
  }

  if (New_LEEFT == New_LEEFT_2) {
    cout << "New_LEEFT equal to New_LEEFT_2" << endl;
  }

  if (New_LEEFT.Gauge_Groups()<New_LEEFT_2.Gauge_Groups()) {
    cout << "New_LEEFT gauge groups < New_LEEFT_2 gauge groups." << endl;
  } else if (New_LEEFT.Gauge_Groups() > New_LEEFT_2.Gauge_Groups()) {
    cout << "New_LEEFT gauge groups > New_LEEFT_2 gauge groups." << endl;
  } else {
    cout << "New_LEEFT gauge groups == New_LEEFT_2 gauge groups." << endl;
  }

  if (New_LEEFT.ST_SUSYs() < New_LEEFT_2.ST_SUSYs()) {
    cout << "New_LEEFT ST SUSY < New_LEEFT_2 ST SUSY." << endl;
  } else if (New_LEEFT.ST_SUSYs() > New_LEEFT_2.ST_SUSYs()) {
    cout << "New_LEEFT ST SUSY > New_LEEFT_2 ST SUSY." << endl;
  } else {
    cout << "New_LEEFT ST SUSY == New_LEEFT_2 ST SUSY." << endl;
  }

  if (New_LEEFT.U1_Factors()<New_LEEFT_2.U1_Factors()) {
    cout << "New_LEEFT U1 Factors < New_LEEFT_2 U1 Factors." << endl;
  } else if (New_LEEFT.U1_Factors() > New_LEEFT_2.U1_Factors()) {
    cout << "New_LEEFT U1 Factors > New_LEEFT_2 U1_Factors." << endl;
  } else {
    cout << "New_LEEFT U1 factors == New_LEEFT_2 U1 factors." << endl;
  }

  cout << "New_LEEFT Total_NA_Reps: " << New_LEEFT.Total_NA_Reps() << endl;
  cout << "New_LEEFT_2 Total_NA_Reps: " << New_LEEFT_2.Total_NA_Reps() << endl;

  bool Equivalent_Matter_Rep_Classes = New_LEEFT.Equal_Matter_Rep_Classes(New_LEEFT.Matter_Rep_Classes(), New_LEEFT_2.Matter_Rep_Classes());
  if (Equivalent_Matter_Rep_Classes) {
    cout << "matter rep classes are equivalent." << endl;
  } else {
    cout << "You screwed up the matter rep equivalence classes." << endl;
  }
  return 0;
}
