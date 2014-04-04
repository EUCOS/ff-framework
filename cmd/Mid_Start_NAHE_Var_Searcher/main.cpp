#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <cstdlib>

#include "FF_Mid_Start_Systematic_Builder.hh"
#include "FF_Model_Builder.hh"

using namespace std;

string Make_Output_File_Name(string Base_Output_File_Name,
			     vector<int> Extension_Orders, vector<bool>
			     Simply_Paired_Layers);
vector<int> Make_Vector(int BV_Array[], int BV_Size);

int main(int argc, char* argv[])
{
  int Large_ST_Dimensions = 4;
  string Kodiak = "/data/rennert/";
  string Base_Output_File_Name = "NAHE_Var_Mid";
  //Base_Output_File_Name = Kodiak+Base_Output_File_Name;

  vector<int> Extension_Orders;
  vector<bool> Simply_Paired_Layers;



  int BV_Array_L1[64] = {1,1,1,0,0,1,1,1,//ST 1,2 xyw 1,2
			 1,0,0,1,0,0,//xyw 3,4
			 1,1,1,1,0,0,//xyw 5,6
			 0,0,0,0,1,1,1,1,1,1,//bar-psi 1-10
			 2,2,1,1,0,0,//bar-eta 1-6
			 0,2,1,1,0,2,//bar-y 1-6
			 2,2,2,2,2,2,//barw 1-6
			 0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2};//bar-phi 1-16

  Extension_Orders.push_back(4);
  Simply_Paired_Layers.push_back(false);

  string Output_File_Name = Make_Output_File_Name(Base_Output_File_Name,
						  Extension_Orders,
						  Simply_Paired_Layers);

  cout<<"Output_File_Name: "<<Output_File_Name<<endl;

  vector<Basis_Vector> Starting_Basis_Vectors;
  Starting_Basis_Vectors.push_back(Basis_Vector(Make_Vector(BV_Array_L1, 64), 
						Extension_Orders.at(0)));

  Model_Builder Initial_Conditions(Large_ST_Dimensions);
  Initial_Conditions.Load_S_Vector(Large_ST_Dimensions);
  Initial_Conditions.Load_NAHE_Variation();

  GSO_Coefficient_Matrix Initial_k_ij = Initial_Conditions.FFHS_Model().k_ij();

  if(!Initial_Conditions.Check_Model_Consistency())
    {
      cout<<"Bad initial list."<<endl;
      return 0;
    }//Close if statement on model consistency.

  Initial_Conditions.rFFHS_Model().Set_k_ij(Initial_k_ij);

  Mid_Start_Systematic_Builder Bob_The_Builder(Initial_Conditions, 
					       Large_ST_Dimensions, 
					       Extension_Orders,
					       Output_File_Name,
					       Simply_Paired_Layers,
					       Starting_Basis_Vectors);
  Bob_The_Builder.Perform_Search();

  cout<<"Search finished."<<endl;
  return 0;
}//Close main.

string Make_Output_File_Name(string Base_Output_File_Name,
			     vector<int> Extension_Orders,
			     vector<bool> Simply_Paired_Layers)
{
  stringstream ss;
  ss<<Base_Output_File_Name;
  for(int a=0; a<Extension_Orders.size(); ++a)
    {
      ss<<"_O"<<Extension_Orders.at(a);
      ss<<"_L"<<(a+1);
      if(!Simply_Paired_Layers.at(a))
	ss<<"_NSP";
    }
  ss<<".txt";
  return ss.str();
}//Close Make_Output_File_Name.

vector<int> Make_Vector(int BV_Array[], int BV_Size)
{
  vector<int> New_Vector;
  for(int a=0; a<BV_Size; ++a)
    New_Vector.push_back(BV_Array[a]);
  return New_Vector;
}//Close Make_Vector.
