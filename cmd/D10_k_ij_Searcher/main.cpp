#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <cstdlib>

#include "FF_Systematic_k_ij_Builder.hh"
#include "FF_Model_Builder.hh"

using namespace std;

string Make_Output_File_Name(string Base_Output_File_Name, 
			     vector<int> Extension_Orders);
vector<int> Make_Vector(int Array[], int BV_Size);

int main() 
{
  int Large_ST_Dimensions = 10;
  string Kodiak = "/data/rennert/";
  string Base_Output_File_Name = "D10_k_ij";
  //Base_Output_File_Name = Kodiak + Base_Output_File_Name;


  vector<int> Extension_Orders;
  vector<Basis_Vector> Basis_Vectors;
  //Set the inputs.
  int L1_BV_Array[40] = {1,1,1,1,1,1,1,1,
			 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			 0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1};
  Extension_Orders.push_back(3);
  Basis_Vectors.push_back(Basis_Vector(Make_Vector(L1_BV_Array, 40), 
				       Extension_Orders.at(0)));

  int L2_BV_Array[40] = {1,1,1,1,1,1,1,1,
			 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,
			 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  // Extension_Orders.push_back(2);
  //Basis_Vectors.push_back(Basis_Vector(Make_Vector(L2_BV_Array, 40),
  //			       Extension_Orders.at(1)));

  /*int L3_BV_Array[40] = {1,1,1,1,1,1,1,1,
			 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,
			 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1};
  Extension_Orders.push_back(2);
  Basis_Vectors.push_back(Basis_Vector(Make_Vector(L3_BV_Array, 40),
  Extension_Orders.at(2)));//*/

  string Output_File_Name = Make_Output_File_Name(Base_Output_File_Name, 
						  Extension_Orders);

  cout<<"Output_File_Name: "<<Output_File_Name<<endl;

  Model_Builder Initial_Conditions(Large_ST_Dimensions);
  //Initial_Conditions.Load_Default_k_ij();
  vector<int> k_ij_Row;
  k_ij_Row.push_back(1);
  Initial_Conditions.Load_k_ij_Row(k_ij_Row);
  GSO_Coefficient_Matrix Initial_k_ij = Initial_Conditions.FFHS_Model().k_ij();

  if(!Initial_Conditions.Check_Model_Consistency())
    {
      cout<<"Bad initial list."<<endl;
      return 0;
    }//Close if statement on model consistency.

  Initial_Conditions.rFFHS_Model().Set_k_ij(Initial_k_ij);

  Systematic_k_ij_Builder Bob_The_Builder(Initial_Conditions, Large_ST_Dimensions,
					  Extension_Orders, Output_File_Name, 
					  Basis_Vectors);

  //Bob_The_Builder.Display_k_ij_Extensions();
  Bob_The_Builder.Perform_Search();
  cout<<"Search finished."<<endl;
  return 0;
}//Close main.

string Make_Output_File_Name(string Base_Output_File_Name, 
			     vector<int> Extension_Orders)
{
  stringstream ss;
  ss<<Base_Output_File_Name;
  for(int a=0; a<Extension_Orders.size(); ++a)
    {
      ss<<"_O"<<Extension_Orders.at(a);
      ss<<"_L"<<(a+1);
    }
  ss<<".txt";
  return ss.str();
}//Close Make_Output_File_Name.

vector<int> Make_Vector(int Array[], int BV_Size)
{
  vector<int> New_Vector;
  for(int a=0; a<BV_Size; ++a)
    New_Vector.push_back(Array[a]);
  return New_Vector;
}
