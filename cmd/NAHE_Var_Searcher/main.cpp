#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <cstdlib>

#include "FF_Systematic_Builder.hh"
#include "FF_Model_Builder.hh"

using namespace std;

string Make_Output_File_Name(string Base_Output_File_Name, 
			     vector<int> Extension_Orders, vector<bool>
			     Simply_Paired_Layers);

int main(int argc, char* argv[])
{
  if(argc == 1)
    {
      cout<<"Add an order or orders to specify search parameters."<<endl;
      return 0;
    }//Close if statement on proper number of arguments.

  int Large_ST_Dimensions = 4;
  string Kodiak = "/data/rennert/";
  string Base_Output_File_Name = "NAHE_Var_NSV_";
  Base_Output_File_Name = Kodiak+Base_Output_File_Name;
  vector<int> Extension_Orders;
  vector<bool> Simply_Paired_Layers;


  for(int a=1; a<argc; ++a)
    Extension_Orders.push_back(atoi(argv[a]));

  for(int a=1; a<argc; ++a)
    {
      if(Extension_Orders.at(a-1)<10&&argv[a][1] == 'n')
	Simply_Paired_Layers.push_back(false);
      else if(Extension_Orders.at(a-1)>=10&&argv[a][2] == 'n')
	Simply_Paired_Layers.push_back(false);
      else
	Simply_Paired_Layers.push_back(true);
    }

  string Output_File_Name = Make_Output_File_Name(Base_Output_File_Name, 
						  Extension_Orders, 
						  Simply_Paired_Layers);

  cout<<"Output_File_Name: "<<Output_File_Name<<endl;

  Model_Builder Initial_Conditions(Large_ST_Dimensions);
  //Initial_Conditions.Load_S_Vector(Large_ST_Dimensions);
  Initial_Conditions.Load_NAHE_Variation();

  GSO_Coefficient_Matrix Initial_k_ij = Initial_Conditions.FFHS_Model().k_ij();
  if(!Initial_Conditions.Check_Model_Consistency())
    {
      cout<<"Bad initial list."<<endl;
      return 0;
    }//Close if statement on model consistency.

  Initial_Conditions.rFFHS_Model().Set_k_ij(Initial_k_ij);

  Systematic_Builder Bob_The_Builder(Initial_Conditions, Large_ST_Dimensions,
				     Extension_Orders, Output_File_Name,
				     Simply_Paired_Layers);
  Bob_The_Builder.Perform_Systematic_Search();
  //Bob_The_Builder.Display_k_ij_Extensions();

  cout<<"Search finished."<<endl;
  return 0;
}//Close main.

string Make_Output_File_Name(string Base_Output_File_Name, 
			     vector<int> Extension_Orders, vector<bool>
			     Simply_Paired_Layers)
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
