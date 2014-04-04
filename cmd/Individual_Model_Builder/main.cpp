/*
Author: Timothy Renner (renner.timothy@gmail.com)
Date: 10/27/2010
Baylor University

This program builds individual models and outputs the basis vectors,
k_ij matrix, and particle content to a LaTeX file.

Uses the FF Framework.
*/

#include <vector>
#include <iostream>
#include <fstream>

#include "FF_Basis_Vector.hh"
#include "FF_Model_Builder.hh"
#include "FF_Latex_Writer.hh"
#include "FF_LEEFT.hh"

using namespace std;

vector<int> Make_Vector(int Array[], int BV_Size);

int main()
{
  int Large_ST_Dimensions = 4;
  int Order = 3;
  int BV_Size = 64;
  int BV_Array[64] = {1,1, //ST1, ST2
		      1,0,0,1,0,0,//x12,y12,w12
		      0,0,1,0,0,1,//x34,y34,w34
		      0,0,1,0,0,1,//x56,y56,w56
		      1,1,1,1,1,1,1,1,1,1,//bar-psi1-10
		      1,1,1,1,1,1,//bar-eta1-6
		      0,0,0,0,1,1,//bar-y1-6
		      1,1,1,1,1,1,//bar-w1-6
		      0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1};//bar-phi1-16
  //*/
  /*int BV_Array[40] = {1,1,1,1,1,1,1,1,
		      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  int L2_BV_Array[40] = {1,1,1,1,1,1,1,1,
		      0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,
		      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  int L2_Order = 2; // */


  Model_Builder NAHE_Builder(Large_ST_Dimensions);
  NAHE_Builder.Load_S_Vector(Large_ST_Dimensions);
  NAHE_Builder.Load_NAHE_Set();
  //NAHE_Builder.Load_NAHE_Variation();
  vector<int> k_ij_Row(5,0);
	k_ij_Row.at(1) = 1;
	//k_ij_Row.at(2) = 1;
	//k_ij_Row.at(3) = 1;
	//k_ij_Row.at(4) = 1;
  //k_ij_Row.push_back(1);
  // k_ij_Row.push_back(0);
  //k_ij_Row.push_back(0);
  //k_ij_Row.push_back(Order);
  //k_ij_Row.push_back(0);

 NAHE_Builder.Load_Basis_Vector
	 (Basis_Vector(Make_Vector(BV_Array, BV_Size), Order));
  //NAHE_Builder.Load_Basis_Vector(Basis_Vector(Make_Vector(BV_Array, BV_Size),
		//			      L2_Order));
  NAHE_Builder.Load_k_ij_Row(k_ij_Row);
  //NAHE_Builder.Load_k_ij_Row(vector<int>(1,1));
  //vector<int> k_ij_Row (1,1);
  //k_ij_Row.push_back(0);
  //NAHE_Builder.Load_k_ij_Row(k_ij_Row);

  if(NAHE_Builder.Check_Linear_Independence())
    cout<<"Linearly independent."<<endl;
  else
    cout<<"Linearly dependent - not consistent."<<endl;

  if(NAHE_Builder.Check_k_ij_Consistency())
    cout<<"k_ij good."<<endl;
  else
    cout<<"k_ij inconsistent - inconsistent model."<<endl;

  if(NAHE_Builder.Check_Modular_Invariance())
    NAHE_Builder.Build_Model();
  else
    {
      cout<<"Not modular invariant - inconsistent model."<<endl;
      return 0;
    }
  NAHE_Builder.FFHS_Model().Display_Particle_Content();
  cout<<endl;
  NAHE_Builder.FFHS_Model().Display_k_ij();
  cout<<endl;
  ofstream Latex_Out("/Users/timothyrenner/FF_Framework/tex/FFHS_Model.tex");
  Latex_Writer Stephen_King;
  Stephen_King.Print_Single_Model_Document(Latex_Out, NAHE_Builder.FFHS_Model());

  Latex_Out.close();

	LEEFT New_LEEFT(NAHE_Builder.FFHS_Model());
	cout<<"New_LEEFT matter representation classes: "<<endl;
	New_LEEFT.Display_Matter_Rep_Classes();
  return 0;
}//Close main.

vector<int> Make_Vector(int Array[], int BV_Size)
{
  vector<int> New_Vector;
  for(int a=0; a<BV_Size; a++)
    New_Vector.push_back(Array[a]);
  return New_Vector;
}//Close Make_Vector.
