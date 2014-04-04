/*
Author: Timothy Renner (renner.timothy@gmail.com)
Date: 11/10/2010
Baylor University

This program builds a model while outputting the intermediate data (sectors,
states from that sector, positive roots for gauge groups, as well as 
the highest weight charge vectors) to a file.

Uses the FF Framework.
*/

//C++ includes.
#include <vector>
#include <set>
#include <list>
#include <string>
#include <iostream>
#include <fstream>

//FF Framework includes.
#include "FF_Output_Writer.hh"
#include "FF_Basis_Vector.hh"
#include "FF_Model.hh"
#include "FF_Basis_Alpha_Builder.hh"
#include "FF_Basis_Alpha.hh"
#include "FF_Modular_Invariance_Checker.hh"
#include "FF_Alpha_Builder.hh"
#include "FF_Alpha.hh"
#include "FF_Alpha_Boson.hh"
#include "FF_Alpha_SUSY.hh"
#include "FF_Alpha_Fermion.hh"
#include "FF_Fermion_Mode_Map_Builder.hh"
#include "FF_GSO_Coefficient_Matrix.hh"
#include "FF_GSO_Coefficient_Matrix_Builder.hh"
#include "FF_State_Builder.hh"
#include "FF_State.hh"
#include "FF_Gauge_Group_Identifier.hh"
#include "FF_Gauge_Group.hh"
#include "FF_Gauge_Group_Name.hh"
#include "FF_Matter_State.hh"
#include "FF_Output_Writer.hh"
#include "FF_NAHE_Variation_Loader.hh"
#include "FF_NAHE_Set_Loader.hh"

using namespace std;

vector<int> Make_Vector(int BV_Array[], int BV_Size);
int Gauge_Dot(const State& State1, const State& State2, 
	      const map<int, int>& Fermion_Mode_Map);
bool Is_SUSY_Partner(const State& New_Matter_State, 
		     const list<State>& SUSY_States);
int Find_Rank_Cuts(const Model& Heidi_Klum, const map<int, int>& Fermion_Mode_Map);

int main()
{
  cout<<"Begin program."<<endl;
  string Output_File_Name = "Model_Output.txt";
  Output_Writer Stephen_King;
  ofstream File_Out(Output_File_Name.c_str());

  //First, set the inputs.
  int Large_ST_Dimensions = 4;
  int Order = 3;
  /*int BV_L1_Array[56] = {1,1,1,1,
			 1,0,0,1,0,0,
			 1,0,0,1,0,0,
			 0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,
			 1,1,1,1,
			 1,1,1,1,
			 1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2};//*/

  int BV_L1_Array[64] = {1,1,//ST1,2
			 1,0,0,1,0,0,//xyw,1,2
			 0,1,0,0,1,0,//xyw3,4
			 0,1,0,0,1,0,//xyw5,6
			 0,0,0,0,0,0,0,0,0,0,//bar-psi
			 1,1,1,1,1,1,//bar-eta
			 0,0,0,0,0,0,//bar-y
			 0,0,0,0,0,0,//bar-w
			 0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1};//bar-phi*/

  /* int BV_L2_Array[64] = {0,0,
			 0,0,0,0,0,0,
			 0,0,0,0,0,0,
			 0,0,0,0,0,0,
			 3,3,3,3,0,0,0,0,1,1,
			 1,1,1,1,3,3,
			 3,3,3,3,0,0,
			 0,0,0,0,0,0,
			 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};//*/
  //int L2_Order = 4;

  /*int BV_L1_Array[40] = {1,1,1,1,1,1,1,1,//ST1-8
    		 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,//Observable
  		 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};//Hidden 
  int BV_L2_Array[40] = {1,1,1,1,1,1,1,1,
			 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,
			 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};//
  int L2_Order = 2;
  int BV_L3_Array[40] = {1,1,1,1,1,1,1,1,//ST1-8
			 0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,
			 0,0,1,1,1,1,1,1,0,0,1,1,0,0,1,1,};//
			 int L3_Order = 2;//*/

  //Load the inputs.
  Model Heidi_Klum(Large_ST_Dimensions);
  Heidi_Klum.Load_S_Vector(Large_ST_Dimensions);

 // NAHE_Set_Loader The_Loader;
  //Heidi_Klum = The_Loader.Load_NAHE_Set(Heidi_Klum, true);
  NAHE_Variation_Loader The_Loader;
  Heidi_Klum = The_Loader.Load_NAHE_Variation(Heidi_Klum, true);
  //Heidi_Klum.Load_k_ij_Row(vector<int>(5,0));
  Heidi_Klum.Load_BV_Set(Basis_Vector(Make_Vector(BV_L1_Array, 64), Order));
  //Heidi_Klum.Load_BV_Set(Basis_Vector(Make_Vector(BV_L2_Array, 40), L2_Order));
  //  Heidi_Klum.Load_BV_Set(Basis_Vector(Make_Vector(BV_L3_Array, 40), L3_Order));

  vector<int> k_ij_Row(5,0);
	k_ij_Row.at(3) = 1;
	k_ij_Row.at(4) = 1;
  //Heidi_Klum.Load_k_ij_Row(k_ij_Row);
	//Heidi_Klum.Load_k_ij_Row(k_ij_Row);
  //k_ij_Row.at(0) = 1;
  //k_ij_Row.push_back(0);
  Heidi_Klum.Load_k_ij_Row(k_ij_Row);
  // k_ij_Row.push_back(0);
  //Heidi_Klum.Load_k_ij_Row(k_ij_Row);*/

  //Write the basis vectors.
  cout<<"Basis Vectors: "<<Heidi_Klum.BV_Set().size()<<endl;
  File_Out<<"Basis Vectors: "<<Heidi_Klum.BV_Set().size()<<endl;
  for(int a=0; a<Heidi_Klum.BV_Set().size(); a++)
    {
      cout<<Heidi_Klum.BV_Set().at(a).Order()<<": ";
      File_Out<<Heidi_Klum.BV_Set().at(a).Order()<<": ";
      for(int b=0; b<Heidi_Klum.BV_Set().at(a).LM_Size(); b++)
	{
	  cout<<Heidi_Klum.BV_Set().at(a).BV().at(b)<<" ";
	  File_Out<<Heidi_Klum.BV_Set().at(a).BV().at(b)<<" ";
	}//Close for loop on actual basis vector output.
      cout<<"|| ";
      File_Out<<"|| ";
      for(int b=Heidi_Klum.BV_Set().at(a).LM_Size(); 
	  b<Heidi_Klum.BV_Set().at(a).BV().size(); b++)
	{
	  cout<<Heidi_Klum.BV_Set().at(a).BV().at(b)<<" ";
	  File_Out<<Heidi_Klum.BV_Set().at(a).BV().at(b)<<" ";
	}//Close for loop on RM basis vector output.
      cout<<endl;
      File_Out<<endl;
    }//Close for loop on basis vector selection.
  cout<<endl;
  File_Out<<endl;

  //Build the basis alphas.
  Basis_Alpha_Builder BA_Baracus(Heidi_Klum);
  BA_Baracus.Build_Basis_Alphas();
  BA_Baracus.Build_Common_Basis_Alphas();

  vector<Basis_Alpha> Basis_Alphas = BA_Baracus.Basis_Alphas();
  vector<Basis_Alpha> Common_Basis_Alphas = BA_Baracus.Common_Basis_Alphas();

  //Display the basis and common basis alphas.
  BA_Baracus.Display_Basis_Alphas();
  BA_Baracus.Display_Common_Basis_Alphas();

  //Write the basis and common basis alphas.
  File_Out<<"Basis Alphas (fractional basis vector form,"<<
    " first number is denominator): "<<Basis_Alphas.size()<<endl;
  for(int a=0; a<Basis_Alphas.size(); a++)
    {
      File_Out<<Basis_Alphas.at(a).Denominator()<<": ";
      for(int b=0; b<Basis_Alphas.at(a).LM_Size(); b++)
	File_Out<<Basis_Alphas.at(a).Numerator().at(b)<<" ";
      File_Out<<"|| ";
      for(int b=Basis_Alphas.at(a).LM_Size(); 
	  b<Basis_Alphas.at(a).Numerator().size();
	  b++)
	File_Out<<Basis_Alphas.at(a).Numerator().at(b)<<" ";
      File_Out<<endl;
    }//Close for loop on Basis_Alphas.
  File_Out<<endl;

  File_Out<<"Common Basis Alphas (basis alphas with common denominator): "<<
    Common_Basis_Alphas.size()<<endl;
  for(int a=0; a<Common_Basis_Alphas.size(); a++)
    {
      File_Out<<Common_Basis_Alphas.at(a).Denominator()<<": ";
      for(int b=0; b<Common_Basis_Alphas.at(a).LM_Size(); b++)
	File_Out<<Common_Basis_Alphas.at(a).Numerator().at(b)<<" ";
      File_Out<<"|| ";
      for(int b=Common_Basis_Alphas.at(a).LM_Size(); 
	  b<Common_Basis_Alphas.at(a).Numerator().size();
	  b++)
	File_Out<<Common_Basis_Alphas.at(a).Numerator().at(b)<<" ";
      File_Out<<endl;
    }//Close for loop on Common_Basis_Alphas.
  File_Out<<endl;

  //Check for modular invariance.
  Modular_Invariance_Checker Cody_Rogers;
  if(!Cody_Rogers.Test_Modular_Invariance(Common_Basis_Alphas))
    {
      cout<<"Model is not modular invariant."<<endl;
      return 0;
    }
  
  //Build the alphas.
  Alpha_Builder Hannibal(Common_Basis_Alphas, Basis_Alphas);
  Hannibal.Build_Alphas();
  set<Alpha_Boson> Alpha_Bosons = Hannibal.Alpha_Bosons();
  set<Alpha_Fermion> Alpha_Fermions = Hannibal.Alpha_Fermions();
  set<Alpha_SUSY> Alpha_SUSYs = Hannibal.Alpha_SUSYs();

  //Display the alphas.
  Hannibal.Display_All_Alphas();

  //Write the Alpha_Bosons.
  File_Out<<"Boson sectors (first number is denominator): "
	  <<Alpha_Bosons.size()<<endl;
  set<Alpha_Boson>::iterator itAlpha_Bosons = Alpha_Bosons.begin();
  for(itAlpha_Bosons; itAlpha_Bosons != Alpha_Bosons.end(); ++itAlpha_Bosons)
    {
      File_Out<<itAlpha_Bosons->Denominator()<<": ";
      for(int a=0; a<itAlpha_Bosons->LM_Size(); a++)
	File_Out<<(itAlpha_Bosons->Numerator()).at(a)<<" ";
      File_Out<<"|| ";
      for(int a=itAlpha_Bosons->LM_Size(); 
	  a<(itAlpha_Bosons->Numerator()).size(); a++)
	File_Out<<(itAlpha_Bosons->Numerator()).at(a)<<" ";
      File_Out<<endl;
    }//Close for loop in boson sectors.
  File_Out<<endl;

  //Write the Alpha_Fermions.
  File_Out<<"Fermion sectors (first number is denominator): "
	  <<Alpha_Fermions.size()<<endl;
  set<Alpha_Fermion>::iterator itAlpha_Fermions = Alpha_Fermions.begin();
  for(itAlpha_Fermions; itAlpha_Fermions != Alpha_Fermions.end(); ++itAlpha_Fermions)
    {
      File_Out<<itAlpha_Fermions->Denominator()<<": ";
      for(int a=0; a<itAlpha_Fermions->LM_Size(); a++)
	File_Out<<(itAlpha_Fermions->Numerator()).at(a)<<" ";
      File_Out<<"|| ";
      for(int a=itAlpha_Fermions->LM_Size();
	  a<(itAlpha_Fermions->Numerator()).size(); a++)
	File_Out<<(itAlpha_Fermions->Numerator()).at(a)<<" ";
      File_Out<<endl;
    }//Close for loop on boson sectors.
  File_Out<<endl;

  //Write the Alpha SUSYs.
  File_Out<<"SUSY sectors (first number is denominator): "
	  <<Alpha_SUSYs.size()<<endl;
  set<Alpha_SUSY>::iterator itAlpha_SUSYs = Alpha_SUSYs.begin();
  for(itAlpha_SUSYs; itAlpha_SUSYs != Alpha_SUSYs.end(); ++itAlpha_SUSYs)
    {
      File_Out<<itAlpha_SUSYs->Denominator()<<": ";
      for(int a=0; a<itAlpha_SUSYs->LM_Size(); a++)
	File_Out<<(itAlpha_SUSYs->Numerator()).at(a)<<" ";
      File_Out<<"|| ";
      for(int a=itAlpha_SUSYs->LM_Size();
	  a<(itAlpha_SUSYs->Numerator()).size(); a++)
	File_Out<<(itAlpha_SUSYs->Numerator()).at(a)<<" ";
      File_Out<<endl;
    }//Close for loop on SUSY sectors.
  File_Out<<endl;

  //Check for linear independence.
 if(!Hannibal.Linearly_Independent_Alphas())
    {
      cout<<"Basis vectors are not linearly independent."<<endl;
      //return 0;
      }

  //Build the fermion mode map.
  Fermion_Mode_Map_Builder Face;
  Face.Build_Fermion_Mode_Map(Common_Basis_Alphas);
  map<int, int> Fermion_Mode_Map = Face.Fermion_Mode_Map();

  //Display the fermion mode map.
  Face.Display_Fermion_Mode_Map();

  //Write the fermion mode map.
  File_Out<<"Fermion mode map:"<<endl;
  map<int, int>::iterator itFermion_Mode_Map = Fermion_Mode_Map.begin();
  for(itFermion_Mode_Map; itFermion_Mode_Map != Fermion_Mode_Map.end();
      ++itFermion_Mode_Map)
    {
      File_Out<<itFermion_Mode_Map->first<<"-> "<<itFermion_Mode_Map->second<<endl;
    }//Close for loop on Fermion_Mode_Map.
  File_Out<<endl;

  //Complete the GSO coefficient matrix.
  GSO_Coefficient_Matrix_Builder HM_Murdock(Heidi_Klum, Common_Basis_Alphas);
  HM_Murdock.Build_Complete_GSO_Matrix();
  Heidi_Klum.Set_k_ij(HM_Murdock.Complete_GSO_Matrix());

  //Display the GSO coefficient matrix.
  Heidi_Klum.k_ij().Display();

  //Write the GSO coefficient matrix.
  File_Out<<"GSO coefficient matrix (first column are the denominators of the "<<
    "row): "<<endl;
  for(int a=0; a<Heidi_Klum.k_ij().Numerators().size(); a++)
    {
      File_Out<<Heidi_Klum.k_ij().Denominators().at(a)<<" | ";
      for(int b=0; b<Heidi_Klum.k_ij().Numerators().at(a).size(); b++)
	File_Out<<Heidi_Klum.k_ij().Numerators().at(a).at(b)<<" ";
      File_Out<<endl;
    }//Close for loop on Heidi_Klum.k_ij().Numerators().
  File_Out<<endl;

  //Check the GSO Coefficient matrix.
  if(!HM_Murdock.Consistent_GSO_Matrix())
    {
      File_Out<<"Inconsistent GSO Matrix."<<endl;
      cout<<"Inconsistent GSO Matrix."<<endl;
      return 0;
    }

  //Build the boson states.
  set<Alpha_Boson>::iterator itBoson_Sectors = Alpha_Bosons.begin();
  list<State> All_Boson_States;//States from all sectors.
  for(itBoson_Sectors; itBoson_Sectors != Alpha_Bosons.end(); ++itBoson_Sectors)
    {
      State_Builder Hulk_Hogan(*itBoson_Sectors, Fermion_Mode_Map, 
			       Common_Basis_Alphas, Heidi_Klum.k_ij());
      Hulk_Hogan.Build_States();

      //Display the sector.
      cout<<"Sector: "<<endl;
      itBoson_Sectors->Display();
      cout<<"Coefficients producing this sector: "<<endl;
      itBoson_Sectors->Display_Coefficients();
      cout<<endl;

      //Write the sector.
      File_Out<<"Sector (first number is denominator): "<<endl;
      for(int a=0; a<itBoson_Sectors->LM_Size(); a++)
	File_Out<<(itBoson_Sectors->Numerator()).at(a)<<" ";
      File_Out<<"|| ";
      for(int a=itBoson_Sectors->LM_Size(); 
	  a<(itBoson_Sectors->Numerator()).size(); a++)
	File_Out<<(itBoson_Sectors->Numerator()).at(a)<<" ";
      File_Out<<endl;
      File_Out<<"Coefficients producing this sector: "<<endl;
      for(int a=0; a<(itBoson_Sectors->Coefficients()).size(); a++)
	File_Out<<(itBoson_Sectors->Coefficients()).at(a)<<" ";
      File_Out<<endl<<endl;

      list<State> Boson_States = Hulk_Hogan.States();//States from this sector.
      list<State>::iterator itBoson_States = Boson_States.begin();

      //Write and display the states from this sector.
      File_Out<<"States from this sector surviving GSOP "<<
	"(first number is denominator): "<<Boson_States.size()<<endl;
      cout<<"States from this sector: "<<Boson_States.size()<<endl;
      for(itBoson_States; itBoson_States != Boson_States.end(); ++itBoson_States)
	{
	  itBoson_States->Display();

	  File_Out<<itBoson_States->Denominator()<<": ";
	  for(int a=0; a<itBoson_States->LM_Size(); a++)
	    File_Out<<(itBoson_States->Numerator()).at(a)<<" ";
	  File_Out<<"|| ";
	  for(int a=itBoson_States->LM_Size(); 
	      a<(itBoson_States->Numerator()).size(); a++)
	    File_Out<<(itBoson_States->Numerator()).at(a)<<" ";
	  File_Out<<endl;
	}//Close for loop on Boson_States.
      File_Out<<endl;
      cout<<endl;

      itBoson_States = Boson_States.begin();
      for(itBoson_States; itBoson_States != Boson_States.end(); ++itBoson_States)
	{
	  if(itBoson_States->Is_Positive(Fermion_Mode_Map))
	    {
	      State Boson_State = *itBoson_States;
	      Boson_State.Calculate_Length_Squared(Fermion_Mode_Map);
	      All_Boson_States.push_back(Boson_State);
	    }//Close if statement on positive states.
	}//Close for loop on Boson_States.
    }//Close for loop on Boson_Sectors.

  //Assemble the gauge groups.
  while(All_Boson_States.size() != 0)
    {
      list<State>::iterator itAll_Boson_States = All_Boson_States.begin();
      list<State> Positive_Roots;
      Positive_Roots.push_back(*itAll_Boson_States);
      All_Boson_States.erase(itAll_Boson_States);

      list<State>::iterator itPositive_Roots = Positive_Roots.begin();
      for(itPositive_Roots; itPositive_Roots != Positive_Roots.end(); 
	  ++itPositive_Roots)
	{
	  for(itAll_Boson_States = All_Boson_States.begin();
	      itAll_Boson_States != All_Boson_States.end(); )
	    {
	      if(Gauge_Dot(*itPositive_Roots, *itAll_Boson_States, 
			   Fermion_Mode_Map) != 0)
		{
		  Positive_Roots.push_back(*itAll_Boson_States);
		  itAll_Boson_States = All_Boson_States.erase(itAll_Boson_States);
		}else
		++itAll_Boson_States;
	    }//Close for loop on All_Boson_States.
	}//Close for loop on Positive_Roots.
      //Identify the gauge group.
      Gauge_Group_Identifier Dangle(Positive_Roots, Fermion_Mode_Map);
      Heidi_Klum.Add_Gauge_Group(Dangle.Get_Group());

			//DEBUG.
			Dangle.Get_Group().Display();
			Dangle.Get_Group().Display_Positive_Roots();
			Dangle.Display_Simple_Roots();
			//END DEBUG.

    }//Close while loop on All_Boson_States.size().
  Heidi_Klum.Sort_Gauge_Groups();

  //Build the SUSY states.
   itAlpha_SUSYs = Alpha_SUSYs.begin();
  list<State> All_SUSY_States;
  for(itAlpha_SUSYs; itAlpha_SUSYs != Alpha_SUSYs.end(); ++itAlpha_SUSYs)
    {
      State_Builder Hulk_Hogan(*itAlpha_SUSYs, Fermion_Mode_Map, Common_Basis_Alphas,
			       Heidi_Klum.k_ij());
      Hulk_Hogan.Build_States();
      
      //Display the sector.
      cout<<"SUSY sector: "<<endl;
      itAlpha_SUSYs->Display();
      cout<<"Coefficients producing this sector: "<<endl;
      itAlpha_SUSYs->Display_Coefficients();
      cout<<endl;
   

      //Write the sector.
      File_Out<<"SUSY sector (first number is denominator): "<<endl;
      for(int a=0; a<itAlpha_SUSYs->LM_Size(); a++)
	File_Out<<(itAlpha_SUSYs->Numerator()).at(a)<<" ";
      File_Out<<"|| ";
      for(int a=itAlpha_SUSYs->LM_Size();
	  a<(itAlpha_SUSYs->Numerator()).size(); a++)
	File_Out<<(itAlpha_SUSYs->Numerator()).at(a)<<" ";
      File_Out<<endl;
      File_Out<<"Coefficients producing this sector: "<<endl;
      for(int a=0; a<(itAlpha_SUSYs->Coefficients()).size(); a++)
	File_Out<<(itAlpha_SUSYs->Coefficients()).at(a)<<" ";
      File_Out<<endl<<endl;

      list<State> SUSY_States = Hulk_Hogan.States();
      list<State>::iterator itSUSY_States = SUSY_States.begin();

      //Write and display the SUSY states.
      File_Out<<"States from this sector surviving GSOP "<<
	"(first number is denominator): "<<SUSY_States.size()<<endl;
      cout<<"States from this sector: "<<SUSY_States.size()<<endl;
      for(itSUSY_States; itSUSY_States != SUSY_States.end(); ++itSUSY_States)
	{
	  itSUSY_States->Display();

	  for(int a=0; a<itSUSY_States->LM_Size(); a++)
	    File_Out<<(itSUSY_States->Numerator()).at(a)<<" ";
	  File_Out<<"|| ";
	  for(int a=itSUSY_States->LM_Size(); a<(itSUSY_States->Numerator()).size();
	      a++)
	    File_Out<<(itSUSY_States->Numerator()).at(a)<<" ";
	  File_Out<<endl;
	}//Close for loop on SUSY_States.
      File_Out<<endl;
      cout<<endl;

      itSUSY_States = SUSY_States.begin();
      for(itSUSY_States; itSUSY_States != SUSY_States.end(); ++itSUSY_States)
	All_SUSY_States.push_back(*itSUSY_States);
    }//Close for loop on SUSY sectors.
  Heidi_Klum.Set_SUSY_States(All_SUSY_States);

  //Build the fermion states.
  list<State> All_Matter_States;
  itAlpha_Fermions = Alpha_Fermions.begin();
  for(itAlpha_Fermions; itAlpha_Fermions != Alpha_Fermions.end(); ++itAlpha_Fermions)
    {
      State_Builder Hulk_Hogan(*itAlpha_Fermions, Fermion_Mode_Map, 
			       Common_Basis_Alphas, Heidi_Klum.k_ij());
      Hulk_Hogan.Build_States();

      //Display the sector.
      cout<<"Matter sector: "<<endl;
      itAlpha_Fermions->Display();
      cout<<"Coefficients producing this sector: "<<endl;
      itAlpha_Fermions->Display_Coefficients();
      cout<<endl;

      //Write the sector.
      File_Out<<"Matter sector (first number is the denominator): "<<endl;
      File_Out<<itAlpha_Fermions->Denominator()<<": ";
      for(int a=0; a<itAlpha_Fermions->LM_Size(); a++)
	File_Out<<(itAlpha_Fermions->Numerator()).at(a)<<" ";
      File_Out<<"|| ";
      for(int a=itAlpha_Fermions->LM_Size();
	  a<(itAlpha_Fermions->Numerator()).size(); a++)
	File_Out<<(itAlpha_Fermions->Numerator()).at(a)<<" ";
      File_Out<<endl;
      File_Out<<"Coefficients producing this sector: "<<endl;
      for(int a=0; a<(itAlpha_Fermions->Coefficients()).size(); a++)
	File_Out<<(itAlpha_Fermions->Coefficients()).at(a)<<" ";
      File_Out<<endl<<endl;
      
      list<State> Matter_States = Hulk_Hogan.States();
      list<State>::iterator itMatter_States = Matter_States.begin();

      //Write and display the states.
      File_Out<<"States from this sector surviving GSOP "<<
	"(first number is denominator): "<<Matter_States.size()<<endl;
      cout<<"States from this sector."<<endl;
      for(itMatter_States; itMatter_States != Matter_States.end(); ++itMatter_States)
	{
	  itMatter_States->Display();
	  File_Out<<itMatter_States->Denominator()<<": ";
	  for(int a=0; a<itMatter_States->LM_Size(); a++)
	    File_Out<<(itMatter_States->Numerator()).at(a)<<" ";
	  File_Out<<"|| ";
	  for(int a=itMatter_States->LM_Size();
	      a<(itMatter_States->Numerator()).size(); a++)
	    File_Out<<(itMatter_States->Numerator()).at(a)<<" ";
	  File_Out<<endl;

	}//Close for loop on Matter_States.
      File_Out<<endl;
      cout<<endl;

			itMatter_States = Matter_States.begin();
			for(itMatter_States; itMatter_States != Matter_States.end(); 
					++itMatter_States)
			{
				if(!Is_SUSY_Partner(*itMatter_States, All_SUSY_States))
				{
					vector<Group_Representation> Representation_Dimensions;
					for(int a=0; a<Heidi_Klum.Gauge_Groups().size(); a++)
					{
						Group_Representation Representation = 
							Heidi_Klum.rGauge_Groups().at(a).
							Compute_Rep_Dimension(*itMatter_States, Fermion_Mode_Map);

						if(Representation.Dimension() == 0)
							break;

						Representation_Dimensions.push_back(Representation);
					}//Close for loop on gauge groups.
					if(Representation_Dimensions.size() == 
							Heidi_Klum.Gauge_Groups().size())
						Heidi_Klum.Add_Matter_State(Matter_State
								(itMatter_States->Numerator(),
								 itMatter_States->Denominator(),
								 itMatter_States->LM_Size(),
								 Representation_Dimensions));
				}//Close if statement on SUSY partner.
			}//Close for loop on Matter_States.
		}//Close for loop on Alpha_Fermions.
	Heidi_Klum.Build_Matter_Representations();
	if(Heidi_Klum.Matter_Representations().size()!= 0)
		Heidi_Klum.Set_Matter_Rep_Sign_Convention();

	//Display the gauge groups and the positive roots.
	cout<<"Gauge groups: "<<Heidi_Klum.Gauge_Groups().size()<<endl;
	Heidi_Klum.Display_Gauge_Groups();
	cout<<endl;

  for(int a=0; a<Heidi_Klum.Gauge_Groups().size(); a++)
    {
      cout<<"Gauge group is: "<<endl;
      Heidi_Klum.Gauge_Groups().at(a).Display();
      Heidi_Klum.Gauge_Groups().at(a).Display_Positive_Roots();
      cout<<endl;
    }//Close for loop on displaying the gauge groups.
  cout<<endl;

  //Write the gauge groups and the positive roots.
  File_Out<<"Gauge groups: "<<Heidi_Klum.Gauge_Groups().size()<<endl;
  for(int a=0; a<Heidi_Klum.Gauge_Groups().size(); a++)
    {
      Gauge_Group_Name Bill = Heidi_Klum.Gauge_Groups().at(a).Name();
      File_Out<<Bill.Class()<<" "<<Bill.Rank()<<" "<<Bill.KM_Level()<<" ";
    }//Close for loop on gauge groups.
  File_Out<<endl<<endl;

  for(int a=0; a<Heidi_Klum.Gauge_Groups().size(); a++)
    {
      Gauge_Group_Name Frank = Heidi_Klum.Gauge_Groups().at(a).Name();
      File_Out<<"Gauge group is: "<<Frank.Class()<<" "<<Frank.Rank()<<" "
	      <<Frank.KM_Level()<<endl;
      //Display the nonzero positive roots.
      list<State> Positive_Roots = Heidi_Klum.Gauge_Groups().at(a).Positive_Roots();
      list<State>::iterator itPositive_Roots = Positive_Roots.begin();
      File_Out<<"Nonzero positive roots for this gauge group are :"<<
	Positive_Roots.size()<<endl;
      for(itPositive_Roots; itPositive_Roots != Positive_Roots.end(); 
	  ++itPositive_Roots)
	{
	  File_Out<<itPositive_Roots->Denominator()<<": ";
	  for(int a=0; a<itPositive_Roots->LM_Size(); a++)
	    File_Out<<(itPositive_Roots->Numerator()).at(a)<<" ";
	  File_Out<<"|| ";
	  for(int a=itPositive_Roots->LM_Size();
	      a<(itPositive_Roots->Numerator()).size(); a++)
	    File_Out<<(itPositive_Roots->Numerator()).at(a)<<" ";
	  File_Out<<endl;
	}//Close for loop on Positive_Roots.
      File_Out<<endl;

      //Display the simple roots.
      list<State> Simple_Roots = Heidi_Klum.Gauge_Groups().at(a).Simple_Roots();
      list<State>::iterator itSimple_Roots = Simple_Roots.begin();
      File_Out<<"Simple roots: "<<Simple_Roots.size()<<endl;
      for(itSimple_Roots; itSimple_Roots != Simple_Roots.end(); ++itSimple_Roots)
	{
	  File_Out<<itSimple_Roots->Denominator()<<": ";
	  for(int a=0; a<itSimple_Roots->LM_Size(); a++)
	    File_Out<<(itSimple_Roots->Numerator()).at(a)<<" ";
	  File_Out<<"|| ";
	  for(int a=(itSimple_Roots->LM_Size()); 
	      a<(itSimple_Roots->Numerator()).size();
	      a++)
	    File_Out<<(itSimple_Roots->Numerator()).at(a)<<" ";
	  File_Out<<endl;
	}//Close for loop on Simple_Roots.
      File_Out<<endl;

			//Display the Dynkin map.
			File_Out<<"Dynkin labels: "<<endl;
			map<vector<int>, Group_Representation> Dynkin_Labels = 
				Heidi_Klum.Gauge_Groups().at(a).Dynkin_Labels();
			map<vector<int>, Group_Representation>::iterator itDynkin_Labels = 
				Dynkin_Labels.begin();
			for(itDynkin_Labels; itDynkin_Labels != 
					Dynkin_Labels.end(); ++itDynkin_Labels)
			{
				vector<int> Dynkin_Label = itDynkin_Labels->first;
				for(int a=0; a<Dynkin_Label.size(); a++)
					File_Out<<Dynkin_Label.at(a)<<" ";
				File_Out<<"-> "<<itDynkin_Labels->second.Dimension();
				File_Out<<itDynkin_Labels->second.Triality()<<endl;
			}//Close for loop on Dynkin_Labels.
			File_Out<<endl;
		}//Close for loop on gauge groups.
	File_Out<<endl;
      
  //Display the highest weight states and their charges under the gauge groups.
  cout<<"Gauge groups: "<<endl;
  Heidi_Klum.Display_Gauge_Groups();
  cout<<endl;

  list<Matter_State> HW_States = Heidi_Klum.Matter_States();
  list<Matter_State>::const_iterator itHW_States = HW_States.begin();
  cout<<"Highest weight matter states: "<<HW_States.size()<<endl;
  for(itHW_States; itHW_States != HW_States.end(); ++itHW_States)
    {
      itHW_States->Display();
      itHW_States->Display_Representations();
    }//Close for loop on highest weight states.

  //Write the highest weight states and their charges under the gauge groups.
	File_Out<<"Gauge groups: "<<endl;
	for(int a=0; a<Heidi_Klum.Gauge_Groups().size(); a++)
	{
		Gauge_Group_Name Dick = Heidi_Klum.Gauge_Groups().at(a).Name();
		File_Out<<Dick.Class()<<" "<<Dick.Rank()<<" "<<Dick.KM_Level()<<" ";
	}//Close for loop on gauge groups.
	File_Out<<endl;

	itHW_States = HW_States.begin();
	File_Out<<"Highest weight matter states (first number is denominator): "
		<<HW_States.size()<<endl;
	for(itHW_States; itHW_States != HW_States.end(); ++itHW_States)
	{
		File_Out<<itHW_States->Denominator()<<": "<<endl;
		for(int a=0; a<itHW_States->LM_Size(); a++)
			File_Out<<(itHW_States->Numerator()).at(a)<<" ";
		File_Out<<"|| ";
		for(int a=itHW_States->LM_Size();
				a<(itHW_States->Numerator()).size(); a++)
			File_Out<<(itHW_States->Numerator()).at(a)<<" ";
		File_Out<<endl;

		for(int a=0; a<(itHW_States->Representations()).size(); a++)
		{
			File_Out<<(itHW_States->Representations()).at(a).Dimension();
			File_Out<<(itHW_States->Representations()).at(a).Triality()<<" ";
		}
		File_Out<<endl<<endl;
	}//Close for loop on highest weight states.
	File_Out<<endl;

	//Compute the U(1)'s.
	int Total_Rank = Heidi_Klum.BV_Set().at(0).BV().size()/4 + 6;
	Total_Rank -= Find_Rank_Cuts(Heidi_Klum,Fermion_Mode_Map);
	int Total_NA_Rank = 0;
	for(int a=0; a<Heidi_Klum.Gauge_Groups().size(); a++)
		Total_NA_Rank += Heidi_Klum.Gauge_Groups().at(a).Name().Rank();
	Heidi_Klum.Set_U1_Factors(Total_Rank - Total_NA_Rank);


  cout<<"Model summary."<<endl;
  Heidi_Klum.Display_Particle_Content();  

  File_Out<<"Model summary: "<<endl;
  Stephen_King.Write_Model_Particle_Content(File_Out, Heidi_Klum);

  File_Out<<"Program finished normally."<<endl;
  File_Out.close();
  cout<<"End program."<<endl;
  return 0;
}//Close main.

vector<int> Make_Vector(int BV_Array[], int BV_Size)
{
  vector<int> BV_Vector;
  for(int a=0; a<BV_Size; a++)
    BV_Vector.push_back(BV_Array[a]);
  return BV_Vector;
}//Close Make_Vector.

int Gauge_Dot(const State& State1, const State& State2, 
	      const map<int, int>& Fermion_Mode_Map)
{
  int Dot = 0;
  map<int, int>::const_iterator itFermion_Mode_Map = Fermion_Mode_Map.
    find(State1.LM_Size());
  for(itFermion_Mode_Map; itFermion_Mode_Map != Fermion_Mode_Map.end();
      ++itFermion_Mode_Map)
    Dot += (State1.Numerator().at(itFermion_Mode_Map->first) * 
	    State2.Numerator().at(itFermion_Mode_Map->first));
  return Dot;
}//Close Gauge_Dot.

bool Is_SUSY_Partner(const State& New_Matter_State, const list<State>& SUSY_States)
{
  if(SUSY_States.size() == 0)
    return false;

  list<State>::const_iterator itSUSY_States = SUSY_States.begin();
  for(itSUSY_States; itSUSY_States != SUSY_States.end(); ++itSUSY_States)
    {
      bool SUSY_State = true;
      for(int a=0; a<itSUSY_States->LM_Size(); a++)
	{
	  if(New_Matter_State.Numerator().at(a) != 
	     (itSUSY_States->Numerator()).at(a))
	    {
	      SUSY_State = false;
	      break;
	    }
	}//Close for loop on LM_Size.
      if(SUSY_State)
	return SUSY_State;
    }//Close for loop on SUSY_States.
  return false;
}//Close Is_SUSY_Partner.

int Find_Rank_Cuts(const Model& Heidi_Klum, const map<int, int>& Fermion_Mode_Map)
{
  int Rank_Cuts = 0;
  int LM_Size = Heidi_Klum.BV_Set().at(0).LM_Size();
  map<int, int>::const_iterator itFermion_Mode_Map = Fermion_Mode_Map.begin();
  map<int, int>::const_iterator itMap_End = Fermion_Mode_Map.end();
  for(itFermion_Mode_Map; itFermion_Mode_Map != itMap_End; ++itFermion_Mode_Map)
    {
      if((itFermion_Mode_Map->first)<LM_Size && 
	 (itFermion_Mode_Map->second) >=LM_Size)
	Rank_Cuts++;
    }
  return Rank_Cuts/2;
}//Close Find_Rank_Cuts.
