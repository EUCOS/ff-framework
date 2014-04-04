/* 
Authors: Timothy Renner (renner.timothy@gmail.com)
				 Doug Moore (Douglas_Moore1@baylor.edu)
				 Gerald Cleaver, Ph.D. (Gerald_Cleaver@baylor.edu)

Institute: Baylor University
Created: 6/3/2011

This program searches any data set which may contain redundancies due to
improper ordering of the gauge groups.
It rebuilds the models it reads in and outputs them to a file.

Uses the FF Framework

*/

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "FF_LEEFT.hh"
#include "FF_Output_Writer.hh"
#include "FF_LEEFT_Tools.hh"
#include "FF_Model_Builder.hh"
#include "FF_Basis_Vector.hh"

using namespace std;

Basis_Vector Make_Basis_Vector(string BV, int Order, 
		int Large_ST_Dimensions);
vector<int> Make_Vector(string The_String);

int main()
{
	cout<<"Begin program."<<endl;
	int Layers = 1; //For LEEFT Tools functions.
	//Set up the initial conditions.
	vector<int> Orders;
	Orders.push_back(2);
	int Large_ST_Dimensions = 4;

	Model_Builder Initial_Conditions(Large_ST_Dimensions);
	//Initial_Conditions.Load_S_Vector(Large_ST_Dimensions);
	Initial_Conditions.Load_NAHE_Set();
	//Initial_Conditions.Load_NAHE_Variation();

	//Set up the input files.
	//string FF_Data = 
		//"/Users/timothyrenner/FF_Framework/fail/";
	string FF_Data = "/Users/timothyrenner/FF_Framework/data/";
	vector<string> Input_File_Names;
	Input_File_Names.push_back(FF_Data + 
			"NAHE_Extensions/One_Layer/NAHE_NSV_O2_L1_NSP.txt");

	//Set up the output file.
	ofstream out("New_Nonredundant_Data_Set.txt");
	Output_Writer Doug_Moore;

	cout<<"Reading in old basis vectors, rebuilding models."<<endl;
	set<LEEFT> Unique_LEEFTs;

	for(int a=0; a<Input_File_Names.size(); ++a)
	{
		ifstream Model_In(Input_File_Names.at(a).c_str());
		if(Model_In.fail())
		{
			cout<<Input_File_Names.at(a)<<" did not open.";
			cout<<"Terminating program."<<endl;
			return 0;
		}//Close if statement on the input file failing.

		LEEFT_Tools::Read_Front_Info(Model_In);
		while(Model_In.peek() != 'T')
		{
			vector<string> LEEFT_Data = LEEFT_Tools::Read_LEEFT_Data(Model_In);
			string blank;
			getline(Model_In, blank);
			vector<string> Model_Inputs = 
				LEEFT_Tools::Read_Model_Inputs(Model_In, Layers);

			getline(Model_In,blank);

			vector<Basis_Vector> BV_Set;
			for(int b=0; b<Model_Inputs.size()/2; ++b)
				BV_Set.push_back(Make_Basis_Vector(Model_Inputs.at(b),
							Orders.at(b), Large_ST_Dimensions));
			Model_Builder The_Builder(Large_ST_Dimensions);
			for(int b=1; b<Initial_Conditions.FFHS_Model().BV_Set().size(); 
					++b)
				The_Builder.rFFHS_Model().Load_BV_Set(Initial_Conditions.
						FFHS_Model().BV_Set().at(b));
			for(int b=0; b<BV_Set.size(); ++b)
				The_Builder.rFFHS_Model().Load_BV_Set(BV_Set.at(b));
			for(int b=0; b<Initial_Conditions.FFHS_Model().k_ij().
					Numerators().size(); ++b)
				The_Builder.rFFHS_Model().
					Load_k_ij_Row(Initial_Conditions.FFHS_Model().k_ij().
						Numerators().at(b));
			for(int b=Model_Inputs.size()/2; b<Model_Inputs.size(); ++b)
			{
				vector<int> k_ij_Row = Make_Vector(Model_Inputs.at(b));
				The_Builder.rFFHS_Model().
					Load_k_ij_Row(k_ij_Row);
			}

			
			The_Builder.Build_Model();
			int Unique_Model_Size = Unique_LEEFTs.size();
			LEEFT New_LEEFT(The_Builder.FFHS_Model());
			Unique_LEEFTs.insert(New_LEEFT);
			if(Unique_Model_Size < Unique_LEEFTs.size())
			{
				Doug_Moore.Write_Model_Particle_Content(out, 
						The_Builder.FFHS_Model()); 
				Doug_Moore.Write_Model_Extension_BVs(out, 
						The_Builder.FFHS_Model(), Layers);
				for(int b=Model_Inputs.size()/2; b<Model_Inputs.size(); ++b)
					out<<Model_Inputs.at(b)<<endl;

				out<<endl;
			}//CLose if statement on unqiue models.
		}//CLose while loop on the file.
		
	}//Close for loop on Input_File_Names.
	out<<"Total unique models: "<<Unique_LEEFTs.size()<<endl;	

	ofstream LEEFT_Out("Unique_Ordered_LEEFTs.txt");
	set<LEEFT>::iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
	for(itUnique_LEEFTs; itUnique_LEEFTs != Unique_LEEFTs.end(); 
			++itUnique_LEEFTs)
		Doug_Moore.Write_LEEFT_Particle_Content(LEEFT_Out, *itUnique_LEEFTs);
	LEEFT_Out<<"Total unique LEEFTs: "<<Unique_LEEFTs.size();

	cout<<"End program. Rename new data set and move old one."<<endl;
	return 0;
}//Close main.
	
Basis_Vector Make_Basis_Vector(string BV, int Order, 
		int Large_ST_Dimensions)
{
	vector<int> BV_Elements;
	stringstream ss(BV);
	int BV_Size = 28 - 2*Large_ST_Dimensions + //Left mover.
		16 + //psi bar and eta bar
		2*(10 - Large_ST_Dimensions) + //ybar and wbar
		16; //phi-bar.
	for(int a=0; a<BV_Size; ++a)
	{
		int BV_Element = -1;
		ss>>BV_Element;
		BV_Elements.push_back(BV_Element);
	}//Close for loop on creating the basis vector from the stringstream.
	return Basis_Vector(BV_Elements, Order, Large_ST_Dimensions);
}
vector<int> Make_Vector(string The_String)
{
	stringstream ss(The_String);
	vector<int> New_Vector;
	//Divide by two to account for the spaces in the string.
	for(int a=0; a<The_String.size()/2; ++a)
	{
		int k_ij_element;
		ss>>k_ij_element;
		New_Vector.push_back(k_ij_element);
	}
	return New_Vector;
}//CLose Make_Vector.
