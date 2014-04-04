/*
Authors: Timothy Renner (renner.timothy@gmail.com)
				 Doug Moore (Douglas_Moore1@baylor.edu)
				 Gerald Cleaver (Gerald_Cleaver@baylor.edu)

Institute: Baylor University
Created: 5/18/2011

This program searches through a list of files
and prints those with three chiral matter
generations for a MSSM gauge group.
The representations which form the chiral
generations are specified explicitly.

 Uses the FF Framework.
 */

#include <vector>
#include <set>
#include <string>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "FF_LEEFT.hh"
#include "FF_Output_Writer.hh"
#include "FF_Gauge_Group_Name.hh"
#include "FF_LEEFT_Tools.hh"
#include "FF_Latex_Writer.hh"
#include "FF_Statistics.hh"
#include "FF_Statistical_Report_Generator.hh"

using namespace std;

bool Has_MSSM_Group(const LEEFT& LEEFT_Input);

LEEFT Find_All_MSSM_Observable_Generations(LEEFT LEEFT_Input);
Observable_Sector Find_MSSM_Observable_Sector_States
	(const LEEFT& LEEFT_Input, int SU3_Index, int SU2_Index);
bool Has_Hidden_Sector_Charge(const vector<int>& Representation, 
		int SU3_Index, int SU2_Index);


int main()
{
	cout<<"Begin program."<<endl;
	int Layers = 1;

	//Set up the input files.
	string FF_Data = "/Users/timothyrenner/FF_Framework/data/";
	vector<string> Input_File_Names;
	Input_File_Names.push_back(FF_Data + 
			"/NAHE_Extensions/One_Layer/NAHE_O3_L1.txt");

	//Set up the output files.
	ofstream ThreeG_Out("3G_MSSM_Models.txt");
	ofstream ThreeAG_Out("3AG_MSSM_Models.txt");
	ofstream ThreeG_ThreeAG_Out("3G_3AG_MSSM_Models.txt");
	ofstream Latex_Out("./tex/3G_MSSM_Stats.tex");
	Output_Writer Doug_Moore;
	Latex_Writer Jared_Greenwald;

	//Read in and analyze the LEEFTs.
	cout<<"Reading in LEEFTs with MSSM gauge groups."<<endl;
	set<LEEFT> Unique_3G_MSSM_LEEFTs;
	set<LEEFT> Unique_3AG_MSSM_LEEFTs;
	set<LEEFT> Unique_3G_3AG_MSSM_LEEFTs;
	for(int a=0; a<Input_File_Names.size(); ++a)
	{
		ifstream LEEFT_In(Input_File_Names.at(a).c_str());
		if(LEEFT_In.fail())
		{
			cout<<Input_File_Names.at(a)<<" did not open."<<endl;
			cout<<"Terminating program."<<endl;
			return 0;
		}//Close if statement on the input file failing.
	
		LEEFT_Tools::Read_Front_Info(LEEFT_In);
		while(LEEFT_In.peek() != 'T')
		{
			vector<string> LEEFT_Data = LEEFT_Tools::Read_LEEFT_Data(LEEFT_In);
			string blank;
			getline(LEEFT_In, blank);
			vector<string> Model_Inputs = 
				LEEFT_Tools::Read_Model_Inputs(LEEFT_In, Layers);
			getline(LEEFT_In, blank);
			
			LEEFT New_LEEFT(LEEFT_Data);
			if(Has_MSSM_Group(New_LEEFT))
			{
				New_LEEFT = Find_All_MSSM_Observable_Generations(New_LEEFT);
				bool ThreeG_Model_Written = false;
				bool ThreeAG_Model_Written = false;
				bool ThreeG_ThreeAG_Model_Written = false;
				for(int a=0; a<New_LEEFT.Observable_Sectors().size(); ++a)
				{
					int Chiral_Generations = New_LEEFT.
						Observable_Sectors().at(a).Count_Chiral_Generations();
					int Chiral_Anti_Generations = New_LEEFT.
						Observable_Sectors().at(a).Count_Chiral_Anti_Generations();

					if(!ThreeG_Model_Written && 
							(Chiral_Generations == 3))
					{
						Doug_Moore.Write_LEEFT_Particle_Content(ThreeG_Out, New_LEEFT);
						for(int b=0; b<Model_Inputs.size(); ++b)
							ThreeG_Out<<Model_Inputs.at(b)<<endl;
						ThreeG_Out<<endl;
						Unique_3G_MSSM_LEEFTs.insert(New_LEEFT);
						ThreeG_Model_Written = true;
					}//Close if statement on three generations.

					if(!ThreeAG_Model_Written && 
							(Chiral_Anti_Generations == 6))
					{
						Doug_Moore.Write_LEEFT_Particle_Content
							(ThreeAG_Out, New_LEEFT);
						for(int b=0; b<Model_Inputs.size(); ++b)
							ThreeAG_Out<<Model_Inputs.at(b)<<endl;
						ThreeAG_Out<<endl;
						Unique_3AG_MSSM_LEEFTs.insert(New_LEEFT);
						ThreeAG_Model_Written = true;
					}//Close if statement on three anti-generations.

					if(!ThreeG_ThreeAG_Model_Written && 
							(Chiral_Generations == 3) && 
							(Chiral_Anti_Generations == 6))
					{
						Doug_Moore.Write_LEEFT_Particle_Content
							(ThreeG_ThreeAG_Out, New_LEEFT);
						for(int b=0; b<Model_Inputs.size(); ++b)
							ThreeG_ThreeAG_Out<<Model_Inputs.at(b)<<endl;
						ThreeG_ThreeAG_Out<<endl;
						Unique_3G_3AG_MSSM_LEEFTs.insert(New_LEEFT);
						ThreeG_ThreeAG_Model_Written = true;
					}//Close if statement on three generations and anti generations.
					if(ThreeG_Model_Written && 
							ThreeAG_Model_Written && 
							ThreeG_ThreeAG_Model_Written)
						break;
				}//Close for loop on Observable_Sectors.
			}//Close if statement on Has_MSSM_Group.
		}//Close while loop on LEEFT_In.
	}//Close for loop on Input_File_Names.
	Statistical_Report_Generator ThreeG_Stat_Reporter
		(Unique_3G_MSSM_LEEFTs);

	string Doc_Title = 
		"Statistics for Three Generation MSSM";
	Doc_Title += " NAHE Based Models, Order 3, Layer 1";
	Jared_Greenwald.Begin_Latex_Document(Latex_Out, Doc_Title);
	//Three Generation models.
	Latex_Out<<"\\textbf{Unique Three Generation Models: ";
	Latex_Out<<Unique_3G_MSSM_LEEFTs.size()<<"}\\\\"<<endl;
	Latex_Out<<"\\vspace{3 mm}"<<endl;
	ThreeG_Stat_Reporter.Generate_Statistical_Report(Latex_Out, 4);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;

	//Three Anti_Generation models.
	Statistical_Report_Generator ThreeAG_Stat_Reporter
		(Unique_3AG_MSSM_LEEFTs);
	Latex_Out<<"\\textbf{Unique Three Anti-Generation Models: ";
	Latex_Out<<Unique_3AG_MSSM_LEEFTs.size()<<"}\\\\"<<endl;
	Latex_Out<<"\\vspace{3 mm}"<<endl;
	ThreeAG_Stat_Reporter.Generate_Statistical_Report(Latex_Out, 4);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;

	//Three Generation and Anti-Generation models.
	Statistical_Report_Generator ThreeG_ThreeAG_Stat_Reporter
		(Unique_3G_3AG_MSSM_LEEFTs);
	Latex_Out<<"\\textbf{Unique Three Generation and ";
	Latex_Out<<"Anti-Generation Models: ";
	Latex_Out<<Unique_3G_3AG_MSSM_LEEFTs.size()<<"}\\\\"<<endl;
	Latex_Out<<"\\vspace{3 mm}"<<endl;
	ThreeG_ThreeAG_Stat_Reporter.Generate_Statistical_Report(Latex_Out, 4);

	//Finish the Latex document.
	Jared_Greenwald.Finish_Latex_Document(Latex_Out);
	cout<<"Program finished."<<endl;
	return 0;
}//Close main.

bool Has_MSSM_Group(const LEEFT& LEEFT_Input)
{
	bool Has_SU3 = false;
	bool Has_SU2 = false;
	
	if(LEEFT_Input.U1_Factors() == 0)
		return false;

	for(int a=0; a<LEEFT_Input.Gauge_Groups().size(); ++a)
	{
		if(LEEFT_Input.Gauge_Groups().at(a).Class() == 'A')
		{
			if(LEEFT_Input.Gauge_Groups().at(a).Rank() == 2)
				Has_SU3 = true;
			else if(LEEFT_Input.Gauge_Groups().at(a).Rank() == 1)
				Has_SU2 = true;

			if(Has_SU3 && Has_SU2)
				return true;
		}//Close if statement on A class gauge groups.
	}//Close for loop on Gauge_Groups.

	return false;
}//Close Has_MSSM_Group.

LEEFT Find_All_MSSM_Observable_Generations(LEEFT LEEFT_Input)
{
	//Get the indices of the left-right symmetric gauge groups.
	vector<int> SU3_Indices;
	vector<int> SU2_Indices;
	for(int a=0; a<LEEFT_Input.Gauge_Groups().size(); ++a)
	{
		if(LEEFT_Input.Gauge_Groups().at(a).Class() == 'A')
		{
			if(LEEFT_Input.Gauge_Groups().at(a).Rank() == 2)
				SU3_Indices.push_back(a);
			else if(LEEFT_Input.Gauge_Groups().at(a).Rank() == 1)
				SU2_Indices.push_back(a);
		}//Close if statement on SU - type gauge group.
	}//Close for loop over the gauge groups.

	//Now find the observable sector generations for each observable sector.
	vector<Observable_Sector> MSSM_Observable_Sectors;
	for(int a=0; a<SU3_Indices.size(); ++a)
	{
		for(int b=0; b<SU2_Indices.size(); ++b)
		{
			int SU3_Index = SU3_Indices.at(a);
			int SU2_Index = SU2_Indices.at(b);
			Observable_Sector New_MSSM_Observable_Sector = 
				Find_MSSM_Observable_Sector_States(LEEFT_Input, SU3_Index, 
						SU2_Index);
			//Make sure it isn't empty.
			if(New_MSSM_Observable_Sector.Generations().size() != 0)
				MSSM_Observable_Sectors.push_back(New_MSSM_Observable_Sector);
		}//Close for loop over SU2_Indices.
	}//Close for loop over SU3_Indices.

	LEEFT_Input.Set_Observable_Sectors(MSSM_Observable_Sectors);
	return LEEFT_Input;
}//Close Find_All_MSSM_Observable_Generations.

Observable_Sector Find_MSSM_Observable_Sector_States(const LEEFT&
		LEEFT_Input, int SU3_Index, int SU2_Index)
{
	//First, find the generations corresponding to the gauge groups.
	vector<vector<int> > Generations;
	vector<int> Generation_Qty;
	vector<vector<int> > Barred_Generations;
	vector<int> Barred_Generation_Qty;
	vector<vector<int> > Anti_Generations;
	vector<int> Anti_Generation_Qty;
	vector<vector<int> > Anti_Barred_Generations;
	vector<int> Anti_Barred_Generation_Qty;

	for(int a=0; a<LEEFT_Input.Representations().size(); ++a)
	{
		if(!Has_Hidden_Sector_Charge(LEEFT_Input.Representations().at(a),
					SU3_Index, SU2_Index))
		{
			if((LEEFT_Input.Representations().at(a).at(SU3_Index) == 3) &&
					(LEEFT_Input.Representations().at(a).at(SU2_Index) == 2))
			{
				Generations.push_back(LEEFT_Input.Representations().at(a));
				Generation_Qty.push_back
					(LEEFT_Input.Representation_Count().at(a));
			}else if((LEEFT_Input.Representations().at(a).at(SU3_Index) == -3)
					&&(LEEFT_Input.Representations().at(a).at(SU2_Index) == 2))
			{
				Barred_Generations.push_back(LEEFT_Input.Representations().at(a));
				Barred_Generation_Qty.push_back
					(LEEFT_Input.Representation_Count().at(a));
			}else if((LEEFT_Input.Representations().at(a).at(SU3_Index) == 3)
					&& (LEEFT_Input.Representations().at(a).at(SU2_Index) == 1))
			{
				Anti_Generations.push_back(LEEFT_Input.Representations().at(a));
				Anti_Generation_Qty.push_back
					(LEEFT_Input.Representation_Count().at(a));
			}else if((LEEFT_Input.Representations().at(a).at(SU3_Index) == -3)
					&& (LEEFT_Input.Representations().at(a).at(SU2_Index) == 1))
			{
				Anti_Barred_Generations.push_back
					(LEEFT_Input.Representations().at(a));
				Anti_Barred_Generation_Qty.push_back
					(LEEFT_Input.Representation_Count().at(a));
			}//Close if/else chain for finding the types of generations.
		}//Close if statement on hidden sector charge.
	}//Close for loop on Representations.

	//Now fill in the blanks for any missing representations.
	vector<int> Generation_Rep
		(LEEFT_Input.Representations().front().size(), 1);
	if(Generations.size() == 0)
	{
		Generation_Rep.at(SU3_Index) = 3;
		Generation_Rep.at(SU2_Index) = 2;
		Generations.push_back(Generation_Rep);
		Generation_Qty.push_back(0);
	}
	if(Barred_Generations.size() == 0)
	{
		Generation_Rep.at(SU3_Index) = -3;
		Generation_Rep.at(SU2_Index) = 2;
		Barred_Generations.push_back(Generation_Rep);
		Barred_Generation_Qty.push_back(0);
	}
	if(Anti_Generations.size() == 0)
	{
		Generation_Rep.at(SU3_Index) = 3;
		Generation_Rep.at(SU2_Index) = 1;
		Anti_Generations.push_back(Generation_Rep);
		Anti_Generation_Qty.push_back(0);
	}
	if(Anti_Barred_Generations.size() == 0)
	{
		Generation_Rep.at(SU3_Index) = -3;
		Generation_Rep.at(SU2_Index) = 1;
		Anti_Barred_Generations.push_back(Generation_Rep);
		Anti_Barred_Generation_Qty.push_back(0);
	}

	//Create a new Observable_Sector object and return it.
	Observable_Sector New_Observable_Sector;
	New_Observable_Sector.Set_Generations(Generations);
	New_Observable_Sector.Set_Generation_Qty(Generation_Qty);
	New_Observable_Sector.Set_Barred_Generations(Barred_Generations);
	New_Observable_Sector.Set_Barred_Generation_Qty(Barred_Generation_Qty);
	New_Observable_Sector.Set_Anti_Generations(Anti_Generations);
	New_Observable_Sector.Set_Anti_Generation_Qty(Anti_Generation_Qty);
	New_Observable_Sector.Set_Anti_Barred_Generations
		(Anti_Barred_Generations);
	New_Observable_Sector.Set_Anti_Barred_Generation_Qty
		(Anti_Barred_Generation_Qty);

	return New_Observable_Sector;
}//Close Find_MSSM_Observable_Sector_States.

bool Has_Hidden_Sector_Charge(const vector<int>& Representation,
		int SU3_Index, int SU2_Index)
{
  for(int a=0; a<Representation.size(); ++a)
	{
		if((a != SU3_Index) && (a != SU2_Index) &&
				(Representation.at(a) != 1))
			return true;
	}//Close for loop on Representation.
	return false;
}//Close Has_Hidden_Sector_Charge.
