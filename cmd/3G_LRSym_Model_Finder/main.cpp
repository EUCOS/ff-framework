/*
Authors: Timothy Renner (renner.timothy@gmail.com)
				 Doug Moore (Douglas_Moore1@baylor.edu)
				 Gerald Cleaver (Gerald_Cleaver@baylor.edu)

Institute: Baylor University
Created: 5/18/2011

This program searches through a list of files
and prints those with three chiral matter
generations for a left-right symmetric gauge group.
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

bool Has_LRSym_Group(const LEEFT& LEEFT_Input);

LEEFT Find_All_LRSym_Observable_Generations(LEEFT LEEFT_Input);
Observable_Sector Find_LRSym_Observable_Sector_States
	(const LEEFT& LEEFT_Input, int SU3_Index, int SU2A_Index,
	 int SU2B_Index);
bool Has_Hidden_Sector_Charge(const vector<int>& Representation, 
		int SU3_Index, int SU2A_Index, int SU2B_Index);


int Count_OS_Charged_Exotics(const LEEFT& LEEFT_Input,
		int SU3_Index, int SU2A_Index, int SU2B_Index);
map<int, int> Find_OS_Charged_Exotics(const set<LEEFT>& Unique_LEEFTs);

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
	ofstream ThreeG_Out("3G_LRSym_Models.txt");
	ofstream ThreeG_ThreeAG_Out("3G_3AG_LRSym_Models.txt");
	ofstream Latex_Out("./tex/3G_LRSym_Stats.tex");
	Output_Writer Doug_Moore;
	Latex_Writer Jared_Greenwald;

	//Read in and analyze the LEEFTs.
	cout<<"Reading in LEEFTs with left-right symmetric gauge groups."<<endl;
	set<LEEFT> Unique_3G_LRSym_LEEFTs;
	set<vector<string> >Unique_3G_LRSym_BVs;
	set<LEEFT> Unique_3G_3AG_LRSym_LEEFTs;
	set<vector<string> > Unique_3G_3AG_LRSym_BVs;
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
			if(Has_LRSym_Group(New_LEEFT))
			{
				New_LEEFT = Find_All_LRSym_Observable_Generations(New_LEEFT);
				bool ThreeG_Model_Written = false;
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
						vector<string> Unique_BV_Loader;
						for(int a=0; a<Model_Inputs.size()/2; ++a)
						{
							Unique_BV_Loader.push_back(Model_Inputs.at(a));
						}
						Unique_3G_LRSym_BVs.insert(Unique_BV_Loader);
						Doug_Moore.Write_LEEFT_Particle_Content(ThreeG_Out, New_LEEFT);
						for(int b=0; b<Model_Inputs.size(); ++b)
							ThreeG_Out<<Model_Inputs.at(b)<<endl;
						ThreeG_Out<<endl;
						Unique_3G_LRSym_LEEFTs.insert(New_LEEFT);
						ThreeG_Model_Written = true;
					}//Close if statement on three generations.

					if(!ThreeG_ThreeAG_Model_Written && 
							(Chiral_Generations == 3) && 
							(Chiral_Anti_Generations == 3))
					{
						vector<string> Unique_BV_Loader;
						for(int a=0; a<Model_Inputs.size()/2; ++a)
						{
							Unique_BV_Loader.push_back(Model_Inputs.at(a));
						}
						Unique_3G_3AG_LRSym_BVs.insert(Unique_BV_Loader);
						Doug_Moore.Write_LEEFT_Particle_Content
							(ThreeG_ThreeAG_Out, New_LEEFT);
						for(int b=0; b<Model_Inputs.size(); ++b)
							ThreeG_ThreeAG_Out<<Model_Inputs.at(b)<<endl;
						ThreeG_ThreeAG_Out<<endl;
						Unique_3G_3AG_LRSym_LEEFTs.insert(New_LEEFT);
						ThreeG_ThreeAG_Model_Written = true;
					}//Close if statement on three generations and anti generations.
					if(ThreeG_Model_Written && 
							ThreeG_ThreeAG_Model_Written)
						break;
				}//Close for loop on Observable_Sectors.
			}//Close if statement on Has_LRSym_Group.
		}//Close while loop on LEEFT_In.
	}//Close for loop on Input_File_Names.

	//Write the unique BVs producing three generation models.
	set<vector<string> >::iterator itUnique_BVs = 
		Unique_3G_LRSym_BVs.begin();
	ThreeG_Out<<"Unique Basis vectors producing these models:"
		<<Unique_3G_LRSym_BVs.size()<<endl;
	for(itUnique_BVs; itUnique_BVs != Unique_3G_LRSym_BVs.end(); 
			++itUnique_BVs)
	{
		for(int a=0; a< itUnique_BVs->size(); ++a)
			ThreeG_Out<<itUnique_BVs->at(a)<<endl;
	}

	//Write the unique BVs producing three generation and anti-generation
	//models.
	itUnique_BVs = Unique_3G_3AG_LRSym_BVs.begin();
	ThreeG_ThreeAG_Out<<"Unique Basis vectors producing these models:"
		<<Unique_3G_3AG_LRSym_BVs.size()<<endl;
	for(itUnique_BVs; itUnique_BVs != Unique_3G_3AG_LRSym_BVs.end(); 
			++itUnique_BVs)
	{
		for(int a=0; a< itUnique_BVs->size(); ++a)
			ThreeG_ThreeAG_Out<<itUnique_BVs->at(a)<<endl;
	}

	Statistical_Report_Generator ThreeG_Stat_Reporter
		(Unique_3G_LRSym_LEEFTs);

	string Doc_Title = 
		"Statistics for Three Generation Left-Right Symmetric";
	Doc_Title += " NAHE Based Models, Order 3, Layer 1";
	Jared_Greenwald.Begin_Latex_Document(Latex_Out, Doc_Title);
	//Three Generation models.
	Latex_Out<<"\\textbf{Unique Three Generation Models: ";
	Latex_Out<<Unique_3G_LRSym_LEEFTs.size()<<"}\\\\"<<endl;
	Latex_Out<<"\\vspace{3 mm}"<<endl;
	ThreeG_Stat_Reporter.Generate_Statistical_Report(Latex_Out, 4);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;

	//Three Generation and Anti-Generation models.
	Statistical_Report_Generator ThreeG_ThreeAG_Stat_Reporter
		(Unique_3G_3AG_LRSym_LEEFTs);
	Latex_Out<<"\\textbf{Unique Three Generation and ";
	Latex_Out<<"Anti-Generation Models: ";
	Latex_Out<<Unique_3G_3AG_LRSym_LEEFTs.size()<<"}\\\\"<<endl;
	Latex_Out<<"\\vspace{3 mm}"<<endl;
	ThreeG_ThreeAG_Stat_Reporter.Generate_Statistical_Report(Latex_Out, 4);
	map<int, int> Charged_3G_3AG_Exotics = 
		Find_OS_Charged_Exotics(Unique_3G_3AG_LRSym_LEEFTs);
	Jared_Greenwald.Write_OS_Charged_Exotics(Latex_Out, Charged_3G_3AG_Exotics);
	Latex_Out<<"\\hspace{\\fill}\\\\vspace{3 mm}"<<endl;
	Jared_Greenwald.Plot_OS_Charged_Exotics(Latex_Out, Charged_3G_3AG_Exotics);

	//Finish the Latex document.
	Jared_Greenwald.Finish_Latex_Document(Latex_Out);
	cout<<"Program finished."<<endl;
	return 0;
}//Close main.

bool Has_LRSym_Group(const LEEFT& LEEFT_Input)
{
	bool Has_SU3 = false;
	bool Has_2_SU2s = false;
	int SU2_Count = 0;

	for(int a=0; a<LEEFT_Input.Gauge_Groups().size(); ++a)
	{
		if(LEEFT_Input.Gauge_Groups().at(a).Class() == 'A')
		{
			if(LEEFT_Input.Gauge_Groups().at(a).Rank() == 2)
				Has_SU3 = true;
			else if((LEEFT_Input.Gauge_Groups().at(a).Rank() == 1)
					&& !Has_2_SU2s)
				SU2_Count++;

			if(SU2_Count == 2)
				Has_2_SU2s = true;
			if(Has_SU3 && Has_2_SU2s)
				return true;
		}//Close if statement on A class gauge groups.
	}//Close for loop on Gauge_Groups.

	return false;
}//Close Has_LRSym_Group.

LEEFT Find_All_LRSym_Observable_Generations(LEEFT LEEFT_Input)
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
	vector<Observable_Sector> LRSym_Observable_Sectors;
	for(int a=0; a<SU3_Indices.size(); ++a)
	{
		for(int b=0; b<SU2_Indices.size(); ++b)
		{
			for(int c=0; c<SU2_Indices.size(); ++c)
			{
				if(c != b)
				{
					int SU3_Index = SU3_Indices.at(a);
					int SU2A_Index = SU2_Indices.at(b);
					int SU2B_Index = SU2_Indices.at(c);
					Observable_Sector New_LRSym_Observable_Sector = 
						Find_LRSym_Observable_Sector_States(LEEFT_Input, SU3_Index, 
								SU2A_Index, SU2B_Index);
					//Make sure it isn't empty.
					if(New_LRSym_Observable_Sector.Generations().size() != 0)
						LRSym_Observable_Sectors.push_back(New_LRSym_Observable_Sector);
				}//Close if statement on c and b.
			}//Close for loop over SU2_Indices.
		}//Close for loop over SU2_Indices.
	}//Close for loop over SU3_Indices.
	
	LEEFT_Input.Set_Observable_Sectors(LRSym_Observable_Sectors);
	return LEEFT_Input;
}//Close Find_All_LRSym_Observable_Generations.

Observable_Sector Find_LRSym_Observable_Sector_States(const LEEFT&
		LEEFT_Input, int SU3_Index, int SU2A_Index, int SU2B_Index)
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
					SU3_Index, SU2A_Index, SU2B_Index))
		{
			if((LEEFT_Input.Representations().at(a).at(SU3_Index) == 3) &&
					(LEEFT_Input.Representations().at(a).at(SU2A_Index) == 2)
					&&(LEEFT_Input.Representations().at(a).at(SU2B_Index)==1))
			{
				Generations.push_back(LEEFT_Input.Representations().at(a));
				Generation_Qty.push_back
					(LEEFT_Input.Representation_Count().at(a));
			}else if((LEEFT_Input.Representations().at(a).at(SU3_Index) == -3)
					&&(LEEFT_Input.Representations().at(a).at(SU2A_Index) == 2)
					&&(LEEFT_Input.Representations().at(a).at(SU2B_Index)==1))
			{
				Barred_Generations.push_back(LEEFT_Input.Representations().at(a));
				Barred_Generation_Qty.push_back
					(LEEFT_Input.Representation_Count().at(a));
			}else if((LEEFT_Input.Representations().at(a).at(SU3_Index) == 3)
					&& (LEEFT_Input.Representations().at(a).at(SU2A_Index) == 1)
					&& (LEEFT_Input.Representations().at(a).at(SU2B_Index)==2))
			{
				Anti_Generations.push_back(LEEFT_Input.Representations().at(a));
				Anti_Generation_Qty.push_back
					(LEEFT_Input.Representation_Count().at(a));
			}else if((LEEFT_Input.Representations().at(a).at(SU3_Index) == -3)
					&& (LEEFT_Input.Representations().at(a).at(SU2A_Index) == 1)
					&& (LEEFT_Input.Representations().at(a).at(SU2B_Index)==2))
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
		Generation_Rep.at(SU2A_Index) = 2;
		Generation_Rep.at(SU2B_Index) = 1;
		Generations.push_back(Generation_Rep);
		Generation_Qty.push_back(0);
	}
	if(Barred_Generations.size() == 0)
	{
		Generation_Rep.at(SU3_Index) = -3;
		Generation_Rep.at(SU2A_Index) = 2;
		Generation_Rep.at(SU2B_Index) = 1;
		Barred_Generations.push_back(Generation_Rep);
		Barred_Generation_Qty.push_back(0);
	}
	if(Anti_Generations.size() == 0)
	{
		Generation_Rep.at(SU3_Index) = 3;
		Generation_Rep.at(SU2A_Index) = 1;
		Generation_Rep.at(SU2B_Index) = 2;
		Anti_Generations.push_back(Generation_Rep);
		Anti_Generation_Qty.push_back(0);
	}
	if(Anti_Barred_Generations.size() == 0)
	{
		Generation_Rep.at(SU3_Index) = -3;
		Generation_Rep.at(SU2A_Index) = 1;
		Generation_Rep.at(SU2B_Index) = 2;
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
}//Close Find_LRSym_Observable_Sector_States.

int Count_OS_Charged_Exotics(const LEEFT& LEEFT_Input, int SU3_Index,
		int SU2A_Index, int SU2B_Index)
{
	int OS_Charged_Exotics = 0;
	for(int a=0; a<LEEFT_Input.Representations().size(); ++a)
	{
		if(((LEEFT_Input.Representations().at(a).at(SU3_Index) != 1) ||
				(LEEFT_Input.Representations().at(a).at(SU2A_Index) != 1) ||
				(LEEFT_Input.Representations().at(a).at(SU2B_Index) != 1)) && 
				(Has_Hidden_Sector_Charge(LEEFT_Input.Representations().at(a),
																	SU3_Index, SU2A_Index, SU2B_Index)))
			OS_Charged_Exotics += LEEFT_Input.Representation_Count().at(a);
	}//Close for loop over the Representations.
	return OS_Charged_Exotics;
}// Close Count_OS_Charged_Exotics.

map<int, int> Find_OS_Charged_Exotics(const set<LEEFT>& Unique_LEEFTs)
{
	map<int, int> OS_Charged_Exotics;
	set<LEEFT>::const_iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
	for(itUnique_LEEFTs; itUnique_LEEFTs != Unique_LEEFTs.end(); 
			++itUnique_LEEFTs)
	{
		vector<int> SU3_Indices;
		for(int a=0; a<itUnique_LEEFTs->Gauge_Groups().size(); ++a)
		{
			if((itUnique_LEEFTs->Gauge_Groups().at(a).Class() == 'A') && 
					(itUnique_LEEFTs->Gauge_Groups().at(a).Rank() == 2))
				SU3_Indices.push_back(a);
		}//Close for loop for finding the SU(3) indices.

		vector<int> SU2_Indices;
		for(int a=0; a<itUnique_LEEFTs->Gauge_Groups().size(); ++a)
		{
			if((itUnique_LEEFTs->Gauge_Groups().at(a).Class() == 'A') && 
					(itUnique_LEEFTs->Gauge_Groups().at(a).Rank() == 1))
				SU2_Indices.push_back(a);
		}//Close for loop for finding the SU(2) indices.

		//Loop over the three groups, count the generations,
		//and add them to the map.
		for(int a=0; a<SU3_Indices.size(); ++a)
		{
			for(int b=0; b<SU2_Indices.size(); ++b)
			{
				for(int c=0; c<SU2_Indices.size(); ++c)
				{
					if(c != b)
					{
						int SU3_Index = SU3_Indices.at(a);
						int SU2A_Index = SU2_Indices.at(b);
						int SU2B_Index = SU2_Indices.at(c);
						int Number_Charged_Exotics = 
							Count_OS_Charged_Exotics(*itUnique_LEEFTs,
									SU3_Index, SU2A_Index, SU2B_Index);
						if(OS_Charged_Exotics.find(Number_Charged_Exotics) !=
								OS_Charged_Exotics.end())
							OS_Charged_Exotics[Number_Charged_Exotics]++;
						else
							OS_Charged_Exotics[Number_Charged_Exotics] = 1;
					}
				}//Close for loop over the second SU(2) indices.
			}//Close for loop over the first SU(2) indices.
		}//Close for loop over the SU(4) indices.
	}//Close for loop on Unique_LEEFTs.
	return OS_Charged_Exotics;
}//Close Find_OS_Charged_Exotics.
bool Has_Hidden_Sector_Charge(const vector<int>& Representation,
		int SU3_Index, int SU2A_Index, int SU2B_Index)
{
  for(int a=0; a<Representation.size(); ++a)
	{
		if((a != SU3_Index) && (a != SU2A_Index) &&
				(a != SU2B_Index) && 
				(Representation.at(a) != 1))
			return true;
	}//Close for loop on Representation.
	return false;
}//Close Has_Hidden_Sector_Charge.
