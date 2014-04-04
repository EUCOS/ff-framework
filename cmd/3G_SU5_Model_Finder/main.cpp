/*
Authors: Timothy Renner (renner.timothy@gmail.com)
				 Doug Moore (Douglas_Moore1@baylor.edu)
				 Gerald Cleaver (Gerald_Cleaver@baylor.edu)

Institute: Baylor University
Created: 5/18/2011

This program searches through a list of files
and prints those with three chiral matter
generations for an SU5 gauge group.
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
#include "FF_Statistics.hh"
#include "FF_Latex_Writer.hh"
#include "FF_Output_Writer.hh"
#include "FF_Gauge_Group_Name.hh"
#include "FF_LEEFT_Tools.hh"
#include "FF_Statistical_Report_Generator.hh"

using namespace std;

bool Has_SU5_Group(const LEEFT& LEEFT_Input);

LEEFT Find_All_SU5_Observable_Generations(LEEFT LEEFT_Input);
Observable_Sector Find_SU5_Observable_Sector_States
	(const LEEFT& LEEFT_Input, int SU5_Index);
int Find_Hidden_Sector_Multiplicity(const vector<int>& Representation, 
		int SU5_Index);

int Count_OS_Charged_Exotics(const LEEFT& LEEFT_Input,
		int SU5_Index);
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
	ofstream ThreeG_Out("3G_SU5_Models.txt");
	ofstream Latex_Out("./tex/3G_SU5_Model_Stats.tex");
	Output_Writer Doug_Moore;
	Latex_Writer Jared_Greenwald;

	//Read in and analyze the LEEFTs.
	cout<<"Reading in LEEFTs with SU(5) gauge groups."<<endl;
	set<LEEFT> Unique_3G_SU5_LEEFTs;
	set< vector<string> > Unique_BVs;
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
			if(Has_SU5_Group(New_LEEFT))
			{
				New_LEEFT = Find_All_SU5_Observable_Generations(New_LEEFT);
				bool ThreeG_Model_Written = false;
				for(int a=0; a<New_LEEFT.Observable_Sectors().size(); ++a)
				{
					int Chiral_Generations = New_LEEFT.
						Observable_Sectors().at(a).Count_Chiral_Generations();

					if(!ThreeG_Model_Written && 
							(Chiral_Generations == 3))
					{
						vector<string> Unique_BV_Loader;
						for(int b=0; b<Model_Inputs.size()/2; ++b)
						{
							Unique_BV_Loader.push_back(Model_Inputs.at(b));
						}
						Unique_BVs.insert(Unique_BV_Loader);
						Doug_Moore.Write_LEEFT_Particle_Content(ThreeG_Out, New_LEEFT);
						for(int b=0; b<Model_Inputs.size(); ++b)
							ThreeG_Out<<Model_Inputs.at(b)<<endl;
						ThreeG_Out<<endl;
						Unique_3G_SU5_LEEFTs.insert(New_LEEFT);
						ThreeG_Model_Written = true;
					}//Close if statement on three generations.

					if(ThreeG_Model_Written)
						break;
				}//Close for loop on Observable_Sectors.
			}//Close if statement on Has_SU5_Group.
		}//Close while loop on LEEFT_In.
	}//Close for loop on Input_File_Names.

	set<vector<string> >::iterator itUnique_BVs = Unique_BVs.begin();
	ThreeG_Out<<"Unique Basis vectors producing these models:"
		<<Unique_BVs.size()<<endl;
	for(itUnique_BVs; itUnique_BVs != Unique_BVs.end(); ++itUnique_BVs)
	{
		for(int a=0; a< itUnique_BVs->size(); ++a)
			ThreeG_Out<<itUnique_BVs->at(a)<<endl;
	}

	//Get OS Charged exotics.
	map<int, int> Exotics = Find_OS_Charged_Exotics(Unique_3G_SU5_LEEFTs);

	//Write the statistics to a Latex file.
	Statistical_Report_Generator ThreeG_Stat_Reporter(Unique_3G_SU5_LEEFTs);
	string Doc_Title = "Statistics for Three Generation $SU(5)$ ";
	Doc_Title += "NAHE Based Models, Order 3, Layer 1";
	Jared_Greenwald.Begin_Latex_Document(Latex_Out, Doc_Title);
	Latex_Out<<"\\textbf{Unique Three Generation Models: ";
	Latex_Out<<Unique_3G_SU5_LEEFTs.size()<<"}\\\\"<<endl;
	Latex_Out<<"\\vspace{3 mm}"<<endl;
	ThreeG_Stat_Reporter.Generate_Statistical_Report(Latex_Out, 4);
	Jared_Greenwald.Write_OS_Charged_Exotics(Latex_Out, Exotics);
	Latex_Out<<"\\hspace{\\fill}\\vspace{3 mm}\\\\"<<endl;
	Jared_Greenwald.Plot_OS_Charged_Exotics(Latex_Out, Exotics);
	Jared_Greenwald.Finish_Latex_Document(Latex_Out);


	cout<<"Program finished."<<endl;
	return 0;
}//Close main.

bool Has_SU5_Group(const LEEFT& LEEFT_Input)
{
	bool Has_SU5 = false;

	if(LEEFT_Input.U1_Factors() == 0)
		return false;

	for(int a=0; a<LEEFT_Input.Gauge_Groups().size(); ++a)
	{
		if((LEEFT_Input.Gauge_Groups().at(a).Class() == 'A')
				&&(LEEFT_Input.Gauge_Groups().at(a).Rank() == 4))
			return true;
	}//Close for loop on Gauge_Groups.

	return false;
}//Close Has_SU5_Group.

LEEFT Find_All_SU5_Observable_Generations(LEEFT LEEFT_Input)
{
	//Get the indices of the SU5 gauge groups.
	vector<int> SU5_Indices;
	for(int a=0; a<LEEFT_Input.Gauge_Groups().size(); ++a)
	{
		if((LEEFT_Input.Gauge_Groups().at(a).Class() == 'A')
				&&(LEEFT_Input.Gauge_Groups().at(a).Rank()==4))
			SU5_Indices.push_back(a);
	}//Close for loop over the gauge groups.

	//Now find the observable sector generations for each observable sector.
	vector<Observable_Sector> SU5_Observable_Sectors;
	for(int a=0; a<SU5_Indices.size(); ++a)
	{
		int SU5_Index = SU5_Indices.at(a);
		Observable_Sector New_SU5_Observable_Sector = 
			Find_SU5_Observable_Sector_States(LEEFT_Input, SU5_Index);
		//Make sure it isn't empty.
	if(New_SU5_Observable_Sector.Generations().size() != 0)
		SU5_Observable_Sectors.push_back(New_SU5_Observable_Sector);
}//Close for loop over SU3_Indices.

LEEFT_Input.Set_Observable_Sectors(SU5_Observable_Sectors);
	return LEEFT_Input;
}//Close Find_All_SU5_Observable_Generations.

Observable_Sector Find_SU5_Observable_Sector_States(const LEEFT&
		LEEFT_Input, int SU5_Index)
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
		int Hidden_Sector_Multiplicity = 
			Find_Hidden_Sector_Multiplicity(LEEFT_Input.Representations().at(a), 
					SU5_Index);
			if(((LEEFT_Input.Representations().at(a).at(SU5_Index) == 10) ||
					(LEEFT_Input.Representations().at(a).at(SU5_Index) == -5))&& 
					(Hidden_Sector_Multiplicity == 1))
			{
				Generations.push_back(LEEFT_Input.Representations().at(a));
				Generation_Qty.push_back
					(LEEFT_Input.Representation_Count().at(a)*Hidden_Sector_Multiplicity);
			}else if(((LEEFT_Input.Representations().at(a).at(SU5_Index) == -10)
					||(LEEFT_Input.Representations().at(a).at(SU5_Index) == 5)) && 
					(Hidden_Sector_Multiplicity == 1))
			{
				Barred_Generations.push_back(LEEFT_Input.Representations().at(a));
				Barred_Generation_Qty.push_back
					(LEEFT_Input.Representation_Count().at(a)*Hidden_Sector_Multiplicity);
			}//Close if/else chain for finding the types of generations.
	}//Close for loop on Representations.

	//Now fill in the blanks for any missing representations.
	vector<int> Generation_Rep
		(LEEFT_Input.Representations().front().size(), 1);
	if(Generations.size() == 0)
	{
		Generation_Rep.at(SU5_Index) = 10;
		Generations.push_back(Generation_Rep);
		Generation_Qty.push_back(0);
		Generation_Rep.at(SU5_Index) = -5;
		Generations.push_back(Generation_Rep);
		Generation_Qty.push_back(0);
	}else if(Generations.size() == 1)
	{
		if(Generations.at(0).at(SU5_Index) == 10)
			Generation_Rep.at(SU5_Index) = -5;
		else if(Generations.at(0).at(SU5_Index) == -5)
			Generation_Rep.at(SU5_Index)=10;
		Generations.push_back(Generation_Rep);
		Generation_Qty.push_back(0);
	}//Close if/else on incomplete generations.
	Generation_Rep = vector<int> 
		(LEEFT_Input.Representations().front().size(),1);
	if(Barred_Generations.size() == 0)
	{
		Generation_Rep.at(SU5_Index) = -10;
		Barred_Generations.push_back(Generation_Rep);
		Barred_Generation_Qty.push_back(0);
		Generation_Rep.at(SU5_Index) = 5;
		Barred_Generations.push_back(Generation_Rep);
		Barred_Generation_Qty.push_back(0);
	}else if(Barred_Generations.size()==1)
	{
		if(Barred_Generations.at(0).at(SU5_Index) == -10)
			Generation_Rep.at(SU5_Index)=5;
		else if(Barred_Generations.at(0).at(SU5_Index)==5)
			Generation_Rep.at(SU5_Index)=-10;
		Barred_Generations.push_back(Generation_Rep);
		Barred_Generation_Qty.push_back(0);
	}//Close if/else on incomplete barred generations.
	Generation_Rep.clear();

	//Create a new Observable_Sector object and return it.
	Observable_Sector New_Observable_Sector;
	New_Observable_Sector.Set_Generations(Generations);
	New_Observable_Sector.Set_Generation_Qty(Generation_Qty);
	New_Observable_Sector.Set_Barred_Generations(Barred_Generations);
	New_Observable_Sector.Set_Barred_Generation_Qty(Barred_Generation_Qty);

	return New_Observable_Sector;
}//Close Find_SU5_Observable_Sector_States.

int Find_Hidden_Sector_Multiplicity(const vector<int>& Representation,
		int SU5_Index)
{
	int Hidden_Sector_Multiplicity = 1;
  for(int a=0; a<Representation.size(); ++a)
	{
		if(a != SU5_Index)
			Hidden_Sector_Multiplicity *= Representation.at(a);
	}//Close for loop on Representation.
	return Hidden_Sector_Multiplicity;
}//Close Has_Hidden_Sector_Charge.

int Count_OS_Charged_Exotics(const LEEFT& LEEFT_Input, int SU5_Index)
{
	int OS_Charged_Exotics = 0;
	for(int a=0; a<LEEFT_Input.Representations().size(); ++a)
	{
		if((LEEFT_Input.Representations().at(a).at(SU5_Index) != 1) &&
				(Find_Hidden_Sector_Multiplicity(LEEFT_Input.Representations().at(a),
																	SU5_Index)>1))
			OS_Charged_Exotics += LEEFT_Input.Representation_Count().at(a);
	}//Close for loop over the Representations.
	return OS_Charged_Exotics;
}//Close Count_OS_Charged_Exotics.

map<int, int> Find_OS_Charged_Exotics(const set<LEEFT>& Unique_LEEFTs)
{
	map<int, int> OS_Charged_Exotics;
	set<LEEFT>::const_iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
	for(itUnique_LEEFTs; itUnique_LEEFTs != Unique_LEEFTs.end(); 
			++itUnique_LEEFTs)
	{
		vector<int> SU5_Indices;
		for(int a=0; a<itUnique_LEEFTs->Gauge_Groups().size(); ++a)
		{
			if((itUnique_LEEFTs->Gauge_Groups().at(a).Class() == 'A') &&
					(itUnique_LEEFTs->Gauge_Groups().at(a).Rank() == 4))
				SU5_Indices.push_back(a);
	}//Close for loop for finding all of the SU(5) indices.
	
		//Loop over the SU(5) indices, count the number of observable sector
		//charged exotics, and add them to the map.
		for(int a=0; a<SU5_Indices.size(); ++a)
		{
			int SU5_Index = SU5_Indices.at(a);
			int Number_OS_Charged_Exotics = 
				Count_OS_Charged_Exotics(*itUnique_LEEFTs, SU5_Index);
			if(OS_Charged_Exotics.find(Number_OS_Charged_Exotics) !=
					OS_Charged_Exotics.end())
				OS_Charged_Exotics[Number_OS_Charged_Exotics]++;
			else
				OS_Charged_Exotics[Number_OS_Charged_Exotics] = 1;
		}//Close for loop over the SO(10) indices.
	}//Close for loop over the unique LEEFTs.
	return OS_Charged_Exotics;
}//Close Find_OS_Charged_Exotics.
