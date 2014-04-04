/*
Authors: Timothy Renner (renner.timothy@gmail.com)
				 Doug Moore (Douglas_Moore1@baylor.edu)
				 Gerald Cleaver (Gerald_Cleaver@baylor.edu)

Institute: Baylor University
Created: 6/2/2011

This program searches through a list of files and 
prints those with mirrored gauge groups only, and
those with mirrored matter representations as well.
Statistics are prepared for each case.

Uses the FF Framework
*/


#include <vector>
#include <algorithm>
#include<numeric>
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

bool Has_Mirrored_Gauge_Groups(const LEEFT& LEEFT_Input);
bool Has_Mirrored_Matter_Reps(const LEEFT& LEEFT_Input);

int main()
{
	cout<<"Begin program."<<endl;
	int Layers = 1;
	
	//Set up the input files.
	string FF_Data = "/Users/timothyrenner/FF_Framework/data/";
	vector<string> Input_File_Names;
	Input_File_Names.push_back(FF_Data + 
			"/NAHE_Var_Extensions/One_Layer/NAHE_Var_O3_L1.txt");
	//Input_File_Names.push_back(FF_Data + 
		//	"/NAHE_Var_Extensions/One_Layer/NAHE_Var_O2_L1_NSP.txt");
	//Set up the output files.
	ofstream Mirrored_Gauge_Groups_Out("Mirrored_Gauge_Models.txt");
	ofstream Mirrored_Matter_Out("Mirrored_Matter_Models.txt");
	ofstream Latex_Out("./tex/Mirrored_Model_Stats.tex");
	Output_Writer Doug_Moore;
	Latex_Writer Jared_Greenwald;

	//Read in and analyze the LEEFTs.
	cout<<"Reading in the LEEFTs with mirroring."<<endl;
	set<LEEFT> Unique_Mirrored_Gauge_LEEFTs;
	set<vector<string> > Unique_Mirrored_Gauge_BVs;
	set<LEEFT> Unique_Mirrored_Matter_LEEFTs;
	set<vector<string> > Unique_Mirrored_Matter_BVs;

	
	for(int a=0; a<Input_File_Names.size(); ++a)
	{
		ifstream LEEFT_In(Input_File_Names.at(a).c_str());
		if(LEEFT_In.fail())
		{
			cout<<Input_File_Names.at(a)<<" did not open."
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

			if(Has_Mirrored_Gauge_Groups(New_LEEFT))
			{
				int Old_Gauge_Size = Unique_Mirrored_Gauge_LEEFTs.size();
				Unique_Mirrored_Gauge_LEEFTs.insert(New_LEEFT);
				if(Unique_Mirrored_Gauge_LEEFTs.size() > Old_Gauge_Size)
				{
					Doug_Moore.Write_LEEFT_Particle_Content(Mirrored_Gauge_Groups_Out,
							New_LEEFT);
					for(int b=0; b<Model_Inputs.size(); ++b)
						Mirrored_Gauge_Groups_Out<<Model_Inputs.at(b)<<endl;
					Mirrored_Gauge_Groups_Out<<endl;
					vector<string> Unique_Mirrored_Gauge_BVs_Loader;
					for(int b=0; b<Model_Inputs.size()/2; ++b)
						Unique_Mirrored_Gauge_BVs_Loader.push_back(Model_Inputs.at(b));
					Unique_Mirrored_Gauge_BVs.insert
						(Unique_Mirrored_Gauge_BVs_Loader);
				}

				if(Has_Mirrored_Matter_Reps(New_LEEFT))
				{
					int Old_Matter_Size = Unique_Mirrored_Matter_LEEFTs.size();
					Unique_Mirrored_Matter_LEEFTs.insert(New_LEEFT);
					if(Unique_Mirrored_Matter_LEEFTs.size() > Old_Matter_Size)
					{
						Doug_Moore.Write_LEEFT_Particle_Content(Mirrored_Matter_Out,
								New_LEEFT);
						for(int b=0; b<Model_Inputs.size(); ++b)
							Mirrored_Matter_Out<<Model_Inputs.at(b)<<endl;
						Mirrored_Matter_Out<<endl;
						vector<string> Unique_Mirrored_Matter_BVs_Loader;
						for(int b=0; b<Model_Inputs.size()/2; ++b)
							Unique_Mirrored_Matter_BVs_Loader.push_back
								(Model_Inputs.at(b));
						Unique_Mirrored_Matter_BVs.insert
							(Unique_Mirrored_Matter_BVs_Loader);
					}
				}//Close if statement on mirrored matter.
			}//Close if statement on classifying LEEFTs.
		}//Close while loop on input file.
	}//Close for loop on input files.

	//Write the Mirrored LEEFT metadata to a file.`
	Mirrored_Gauge_Groups_Out<<
		"Total unique LEEFTs with mirrored gauge groups: "
		<<Unique_Mirrored_Gauge_LEEFTs.size()<<endl;
	Mirrored_Gauge_Groups_Out<<
		"Total unique BVs producing mirrored gauge groups: "
		<<Unique_Mirrored_Gauge_BVs.size()<<endl;
	set<vector<string> >::iterator itUnique_BVs = 
		Unique_Mirrored_Gauge_BVs.begin();
	for(itUnique_BVs; itUnique_BVs != Unique_Mirrored_Gauge_BVs.end();
			++itUnique_BVs)
	{
		for(int a=0; a<itUnique_BVs->size(); ++a)
		Mirrored_Gauge_Groups_Out<<itUnique_BVs->at(a)<<endl;
	}

	Mirrored_Matter_Out<<"Total unique LEEFTs with mirrored matter: "
		<<Unique_Mirrored_Matter_LEEFTs.size()<<endl;
	Mirrored_Matter_Out<<"Total unique BVs producing mirrored gauge groups: "
		<<Unique_Mirrored_Matter_BVs.size()<<endl;
	itUnique_BVs = Unique_Mirrored_Matter_BVs.begin();
	for(itUnique_BVs; itUnique_BVs != Unique_Mirrored_Matter_BVs.end();
			++itUnique_BVs)
	{
		for(int a=0; a<itUnique_BVs->size(); ++a)
		Mirrored_Matter_Out<<itUnique_BVs->at(a)<<endl;
	}

	//Generate and write the statistics for the mirrored gauge groups.
	string Doc_Title = "Statistics for mirrored models";
	Doc_Title += " from the NAHE Variation, Order 2, Layer 1";
	Jared_Greenwald.Begin_Latex_Document(Latex_Out, Doc_Title);
	Statistical_Report_Generator Satheeshkumar(Unique_Mirrored_Gauge_LEEFTs);
  Latex_Out<<"\\textbf{Unique Mirrored Gauge Group Models: "<<
		Unique_Mirrored_Gauge_LEEFTs.size()<<"}\\\\"<<endl;
	Latex_Out<<"\\vspace{3 mm}"<<endl;
	Satheeshkumar.Generate_Statistical_Report(Latex_Out, 4);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;

	//Generate and write the statistics for the mirrored matter reps.
	Statistical_Report_Generator Jonathan_Perry
		(Unique_Mirrored_Matter_LEEFTs);
	Latex_Out<<"\\textbf{Unique Mirrored Matter Rep Models: "<<
		Unique_Mirrored_Matter_LEEFTs.size()<<"}\\\\"<<endl;
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;
	Jonathan_Perry.Generate_Statistical_Report(Latex_Out, 4);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;

	//Finish the Latex Document
	Jared_Greenwald.Finish_Latex_Document(Latex_Out);

	cout<<"Program finished."<<endl;
	return 0;
}

bool Has_Mirrored_Gauge_Groups(const LEEFT& LEEFT_Input)
{
	map<Gauge_Group_Name, int> Gauge_Group_Counts;
	for(int a=0; a<LEEFT_Input.Gauge_Groups().size(); ++a)
	{
		if(Gauge_Group_Counts.find(LEEFT_Input.Gauge_Groups().at(a)) !=
				Gauge_Group_Counts.end())
			Gauge_Group_Counts[LEEFT_Input.Gauge_Groups().at(a)]++;
		else
			Gauge_Group_Counts[LEEFT_Input.Gauge_Groups().at(a)] = 1;
	}//Close for loop over the gauge group counting.
	
	map<Gauge_Group_Name,int>::iterator itGauge_Group_Counts = 
		Gauge_Group_Counts.begin();
	for(itGauge_Group_Counts; itGauge_Group_Counts != 
			Gauge_Group_Counts.end(); ++itGauge_Group_Counts)
	{
		if((itGauge_Group_Counts->second)%2!=0)
			return false;
	}//Close for loop over the map.
	return true;
}//Close Has_Mirrored_Gauge_Groups.

bool Has_Mirrored_Matter_Reps(const LEEFT& LEEFT_Input)
{
	//Changing the 1's to 0's in the representations allows for
	//them to be formed into mutually orthogonal sets of states.
	//This will be used to deal with the mirroring.

	//First, create a master list of states, and convert the 1's to 0's
	list<vector<int> > Zeroed_Reps;
	for(int a=0; a<LEEFT_Input.Representations().size(); ++a)
	{
		vector<int> New_Zeroed_Representation = 
			LEEFT_Input.Representations().at(a);
		//Convert the 1's to 0's.
		for(int b=0; b<New_Zeroed_Representation.size(); ++b)
		{
			if(New_Zeroed_Representation.at(b) == 1)
				New_Zeroed_Representation.at(b) = 0;
		}//Close for loop for converting the 1's to 0's.
		for(int b=0; b<LEEFT_Input.Representation_Count().at(a); ++b)
			Zeroed_Reps.push_back(New_Zeroed_Representation);
	}//Close for loop on Representations.

	//Arrange the newly cast representations into mutually orthogonal sets.
	vector<int> Rep_Set_Sizes;
	while(Zeroed_Reps.size() != 0)
	{
		list<vector<int> >::iterator itZeroed_Reps = Zeroed_Reps.begin();
		list<vector<int> > Zeroed_Rep_Subset;
		Zeroed_Rep_Subset.push_back(*itZeroed_Reps);
		Zeroed_Reps.erase(itZeroed_Reps);

		list<vector<int> >::iterator itSubset = Zeroed_Rep_Subset.begin();
		for(itSubset; itSubset != Zeroed_Rep_Subset.end(); ++itSubset)
		{
			//Compute the dot products with the member of the main list.
			for(itZeroed_Reps = Zeroed_Reps.begin(); 
					itZeroed_Reps != Zeroed_Reps.end(); )
			{
				if(inner_product(itSubset->begin(), itSubset->end(),
							itZeroed_Reps->begin(),0)!=0)
				{
					Zeroed_Rep_Subset.push_back(*itZeroed_Reps);
					itZeroed_Reps = Zeroed_Reps.erase(itZeroed_Reps);
				}else
					++itZeroed_Reps;
			}
		}//Close for loop on the subset.
		Rep_Set_Sizes.push_back(Zeroed_Rep_Subset.size());
	}//Close while loop.
	
	//Now determine if there are an even number of equally sized reps.
	map<int, int> Rep_Subset_Counts;
	for(int a=0; a<Rep_Set_Sizes.size(); ++a)
	{
		if(Rep_Subset_Counts.find(Rep_Set_Sizes.at(a)) != 
				Rep_Subset_Counts.end())
			Rep_Subset_Counts[Rep_Set_Sizes.at(a)]++;
		else
			Rep_Subset_Counts[Rep_Set_Sizes.at(a)] = 1;
	}
	map<int, int>::iterator itRep_Subset_Counts = Rep_Subset_Counts.begin();
	for(itRep_Subset_Counts; itRep_Subset_Counts != Rep_Subset_Counts.end();
			++itRep_Subset_Counts)
	{
		if((itRep_Subset_Counts->second)%2!=0)
			return false;
	}
	return true;
}//Close Has_Mirrored_Matter_Reps.
