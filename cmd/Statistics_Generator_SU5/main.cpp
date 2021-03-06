/*
	Authors: Timothy Renner (renner.timothy@gmail.com)
				   Doug Moore (Douglas_Moore1@baylor.edu)
					 Gerald Cleaver, Ph.D. (Gerald_Cleaver@baylor.edu)
	Date: 5/17/2011
	Baylor University

	This program gathers relevant statistical information for models
	containing at least one SU(5)xU(1) group. Tabulates gauge group factors,
	number of ST SUSYs, number of chiral matter representations, etc.

	Uses the FF Framework
*/

#include <string>
#include <algorithm>
#include <set>
#include <map>
#include <assert.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "FF_LEEFT.hh"
#include "FF_Gauge_Group_Name.hh"
#include "FF_Statistics.hh"
#include "FF_Latex_Writer.hh"

using namespace std;

void Read_Front_Info(ifstream& LEEFT_In);
vector<string> Read_LEEFT_Data(ifstream& LEEFT_In);
void Skip_Model_Inputs(ifstream& LEEFT_In, int Skipped_Lines);

bool Has_SU5_Group(const LEEFT& LEEFT_Input);

int Count_HS_Charged_Generations(const LEEFT& LEEFT_Input, 
		int SU5_Index);
map<int, int> Find_HS_Charged_Generations
	(const set<LEEFT>& Unique_LEEFTs);

int Count_HS_Uncharged_Generations(const LEEFT& LEEFT_Input, 
		int SU5_Index);
map<int, int> Find_HS_Uncharged_Generations
	(const set<LEEFT>& Unique_LEEFTs);

int Count_OS_Charged_Exotics(const LEEFT& LEEFT_Input,
		int SU5_Index);
map<int, int> Find_OS_Charged_Exotics(const set<LEEFT>& Unique_LEEFTs);

bool Has_Hidden_Sector_Charge(const vector<Group_Representation>& 
		Representations, int SU5_Index);
int Get_HS_Multiplicity(const vector<Group_Representation>& 
		Representations, int SU5_Index, 
		const vector<Gauge_Group_Name>& Gauge_Groups);

int main()
{
	string FF_Data = "/Users/timothyrenner/FF_Framework/data/";
	int Large_ST_Dimensions = 4;
	//Set the output file name.
	string Output_File_Name = "./tex/SU5_Model_Statistics.tex";

	//Set the Latex document title.
	string Doc_Title = "Statistical Report for NAHE Based $SU(5)\\otimes U(1)$ Models, ";
		Doc_Title += "Order 3, Layer 1";

	//Load the file names.
	vector<string> Input_File_Names;
	Input_File_Names.push_back
		(FF_Data + "NAHE_Extensions/One_Layer/NAHE_O3_L1.txt");
	//Input_File_Names.push_back
		//(FF_Data + "NAHE_Extensions/One_Layer/NR_NAHE_O2_L1_NSP.txt");

	//int Skipped_Lines = 2;//For no layers.
	//int Skipped_Lines = 3; //For NR NAHE.
	int Skipped_Lines = 4;//For one layer.

	cout<<"Reading unique LEEFTs containing SU(5) gauge groups."<<endl;
	set<LEEFT> Unique_LEEFTs;
	for(int a=0; a<Input_File_Names.size(); ++a)
	{
		ifstream LEEFT_In(Input_File_Names.at(a).c_str());
		if(LEEFT_In.fail())
		{
			cout<<Input_File_Names.at(a)<<" did not open."<<endl;
			assert(!LEEFT_In.fail());
		}//Close assertion.

		Read_Front_Info(LEEFT_In);
		while(!LEEFT_In.eof())
		{
			if(LEEFT_In.peek()=='T')
				break;

			vector<string> LEEFT_Data = Read_LEEFT_Data(LEEFT_In);
			LEEFT New_LEEFT(LEEFT_Data);
			if(Has_SU5_Group(New_LEEFT))
				Unique_LEEFTs.insert(New_LEEFT);

			Skip_Model_Inputs(LEEFT_In, Skipped_Lines);
		}//Close while loop on LEEFT_In.
	}//Close for loop on the data files.

	//Get statistics.
	cout<<"Computing statistics."<<endl;
	Statistics Doug_Moore;
	//Gauge group factors.
	map<Gauge_Group_Name, int> Gauge_Groups = 
		Doug_Moore.Count_Gauge_Groups(Unique_LEEFTs);
	//ST SUSYs.
	map<int, int> ST_SUSYs = 
		Doug_Moore.Count_ST_SUSYs(Unique_LEEFTs);
	//Gauge group combinations.
	map<vector<Gauge_Group_Name>, int> Gauge_Group_Combinations = 
		Doug_Moore.Count_Gauge_Group_Combinations(Unique_LEEFTs);
	//U(1) factors.
	map<int, int> U1_Factors = 
		Doug_Moore.Count_U1_Factors(Unique_LEEFTs);
	//Gauge group factors.
	map<int, int> Gauge_Group_Factors = 
		Doug_Moore.Count_Gauge_Group_Factors(Unique_LEEFTs);
	//Chiral matter generations.
	map<int, int> HS_Charged_Generations = 
		Find_HS_Charged_Generations(Unique_LEEFTs);
	map<int, int> HS_Uncharged_Generations = 
		Find_HS_Uncharged_Generations(Unique_LEEFTs);
	//NA Singlets.
	map<int, int> NA_Singlets = 
		Doug_Moore.Count_NA_Singlets(Unique_LEEFTs);
	//OS Charged exotics
	map<int, int> OS_Charged_Exotics = 
		Find_OS_Charged_Exotics(Unique_LEEFTs);


	//Write statistics to a Latex document.
	cout<<"Writing statistics."<<endl;
	Latex_Writer Jared_Greenwald;
	ofstream Latex_Out(Output_File_Name.c_str());
	Jared_Greenwald.Begin_Latex_Document(Latex_Out, Doc_Title);
	Latex_Out<<"\\textbf{Unique Models: "<<Unique_LEEFTs.size()<<"}\\\\";
	Latex_Out<<endl;

	//Gauge group factors.
	Jared_Greenwald.Write_Gauge_Group_Count(Latex_Out, Gauge_Groups, 
			Unique_LEEFTs.size());
	Latex_Out<<"\\vspace{3 mm}"<<endl;

	//Gauge group combinations.
	Jared_Greenwald.Write_Gauge_Group_Product_Count(Latex_Out, 
			Gauge_Group_Combinations, Unique_LEEFTs.size());
	Latex_Out<<"\\vspace{3 mm}"<<endl;

	//ST SUSYs.
	Jared_Greenwald.Write_ST_SUSY_Count(Latex_Out, ST_SUSYs,
			Unique_LEEFTs.size());
	Latex_Out<<"\\vspace{3 mm}"<<endl;
	Jared_Greenwald.Plot_ST_SUSY_Count(Latex_Out, ST_SUSYs, 
			Large_ST_Dimensions);
	Latex_Out<<"\\vspace{3 mm}"<<endl;

	//U(1) factors.
	Jared_Greenwald.Write_U1_Factor_Count(Latex_Out, U1_Factors,
			Unique_LEEFTs.size());
	Latex_Out<<"\\vspace{3 mm}"<<endl;
	Jared_Greenwald.Plot_U1_Factor_Count(Latex_Out, U1_Factors);
	Latex_Out<<"\\vspace{3 mm}"<<endl;

	//Number of gauge group factors.
	Jared_Greenwald.Write_Gauge_Group_Factor_Count(Latex_Out,
			Gauge_Group_Factors, Unique_LEEFTs.size());
	Latex_Out<<"\\vspace{3 mm}"<<endl;
	Jared_Greenwald.Plot_Gauge_Group_Factor_Count(Latex_Out,
			Gauge_Group_Factors);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;
	
	//Chiral matter generations.
	Latex_Out<<"\\textbf{Hidden sector charged generations.}"<<endl;
	Latex_Out<<"\\vspace{3 mm}\\\\"<<std::endl;
	Jared_Greenwald.Write_Matter_Generations(Latex_Out, 
			HS_Charged_Generations);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;
	Jared_Greenwald.Plot_Matter_Generations(Latex_Out,
			HS_Charged_Generations);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;
	Latex_Out<<"\\textbf{Hidden sector uncharged generations.}"<<endl;
	Latex_Out<<"\\vspace{3 mm}\\\\"<<std::endl;
	Jared_Greenwald.Write_Matter_Generations(Latex_Out,
			HS_Uncharged_Generations);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;
	Jared_Greenwald.Plot_Matter_Generations(Latex_Out, 
			HS_Uncharged_Generations);

	//NA Singlets.
	Jared_Greenwald.Write_NA_Singlets(Latex_Out, NA_Singlets);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<std::endl;
	Jared_Greenwald.Plot_NA_Singlets(Latex_Out, NA_Singlets);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<std::endl;

	//Observable sector charged exotics.
	Latex_Out<<"\\textbf{Observable sector charged exotics.}"<<endl;
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;
	Jared_Greenwald.Write_OS_Charged_Exotics(Latex_Out, 
			OS_Charged_Exotics);
	Latex_Out<<"\\hspace{\\fill}\\vspace{3 mm}\\\\"<<endl;
	Jared_Greenwald.Plot_OS_Charged_Exotics(Latex_Out,
			OS_Charged_Exotics);

	//Finish document.
	Jared_Greenwald.Finish_Latex_Document(Latex_Out);

	cout<<"Statistics written for "<<Doc_Title<<endl;
	cout<<"Statistics written to "<<Output_File_Name<<endl;
	cout<<"Program finished."<<endl;
	return 0;
}//Close main.

void Read_Front_Info(ifstream& LEEFT_In)
{
	string Front_Info_Line;
	getline(LEEFT_In, Front_Info_Line);
	while(Front_Info_Line.size() == 0 || Front_Info_Line.at(0) != '-')
		getline(LEEFT_In, Front_Info_Line);
}//Close Read_Front_Info.

vector<string> Read_LEEFT_Data(ifstream& LEEFT_In)
{
	vector<string> LEEFT_Data;
	string LEEFT_Line;
	getline(LEEFT_In, LEEFT_Line);
	LEEFT_Data.push_back(LEEFT_Line);
	while(LEEFT_Line.at(0) != 'S')
	{
		getline(LEEFT_In, LEEFT_Line);
		LEEFT_Data.push_back(LEEFT_Line);
	}//Close while loop.
	return LEEFT_Data;
}//Close Read_LEEFT_Data.

void Skip_Model_Inputs(ifstream& LEEFT_In, int Skipped_Lines)
{
	string junk;
	for(int a=0; a<Skipped_Lines; ++a)
		getline(LEEFT_In, junk);
}//Close Skip_Model_Inputs.

bool Has_SU5_Group(const LEEFT& LEEFT_Input)
{
	//Check the U(1) factors first.
	if(LEEFT_Input.U1_Factors() == 0)
		return false;

	for(int a=0; a<LEEFT_Input.Gauge_Groups().size(); ++a)
	{
		if((LEEFT_Input.Gauge_Groups().at(a).Class() == 'A') &&
				(LEEFT_Input.Gauge_Groups().at(a).Rank() == 4))
			return true;
	}
	return false;
}//Close Has_SU5_Group.

int Count_HS_Charged_Generations(const LEEFT& LEEFT_Input, int SU5_Index)
{
	//Count the number of 10s, with hidden sector charges less than 
	//the rank.
	int Number_of_10s = 0;//For unbarred generations.
	int Number_of_b5s = 0;//For unbarred generations.
	int Number_of_b10s = 0;//For barred generaitons.
	int Number_of_5s = 0;//For bared generations.
	for(int a=0; a<LEEFT_Input.Matter_Representations().size(); ++a)
	{
		int HS_Multiplicity = Get_HS_Multiplicity
			(LEEFT_Input.Matter_Representations().at(a).Rep_Dimension(),	
			 SU5_Index, LEEFT_Input.Gauge_Groups());

		switch(LEEFT_Input.Matter_Representations().at(a).
				Rep_Dimension().at(SU5_Index).Dimension())
		{
			case 10:
				Number_of_10s += LEEFT_Input.Matter_Representations().at(a).
					Duplicates()*HS_Multiplicity;
				break;
			case -10:
				Number_of_b10s += LEEFT_Input.Matter_Representations().at(a).
					Duplicates()*HS_Multiplicity;
				break;
			case 5:
				Number_of_5s += LEEFT_Input.Matter_Representations().at(a).
					Duplicates()*HS_Multiplicity;
				break;
			case -5:
				Number_of_b5s += LEEFT_Input.Matter_Representations().at(a).
					Duplicates()*HS_Multiplicity;
				break;
		};
	}//Close for loop on representations.
	int Unbarred_Generations = min(Number_of_10s, Number_of_b5s);
	int Barred_Generations = min(Number_of_b10s, Number_of_5s);

	return abs(Unbarred_Generations - Barred_Generations);
}//Close Count_Generations.

map<int, int> Find_HS_Charged_Generations(const set<LEEFT>& Unique_LEEFTs)
{
	map<int, int> HS_Charged_Generations;
	set<LEEFT>::iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
	for(itUnique_LEEFTs; itUnique_LEEFTs != Unique_LEEFTs.end(); 
			++itUnique_LEEFTs)
	{
		//First, find all of the indices of the SU(5) groups.
		vector<int> SU5_Indices;
		for(int a=0; a<itUnique_LEEFTs->Gauge_Groups().size(); ++a)
		{
			if((itUnique_LEEFTs->Gauge_Groups().at(a).Class() == 'A') && 
					(itUnique_LEEFTs->Gauge_Groups().at(a).Rank() == 4))
				SU5_Indices.push_back(a);
		}//Close for loop on the gauge groups of *itUnique_LEEFTs.

		//Loop over the SU5() groups, count the generations, 
		//and add them to the map.
		for(int a=0; a<SU5_Indices.size(); ++a)
		{
			int Number_HS_Charged_Generations = Count_HS_Charged_Generations
				(*itUnique_LEEFTs, SU5_Indices.at(a));
			if(HS_Charged_Generations.find(Number_HS_Charged_Generations) !=
					HS_Charged_Generations.end())
				HS_Charged_Generations[Number_HS_Charged_Generations]++;
			else
				HS_Charged_Generations[Number_HS_Charged_Generations] = 1;
			if(Number_HS_Charged_Generations == 3)
				itUnique_LEEFTs->Display();
		}//Close for loop on the different SU(5) indices.
	}//Close for loop on Unique_LEEFTs.
	return HS_Charged_Generations;
}//Close HS_Charged_Generations.

int Count_HS_Uncharged_Generations(const LEEFT& LEEFT_Input, int SU5_Index)
{
	int Number_of_10s = 0;//For unbarred generations.
	int Number_of_b5s = 0;//For unbarred generations.
	int Number_of_b10s = 0;//For barred generations.
	int Number_of_5s = 0;//For barred generations.
	for(int a=0; a<LEEFT_Input.Matter_Representations().size(); ++a)
	{
		if(!Has_Hidden_Sector_Charge
				(LEEFT_Input.Matter_Representations().at(a).Rep_Dimension(), 
				 SU5_Index))
		{
			switch(LEEFT_Input.Matter_Representations().at(a).Rep_Dimension().
					at(SU5_Index).Dimension())
			{
				case 10:
					Number_of_10s += LEEFT_Input.Matter_Representations().at(a).
						Duplicates();
					break;
				case -10:
					Number_of_b10s += LEEFT_Input.Matter_Representations().at(a).
						Duplicates();
					break;
				case 5:
					Number_of_5s += LEEFT_Input.Matter_Representations().at(a).
						Duplicates();
					break;
				case -5: 
					Number_of_b5s += LEEFT_Input.Matter_Representations().at(a).
						Duplicates();
					break;
			};//Close switch.
		}//Close if statement on hidden sector charges.
	}//Close for loop on Representations.
	int Unbarred_Generations = min(Number_of_10s, Number_of_b5s);
	int Barred_Generations = min(Number_of_b10s, Number_of_5s);
	return abs(Unbarred_Generations - Barred_Generations);
}//Close Count_HS_Uncharged_Generations.

map<int, int> Find_HS_Uncharged_Generations
	(const set<LEEFT>& Unique_LEEFTs)
{
	map<int, int> HS_Uncharged_Generations;
	set<LEEFT>::const_iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
	for(itUnique_LEEFTs; itUnique_LEEFTs != Unique_LEEFTs.end();
			++itUnique_LEEFTs)
	{
		//First, find all of the indices of the SO(10) groups.
		vector<int> SU5_Indices;
		for(int a=0; a<itUnique_LEEFTs->Gauge_Groups().size(); ++a)
		{
			if((itUnique_LEEFTs->Gauge_Groups().at(a).Class() == 'A') &&
					(itUnique_LEEFTs->Gauge_Groups().at(a).Rank() == 4))
				SU5_Indices.push_back(a);
		}//Close for loop for finding all of the SU(5) indices.

		//Loop over the SU(5) groups, count the generations, and 
		//add them to the map.
		for(int a=0; a<SU5_Indices.size(); ++a)
		{
			int Number_HS_Uncharged_Generations = 
				Count_HS_Uncharged_Generations
					(*itUnique_LEEFTs, SU5_Indices.at(a));
			if(HS_Uncharged_Generations.find(Number_HS_Uncharged_Generations) !=
					HS_Uncharged_Generations.end())
				HS_Uncharged_Generations[Number_HS_Uncharged_Generations]++;
			else
				HS_Uncharged_Generations[Number_HS_Uncharged_Generations] = 1;
		}//Close for loop on the different SU(5) indices.
	}//Close for loop on Unique_LEEFTs.
	return HS_Uncharged_Generations;
}//Close HS_Uncharged_Generations.

bool Has_Hidden_Sector_Charge(const vector<Group_Representation>& 
		Representations, int SU5_Index)
{
	for(int a=0; a<Representations.size(); ++a)
	{
		if((a != SU5_Index) && (Representations.at(a).Dimension() != 1))
			return true;
	}
	return false;
}//Close Has_Hidden_Sector_Charge.

int Count_OS_Charged_Exotics(const LEEFT& LEEFT_Input, int SU5_Index)
{
	int OS_Charged_Exotics = 0;
	for(int a=0; a<LEEFT_Input.Matter_Representations().size(); ++a)
	{
		if((LEEFT_Input.Matter_Representations().at(a).Rep_Dimension().
					at(SU5_Index).Dimension() != 1) &&
				(Has_Hidden_Sector_Charge(LEEFT_Input.Matter_Representations().
																	at(a).Rep_Dimension(), SU5_Index)))
			OS_Charged_Exotics += LEEFT_Input.Matter_Representations().at(a).
				Duplicates();
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

int Get_HS_Multiplicity(const vector<Group_Representation>& 
		Representations, int SU5_Index, 
		const vector<Gauge_Group_Name>& Gauge_Groups)
{
	int HS_Multiplicity = 1;
	for(int a=0; a<Representations.size(); ++a)
	{
		if(a != SU5_Index)
		{
			HS_Multiplicity*=abs(Representations.at(a).Dimension());
		}//Close if statement on SU5 index.
	}//Close for loop on Representations.
	return HS_Multiplicity;
}//Close Get_HS_Multiplicity.
