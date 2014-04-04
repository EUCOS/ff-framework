/*
	Authors: Timothy Renner (renner.timothy@gmail.com)
				   Doug Moore (Douglas_Moore1@baylor.edu)
					 Gerald Cleaver, Ph.D. (Gerald_Cleaver@baylor.edu)
	Date: 5/17/2011
	Baylor University

	This program gathers relevant statistical information for models
	containing at least one Left-Right Symmetric group. 
	Tabulates gauge group factors,
	number of ST SUSYs, number of chiral matter representations, etc.

	Uses the FF Framework
*/

#include <string>
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

bool Has_LRSym_Group(const LEEFT& LEEFT_Input);

int Count_Generations(const LEEFT& LEEFT_Input, 
		int SU3_Index, int SU2A_Index, int SU2B_Index);
map<int, int> Find_Generations
	(const set<LEEFT>& Unique_LEEFTs);

int Count_OS_Charged_Exotics(const LEEFT& LEEFT_Input,
		int SU3_Index, int SU2A_Index, int SU2B_Index);
map<int, int> Find_OS_Charged_Exotics(const set<LEEFT>& Unique_LEEFTs);

bool Has_Hidden_Sector_Charge(const vector<Group_Representation>& Representations,
		int SU3_Index, int SU2A_Index, int SU2B_Index);

int Get_HS_Multiplicity(const vector<Group_Representation>& 
		Representations, int SU3_Index, int SU2A_Index, int SU2B_Index,
		const vector<Gauge_Group_Name>& Gauge_Groups);
int main()
{
	for(int a=0; a<10; ++a)
		cout<<"MODIFIED TO COMPUTE HIDDEN SECTOR MULTIPLICITIES."<<endl;
	string FF_Data = "/Users/timothyrenner/FF_Framework/data/";
	int Large_ST_Dimensions = 4;
	//Set the output file name.
	string Output_File_Name = "./tex/LRSym_Model_Statistics.tex";

	//Set the Latex document title.
	string Doc_Title = 
		"Statistical Report for NAHE Based Left-Right Symmetric Models," ;
	Doc_Title+=" Without S Vector, ";
	Doc_Title += "Order 3, Layer 1";

	//Load the file names.
	vector<string> Input_File_Names;
	Input_File_Names.push_back
		(FF_Data + "NAHE_Extensions/One_Layer/NAHE_NSV_O3_L1.txt");
	//Input_File_Names.push_back
		//(FF_Data + "NAHE_Extensions/One_Layer/NR_NAHE_O2_L1_NSP.txt");

	//int Skipped_Lines = 2;//For no layers.
	//int Skipped_Lines = 3; //For the NR_NAHE data set.
	int Skipped_Lines = 4;//For one layer.

	cout<<"Reading unique LEEFTs containing Left_Right Symmetric gauge groups."<<endl;
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
			if(Has_LRSym_Group(New_LEEFT))
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
	map<int, int> Generations = 
		Find_Generations(Unique_LEEFTs);
	//NA Singlets.
	map<int, int> NA_Singlets = 
		Doug_Moore.Count_NA_Singlets(Unique_LEEFTs);
	//Observable sector charged exotics.
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
	Latex_Out<<"\\textbf{(3,2,1) generations.}"<<endl;
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;
	Jared_Greenwald.Write_Matter_Generations(Latex_Out,
			Generations);
	Latex_Out<<"\\vspace{3 mm}"<<endl;
	Jared_Greenwald.Plot_Matter_Generations(Latex_Out, 
			Generations);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<endl;

	//NA Singlets.
	Jared_Greenwald.Write_NA_Singlets(Latex_Out, NA_Singlets);
	Latex_Out<<"\\vspace{3 mm}\\\\"<<std::endl;
	Jared_Greenwald.Plot_NA_Singlets(Latex_Out, NA_Singlets);

	//Observable sector charged exotics.
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

bool Has_LRSym_Group(const LEEFT& LEEFT_Input)
{
	bool Has_SU3 = false;
	bool Has_2_SU2 = false;
	int SU2_Count = 0;
	for(int a=0; a<LEEFT_Input.Gauge_Groups().size(); ++a)
	{
		if(!Has_SU3 && (LEEFT_Input.Gauge_Groups().at(a).Class() == 'A')
				&& (LEEFT_Input.Gauge_Groups().at(a).Rank() == 2))
			Has_SU3 = true;

		if(!Has_2_SU2 && (LEEFT_Input.Gauge_Groups().at(a).Class() == 'A')
				&& (LEEFT_Input.Gauge_Groups().at(a).Rank() == 1))
			SU2_Count++;

		if(SU2_Count == 2)
			Has_2_SU2 = true;

		if(Has_SU3 && Has_2_SU2)
			return true;
	}
	return false;
}//Close Has_LRSym_Group.

int Count_Generations(const LEEFT& LEEFT_Input, int SU3_Index, 
		int SU2A_Index, int SU2B_Index)
{
	int Unbarred_Generations = 0;
	int Barred_Generations = 0;
	for(int a=0; a<LEEFT_Input.Matter_Representations().size(); ++a)
	{
		//if(!Has_Hidden_Sector_Charge(LEEFT_Input.Matter_Representations().
			//		at(a).Rep_Dimension(), SU3_Index, SU2A_Index, SU2B_Index))
		//{
			if((LEEFT_Input.Matter_Representations().at(a).Rep_Dimension().
						at(SU3_Index).Dimension() == -3) && 
					(LEEFT_Input.Matter_Representations().at(a).Rep_Dimension().
					 at(SU2A_Index).Dimension() == 2) &&
					(LEEFT_Input.Matter_Representations().at(a).Rep_Dimension().
					 at(SU2B_Index).Dimension() == 1))
			{
				int Current_Barred_Generations = 0;
				Current_Barred_Generations += 
					LEEFT_Input.Matter_Representations().at(a).Duplicates()*
					Get_HS_Multiplicity(LEEFT_Input.Matter_Representations().at(a).
							Rep_Dimension(), SU3_Index, SU2A_Index, SU2B_Index,
							LEEFT_Input.Gauge_Groups());

				Barred_Generations += Current_Barred_Generations;

			}else if((LEEFT_Input.Matter_Representations().at(a).Rep_Dimension().
						at(SU3_Index).Dimension() == 3) &&
					(LEEFT_Input.Matter_Representations().at(a).Rep_Dimension().
					 at(SU2A_Index).Dimension() == 2) &&
					(LEEFT_Input.Matter_Representations().at(a).Rep_Dimension().
					 at(SU2B_Index).Dimension() == 1))
			{
				int Current_Unbarred_Generations = 0;
				Current_Unbarred_Generations +=
					LEEFT_Input.Matter_Representations().at(a).Duplicates()*
					Get_HS_Multiplicity(LEEFT_Input.Matter_Representations().at(a).
							Rep_Dimension(), SU3_Index, SU2A_Index, SU2B_Index,
							LEEFT_Input.Gauge_Groups());
				Unbarred_Generations += Current_Unbarred_Generations;

			}//Close if/else for counting the barred and unbarred
			//representations.
		//}//Close if statement on hidden sector charge.
	}//Close for loop on Representations.
	return abs(Unbarred_Generations - Barred_Generations);
}//Close Count_Generations.

map<int, int> Find_Generations(const set<LEEFT>& Unique_LEEFTs)
{
	map<int, int> Generations;
	set<LEEFT>::iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
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
						int Number_Generations = Count_Generations
							(*itUnique_LEEFTs, SU3_Indices.at(a), SU2_Indices.at(b),
							 SU2_Indices.at(c));
						if(Number_Generations == 3)
							itUnique_LEEFTs->Display();
						if(Generations.find(Number_Generations) !=
								Generations.end())
							Generations[Number_Generations]++;
						else
							Generations[Number_Generations] = 1;
					}//Close if statement on b != c.
				}//Close for loop on the second SU(2) indices.
			}//Close for loop on the first SU(2) indices.
		}//Close for loop on the different SU(3) indices.
	}//Close for loop on Unique_LEEFTs.
	return Generations;
}//Close Generations.

int Count_OS_Charged_Exotics(const LEEFT& LEEFT_Input, int SU3_Index,
		int SU2A_Index, int SU2B_Index)
{
	int OS_Charged_Exotics = 0;
	for(int a=0; a<LEEFT_Input.Matter_Representations().size(); ++a)
	{
		if(((LEEFT_Input.Matter_Representations().at(a).Rep_Dimension().
						at(SU3_Index).Dimension() != 1) ||
				(LEEFT_Input.Matter_Representations().at(a).Rep_Dimension().
				 at(SU2A_Index).Dimension() != 1) ||
				(LEEFT_Input.Matter_Representations().at(a).Rep_Dimension().
				 at(SU2B_Index).Dimension() != 1)) && 
				(Has_Hidden_Sector_Charge(LEEFT_Input.Matter_Representations().
																	at(a).Rep_Dimension(),
																	SU3_Index, SU2A_Index, SU2B_Index)))
			OS_Charged_Exotics += LEEFT_Input.Matter_Representations().
				at(a).Duplicates();
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

bool Has_Hidden_Sector_Charge(const vector<Group_Representation>& Representations, 
		int SU3_Index, int SU2A_Index, int SU2B_Index)
{
	for(int a=0; a<Representations.size(); ++a)
	{
		if((a != SU3_Index) && (a != SU2A_Index) && 
				(a != SU2B_Index) && (Representations.at(a).Dimension() != 1))
			return true;
	}
	return false;
}//Close Has_Hidden_Sector_Charge.

int Get_HS_Multiplicity(const vector<Group_Representation>& 
		Representations, int SU3_Index, int SU2A_Index, int SU2B_Index, 
		const vector<Gauge_Group_Name>& Gauge_Groups)
{
	int HS_Multiplicity = 1;
	for(int a=0; a<Representations.size(); ++a)
	{
		if((a != SU3_Index)&&(a != SU2A_Index) && (a != SU2B_Index))
		{
			HS_Multiplicity*=abs(Representations.at(a).Dimension());
		}//Close if statement on SU5 index.
	}//Close for loop on Representations.
	return HS_Multiplicity;
}//Close Get_HS_Multiplicity.
