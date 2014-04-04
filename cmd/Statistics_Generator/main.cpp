/*
Author: Timothy Renner (renner.timothy@gmail.com)
				Doug Moore (Douglas_Moore1@baylor.edu)
				Gerald Cleaver, Ph.D. (Gerald_Cleaver@baylor.edu)
Date: 11/4/2010
Baylor University

This program gathers statistics from a file by reading the LEEFTs into an STL set, 
then tabulating gauge group factors, number of ST SUSYs, etc. May be incorporated 
into a class at some point.
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
void Display_LEEFT_Data(const vector<string>& LEEFT_Data);
void Skip_Model_Inputs(ifstream& LEEFT_In, int Skipped_Lines);

int main()
{
  string FF_Data = "/Users/timothyrenner/FF_Framework/data/";
  int Large_ST_Dimensions = 4;
  //Set the output file name.
  string Output_File_Name = "./tex/Statistical_Report.tex";

  //Set the Latex Doc title
  string Doc_Title = 
  "Statistical Report for NAHE Extensions,";
	Doc_Title += " Without S vector, ";
  Doc_Title += "Order 3, Layer 1";
	//Doc_Title += " Simply Paired";

  //Load the files into the vector.
  vector<string> Input_File_Names;
  Input_File_Names.push_back(FF_Data+
  	     "NAHE_Extensions/One_Layer/NAHE_NSV_O3_L1.txt");
	//Input_File_Names.push_back(FF_Data+
		//		"NAHE_Extensions/One_Layer/NR_NAHE_NSV_O2_L1_NSP.txt");

  //Set input skip lines.
	//int Skipped_Lines = 2; //For no basis vectors.
	//int Skipped_Lines = 3; //For NAHE models.
	int Skipped_Lines = 4; //For one layer.
	//int Skipped_Lines = 6; //For two layers/

	//Load the data into Unique_LEEFTs.
	cout<<"Reading unique LEEFTs."<<endl;
	set<LEEFT> Unique_LEEFTs;
	for(int a=0; a<Input_File_Names.size(); a++)
	{
		ifstream LEEFT_In(Input_File_Names.at(a).c_str());
		if(LEEFT_In.fail())
		{
			cout<<Input_File_Names.at(a)<<" did not open."<<endl;
			assert(!LEEFT_In.fail());
		}
		Read_Front_Info(LEEFT_In);
		while(!LEEFT_In.eof())
		{
			if(LEEFT_In.peek() == 'T')
				break;
			vector<string> LEEFT_Data = Read_LEEFT_Data(LEEFT_In);
			LEEFT New_LEEFT(LEEFT_Data);
			//New_LEEFT.Display();
			Unique_LEEFTs.insert(New_LEEFT);
			Skip_Model_Inputs(LEEFT_In, Skipped_Lines);
		}//Close while loop on LEEFT_In.
	}//Close for loop on Input_File_Names.

	//Display the unique LEEFTs.
	/*cout<<"Unique LEEFTs:"<<endl;
	set<LEEFT>::const_iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
	set<LEEFT>::const_iterator itEnd = Unique_LEEFTs.end();
	for(itUnique_LEEFTs; itUnique_LEEFTs != itEnd; ++itUnique_LEEFTs)
		itUnique_LEEFTs->Display();
	cout<<endl;
	cout<<"Total Unique LEEFTs: "<<Unique_LEEFTs.size()<<endl;\\*/

	//Get Statistics.
	cout<<"Getting statistics."<<endl;
	Statistics Doug_Moore;
	//Gauge groups.
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
	//Non-Abelian singlets.
	map<int, int> NA_Singlets = 
		Doug_Moore.Count_NA_Singlets(Unique_LEEFTs);

  //Display Statistics.
	//Gauge groups.
  cout<<"Gauge Groups:"<<endl;
  map<Gauge_Group_Name, int>::const_iterator itGauge_Groups = 
    Gauge_Groups.begin();
  map<Gauge_Group_Name, int>::const_iterator itMapEnd = 
    Gauge_Groups.end();
  for(itGauge_Groups; itGauge_Groups != itMapEnd; 
      ++itGauge_Groups)
    {
      cout<<(itGauge_Groups->first).Class()<<" ";
      cout<<(itGauge_Groups->first).Rank()<<" ";
      cout<<(itGauge_Groups->first).KM_Level()<<" ";
      cout<<"--->"<<itGauge_Groups->second<<endl;
    }
  cout<<endl;
	//Gauge group combinations.
  cout<<"Gauge Group Combinations:"<<endl;
  map<vector<Gauge_Group_Name>, int>::const_iterator itGauge_Group_Combinations = 
    Gauge_Group_Combinations.begin();
  map<vector<Gauge_Group_Name>, int>::const_iterator itCombEnd = 
    Gauge_Group_Combinations.end();
  for(itGauge_Group_Combinations; itGauge_Group_Combinations != itCombEnd;
      ++itGauge_Group_Combinations)
    {
      cout<<(itGauge_Group_Combinations->first).at(0).Class()<<" ";
      cout<<(itGauge_Group_Combinations->first).at(0).Rank()<<" ";
      cout<<(itGauge_Group_Combinations->first).at(0).KM_Level()<<" x ";
      cout<<(itGauge_Group_Combinations->first).at(1).Class()<<" ";
      cout<<(itGauge_Group_Combinations->first).at(1).Rank()<<" ";
      cout<<(itGauge_Group_Combinations->first).at(1).KM_Level()<<" ";
      cout<<"--->"<<itGauge_Group_Combinations->second<<endl;
    }

	//ST SUSYs.
  cout<<"ST SUSYs: "<<endl;
  map<int, int>::const_iterator itST_SUSYs = ST_SUSYs.begin();
  map<int, int>::const_iterator itSSEnd = ST_SUSYs.end();
  for(itST_SUSYs; itST_SUSYs != itSSEnd; ++itST_SUSYs)
    {
      cout<<(itST_SUSYs->first)<<" ";
      cout<<"--->";
      cout<<(itST_SUSYs->second)<<endl;
    }
  cout<<endl;

	//U(1) factors.
  cout<<"U1 Factors: "<<endl;
  map<int, int>::const_iterator itU1_Factors = U1_Factors.begin();
  map<int, int>::const_iterator itU1_End = U1_Factors.end();
  for(itU1_Factors; itU1_Factors != itU1_End; ++itU1_Factors)
    {
      cout<<(itU1_Factors->first)<<" ";
      cout<<"--->";
      cout<<(itU1_Factors->second)<<endl;
    }
  cout<<endl;

	//Gauge group factors.
  cout<<"Gauge Group Factors: "<<endl;
  map<int, int>::const_iterator itGauge_Group_Factors = Gauge_Group_Factors.begin();
  map<int, int>::const_iterator itGGF_End = Gauge_Group_Factors.end();
  for(itGauge_Group_Factors; itGauge_Group_Factors != itGGF_End; 
      ++itGauge_Group_Factors)
    {
      cout<<(itGauge_Group_Factors->first)<<" ";
      cout<<"--->";
      cout<<(itGauge_Group_Factors->second)<<endl;
    }
  cout<<endl;

	//NA Singlets.
	cout<<"Non-Abelian singlets: "<<endl;
	map<int, int>::const_iterator itNA_Singlets = NA_Singlets.begin();
	for(itNA_Singlets; itNA_Singlets != NA_Singlets.end(); 
			++itNA_Singlets)
	{
		cout<<(itNA_Singlets->first)<<" ";
		cout<<"--->";
		cout<<(itNA_Singlets->second)<<endl;
	}//Close for loop.
	cout<<endl;

  //Write Statistics to a Latex document.
  cout<<"Writing statistics."<<endl;
  Latex_Writer Jared_Greenwald;
  ofstream Latex_Out(Output_File_Name.c_str());
	//Document front matter.
  Jared_Greenwald.Begin_Latex_Document(Latex_Out, Doc_Title);
  Latex_Out<<"\\textbf{Unique Models: ";
	Latex_Out<<Unique_LEEFTs.size()<<"}\\\\"<<endl;

	//Gauge group content.
  Jared_Greenwald.Write_Gauge_Group_Count(Latex_Out,
						 Gauge_Groups,
						 Unique_LEEFTs.size());
  Latex_Out<<"\\vspace{3 mm}"<<endl;

	//Gauge group products.
  Jared_Greenwald.Write_Gauge_Group_Product_Count(Latex_Out,
						  Gauge_Group_Combinations,
						  Unique_LEEFTs.size());
  //Latex_Out<<"\\\\"<<endl;
  Latex_Out<<"\\vspace{3 mm}"<<endl;

	//ST SUSYs.
  Jared_Greenwald.Write_ST_SUSY_Count(Latex_Out,
				       ST_SUSYs,
				       Unique_LEEFTs.size());
  Latex_Out<<"\\vspace{3 mm}"<<endl;
  Jared_Greenwald.Plot_ST_SUSY_Count(Latex_Out, ST_SUSYs, Large_ST_Dimensions);
  Latex_Out<<"\\vspace{3 mm}"<<endl;

	//U(1) factors.
  Jared_Greenwald.Write_U1_Factor_Count(Latex_Out, U1_Factors, Unique_LEEFTs.size());
  Latex_Out<<"\\vspace{3 mm}"<<endl;
  Jared_Greenwald.Plot_U1_Factor_Count(Latex_Out, U1_Factors);
  Latex_Out<<"\\vspace{3 mm}"<<endl;

	//Gauge group factors.
  Jared_Greenwald.Write_Gauge_Group_Factor_Count(Latex_Out, Gauge_Group_Factors,
						 Unique_LEEFTs.size());
  Latex_Out<<"\\vspace{3 mm}"<<endl;
  Jared_Greenwald.Plot_Gauge_Group_Factor_Count(Latex_Out, Gauge_Group_Factors);
	Latex_Out<<"\\vspace{3 mm} \\\\ "<<endl;

	//Non-Abelian singlets.
	Jared_Greenwald.Write_NA_Singlets(Latex_Out, NA_Singlets);
	Latex_Out<<"\\vspace{3 mm} \\\\ "<<endl;
	Jared_Greenwald.Plot_NA_Singlets(Latex_Out, NA_Singlets);

	//Finish Latex document.
  Jared_Greenwald.Finish_Latex_Document(Latex_Out);

  cout<<"Statistics written for "<<Doc_Title<<endl;
  cout<<"Statistics written to "<<Output_File_Name<<endl;
  return 0;
}//Close main.

void Read_Front_Info(ifstream& LEEFT_In)
{
  string Front_Info_Line;
  getline(LEEFT_In, Front_Info_Line);
  while(Front_Info_Line.size() == 0 || Front_Info_Line.at(0) != '-')
    {
      getline(LEEFT_In, Front_Info_Line);
    }
}//Close Read_Front_Info.

vector<string> Read_LEEFT_Data(ifstream& LEEFT_In)
{
  vector<string> LEEFT_Data;
  string LEEFT_Line;
  getline(LEEFT_In, LEEFT_Line);//Gauge Groups.
  LEEFT_Data.push_back(LEEFT_Line);
  while(LEEFT_Line.at(0) != 'S')
    {
      getline(LEEFT_In, LEEFT_Line);
      LEEFT_Data.push_back(LEEFT_Line);
    }//Close while loop.
  return LEEFT_Data;
}//Close Read_LEEFT_Data.

void Display_LEEFT_Data(const vector<string>& LEEFT_Data)
{
  for(int a=0; a<LEEFT_Data.size(); a++)
    cout<<LEEFT_Data.at(a)<<endl;
  cout<<endl;
}//Close Display_LEEFT_Data.

void Skip_Model_Inputs(ifstream& LEEFT_In, int Lines_To_Skip)
{
  string junk;
  for(int a=0; a<Lines_To_Skip; ++a)
    getline(LEEFT_In, junk);
}//Close Skip_Model_Inputs.
