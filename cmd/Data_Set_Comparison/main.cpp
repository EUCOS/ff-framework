/*
AUTHORS: Timothy Renner (renner.timothy@gmail.com)
				 Gerald Cleaver, Ph.D. (Gerald_Cleaver@Baylor.edu)
INSTITUTE: Baylor University

Takes two data sets, reads the unique LEEFTs, and compares 
the number of LEEFTs common to both sets.
*/

#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>

#include "FF_LEEFT_Tools.hh"
#include "FF_LEEFT.hh"

using namespace std;

int main()
{
	cout<<"Begin program."<<endl;
	int Layers = 1;

	//Set up the input files.
	string FF_Data = "/Users/timothyrenner/FF_Framework/data/";
	vector<string> Set1_File_Names;
	Set1_File_Names.push_back(FF_Data + 
			"/NAHE_Extensions/One_Layer/NAHE_O3_L1.txt");
	//Set1_File_Names.push_back(FF_Data + 
		//	"/NAHE_Extensions/One_Layer/NR_NAHE_O2_L1_NSP.txt");

	vector<string> Set2_File_Names;
	Set2_File_Names.push_back(FF_Data + 
			"/NAHE_Extensions/One_Layer/NAHE_NSV_O3_L1.txt");
	//Set2_File_Names.push_back(FF_Data + 
		//	"/NAHE_Extensions/One_Layer/NR_NAHE_NSV_O2_L1_NSP.txt");

	//Read in the LEEFTs.
	cout<<"Reading in LEEFTs."<<endl;
	set<LEEFT> Set1_LEEFTs;
	
	for(int a=0; a<Set1_File_Names.size(); ++a)
	{
		ifstream LEEFT_In(Set1_File_Names.at(a).c_str());
		if(LEEFT_In.fail())
		{
			cout<<Set1_File_Names.at(a)<<" did not open."<<endl;
			return 0;
		}//Close if statement on file opening.

		LEEFT_Tools::Read_Front_Info(LEEFT_In);
		while(LEEFT_In.peek() != 'T')
		{
			vector<string> LEEFT_Data = LEEFT_Tools::Read_LEEFT_Data(LEEFT_In);
			string blank;
			getline(LEEFT_In, blank);
			vector<string> Model_Inputs = 
				LEEFT_Tools::Read_Model_Inputs(LEEFT_In, Layers);
			//getline(LEEFT_In,blank); // For the NAHE_NR data sets.
			getline(LEEFT_In, blank);

			//Insert the LEEFT into the set.
			Set1_LEEFTs.insert(LEEFT(LEEFT_Data));
		}//Close while loop.
	}//Close for loop on Set1_File_Names.
	

	set<LEEFT> Set2_LEEFTs;
	
	for(int a=0; a<Set2_File_Names.size(); ++a)
	{
		ifstream LEEFT_In(Set2_File_Names.at(a).c_str());
		if(LEEFT_In.fail())
		{
			cout<<Set2_File_Names.at(a)<<" did not open."<<endl;
			return 0;
		}//Close if statement on file opening.

		LEEFT_Tools::Read_Front_Info(LEEFT_In);
		while(LEEFT_In.peek() != 'T')
		{
			vector<string> LEEFT_Data = LEEFT_Tools::Read_LEEFT_Data(LEEFT_In);
			string blank;
			getline(LEEFT_In, blank);
			vector<string> Model_Inputs = 
				LEEFT_Tools::Read_Model_Inputs(LEEFT_In, Layers);
			//getline(LEEFT_In, blank);//For the NAHE_NR data sets.
			getline(LEEFT_In, blank);

			//Insert the LEEFT into the set.
			Set2_LEEFTs.insert(LEEFT(LEEFT_Data));
		}//Close while loop.
	}//Close for loop.

		set<LEEFT> Intersection; 
			set_intersection(Set1_LEEFTs.begin(), Set1_LEEFTs.end(),
					Set2_LEEFTs.begin(), Set2_LEEFTs.end(),
					inserter(Intersection, Intersection.begin()));
	

	cout<<"Set 1 has: "<<Set1_LEEFTs.size()<<endl;
	cout<<"Set 2 has: "<<Set2_LEEFTs.size()<<endl;
	cout<<"Common to both sets: "<<Intersection.size()<<endl;
	cout<<"End program."<<endl;
	return 0;
}//Close main.
