#include "LEEFTTools.h"

void LEEFTTools::Read_Front_Info(std::ifstream& LEEFT_In)
{
	std::string Front_Info_Line;
	getline(LEEFT_In, Front_Info_Line);
	while(Front_Info_Line.size() == 0 || Front_Info_Line.at(0) != '-')
		getline(LEEFT_In, Front_Info_Line);
}//Close Read_Front_Info.

std::vector<std::string> LEEFTTools::Read_LEEFT_Data
	(std::ifstream& LEEFT_In)
{
	std::vector<std::string> LEEFT_Data;
	std::string LEEFT_Line;
	getline(LEEFT_In, LEEFT_Line); //For the gauge groups.
	LEEFT_Data.push_back(LEEFT_Line);
	while(LEEFT_Line.at(0) != 'S')
	{
		getline(LEEFT_In, LEEFT_Line);
		LEEFT_Data.push_back(LEEFT_Line);
	}//Close while loop on reading in strings of LEEFT data.
	return LEEFT_Data;
}//Close Read_LEEFT_Data.

std::vector<std::string> LEEFTTools::Read_Model_Inputs
	(std::ifstream& LEEFT_In, int Layers)
{
	std::vector<std::string> Model_Inputs;
	std::string Model_Input_Line;
	for(int a=0; a<2*Layers; ++a)
	{
		getline(LEEFT_In, Model_Input_Line);
		Model_Inputs.push_back(Model_Input_Line);
	}//Close for loop on reading in the model inputs.
	return Model_Inputs;
}//Close Read_Model_Inputs.
