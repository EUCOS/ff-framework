#include "StatisticalReportGenerator.h"

StatisticalReportGenerator::StatisticalReportGenerator
(const std::set<LEEFT>& Unique_LEEFTs)
{
	Unique_LEEFTs_ = Unique_LEEFTs;
}//Close constructor.

StatisticalReportGenerator::StatisticalReportGenerator
(const StatisticalReportGenerator& New_StatisticalReportGenerator)
{
	Unique_LEEFTs_ = New_StatisticalReportGenerator.Unique_LEEFTs();
}//Close copy constructor.

//INTERFACE.
void StatisticalReportGenerator::Generate_Statistical_Report
	(std::ofstream& Latex_Out, int Large_ST_Dimensions)
{
	//Gauge groups.
	std::map<GaugeGroupName, int> Gauge_Groups = 
		Stats_.Count_Gauge_Groups(Unique_LEEFTs());
	Writer_.Write_Gauge_Group_Count(Latex_Out, Gauge_Groups, 
			Unique_LEEFTs().size());
	Latex_Out<<"\\vspace{3 mm}"<<std::endl;
	
	//Gauge group products.
	std::map<std::vector<GaugeGroupName>, int> Gauge_Group_Combinations = 
		Stats_.Count_Gauge_Group_Combinations(Unique_LEEFTs());
	Writer_.Write_Gauge_Group_Product_Count(Latex_Out, 
			Gauge_Group_Combinations, Unique_LEEFTs().size());
	Latex_Out<<"\\vspace{3 mm}"<<std::endl;

	//Gauge group factors.
	std::map<int, int> Gauge_Group_Factors = 
		Stats_.Count_Gauge_Group_Factors(Unique_LEEFTs());
	Writer_.Write_Gauge_Group_Factor_Count(Latex_Out, Gauge_Group_Factors,
			Unique_LEEFTs().size());
	Latex_Out<<"\\vspace{3 mm}"<<std::endl;
	Writer_.Plot_Gauge_Group_Factor_Count(Latex_Out, Gauge_Group_Factors);
	Latex_Out<<"\\vspace{3 mm}"<<std::endl;

	//U(1) factors.
	std::map<int, int> U1_Factors = 
		Stats_.Count_U1_Factors(Unique_LEEFTs());
	Writer_.Write_U1_Factor_Count(Latex_Out, U1_Factors,
			Unique_LEEFTs().size());
	Latex_Out<<"\\vspace{3 mm}"<<std::endl;
	Writer_.Plot_U1_Factor_Count(Latex_Out, U1_Factors);
	Latex_Out<<"\\vspace{3 mm}"<<std::endl;

	//ST SUSYs.
	std::map<int, int> ST_SUSYs = 
		Stats_.Count_ST_SUSYs(Unique_LEEFTs());
	Writer_.Write_ST_SUSY_Count(Latex_Out, ST_SUSYs, Unique_LEEFTs().size());
	Latex_Out<<"\\vspace{3 mm}"<<std::endl;
	Writer_.Plot_ST_SUSY_Count(Latex_Out, ST_SUSYs, Large_ST_Dimensions);
	Latex_Out<<"\\vspace{3 mm}"<<std::endl;

	//NA Singlets.
	std::map<int, int> NA_Singlets = 
		Stats_.Count_NA_Singlets(Unique_LEEFTs());
	Writer_.Write_NA_Singlets(Latex_Out, NA_Singlets);
	Latex_Out<<"\\vspace{3 mm}"<<std::endl;
	Writer_.Plot_NA_Singlets(Latex_Out, NA_Singlets);
	Latex_Out<<"\\vspace{3 mm}"<<std::endl;

}//Close Generate_Statistical_Report.
