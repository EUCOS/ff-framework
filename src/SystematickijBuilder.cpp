#include "SystematickijBuilder.h"

SystematickijBuilder::SystematickijBuilder(const ModelBuilder& Initial_Build,
						 int Large_ST_Dimensions,
						 const std::vector<int>& 
						 Extension_Orders,
						 const std::string& 
						 Output_File_Name,
						 const std::vector<BasisVector>&
						 Basis_Vectors):
  SystematicBuilder(Initial_Build, Large_ST_Dimensions, Extension_Orders, 
		     Output_File_Name)
{
  Basis_Vectors_ = Basis_Vectors;
}//Close constructor.

SystematickijBuilder::SystematickijBuilder(const SystematickijBuilder& 
						 New_SystematickijBuilder):
  SystematicBuilder(New_SystematickijBuilder)
{
  Basis_Vectors_ = New_SystematickijBuilder.Basis_Vectors();
}//Close copy constructor.

void SystematickijBuilder::Perform_Search()
{
  Data_Writer_.Write_BV_Search_Front_Info(Model_Out_, Initial_Model(),
					  Large_ST_Dimensions(), Extension_Orders());
  Build_k_ij_Extensions();
  Set_Consistent_BVs(1);
  Set_Total_BVs(Basis_Vectors().size());
  double Start_Time = clock();
  Build_Models(Basis_Vectors());
  double Finish_Time = clock();
  double Total_Time = (Finish_Time - Start_Time)/double(CLOCKS_PER_SEC);
  Write_End_MD(Total_Time);
}//Close Perform_Search.
