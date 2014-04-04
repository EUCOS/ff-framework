#include "SystematicD10Builder.h"

SystematicD10Builder::SystematicD10Builder(const ModelBuilder& Initial_Build,
					       int Large_ST_Dimensions,
					       const std::vector<int>& 
					       Extension_Orders,
					       const std::string& 
					       Output_File_Name):
  SystematicBuilder(Initial_Build, Large_ST_Dimensions, Extension_Orders,
		     Output_File_Name)
{
  ;
}//Close constructor.

SystematicD10Builder::SystematicD10Builder(const SystematicD10Builder&
					       New_SystematicD10Builder):
  SystematicBuilder(New_SystematicD10Builder)
{
  ;
}//Close copy constructor.

//INTERFACE
void SystematicD10Builder::Perform_Search()
{
Data_Writer_.Write_BV_Search_Front_Info(Model_Out_, Initial_Model(), 
					  Large_ST_Dimensions(),
	    				  Extension_Orders());
  if(k_ij_Extensions().size() == 0)//In case the user wants k_ij fixed.
    Build_k_ij_Extensions();

  std::vector<BasisVector> Basis_Vectors(Extension_Orders().size(),
					  BasisVector());
  double Start_Time = clock();
  Build_Extensions(0, Basis_Vectors);
  double Finish_Time = clock();
  double Total_Time = (Finish_Time - Start_Time)/double(CLOCKS_PER_SEC);
  Write_End_MD(Total_Time);
}//Close Perform_Search.

//PRIVATE.
//HELPERS.
void SystematicD10Builder::Build_Extensions(int Layer, std::vector<BasisVector> 
					      Basis_Vectors)
{
 if(Layer<static_cast<int>(Extension_Orders().size()))
    Build_Basis_Vectors(Layer, Basis_Vectors);
  else
    Build_Models(Basis_Vectors);
}//Close Build_Extensions.

void SystematicD10Builder::Build_Basis_Vectors(int Layer, std::vector<BasisVector>
						 Basis_Vectors)
{
 MIBVGenerator BV_Generator(Extension_Orders().at(Layer), Large_ST_Dimensions());
  if(Layer == 0)
    BV_Generator.Build_Full_Chunks(Initial_Common_Basis_Alphas());
  else
    {
      std::vector<BasisVector> All_Basis_Vectors = Initial_Model().BV_Set();
      for(int a=0; a<Layer; ++a)
	All_Basis_Vectors.push_back(Basis_Vectors.at(a));
      BasisAlphaBuilder BA_Builder(All_Basis_Vectors);
      BA_Builder.Build_Basis_Alphas();
      BA_Builder.Build_Common_Basis_Alphas();
      BV_Generator.Build_Full_Chunks(BA_Builder.Common_Basis_Alphas());
    }

  std::list<Chunk> LMs = BV_Generator.SP_LMs();
  std::list<Chunk> Obs = BV_Generator.RMs_Observable();
  std::list<Chunk> Hid = BV_Generator.RMs_Hidden();

  std::list<Chunk>::iterator itLMs = LMs.begin();
  std::list<Chunk>::iterator itObs = Obs.begin();
  std::list<Chunk>::iterator itHid = Hid.begin();

  std::list<Chunk>::iterator itLMs_End = LMs.end();
  std::list<Chunk>::iterator itObs_End = Obs.end();
  std::list<Chunk>::iterator itHid_End = Hid.end();

  Total_BVs_ += LMs.size()*Obs.size()*Hid.size();

  ChunkConsistencyChecker Checker(FF::LCM(2, Extension_Orders().at(Layer)),
				    Initial_Common_Basis_Alphas().
				    at(0).Denominator());

  for(; itLMs != itLMs_End; ++itLMs)
    {
      for(itObs = Obs.begin(); itObs != itObs_End; ++itObs)
	{
	  for(itHid = Hid.begin(); itHid != itHid_End; ++itHid)
	    {
	      if(Checker.Check_D10_Modular_Invariance(*itLMs, *itObs, *itHid))
		{
		  std::vector<int> New_BV;
		  for(int a=0; a<static_cast<int>(itLMs->BV_Chunk().size()); ++a)
		    New_BV.push_back(itLMs->BV_Chunk().at(a));
		  for(int a=0; a<static_cast<int>(itObs->BV_Chunk().size()); ++a)
		    New_BV.push_back(itObs->BV_Chunk().at(a));
		  for(int a=0; a<static_cast<int>(itHid->BV_Chunk().size()); ++a)
		    New_BV.push_back(itHid->BV_Chunk().at(a));

		  Basis_Vectors.at(Layer) = 
		    (BasisVector(New_BV, Extension_Orders().at(Layer), 
				  Large_ST_Dimensions()));
		  if(Layer == static_cast<int>(Extension_Orders().size())-1)
		    Consistent_BVs_++;
		  Build_Extensions(Layer+1, Basis_Vectors);
		}//Close if statement on modular invariance.
	    }//Close for loop on itHid.
	}//Close for loop on itObs.
    }//Close for loop on LMs.
}//Close Build_Basis_Vectors.

