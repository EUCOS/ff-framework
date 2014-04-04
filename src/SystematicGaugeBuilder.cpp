#include "SystematicGaugeBuilder.h"

//CONSTRUCTORS.
SystematicGaugeBuilder::SystematicGaugeBuilder(const ModelBuilder& 
						   Initial_Build, 
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

SystematicGaugeBuilder::SystematicGaugeBuilder(const SystematicGaugeBuilder&
						   New_SystematicGaugeBuilder):
  SystematicBuilder(New_SystematicGaugeBuilder)
{
  ;
}//Close copy constructor.

//INTERFACE.
void SystematicGaugeBuilder::Perform_Search()
{
  Data_Writer_.Write_BV_Search_Front_Info(Model_Out_, Initial_Model(),
					  Large_ST_Dimensions(), Extension_Orders());
  if(k_ij_Extensions().size() == 0)//Check to see if k_ij is fixed.
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
void SystematicGaugeBuilder::Build_Extensions(int Layer, std::vector<BasisVector>
						Basis_Vectors)
{
  if(Layer<static_cast<int>(Extension_Orders().size()))
    Build_Basis_Vectors(Layer, Basis_Vectors);
  else
    {
      Build_Models(Basis_Vectors);
    }
}//Build_Extensions.

void SystematicGaugeBuilder::Build_Basis_Vectors(int Layer, 
						   std::vector<BasisVector> 
						   Basis_Vectors)
{
   MIBVGenerator BV_Generator(Extension_Orders().at(Layer), Large_ST_Dimensions());
 if(Layer == 0)
   BV_Generator.Build_Gauge_Chunks(Initial_Common_Basis_Alphas());
 else
   {
     std::vector<BasisVector> All_Basis_Vectors = Initial_Model().BV_Set();
     for(int a=0; a<Layer; ++a)
       All_Basis_Vectors.push_back(Basis_Vectors.at(a));
     BasisAlphaBuilder BA_Builder(All_Basis_Vectors);
     BA_Builder.Build_Basis_Alphas();
     BA_Builder.Build_Common_Basis_Alphas();
     BV_Generator.Build_Gauge_Chunks(BA_Builder.Common_Basis_Alphas());
   }

 
  std::list<Chunk> LMs;
  LMs = BV_Generator.SP_LMs();

  std::list<Chunk> Comp;
  Comp = BV_Generator.SP_RMs_Compact();
 

  std::list<Chunk> Obs = BV_Generator.RMs_Observable();
  std::list<Chunk> Hid = BV_Generator.RMs_Hidden();

  std::list<Chunk>::iterator itObs = Obs.begin();
  std::list<Chunk>::iterator itComp = Comp.begin();
  std::list<Chunk>::iterator itHid = Hid.begin();

  std::list<Chunk>::iterator itObs_End = Obs.end();
  std::list<Chunk>::iterator itComp_End = Comp.end();
  std::list<Chunk>::iterator itHid_End = Hid.end();

  Total_BVs_ += LMs.size()*Obs.size()*Comp.size()*Hid.size();

  ChunkConsistencyChecker Checker(Extension_Orders().at(Layer),
				    Initial_Common_Basis_Alphas().at(0).
				    Denominator());

  for(itComp = Comp.begin(); itComp != itComp_End; ++itComp)
    {
      for(itObs = Obs.begin(); itObs != itObs_End; ++itObs)
	{

	  for(itHid = Hid.begin(); itHid != itHid_End; ++itHid)
	    {
	      if(Checker.Check_Modular_Invariance(LMs.front(), *itObs, 
						  *itComp, *itHid))
		{
		  std::vector<int> New_BV;
		  for(int a=0; a<static_cast<int>(LMs.front().BV_Chunk().size()); ++a)
		    New_BV.push_back(LMs.front().BV_Chunk().at(a));
		  for(int a=0; a<static_cast<int>(itObs->BV_Chunk().size()); ++a)
		    New_BV.push_back(itObs->BV_Chunk().at(a));
		  for(int a=0; a<static_cast<int>(itComp->BV_Chunk().size()); ++a)
		    New_BV.push_back(itComp->BV_Chunk().at(a));
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
    }//Close for loop on itComp.
}//Close Build_Basis_Vectors.

void SystematicGaugeBuilder::Build_Models(std::vector<BasisVector> Basis_Vectors)
{
  bool Linear_Independence_Checked = false;
  std::list<std::vector<std::vector<int> > >::const_iterator itk_ij_Extensions = 
    k_ij_Extensions_.begin();
  std::list<std::vector<std::vector<int> > >::const_iterator itk_ij_Extensions_End = 
    k_ij_Extensions_.end();
  for(; itk_ij_Extensions != itk_ij_Extensions_End;
      ++itk_ij_Extensions)
    {
      ModelBuilder The_Builder(Large_ST_Dimensions());
      //Load the initial BVs first.
      for(int b=1; b<static_cast<int>(Initial_Model().BV_Set().size()); ++b)
	  The_Builder.Load_Basis_Vector(Initial_Model().BV_Set().at(b));

      //Load the extension BVs.
      for(int b=0; b<static_cast<int>(Basis_Vectors.size()); ++b)
	The_Builder.Load_Basis_Vector(Basis_Vectors.at(b));

    
      //Check for linear independence.
      if(!Linear_Independence_Checked)
	{
	  if(!The_Builder.Check_Linear_Independence())
	    break;
	  else
	    Linear_Independence_Checked = true;
	}

      //Load the initial k_ijs.
      for(int b=0; b<static_cast<int>(Initial_Model().k_ij().Numerators().size()); ++b)
	The_Builder.Load_k_ij_Row(Initial_Model().k_ij().Numerators().at(b));

      //Load the extension k_ijs.
      for(int b=0; b<static_cast<int>(itk_ij_Extensions->size()); ++b)
	The_Builder.Load_k_ij_Row(itk_ij_Extensions->at(b));

      //Check the consistency of the k_ij matrix.
      if(The_Builder.Check_k_ij_Consistency())
	{
	  The_Builder.Build_Gauge_Group_Model();
	  int Previous_Unique_Model_Count = Unique_LEEFTs().size();
	  Unique_LEEFTs_.insert(LEEFT(The_Builder.FFHS_Model()));
	  Consistent_Models_++;
	  if(static_cast<int>(Unique_LEEFTs().size()) > Previous_Unique_Model_Count)
	    {
	      MD_Out_<<Consistent_Models()<<" "<<Unique_LEEFTs().size()<<std::endl;
	      Data_Writer_.Write_Model_Particle_Content(Model_Out_,
						       The_Builder.FFHS_Model());
	      Model_Out_<<std::endl;
	      Data_Writer_.Write_Model_Extension_BVs(Model_Out_, 
						    The_Builder.FFHS_Model(),
						    Extension_Orders().size());
	      Data_Writer_.Write_Model_Extension_k_ij(Model_Out_, 
						      *itk_ij_Extensions);
	      Model_Out_<<std::endl;
	    }//Close if statement on unique LEEFTs.
	}//Close if statement on k_ij consistency.
    }//Close for loop on k_ij extensions.
}//Close Build_Models.
