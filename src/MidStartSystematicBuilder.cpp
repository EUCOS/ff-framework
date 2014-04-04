#include "MidStartSystematicBuilder.h"

MidStartSystematicBuilder::MidStartSystematicBuilder
(const ModelBuilder& Initial_Build, int Large_ST_Dimensions, 
 const std::vector<int>& Extension_Orders,
 const std::string&  Output_File_Name, 
 const std::vector<BasisVector>& Starting_Basis_Vectors):
  SystematicBuilder(Initial_Build, Large_ST_Dimensions, Extension_Orders, 
		     Output_File_Name)
{
  Starting_Basis_Vectors_ = Starting_Basis_Vectors;
  Layer_Built_ = std::vector<bool>(Starting_Basis_Vectors.size(), false);
}//Close contructor.

MidStartSystematicBuilder::MidStartSystematicBuilder
(const ModelBuilder& Initial_Build, int Large_ST_Dimensions,
 const std::vector<int>& Extension_Orders, const std::string& Output_File_Name,
 const std::vector<bool>& Simply_Paired_Layers, const std::vector<BasisVector>&
 Starting_Basis_Vectors):SystematicBuilder(Initial_Build, Large_ST_Dimensions,
				   Extension_Orders, Output_File_Name,
				   Simply_Paired_Layers)
{
  Starting_Basis_Vectors_ = Starting_Basis_Vectors;
  Layer_Built_ = std::vector<bool>(Starting_Basis_Vectors.size(), false);
}//Close second constructor.

MidStartSystematicBuilder::MidStartSystematicBuilder
(const MidStartSystematicBuilder& New_MidStartSystematicBuilder):
  SystematicBuilder(New_MidStartSystematicBuilder)
{
  Starting_Basis_Vectors_ = 
    New_MidStartSystematicBuilder.Starting_Basis_Vectors();
  Layer_Built_ = New_MidStartSystematicBuilder.Layer_Built();
}//Close copy constructor.

//INTERFACE.
void MidStartSystematicBuilder::Perform_Search()
{
  Data_Writer_.Write_BV_Search_Front_Info(Model_Out_, Initial_Model(),
					  Large_ST_Dimensions(), 
					  Extension_Orders());

  if(k_ij_Extensions().size()==0)//Check to see if k_ij is fixed for this search.
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
void MidStartSystematicBuilder::Build_Extensions(int Layer, 
						    std::vector<BasisVector> 
						    Basis_Vectors)
{
  if(Layer<static_cast<int>(Extension_Orders().size()))
    Build_Basis_Vectors(Layer, Basis_Vectors);
  else
    Build_Models(Basis_Vectors);
}//Close Build_Extensions.

void MidStartSystematicBuilder::Build_Basis_Vectors(int Layer, 
						       std::vector<BasisVector>
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
  //DEBUG.
  std::cout<<"Layer: "<<Layer<<std::endl;
  //END DEBUG.
  std::list<Chunk> LMs;
  if(Simply_Paired_Layers().at(Layer))
    LMs = BV_Generator.SP_LMs();
  else
    LMs = BV_Generator.NSP_LMs();

  std::list<Chunk> Comp;
  if(Simply_Paired_Layers().at(Layer))
    Comp = BV_Generator.SP_RMs_Compact();
  else
    Comp = BV_Generator.NSP_RMs_Compact();

  std::list<Chunk> Obs = BV_Generator.RMs_Observable();
  std::list<Chunk> Hid = BV_Generator.RMs_Hidden();

  std::list<Chunk>::const_iterator itLMs = LMs.begin();
  std::list<Chunk>::const_iterator itObs = Obs.begin();
  std::list<Chunk>::const_iterator itComp = Comp.begin();
  std::list<Chunk>::const_iterator itHid = Hid.begin();

  if(!Layer_Built().at(Layer)) //If we haven't already built this layer before.
    {
      itLMs = Set_LM_Start_Point(LMs, Starting_Basis_Vectors().at(Layer));
      itComp = Set_Comp_Start_Point(Comp, 
				    Starting_Basis_Vectors().at(Layer));
      itObs = Set_Obs_Start_Point(Obs, Starting_Basis_Vectors().at(Layer));
      itHid = Set_Hid_Start_Point(Hid, Starting_Basis_Vectors().at(Layer));
      Layer_Built_.at(Layer) = true;
      }//Close if statement.

  std::list<Chunk>::const_iterator itLMs_End = LMs.end();
  std::list<Chunk>::const_iterator itObs_End = Obs.end();
  std::list<Chunk>::const_iterator itComp_End = Comp.end();
  std::list<Chunk>::const_iterator itHid_End = Hid.end();

  Total_BVs_ += LMs.size()*Obs.size()*Comp.size()*Hid.size();

  ChunkConsistencyChecker Checker(FF::LCM(2, Extension_Orders().at(Layer)),
				    Initial_Common_Basis_Alphas().at(0).
				    Denominator());

  for(; itLMs != itLMs_End; ++itLMs)
    {
      for(; itComp != itComp_End; ++itComp)
	{
	  for(; itObs != itObs_End; ++itObs)
	    {
	      //Make sure the pairings are good before proceeding.
	      if(!Simply_Paired_Layers().at(Layer) && 
		 !Checker.Check_Simultaneous_Periodic_Modes(*itLMs, *itComp))
		break;

	      for(; itHid != itHid_End; ++itHid)
		{
		  if(Checker.Check_Modular_Invariance(*itLMs, *itObs, 
						      *itComp, *itHid))
		    {
		      std::vector<int> New_BV;
		      for(int a=0; a<static_cast<int>(itLMs->BV_Chunk().size()); ++a)
			New_BV.push_back(itLMs->BV_Chunk().at(a));
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
	      itHid = Hid.begin();
	    }//Close for loop on itObs.
	  itObs = Obs.begin();
	}//Close for loop on itComp.
      itComp = Comp.begin();
    }//Close for loop on itLMs.
}//Close Build_Basis_Vectors.

std::list<Chunk>::const_iterator MidStartSystematicBuilder::Set_LM_Start_Point
(const std::list<Chunk>& LMs, const BasisVector& Basis_Vector)
{
  std::vector<int> BV_Chunk;
  for(int a=0; a<Basis_Vector.LM_Size(); ++a)
    BV_Chunk.push_back(Basis_Vector.BV().at(a));
  std::list<Chunk>::const_iterator itLMs = LMs.begin();
  std::list<Chunk>::const_iterator itLMs_End = LMs.end();
  for(; itLMs != itLMs_End; ++itLMs)
    {
      if(itLMs->BV_Chunk() == BV_Chunk)
	return itLMs;
    }//Close for loop on the list LMs.
  //DEBUG.
  std::cout<<"Did not find LM Chunk."<<std::endl;
  //END DEBUG.
  return itLMs;
}//Close Set_LM_Start_Point.

std::list<Chunk>::const_iterator MidStartSystematicBuilder::
Set_Comp_Start_Point(const std::list<Chunk>& Comp, const BasisVector& 
		     Basis_Vector)
{
  std::vector<int> BV_Chunk;
  for(int a=(Basis_Vector.LM_Size()+16); a<(Basis_Vector.LM_Size()+16+
					    Basis_Vector.RM_Compact_Size());++a)
    BV_Chunk.push_back(Basis_Vector.BV().at(a));

  std::list<Chunk>::const_iterator itComp = Comp.begin();
  std::list<Chunk>::const_iterator itComp_End = Comp.end();
  for(; itComp != itComp_End; ++itComp)
    {
      if(itComp->BV_Chunk() == BV_Chunk)
	return itComp;
    }//Close for loop on the list Comp.
  //DEBUG.
  std::cout<<"Did not find Comp chunk."<<std::endl;
  //END DEBUG.
  return itComp;
}//Close Set_Comp_Start_Point.

std::list<Chunk>::const_iterator MidStartSystematicBuilder::
Set_Obs_Start_Point(const std::list<Chunk>& Obs, const BasisVector&
		    Basis_Vector)
{
  std::vector<int> BV_Chunk;
  for(int a=(Basis_Vector.LM_Size()); a<(Basis_Vector.LM_Size()+16); ++a)
    BV_Chunk.push_back(Basis_Vector.BV().at(a));

  std::list<Chunk>::const_iterator itObs = Obs.begin();
  std::list<Chunk>::const_iterator itObs_End = Obs.end();
  for(; itObs != itObs_End; ++itObs)
    {
      if(itObs->BV_Chunk() == BV_Chunk)
	return itObs;
    }//Close for loop on the list Obs.

  //DEBUG.
  std::cout<<"Did not find Obs chunk."<<std::endl;
  //END DEBUG.
  return itObs;
}//Close Set_Obs_Start_Point.

std::list<Chunk>::const_iterator MidStartSystematicBuilder::
Set_Hid_Start_Point(const std::list<Chunk>& Hid, const BasisVector&
		    Basis_Vector)
{
  std::vector<int> BV_Chunk;
  for(int a=(Basis_Vector.LM_Size() + 16 + Basis_Vector.RM_Compact_Size());
      a<(Basis_Vector.LM_Size() + 16 + Basis_Vector.RM_Compact_Size() + 16);
      ++a)
    BV_Chunk.push_back(Basis_Vector.BV().at(a));

  std::list<Chunk>::const_iterator itHid = Hid.begin();
  std::list<Chunk>::const_iterator itHid_End = Hid.end();
  for(; itHid != itHid_End; ++itHid)
    {
      if(itHid->BV_Chunk() == BV_Chunk)
	return itHid;
    }//Close for loop on the list Obs.
  //DEBUG.
  std::cout<<"Did not find Hid chunk."<<std::endl;
  //END DEBUG.
  return itHid;
}//Close Set_Hid_Start_Point.
