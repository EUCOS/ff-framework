#include "RealChunkGenerator.h"

//PUBLIC.
RealChunkGenerator::RealChunkGenerator
(std::vector<BasisAlpha> Common_Basis_Alphas, int Order)
/*
The constructor takes as an argument the common basis alphas of the model which is
being extended.In addition, it takes the order of the basis vector being generated.
 It initializes the privates member of the same name. The reason
the common basis alphas are used is that a common denominator must be found
to pick out the matching boundary conditions.
*/
{
  Common_Basis_Alphas_ = Common_Basis_Alphas;
  Order_ = Order;
}//Close constructor.

RealChunkGenerator::RealChunkGenerator(const std::vector<BasisAlpha>& 
					   Common_Basis_Alphas, int Order, 
					   bool First_Complex_Extension)
{
  Common_Basis_Alphas_ = Common_Basis_Alphas;
  Order_ = Order;
  First_Complex_Extension_ = First_Complex_Extension;
}//Close second constructor.

RealChunkGenerator::RealChunkGenerator(const RealChunkGenerator& 
					   New_RealChunkGenerator)
/*
The copy constructor copies the private members of New_Real_Chunk_Builder
into *this.
*/
{
  Common_Basis_Alphas_ = New_RealChunkGenerator.Common_Basis_Alphas();
  Order_ = New_RealChunkGenerator.Order();
  First_Complex_Extension_ = New_RealChunkGenerator.First_Complex_Extension();
  SP_Real_Chunks_ = New_RealChunkGenerator.SP_Real_Chunks();
  NSP_Real_Chunks_ = New_RealChunkGenerator.NSP_Real_Chunks();
}//Close copy constructor.

//INTERFACE.
void RealChunkGenerator::Build_Real_Chunks()
{
  std::vector<std::vector<int> > Matching_BCs = Find_Matching_BCs();

  std::vector<int> Real_Chunks_Loader
    (Common_Basis_Alphas().at(0).RM_Compact_Size(),0);
  Generate_Real_Chunks(Matching_BCs, 0, 0, Real_Chunks_Loader);
}//Close Build_Real_Chunks.

void RealChunkGenerator::Display_Real_Chunks() const
/*
This function serves a debugging purpose. It displays the SP_Real_Chunks_ 
and NSP_Real_Chunks_ members on the screen.
*/
{
  std::cout<<"Simply paired real chunks: "<<SP_Real_Chunks().size()<<std::endl;
  std::set<std::vector<int> >::const_iterator itSP_Real_Chunks = 
    SP_Real_Chunks_.begin();
  for(; itSP_Real_Chunks != SP_Real_Chunks_.end(); 
      ++itSP_Real_Chunks)
    {
      for(int a=0; a<static_cast<int>(itSP_Real_Chunks->size()); a++)
	std::cout<<itSP_Real_Chunks->at(a)<<" ";
      std::cout<<std::endl;
    }//Close for loop on itSP_Real_Chunks.

  std::cout<<"Non-simply paired real chunks: "<<NSP_Real_Chunks().size()<<std::endl;
  std::set<std::vector<int> >::const_iterator itNSP_Real_Chunks = 
    NSP_Real_Chunks_.begin();
  for(; itNSP_Real_Chunks != NSP_Real_Chunks_.end(); 
      ++itNSP_Real_Chunks)
    {
      for(int a=0; a<static_cast<int>(itNSP_Real_Chunks->size()); a++)
	std::cout<<itNSP_Real_Chunks->at(a)<<" ";
      std::cout<<std::endl;
    }//Close for loop on NSP_Real_Chunks.
}//Close Display_Real_Chunks.

//PRIVATE.
std::vector<std::vector<int> > RealChunkGenerator::Find_Matching_BCs()
{
  int Maximum_Element_Value = Common_Basis_Alphas().at(0).Denominator();
  int Minimum_Element_Value = -1*Maximum_Element_Value;

  std::vector<int> Matching_BCs_Loader;
  for(int a=0; a<Common_Basis_Alphas().at(0).RM_Compact_Size(); a++)
    Matching_BCs_Loader.push_back(a + Common_Basis_Alphas().at(0).LM_Size() + 16);

  std::vector<std::vector<int> > Matching_BCs;
  Matching_BCs.push_back(Matching_BCs_Loader);

  Matching_BCs_Loader.clear();
  for(int CBA=0; CBA<static_cast<int>(Common_Basis_Alphas().size()); CBA++)
    {
      std::vector<std::vector<int> > New_Matching_BCs;
      for(int MBC_Row=0; MBC_Row<static_cast<int>(Matching_BCs.size()); MBC_Row++)
	{
	  for(int E_Val=Minimum_Element_Value; E_Val<=Maximum_Element_Value; E_Val++)
	    {
	      for(int MBC_Col=0; MBC_Col<static_cast<int>(Matching_BCs.at(MBC_Row).size()); MBC_Col++)
		{
		  if(Common_Basis_Alphas().at(CBA).Numerator().at
		     (Matching_BCs.at(MBC_Row).at(MBC_Col)) == E_Val)
		    Matching_BCs_Loader.push_back
		      (Matching_BCs.at(MBC_Row).at(MBC_Col));
		}//Close for loop on Matching_BCs.at(MBC_Row).size().
	      if(Matching_BCs_Loader.size() != 0)
		New_Matching_BCs.push_back(Matching_BCs_Loader);
	      Matching_BCs_Loader.clear();
	    }//Close for loop on Maximum_Element_Value.
	}//Close for loop on Matching_BCs.size().
      Matching_BCs.swap(New_Matching_BCs);
      New_Matching_BCs.clear();
    }//Close for loop on Common_Basis_Alphas.
  return Matching_BCs;
}//Close Find_Matching_BCs.

void RealChunkGenerator::Generate_Real_Chunks
(const std::vector<std::vector<int> >& Matching_BCs, int MBC_Row_Tracker, 
 int MBC_Col_Tracker, std::vector<int> Real_Chunks_Loader)
{
  if(MBC_Row_Tracker<static_cast<int>(Matching_BCs.size()))
    {
      if(MBC_Col_Tracker<static_cast<int>(Matching_BCs.at(MBC_Row_Tracker).size()))
	{
	  int Full_BV_Adjust = Common_Basis_Alphas().at(0).LM_Size()+16;
	  //Compensates for the size difference between the Real_Chunks_Loader vector
	  //and the Common_Basis_Alphas elements indexed by Matching_BCs.
	  int E_Val_Start = 0;
	  int E_Val_End = Order()-1;
	  if(First_Complex_Extension())
	    E_Val_End = Order()/2;//Rounds down to an integer by truncation.
	  if(MBC_Col_Tracker != 0)
	    E_Val_Start = Real_Chunks_Loader.at(Matching_BCs.at(MBC_Row_Tracker)
						.at(MBC_Col_Tracker-1) - 
						Full_BV_Adjust);
	  if((Matching_BCs.at(MBC_Row_Tracker).size()%2==1) && 
	     (MBC_Col_Tracker == static_cast<int>(Matching_BCs.at(MBC_Row_Tracker).size())-1))
	    {
	      for(int E_Val = E_Val_Start; E_Val<=E_Val_End; 
		  E_Val+=Common_Basis_Alphas().at(0).Denominator())
		{
		  Real_Chunks_Loader.at(Matching_BCs.at(MBC_Row_Tracker).at
					(MBC_Col_Tracker)-Full_BV_Adjust) = E_Val;
		  Generate_Real_Chunks(Matching_BCs, MBC_Row_Tracker, 
				       MBC_Col_Tracker+1,Real_Chunks_Loader);
		}//Close for loop on E_Val.
	    }else//if/else on real modes.
	    {
	      for(int E_Val = E_Val_Start; E_Val<=E_Val_End; E_Val++)
		{
		  Real_Chunks_Loader.at(Matching_BCs.at(MBC_Row_Tracker).at
					(MBC_Col_Tracker)-Full_BV_Adjust) = E_Val;
		  if(Order()%2 != 0)
		    Real_Chunks_Loader.at(Matching_BCs.at(MBC_Row_Tracker).at
					  (MBC_Col_Tracker+1)-Full_BV_Adjust) = 
		      E_Val;
		  Generate_Real_Chunks(Matching_BCs, MBC_Row_Tracker,
				       MBC_Col_Tracker+1+(Order()%2), 
				       Real_Chunks_Loader);
		}//Close for loop on E_Val.
	    }//Close if/else on real fermion modes.
	}else
	Generate_Real_Chunks(Matching_BCs, MBC_Row_Tracker+1, 0, Real_Chunks_Loader);
    }else//Close if/else on MBC_Row_Tracker.
    Add_Real_Chunk(Matching_BCs, Real_Chunks_Loader);
}//Close Generate_Real_Chunks.

void RealChunkGenerator::Add_Real_Chunk(const std::vector<std::vector<int> >&
					  Matching_BCs, const std::vector<int>&
					  Real_Chunks_Loader)
{
  if(Is_Simply_Paired(Matching_BCs, Real_Chunks_Loader))
    SP_Real_Chunks_.insert(Real_Chunks_Loader);
  else
    NSP_Real_Chunks_.insert(Real_Chunks_Loader);
}//Close Add_Real_Chunk.

bool RealChunkGenerator::Is_Simply_Paired(const std::vector<std::vector<int> >& 
					    Matching_BCs, 
					    const std::vector<int>& 
					    Real_Chunks_Loader)
{
  for(int MBC_Row=0; MBC_Row<static_cast<int>(Matching_BCs.size()); MBC_Row++)
    {
      for(int E_Val = 0; E_Val < Order(); E_Val++)
	{
	  int E_Counter = 0;
	  for(int MBC_Col=0; MBC_Col < static_cast<int>(Matching_BCs.at(MBC_Row).size()); MBC_Col++)
	    {
	      if(Real_Chunks_Loader.at(Matching_BCs.at(MBC_Row).at(MBC_Col)-
				       (Common_Basis_Alphas().at(0).LM_Size()+16)) 
		 == E_Val)
		E_Counter++;
	    }//Close for loop on Matching_BCs.at(MBC_Row).size().
	  if(E_Counter%2 != 0)
	    return false;
	}//Close for loop on Order().
    }//Close for loop on Matching_BCs.size().
  return true;
}//Close Is_Simply_Paired.
