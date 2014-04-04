#include "LMGenerator.h"

//PUBLIC.
LMGenerator::LMGenerator(int Large_ST_Dimensions, const std::vector<BasisAlpha>&
			   Common_Basis_Alphas)
{
  Large_ST_Dimensions_ = Large_ST_Dimensions;
  Build_Compact_LM_Triplets();
  Find_Matching_BCs(Common_Basis_Alphas);

}//Close constructor.

LMGenerator::LMGenerator(const LMGenerator& New_LMGenerator)
{
  Large_ST_Dimensions_ = New_LMGenerator.Large_ST_Dimensions();
  Compact_LM_Triplets_ = New_LMGenerator.Compact_LM_Triplets();
  Matching_BCs_ = New_LMGenerator.Matching_BCs();
  SP_LMs_ = New_LMGenerator.SP_LMs();
  NSP_LMs_ = New_LMGenerator.NSP_LMs();
}//Close copy constructor.

//INTERFACE.
void LMGenerator::Build_LM_Chunks()
{
  Build_LMs(0, std::vector<int> ((28 - 2*Large_ST_Dimensions()),0));
}//Close Build_LMs.

//DEBUG.
void LMGenerator::Display_Large_ST_Dimensions() const
{
  std::cout<<"Large ST Dimensions: "<<Large_ST_Dimensions()<<std::endl;
}//Close Display_Large_ST_Dimensions.

void LMGenerator::Display_Compact_LM_Triplets() const 
{
  std::cout<<"Compact LM Triplets: "<<Compact_LM_Triplets().size()<<std::endl;
  for(int a=0; a<static_cast<int>(Compact_LM_Triplets().size()); a++)
    {
      for(int b=0; b<static_cast<int>(Compact_LM_Triplets().at(a).size()); b++)
	std::cout<<Compact_LM_Triplets().at(a).at(b)<<" ";
      std::cout<<std::endl;
    }
  std::cout<<std::endl;
}//Close Display_Compact_LM_Triplets.

void LMGenerator::Display_Matching_BCs() const
{
  std::cout<<"Sets of matching boundary values: "<<std::endl;
  for(int a=0; a<static_cast<int>(Matching_BCs().size()); a++)
    {
      for(int b=0; b<static_cast<int>(Matching_BCs().at(a).size()); b++)
	std::cout<<Matching_BCs().at(a).at(b)<<" ";
      std::cout<<std::endl;
    }//Close for loop on Matching_BCs().
  std::cout<<std::endl;
}//Close Display_Matching_BCs.

void LMGenerator::Display_LMs() const
{
  std::list<std::vector<int> >::const_iterator itSP_LMs = SP_LMs_.begin();
  std::cout<<"Simply paired LMs: "<<SP_LMs().size()<<std::endl;
  for(; itSP_LMs != SP_LMs_.end(); ++itSP_LMs)
    {
      for(int b=0; b<static_cast<int>(itSP_LMs->size()); b++)
	std::cout<<itSP_LMs->at(b)<<" ";
      std::cout<<std::endl;
    }
  std::list<std::vector<int> >::const_iterator itNSP_LMs = NSP_LMs_.begin();
  std::cout<<"Non-simply paired LMs: "<<NSP_LMs().size()<<std::endl;
  for(; itNSP_LMs != NSP_LMs_.end(); ++itNSP_LMs)
    {
      for(int a=0; a<static_cast<int>(itNSP_LMs->size()); a++)
	std::cout<<itNSP_LMs->at(a)<<" ";
      std::cout<<std::endl;
    }//Close for loop on NSP_LMs.
  std::cout<<"Total LMs: "<<SP_LMs().size()+NSP_LMs().size()<<std::endl;
  std::cout<<std::endl;
}//Close Display_LMs.

//PRIVATE.
void LMGenerator::Build_Compact_LM_Triplets()
{
  std::vector<int> x;//(1 0 0)
  std::vector<int> y;//(0 1 0)
  std::vector<int> w;//(0 0 1)
  std::vector<int> xyw;//(1 1 1)

  x.push_back(1);
  x.push_back(0);
  x.push_back(0);

  y.push_back(0);
  y.push_back(1);
  y.push_back(0);

  w.push_back(0);
  w.push_back(0);
  w.push_back(1);

  xyw.push_back(1);
  xyw.push_back(1);
  xyw.push_back(1);

  Compact_LM_Triplets_.push_back(x);
  Compact_LM_Triplets_.push_back(xyw);
  Compact_LM_Triplets_.push_back(y);
  Compact_LM_Triplets_.push_back(w);
}//Close Build_Compact_LM_Triplets.

void LMGenerator::Build_LMs(int Tracker, std::vector<int> LMs_Loader)
{
  if(Tracker<(Large_ST_Dimensions()-2))
    {
      LMs_Loader.at(Tracker) = 1;
      LMs_Loader.at(Tracker+1) = 1;
      Build_LMs(Tracker+2, LMs_Loader);
    }else if(Tracker<(28 - 2*Large_ST_Dimensions() ))
    {
      for(int a=0; a<static_cast<int>(Compact_LM_Triplets().size()); a++)
	{
	  for(int b=0; b<static_cast<int>(Compact_LM_Triplets().at(a).size()); b++)
	    LMs_Loader.at(Tracker+b) = Compact_LM_Triplets().at(a).at(b);
	  if(a<2)
	    {
	      for(int b=0; b<2; b++)
		{
		  for(int c=0; c<static_cast<int>(Compact_LM_Triplets().at(b).size()); c++)
		    LMs_Loader.at(Tracker+3+c) = Compact_LM_Triplets().at(b).at(c);

		  Build_LMs(Tracker+6, LMs_Loader);    
		}//Close for loop on x=1 triplets.
	    }else//a>=2
	    {
	      for(int b=2; b<static_cast<int>(Compact_LM_Triplets().size()); b++)
		{
		  for(int c=0; c<static_cast<int>(Compact_LM_Triplets().at(b).size()); c++)
		    LMs_Loader.at(Tracker+3+c) = Compact_LM_Triplets().at(b).at(c);

		  Build_LMs(Tracker+6, LMs_Loader);
		}//Close for loop on x=0 triplets.
	    }//Close else statement on a.
	}//Close for loop on a.
    }else//Close else statement on Tracker.
    Add_LM_Chunk(LMs_Loader);
}//Close Build_LMs.

void LMGenerator::Find_Matching_BCs(const std::vector<BasisAlpha>& 
				     Common_Basis_Alphas)
{
  //Firstly, add the initial coordinates to the set (y,w)1,...,6.
  std::vector<int> Initial_Matching_BCs_Loader;
  for(int a=Large_ST_Dimensions(); a<Common_Basis_Alphas.at(0).LM_Size(); a+=3)
    {
      //x is a-2.
      Initial_Matching_BCs_Loader.push_back(a-1);//y
      Initial_Matching_BCs_Loader.push_back(a);//w
    }//Close for loop on Common_Basis_Alphas.at(0).LM_Size().
  std::vector<std::vector<int> > Matching_BCs;
  Matching_BCs.push_back(Initial_Matching_BCs_Loader);
  Initial_Matching_BCs_Loader.clear();

  for(int CBA = 0; CBA<static_cast<int>(Common_Basis_Alphas.size()); CBA++)
    {
      std::vector<std::vector<int> > New_Matching_BCs;
      for(int MBC_Row=0; MBC_Row < static_cast<int>(Matching_BCs.size()); MBC_Row++)
	{
	  for(int E_Val=0; E_Val<=Common_Basis_Alphas.at(0).Denominator(); 
	      E_Val+=Common_Basis_Alphas.at(0).Denominator())
	    {
	      std::vector<int> New_Matching_BCs_Loader;
	      for(int MBC_Col=0; MBC_Col < static_cast<int>(Matching_BCs.at(MBC_Row).size());
		  MBC_Col++)
		{
		  if(Common_Basis_Alphas.at(CBA).Numerator().at
		     (Matching_BCs.at(MBC_Row).at(MBC_Col))==E_Val)
		    New_Matching_BCs_Loader.push_back
		      (Matching_BCs.at(MBC_Row).at(MBC_Col));
		}//Close for loop on Matching_BCs.at(MBC_Row).size()
	      if(New_Matching_BCs_Loader.size() != 0)
		New_Matching_BCs.push_back(New_Matching_BCs_Loader);
	      New_Matching_BCs_Loader.clear();
	    }//Close for loop on E_Val.
	}//Close for loop on Matching_BCs.size().
      Matching_BCs.swap(New_Matching_BCs);
      New_Matching_BCs.clear();
    }//Close for loop on Common_Basis_Alphas.size().
  Matching_BCs_ = Matching_BCs;
}//Close Find_Matching_BCs.

void LMGenerator::Add_LM_Chunk(const std::vector<int>& LMs_Loader)
{
  if(Is_Simply_Paired(LMs_Loader))
    SP_LMs_.push_back(LMs_Loader);
  else
    NSP_LMs_.push_back(LMs_Loader);
}//Close Add_LM_Chunk.

bool LMGenerator::Is_Simply_Paired(const std::vector<int>& LMs_Loader)
{
  for(int MBC_Row=0; MBC_Row<static_cast<int>(Matching_BCs().size()); MBC_Row++)
    {
      for(int E_Val = 0; E_Val <= 1; E_Val++)
	{
	  int E_Counter = 0;
	  for(int MBC_Col=0; MBC_Col < static_cast<int>(Matching_BCs().at(MBC_Row).size()); MBC_Col++)
	    {
	      if(LMs_Loader.at(Matching_BCs().at(MBC_Row).at(MBC_Col))==E_Val)
		E_Counter++;
	    }//Close for loop on Matching_BCs().at(MBC_Row).size().
	  if(E_Counter%2 != 0)
	    return false;
	}//Close for loop on E_Val.
    }//Close for loop on Matching_BCs().size().
  return true;
}//Close Is_Simply_Paired.

