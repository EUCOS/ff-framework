#include "StateLMBuilder.h"
//PUBLIC.
StateLMBuilder::StateLMBuilder(const std::vector<int>& Alpha_LM_Numerator, 
				   int Alpha_LM_Denominator, 
				   int Large_ST_Dimensions,
				   const std::map<int, int>& Fermion_Mode_Map)
{
  Alpha_LM_Numerator_ = Alpha_LM_Numerator;
  Alpha_LM_Denominator_ = Alpha_LM_Denominator;
  Large_ST_Dimensions_ = Large_ST_Dimensions;
  Fermion_Mode_Map_ = Fermion_Mode_Map;

}//Close constructor.

StateLMBuilder::StateLMBuilder(const Alpha& The_Alpha, int Large_ST_Dimensions,
				   const std::map<int, int>& Fermion_Mode_Map)
{
  for(int a=0; a<The_Alpha.LM_Size(); a++)
    Alpha_LM_Numerator_.push_back(The_Alpha.Numerator().at(a));

  Alpha_LM_Denominator_ = The_Alpha.Denominator();
  Large_ST_Dimensions_ = Large_ST_Dimensions;
  Fermion_Mode_Map_ = Fermion_Mode_Map;

}//Close second constructor.

StateLMBuilder::StateLMBuilder(const StateLMBuilder& New_StateLMBuilder)
{
  Large_ST_Dimensions_ = New_StateLMBuilder.Large_ST_Dimensions();
  Alpha_LM_Numerator_ = New_StateLMBuilder.Alpha_LM_Numerator();
  Alpha_LM_Denominator_ = New_StateLMBuilder.Alpha_LM_Denominator();
  Fermion_Mode_Map_ = New_StateLMBuilder.Fermion_Mode_Map();
  ST_LM_Modes_ = New_StateLMBuilder.ST_LM_Modes();
  Compact_LM_Modes_ = New_StateLMBuilder.Compact_LM_Modes();
  Massless_State_LMs_ = New_StateLMBuilder.Massless_State_LMs();
}//Close copy constructor.

//INTERFACE.
void StateLMBuilder::Build_Massless_State_LMs()
{
  std::vector<int> LM_ST;
  std::vector<int> LM_Compact;

  for(int a=0; a<static_cast<int>(Alpha_LM_Numerator().size()); a++)
    {
      if(a<(Large_ST_Dimensions()-2))
	LM_ST.push_back(Alpha_LM_Numerator().at(a));
      else
	LM_Compact.push_back(Alpha_LM_Numerator().at(a));
    }
  Build_ST_LM_Modes();
  Build_Compact_LM_Modes(0, LM_Compact);

  std::list<std::vector<int> >::iterator itST = ST_LM_Modes_.begin();
  std::list<std::vector<int> >::iterator itCompact = Compact_LM_Modes_.begin();
  for(; itST != ST_LM_Modes_.end(); itST++)
    {
      if(Large_ST_Dimensions()<10)
	{
	  for(; itCompact != Compact_LM_Modes_.end(); itCompact++)
	    {
	      std::vector<int> Massless_State_LMs_Loader;
	      for(int c=0; c<static_cast<int>(itST->size()); c++)
		Massless_State_LMs_Loader.push_back(itST->at(c));
	      for(int c=0; c<static_cast<int>(itCompact->size()); c++)
		Massless_State_LMs_Loader.push_back(itCompact->at(c));

	      Massless_State_LMs_.push_back(Massless_State_LMs_Loader);
	    }//Close inner for loop on compact LM modes.
	}else //Close if statement for large st dimensions.
	Massless_State_LMs_.push_back(*itST);
    }//Close outer for loop on ST LM modes.

}//Close Build_Massless_State_LMs.

//DEBUG.
void StateLMBuilder::Display_Massless_State_LMs() const
{
  std::list<std::vector<int> >::const_iterator itLMs = 
    Massless_State_LMs_.begin();
  for(; itLMs != Massless_State_LMs_.end(); itLMs++)
    {
      for(int b=0; b<static_cast<int>(itLMs->size()); b++)
	std::cout<<itLMs->at(b)<<" ";
      std::cout<<std::endl;
    }//Close for loop on Massless_State_LMs.
}//Close Display_Massless_State_LMs.

void StateLMBuilder::Display_Fermion_Mode_Map() const
{
  std::map<int, int>::const_iterator itMap = Fermion_Mode_Map_.begin();
  std::cout<<"Fermion Mode Map:"<<std::endl;
  for(; itMap != Fermion_Mode_Map_.end(); itMap++)
    std::cout<<itMap->first<<" "<<itMap->second<<std::endl;
  std::cout<<std::endl;
}//Close Display_Fermion_Mode_Map.

//PRIVATE.

void StateLMBuilder::Build_ST_LM_Modes()
{
  std::vector<int> ST_LM_Loader;
  //First, build the even parity ST modes.
  int Loop_End = Large_ST_Dimensions()-2;
  for(int a=0; a<Loop_End; a++)
    ST_LM_Loader.push_back(Alpha_LM_Numerator().at(a));
  ST_LM_Modes_.push_back(ST_LM_Loader);
  ST_LM_Loader.clear();

  //Now build the odd parity ST modes (D=10 only).
  if(Large_ST_Dimensions() == 10)
    {
      for(int a=0; a<6; a++)
	ST_LM_Loader.push_back(Alpha_LM_Numerator().at(a));
      for(int a=6; a<8; a++)
	ST_LM_Loader.push_back(-Alpha_LM_Numerator().at(a));
      ST_LM_Modes_.push_back(ST_LM_Loader);
    }
}//Close Build_ST_LM_Modes.

void StateLMBuilder::Build_Compact_LM_Modes(int element, 
					      std::vector<int>& LM_Compact)
{
  if(element < static_cast<int>(LM_Compact.size()))
    {
      if((Alpha_LM_Numerator().at(element+Large_ST_Dimensions()-2) != 0) && 
	 In_Map(element) && !Real_LM_Mode(element))
	{
	  for(int a=-1; a<=0; a++)
	    {
	      int Complex_Fermion_Partner = 
		(Fermion_Mode_Map_[element + Large_ST_Dimensions()-2]
		 - (Large_ST_Dimensions()-2));//Index for the compact vector.

	      LM_Compact.at(element) = 
		Alpha_LM_Numerator().at(element + Large_ST_Dimensions()-2) +
		a*2*Alpha_LM_Denominator();
	      LM_Compact.at(Complex_Fermion_Partner) = 
		Alpha_LM_Numerator().at(Complex_Fermion_Partner + 
					Large_ST_Dimensions()-2)
		+ a*2*Alpha_LM_Denominator();

	      Build_Compact_LM_Modes(element+1, LM_Compact);
	    }//Close for loop for applying the lowering operator.
	}else
	Build_Compact_LM_Modes(element+1, LM_Compact);
    }else
    Compact_LM_Modes_.push_back(LM_Compact);
}//Close Build_Compact_LM_Modes.

bool StateLMBuilder::Real_LM_Mode(int element)//Works only for compact modes.
//ST fermions are always complex.
{
  element += (Large_ST_Dimensions() - 2);
  int Complex_Fermion_Partner = (Fermion_Mode_Map().find(element))->second;
  if(Complex_Fermion_Partner > static_cast<int>(Alpha_LM_Numerator().size()))
    return true;
  else
    return false;
}//Close Real_LM_Mode.

bool StateLMBuilder::In_Map(int element)
{
  element += Large_ST_Dimensions() - 2;
  return Fermion_Mode_Map_.find(element) != Fermion_Mode_Map_.end();
}//Close In_Map.
