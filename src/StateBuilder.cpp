#include "StateBuilder.h"

//PUBLIC.

StateBuilder::StateBuilder(const Alpha& The_Alpha, const std::map<int, int>&
			     Fermion_Mode_Map, 
			     const std::vector<BasisAlpha>& Common_Basis_Alphas,
			     const GSOCoefficientMatrix& k_ij)
{
  The_Alpha_ = The_Alpha;
  Alpha_Type_ = The_Alpha.Type();
  Fermion_Mode_Map_ = Fermion_Mode_Map;

  GSOProjector GSO(Common_Basis_Alphas, The_Alpha, k_ij, Fermion_Mode_Map);
  GSO_ = GSO;
}//Close constructor.

//INTERFACE.
void StateBuilder::Build_States()
{
  switch(Alpha_Type())
    {
    case 'b':
      Build_Boson_States();
      break;

    case 'f':
      Build_Fermion_States();
      break;

    case 's':
      Build_SUSY_States();
      break;
    }//Close switch.
}//Close Build_States.

//DEBUG.
void StateBuilder::Display_The_Alpha() const
{
  std::cout<<"Alpha: "<<std::endl;
  The_Alpha().Display();
  std::cout<<std::endl;
}//Close Display_The_Alpha.

void StateBuilder::Display_States() const
{
  std::cout<<"States: "<<States().size()<<std::endl;
  std::list<State>::const_iterator itStates = States_.begin();
  for(; itStates != States_.end(); itStates++)
    itStates->Display();
  std::cout<<std::endl;
}//Close Display_States.

//PRIVATE.
//HELPERS.
void StateBuilder::Build_Fermion_States()
{
  int Large_ST_Dimensions = Calculate_Large_ST_Dimensions(The_Alpha().LM_Size());

  StateLMBuilder Fermion_LM_Builder(The_Alpha(), Large_ST_Dimensions, 
				      Fermion_Mode_Map());
  StateRMBuilder Fermion_RM_Builder(The_Alpha(), Large_ST_Dimensions,
				      Fermion_Mode_Map());

  Fermion_LM_Builder.Build_Massless_State_LMs();
  Fermion_RM_Builder.Build_Massless_State_RMs();

  std::list<std::vector<int> > Massless_State_LMs = 
    Fermion_LM_Builder.Massless_State_LMs();
  std::list<std::vector<int> > Massless_State_RMs = 
    Fermion_RM_Builder.Massless_State_RMs();

  std::list<std::vector<int> >::iterator itLMs = Massless_State_LMs.begin();
  for(; itLMs != Massless_State_LMs.end(); ++itLMs)
    {
      std::list<std::vector<int> >::iterator itRMs = 
	Massless_State_RMs.begin();
      for(; itRMs != Massless_State_RMs.end(); ++itRMs)
	{
	  std::vector<int> New_State_Numerator = *itLMs;
	 
	  for(int c=0; c<static_cast<int>(itRMs->size()); c++)
	    New_State_Numerator.push_back(itRMs->at(c));

	  int New_State_Denominator = 2*The_Alpha().Denominator();
	  int New_State_LM_Size = The_Alpha().LM_Size();
			       

	  State New_State(New_State_Numerator, New_State_Denominator, 
			  New_State_LM_Size);

	  if(GSO().GSOP(New_State))
	    States_.push_back(New_State);

	}//Close for loop on Massless_State_RMs.
    }//Close for loop on Massless_State_LMs.
}//Close Build_Fermion_States.

void StateBuilder::Build_Boson_States()
{
  int Large_ST_Dimensions = Calculate_Large_ST_Dimensions(The_Alpha().LM_Size());

  StateRMBuilder Boson_RM_Builder(The_Alpha(), Large_ST_Dimensions, 
				    Fermion_Mode_Map());

 
  Boson_RM_Builder.Build_Massless_State_RMs();

  std::list<std::vector<int> > Massless_State_RMs = 
    Boson_RM_Builder.Massless_State_RMs();

  std::list<std::vector<int> >::iterator itRMs = Massless_State_RMs.begin();
  for(; itRMs != Massless_State_RMs.end(); itRMs++)
    {
      int New_State_Denominator = 2*The_Alpha().Denominator();
      std::vector<int> New_State_Numerator;
      for(int b=0; b<2; b++)
	New_State_Numerator.push_back(New_State_Denominator);
      for(int b=2; b<Large_ST_Dimensions-2; b++)
	New_State_Numerator.push_back(0);
      for(int b=Large_ST_Dimensions-2; b<The_Alpha().LM_Size(); b++)
	New_State_Numerator.push_back(0);
      for(int b=0; b<static_cast<int>(itRMs->size()); b++)
	New_State_Numerator.push_back(itRMs->at(b));

      State New_State(New_State_Numerator, New_State_Denominator,
		      The_Alpha().LM_Size());

      if(GSO().GSOP(New_State))
	States_.push_back(New_State);
    }//Close for loop on Massless_State_RMs.
}//Close Build_Boson_States.

void StateBuilder::Build_SUSY_States()
{
  int Large_ST_Dimensions = Calculate_Large_ST_Dimensions(The_Alpha().LM_Size());

  StateLMBuilder SUSY_LM_Builder(The_Alpha(), Large_ST_Dimensions, 
				   Fermion_Mode_Map());

  SUSY_LM_Builder.Build_Massless_State_LMs();

  std::list<std::vector<int> > Massless_State_LMs = 
    SUSY_LM_Builder.Massless_State_LMs();

  std::list<std::vector<int> >::iterator itLMs = Massless_State_LMs.begin();
  for(; itLMs != Massless_State_LMs.end(); itLMs++)
    {
      std::vector<int> New_State_Numerator = *itLMs;

      for(int b=itLMs->size(); b<static_cast<int>(The_Alpha().Numerator().size());
	  b++)
	New_State_Numerator.push_back(0);

      int New_State_Denominator = 2*The_Alpha().Denominator();

      State New_State(New_State_Numerator, New_State_Denominator,
		      The_Alpha().LM_Size());

      if(GSO().GSOP(New_State))
	States_.push_back(New_State);
    }//Close for loop on Massless_State_LMs.
}//Close Build_SUSY_States.

int StateBuilder::Calculate_Large_ST_Dimensions(int LM_Size)
{
  return (28 - LM_Size)/2;
}//Close Calculate_Large_ST_Dimensions.
