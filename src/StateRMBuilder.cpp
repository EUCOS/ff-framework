#include "StateRMBuilder.h"

//PUBLIC.
StateRMBuilder::StateRMBuilder(const std::vector<int>& Alpha_RM_Numerator,
				   int Alpha_RM_Denominator, 
				   int Large_ST_Dimensions,
				   const std::map<int, int>& Fermion_Mode_Map)
{
  Alpha_RM_Numerator_ = Alpha_RM_Numerator;
  Alpha_RM_Denominator_ = Alpha_RM_Denominator;
  Large_ST_Dimensions_ = Large_ST_Dimensions;
  Fermion_Mode_Map_ = Fermion_Mode_Map;

  Build_Reverse_Fermion_Mode_Map();

  LM_Size_ = 28 - 2*Large_ST_Dimensions;
  Mass_Limit_ = 2*4*Alpha_RM_Denominator*Alpha_RM_Denominator;
}//Close constructor.

StateRMBuilder::StateRMBuilder(const Alpha& The_Alpha, int Large_ST_Dimensions,
				   const std::map<int, int>& Fermion_Mode_Map)
{
  for(int a=The_Alpha.LM_Size(); a<static_cast<int>(The_Alpha.Numerator().size()); a++)
    Alpha_RM_Numerator_.push_back(The_Alpha.Numerator().at(a));

  Alpha_RM_Denominator_ = The_Alpha.Denominator();
  Large_ST_Dimensions_ = Large_ST_Dimensions;
  Fermion_Mode_Map_ = Fermion_Mode_Map;
  LM_Size_ = The_Alpha.LM_Size();

  Build_Reverse_Fermion_Mode_Map();

  Mass_Limit_ = 2*4*Alpha_RM_Denominator()*Alpha_RM_Denominator();

}//Close second constructor.

StateRMBuilder::StateRMBuilder(const StateRMBuilder& New_StateRMBuilder)
{
  Large_ST_Dimensions_ = New_StateRMBuilder.Large_ST_Dimensions();
  Alpha_RM_Numerator_ = New_StateRMBuilder.Alpha_RM_Numerator();
  Alpha_RM_Denominator_ = New_StateRMBuilder.Alpha_RM_Denominator();
  Fermion_Mode_Map_ = New_StateRMBuilder.Fermion_Mode_Map();
  Mass_Limit_ = New_StateRMBuilder.Mass_Limit();
  Reverse_Fermion_Mode_Map_ = New_StateRMBuilder.Reverse_Fermion_Mode_Map();
  Massless_State_RMs_ = New_StateRMBuilder.Massless_State_RMs();
}//Close copy constructor.

//INTERFACE.
void StateRMBuilder::Build_Massless_State_RMs()
{
  std::vector<int> RM = Alpha_RM_Numerator();

  Select_F_Operator(0, RM, Compute_Mass(RM));
}//Close Build_Massless_State_RMs.

//DEBUG.
void StateRMBuilder::Display_Massless_State_RMs() const
{
  std::cout<<"Massless state RMs: "<<Massless_State_RMs().size();
  std::list<std::vector<int> >::const_iterator itRM = 
    Massless_State_RMs_.begin();
  for(; itRM != Massless_State_RMs_.end(); itRM++)
    {
      for(int b=0; b<static_cast<int>(itRM->size()); b++)
	  std::cout<<itRM->at(b)<<" ";
      std::cout<<std::endl;
    }//Close for loop on Massless_State_RMs()
  std::cout<<std::endl;
}//Close Display_Massless_State_RMs.

//PRIVATE.
//HELPERS.
void StateRMBuilder::Build_Reverse_Fermion_Mode_Map()
{
  std::map<int, int>::iterator itMap = Fermion_Mode_Map_.begin();
  std::map<int, int>::iterator itMap_End = Fermion_Mode_Map_.end();
  for(; itMap != itMap_End; ++itMap)
    Reverse_Fermion_Mode_Map_[itMap->second] = itMap->first;
}//Close Build_Reverse_Fermion_Mode_Map.

bool StateRMBuilder::In_Map(int element)
{
  element += LM_Size();
  return Fermion_Mode_Map_.find(element) != Fermion_Mode_Map_.end();
}//Close In_Map.

bool StateRMBuilder::Real_RM_Mode(int element)
{
  element += LM_Size();
  if(Reverse_Fermion_Mode_Map_.find(element)
     != Reverse_Fermion_Mode_Map_.end())
    {
      if(Reverse_Fermion_Mode_Map_[element] < LM_Size())
	return true;
      else 
	return false;
    }else
    return false;
}//Close Real_RM_Mode.

void StateRMBuilder::Select_F_Operator(int element, std::vector<int>& RM, int Mass)
{

  if(element < static_cast<int>(RM.size())-1)//Because the last element is the second of a complex
    //pair.
    {
      if(In_Map(element))//Complex element.
	{
	  bool Becomes_Massive_Raised = (Mass_Increase_Raise(element) + Mass > 
					 Mass_Limit());
	  bool Becomes_Massive_Lowered = (Mass_Increase_Lower(element) + Mass > 
					  Mass_Limit());

	  if(Becomes_Massive_Raised && Becomes_Massive_Lowered)
	    Select_F_Operator(element+1, RM, Mass);//Already massive, move on.
	  else if(!Becomes_Massive_Lowered && !Becomes_Massive_Raised)
	    Apply_Complex_F_Operator(element, RM, -1, 1, Mass);//Apply both.
	  else if(!Becomes_Massive_Raised)
	    Apply_Complex_F_Operator(element, RM, 0, 1, Mass);//Apply raising.
	  else 
	    Apply_Complex_F_Operator(element, RM, -1, 0, Mass);//Apply lowering.

	}else if(Real_RM_Mode(element))//Real element.
	{

	  bool Becomes_Massive_Lowered = (Mass_Increase_Lower(element) + Mass >
					  Mass_Limit());

	  if(!Becomes_Massive_Lowered)
	    Apply_Real_F_Operator(element, RM, Mass);
	  else
	    Select_F_Operator(element+1, RM, Mass);
	}else //...or neither. Move on.
	Select_F_Operator(element+1, RM, Mass);

    }else if(Mass== Mass_Limit())
    Massless_State_RMs_.push_back(RM);
}//Close Select_F_Operator.

int StateRMBuilder::Compute_Mass(std::vector<int>& RM)
{
  int Mass = 0;
  //The F operator doubles the mass contribution of a real RM mode due to a 
  //degeneracy. Since these count as half of a complex mode towards the mass of
  //the state, doubling the mass contributions only increases the counting of the 
  //untwisted modes, giving them equal weight in the mass contribution.
  //Twisted periodic modes do not gain mass when lowered, and thus still count
  //as only half that of a complex mode.
  int RM_Size = RM.size();
  for(int a=0; a<RM_Size; ++a)
    {
      bool Is_Real = Real_RM_Mode(a);
      bool InMap = In_Map(a);

      if(InMap)
	Mass += (RM.at(a)*RM.at(a));
      else if(Is_Real && abs(RM.at(a)) == Alpha_RM_Denominator())
	Mass += (RM.at(a)*RM.at(a)/2);
      else if(Is_Real)
	Mass += (RM.at(a)*RM.at(a));
    }

  return Mass;
}//Close Compute_Mass.

int StateRMBuilder::Mass_Increase_Raise(int element)
{
  return 4*Alpha_RM_Denominator()*Alpha_RM_Denominator() + 
    2*Alpha_RM_Numerator().at(element)*2*Alpha_RM_Denominator();
}//Close Mass_Increase_Raise.

int StateRMBuilder::Mass_Increase_Lower(int element)
{
  return 4*Alpha_RM_Denominator()*Alpha_RM_Denominator() - 
    2*Alpha_RM_Numerator().at(element)*2*Alpha_RM_Denominator();
}//Close Mass_Increase_Lower.

void StateRMBuilder::Apply_Real_F_Operator(int element, std::vector<int> RM, 
					     int Mass)
{
  for(int a=-1; a<=0; a++)
    {
      RM.at(element) = Alpha_RM_Numerator().at(element) + 
	2*a*Alpha_RM_Denominator();

      int New_Mass = Mass + abs(a)*4*Alpha_RM_Denominator()*Alpha_RM_Denominator() + 
	a*2*Alpha_RM_Numerator().at(element)*2*Alpha_RM_Denominator();

      Select_F_Operator(element+1, RM, New_Mass);
    }//Close for loop for the real element.
}//Close Apply_Real_F_Operator.

void StateRMBuilder::Apply_Complex_F_Operator(int element, 
						std::vector<int> RM,
						int F_Lower, int F_Raise, int Mass)
{
  for(int a=F_Lower; a<=F_Raise; a++)
    {
      int Complex_Partner_Element = 
	Fermion_Mode_Map_[element+LM_Size()] - LM_Size();

      RM.at(element) = Alpha_RM_Numerator().at(element) + 
	2*a*Alpha_RM_Denominator();
      RM.at(Complex_Partner_Element) = 
	Alpha_RM_Numerator().at(Complex_Partner_Element) + 
	2*a*Alpha_RM_Denominator();

      int New_Mass = Mass + abs(a)*4*Alpha_RM_Denominator()*Alpha_RM_Denominator() + 
	a*2*Alpha_RM_Numerator().at(element)*2*Alpha_RM_Denominator();

      Select_F_Operator(element+1, RM, New_Mass);
    }//Close for loop for the complex element.
}//Close Apply_F_Operator.




