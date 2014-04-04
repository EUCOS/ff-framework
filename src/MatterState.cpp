#include "MatterState.h"

//PUBLIC.
MatterState::MatterState(const std::vector<int>& Numerator, 
			   int Denominator, int LM_Size,
			   const std::vector<GroupRepresentation>& Representations):
  State(Numerator, Denominator, LM_Size)
{
  Representations_ = Representations;
}//Close constructor.

MatterState::MatterState(const MatterState& New_MatterState):
  State(New_MatterState)
{
  Representations_ = New_MatterState.Representations();
}//Close copy constructor.

//DEBUG.
void MatterState::Display_Representations() const
{
  for(int a=0; a<static_cast<int>(Representations().size()); a++)
		Representations().at(a).Display();
  std::cout<<std::endl;
}//Close Display_Representations.
