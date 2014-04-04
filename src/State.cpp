#include "State.h"

//PUBLIC.
State::State(const std::vector<int>& Numerator, int Denominator, int LM_Size)
{
  Numerator_ = Numerator;
  Denominator_ = Denominator;
  LM_Size_ = LM_Size;
}//Close constructor.

State::State(const State& New_State)
{
  Numerator_ = New_State.Numerator();
  Denominator_ = New_State.Denominator();
  LM_Size_ = New_State.LM_Size();
  Length_Squared_Numerator_ = New_State.Length_Squared_Numerator();
  Length_Squared_Denominator_ = New_State.Length_Squared_Denominator();
}//Close copy constructor.

//INTERFACE.
void State::Calculate_Length_Squared(const std::map<int, int>& Fermion_Mode_Map)
{
  int Dot = 0;
  //Next line may cause issues.
  std::map<int, int>::const_iterator itMap = Fermion_Mode_Map.find(LM_Size());
  for(; itMap != Fermion_Mode_Map.end(); itMap++)
      Dot += Numerator().at(itMap->first)*Numerator().at(itMap->first);
    
  Length_Squared_Numerator_ = Dot;
  Length_Squared_Denominator_ = Denominator()*Denominator();
  
}//Close Calculate_Length_Squared.

bool State::Is_Positive(const std::map<int, int>& Fermion_Mode_Map)
{
  std::map<int, int>::const_iterator itMap = Fermion_Mode_Map.find(LM_Size());
  for(; itMap != Fermion_Mode_Map.end(); ++itMap)
    {
      if(Numerator().at(itMap->first) > 0)
	return true;
      else if(Numerator().at(itMap->first) < 0)
	return false;
    }//Close for loop on Fermion_Mode_Map.
  return false;
}//Close Is_Positive.

//DEBUG.
void State::Display() const
{
  std::cout<<Denominator()<<": ";

  for(int a=0; a<LM_Size(); a++)
    std::cout<<Numerator().at(a)<<" ";

  std::cout<<"|| ";

  for(int a=LM_Size(); a<static_cast<int>(Numerator().size()); a++)
    std::cout<<Numerator().at(a)<<" ";

  std::cout<<std::endl;
}//Close Display.
