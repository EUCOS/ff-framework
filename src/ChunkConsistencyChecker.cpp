#include "ChunkConsistencyChecker.h"
//PUBLIC.
ChunkConsistencyChecker::ChunkConsistencyChecker(int BV_Denom,
						     int CBA_Denom)
{
  BV_Denom_ = BV_Denom;
  CBA_Denom_ = CBA_Denom;
}//Close constructor.

ChunkConsistencyChecker::ChunkConsistencyChecker(const ChunkConsistencyChecker&
						     New_ChunkConsistencyChecker)
{
  BV_Denom_ = New_ChunkConsistencyChecker.BV_Denom();
  CBA_Denom_ = New_ChunkConsistencyChecker.CBA_Denom();
}//Close copy constructor.

//INTERFACE.
bool ChunkConsistencyChecker::Check_Modular_Invariance(const Chunk& LM, 
							 const Chunk& Obs,
							 const Chunk& Comp, 
							 const Chunk& Hid)
{
  for(int a=0; a<static_cast<int>(LM.MI_Dot_Products().size()); a++)//Any chunk will do; they
    //should be the same unless D=10, in which case there are no RM Compact chunks.
    {
      int Dot = LM.MI_Dot_Products().at(a) - (Obs.MI_Dot_Products().at(a)+
					      Comp.MI_Dot_Products().at(a)+
					      Hid.MI_Dot_Products().at(a));
      if(a==0&&(BV_Denom()%2==0))
	Dot%=(16*BV_Denom()*BV_Denom());
      else if(a==0)
	Dot %=(8*BV_Denom()*BV_Denom());
      else
	Dot%=(8*BV_Denom()*CBA_Denom());

      if(Dot != 0)
	return false;
    }//Close for loop on the Dot Products.
  return true;
}//Close Check_Modular_Invariance.

bool ChunkConsistencyChecker::Check_D10_Modular_Invariance(const Chunk& LM,
							     const Chunk& Obs,
							     const Chunk& Hid)
{
  //This function is needed because D=10 has no compact chunks, so the regular
  //Check_Modular_Invariance function will throw an exception.
 for(int a=0; a<static_cast<int>(LM.MI_Dot_Products().size()); a++)//Any chunk will do; they
    //should be the same.
    {
      int Dot = LM.MI_Dot_Products().at(a) - (Obs.MI_Dot_Products().at(a)+
					      Hid.MI_Dot_Products().at(a));
      if(a==0&&(BV_Denom()%2==0))
	Dot%=(16*BV_Denom()*BV_Denom());
      else
	Dot%=(8*BV_Denom()*CBA_Denom());

      if(Dot != 0)
	return false;
    }//Close for loop on the Dot Products.
  return true;
}//Close Check_D10_Modular_Invariance.

bool ChunkConsistencyChecker::Check_Simultaneous_Periodic_Modes(const Chunk& LM,
								  const Chunk& Comp)
{
  bool Even_Simultaneous_Periodic_Modes = true;
  for(int a=0; a<static_cast<int>(LM.Simultaneous_Periodic_Modes().size()); a++)//Should 
    //be the same size for both.
    {
      Even_Simultaneous_Periodic_Modes = !(LM.Simultaneous_Periodic_Modes().at(a)-
					   Comp.Simultaneous_Periodic_Modes().at(a));
      //The !(LM-Comp) boolean equation matches the truth table as follows:
      //even+even = even, odd+odd = even, odd+even = odd.
      if(Even_Simultaneous_Periodic_Modes == false)
	return false;
    }//Close for loop on Simultaneous_Periodic_Modes().size().
  return true;
}//Close Check_Simultaneous_Periodic_Modes.

bool ChunkConsistencyChecker::Check_BV_Consistency(const Chunk& LM,
						     const Chunk& Obs, 
						     const Chunk& Comp,
						     const Chunk& Hid)
{
  if(!Check_Modular_Invariance(LM,Obs,Comp,Hid))
    return false;
  else if(!Check_Simultaneous_Periodic_Modes(LM,Comp))
    return false;
 
  return true;
}//Close Check_BV_Consistency.

bool ChunkConsistencyChecker::Check_D10_BV_Consistency(const Chunk& LM,
							 const Chunk& Obs,
							 const Chunk& Hid)
{
  return Check_D10_Modular_Invariance(LM,Obs,Hid);
}//Close Check_D10_BV_Consistency.
