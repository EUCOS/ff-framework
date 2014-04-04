#include "ModularInvarianceChecker.h"

//PUBLIC.

ModularInvarianceChecker::ModularInvarianceChecker
(const ModularInvarianceChecker& New_ModularInvarianceChecker)
{
  ;
}//Close copy constructor.

bool ModularInvarianceChecker::Test_Modular_Invariance
(const std::vector<BasisAlpha>& Basis_Alpha_Set)
{
  return (Check_Dot_Products(Basis_Alpha_Set)&&
	  Check_Simultaneous_Periodic_Modes(Basis_Alpha_Set));
}//Close Test_Modular_Invariance.

//PRIVATE.
//HELPERS.
bool ModularInvarianceChecker::Check_Dot_Products
(const std::vector<BasisAlpha>& Basis_Alpha_Set)
{
  bool Pass = true;
  for(int a=0; a<static_cast<int>(Basis_Alpha_Set.size());a++)
    {
      for(int b=0; b<static_cast<int>(Basis_Alpha_Set.size()); b++)
	{
	  if(Pass == true)
	    {
	      int N = FF::LCM(Basis_Alpha_Set.at(a).Denominator(),
			      Basis_Alpha_Set.at(b).Denominator());
	      int Denom = Basis_Alpha_Set.at(a).Denominator()*
		Basis_Alpha_Set.at(b).Denominator();
	      int Num_Dot = Basis_Alpha_Set.at(a).Lorentz_Dot
		(Basis_Alpha_Set.at(b));

	      //MAKE SURE THIS RULE FOR EVENNESS APPLIES TO THE
	      //COMBINED ORDER RATHER THAN THE RM ORDER.
	      bool Is_Even = ((Basis_Alpha_Set.at(a).Denominator()%2==0)||
			      (Basis_Alpha_Set.at(b).Denominator()%2==0));
	      if((a==b)&&Is_Even)
		Pass = ((N*Num_Dot)%(16*Denom)==0);
	      else 
		Pass = ((N*Num_Dot)%(8*Denom)==0);
	    }//Close if statement.
	  else
	    {
	      return Pass;
	    }
	}//Close inner for loop on Basis_Alpha_Set.
    }//Close for loop on Basis_Alpha_Set.
  return Pass;
}//Close Check_Dot_Products.

bool ModularInvarianceChecker::Check_Simultaneous_Periodic_Modes
(const std::vector<BasisAlpha>& Basis_Alpha_Set)
{
  bool Pass = true;

  for(int a = 0; a<static_cast<int>(Basis_Alpha_Set.size()); a++)
    {
      for(int b = 0; b<static_cast<int>(Basis_Alpha_Set.size()); b++)
	{
	  for(int c = b; c<static_cast<int>(Basis_Alpha_Set.size()); c++)
	    {
	      int Simultaneous_Periodic = 0;

	      if(Pass == true)
		{
		  for(int d=0; d<static_cast<int>(Basis_Alpha_Set.at(a).Numerator().size()); d++)
		    {

		      double sum = Basis_Alpha_Set.at(a).Numerator().at(d) + 
			Basis_Alpha_Set.at(b).Numerator().at(d) +
			Basis_Alpha_Set.at(c).Numerator().at(d);

		      double denom = Basis_Alpha_Set.at(a).Denominator() + 
			Basis_Alpha_Set.at(b).Denominator() + 
			Basis_Alpha_Set.at(c).Denominator();

		      if(sum==denom)
			Simultaneous_Periodic++;
		    }//Close for loop on basis alpha elements.
		  Pass = (Simultaneous_Periodic%2==0);
		}//Close if statement on Pass.
	      else
		{
		  return Pass;
		}
	    }//Close inner for loop on Basis_Alpha_Set.
	}//Close middle for loop on Basis_Alpha_Set.
    }//Close outer for loop on Basis_Alpha_Set.
  return Pass;
}//Close Check_Simultaneous_Periodic_Modes.
