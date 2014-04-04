#include "BasisAlphaBuilder.h"

BasisAlphaBuilder::BasisAlphaBuilder(const Model& FFHS_Model)
{
  Basis_Vectors_ = FFHS_Model.BV_Set();
}//Close constructor.

BasisAlphaBuilder::BasisAlphaBuilder(const std::vector<BasisVector>& 
					 Basis_Vectors)
{
  Basis_Vectors_ = Basis_Vectors;
}//Close constructor.

BasisAlphaBuilder::BasisAlphaBuilder
(const BasisAlphaBuilder& New_BasisAlphaBuilder)
{
  Basis_Vectors_ = New_BasisAlphaBuilder.Basis_Vectors();
  Basis_Alphas_ = New_BasisAlphaBuilder.Basis_Alphas();
  Common_Basis_Alphas_ = New_BasisAlphaBuilder.Common_Basis_Alphas();
}//Close copy constructor.

//INTERFACE.
void BasisAlphaBuilder::Build_Basis_Alphas()
{
  for(int a=0; a<static_cast<int>(Basis_Vectors().size()); a++)
    {
      int New_Order = Basis_Vectors().at(a).Order();
      if(Basis_Vectors().at(a).BV().at(0) == 1)//If there's a left moving ST Fermion.
	New_Order = FF::LCM(2, Basis_Vectors().at(a).Order());
      int Numerator_Conversion = 2*New_Order / Basis_Vectors().at(a).Order();
      //The proper form for a basis alpha is 2m/N, where N is the order and 
      //m is the basis vector's integer code.
      std::vector<int> Basis_Alpha_Num_Loader;
      //Left mover.
      for(int b=0; b<Basis_Vectors().at(a).LM_Size(); b++)
	Basis_Alpha_Num_Loader.push_back(New_Order*Basis_Vectors().at(a).BV().at(b));
      //Right movers.
      for(int b=Basis_Vectors().at(a).LM_Size(); 
	  b<static_cast<int>(Basis_Vectors().at(a).BV().size());  b++)
	{
	  if(Basis_Vectors().at(a).BV().at(b)<=
	     (double(Basis_Vectors().at(a).Order())/double(2)))
	    {
	      Basis_Alpha_Num_Loader.push_back
		(Numerator_Conversion*Basis_Vectors().at(a).BV().at(b));
	    }else
	    {
	      Basis_Alpha_Num_Loader.push_back(Numerator_Conversion*
					       (Basis_Vectors().at(a).BV().at(b) - 
						Basis_Vectors().at(a).Order()));
	    }
	}//Close for loop on Right mover.
      Basis_Alphas_.push_back(BasisAlpha(Basis_Alpha_Num_Loader, New_Order, 
					  Basis_Vectors().at(a)));
    }//Close for loop on BV_Set().size().
}//Close Build_Basis_Alphas.

void BasisAlphaBuilder::Build_Common_Basis_Alphas()
  {
  //First, the common denominator needs to be found for the entire
  //set of basis alphas. The loop starts with the lowest nontrivial order, 2.

  std::vector<int> Common_Basis_Alpha_Num_Loader;

  int Common_Denom = 2;
  for(int a=0; a<static_cast<int>(Basis_Alphas().size()); a++)
    Common_Denom = FF::LCM(Common_Denom, Basis_Alphas().at(a).Denominator());

  //Now the conversion factor needed for the numerators is used to build the 
  //numerators.

  for(int a=0; a<static_cast<int>(Basis_Alphas().size()); a++)
    {
      int Numerator_Conversion = Common_Denom / Basis_Alphas().at(a).Denominator();

      for(int b=0; b<static_cast<int>(Basis_Alphas().at(a).Numerator().size()); b++)
	{
	  int Converted_Element = Basis_Alphas().at(a).Numerator().at(b)*
	    Numerator_Conversion;

	  Common_Basis_Alpha_Num_Loader.push_back(Converted_Element);
	}//Close for loop on the individual vector.

      Common_Basis_Alphas_.push_back(BasisAlpha(Common_Basis_Alpha_Num_Loader,
						 Common_Denom));
      Common_Basis_Alpha_Num_Loader.clear();
    }//Close for loop on Basis_Alpha_Num.
}//Close Build_Common_Basis_Alphas.

//DEBUG.
void BasisAlphaBuilder::Display_Basis_Alphas() const
{
  std::cout<<"Basis Alphas: "<<Basis_Alphas().size()<<std::endl;
  for(int a=0; a<static_cast<int>(Basis_Alphas().size()); a++)
    {
      Basis_Alphas().at(a).Display();
    }//Close for loop on Basis Alphas.
  std::cout<<std::endl;
}//Close Display_Basis_Alphas.

void BasisAlphaBuilder::Display_Common_Basis_Alphas() const
{
  std::cout<<"Common Basis Alphas: "<<Common_Basis_Alphas().size()<<std::endl;
  for(int a=0; a<static_cast<int>(Common_Basis_Alphas().size()); a++)
    {
      Common_Basis_Alphas().at(a).Display();
    }//Close for loop on Common_Basis_Alphas.
  std::cout<<std::endl;
}//Close Display_Common_Basis_Alphas.
