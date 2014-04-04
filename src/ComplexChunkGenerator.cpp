#include "ComplexChunkGenerator.h"

//PUBLIC.
ComplexChunkGenerator::ComplexChunkGenerator
(const std::vector<BasisAlpha>& Common_Basis_Alphas, int Order, 
 const std::vector<int>& Complex_Elements)
{
  Common_Basis_Alphas_ = Common_Basis_Alphas;
  Order_ = Order;
  Complex_Elements_ = Complex_Elements;
}//Close constructor.

ComplexChunkGenerator::ComplexChunkGenerator
(const std::vector<BasisAlpha>& Common_Basis_Alphas, int Order,
 const std::vector<int>& Complex_Elements, bool First_Complex_Extension)
{
  Common_Basis_Alphas_ = Common_Basis_Alphas;
  Order_ = Order;
  Complex_Elements_ = Complex_Elements;
  First_Complex_Extension_ = First_Complex_Extension;
}//Close second extension.

ComplexChunkGenerator::ComplexChunkGenerator(const ComplexChunkGenerator&
						 New_ComplexChunkGenerator)
{
  Common_Basis_Alphas_ = New_ComplexChunkGenerator.Common_Basis_Alphas();
  Order_ = New_ComplexChunkGenerator.Order();
  Complex_Elements_ = New_ComplexChunkGenerator.Complex_Elements();
  First_Complex_Extension_ = New_ComplexChunkGenerator.First_Complex_Extension();
  Complex_Chunks_ = New_ComplexChunkGenerator.Complex_Chunks();
}//Close copy constructor.

//INTERFACE.
void ComplexChunkGenerator::Build_Complex_Chunks()
{
  std::vector<std::vector<int> > Matching_BCs = Find_Matching_BCs();
  std::vector<int> Complex_Chunks_Loader(16,0);
  Generate_Complex_Chunks(Matching_BCs, 0, 0, Complex_Chunks_Loader);
}//Close Build_Complex_Chunks.

//DEBUG.
void ComplexChunkGenerator::Display_Complex_Chunks() const
{
  std::cout<<"Complex chunks: "<<Complex_Chunks().size()<<std::endl;
  std::set<std::vector<int> >::const_iterator itComplex_Chunks = 
    Complex_Chunks_.begin();
  for(; itComplex_Chunks != Complex_Chunks_.end(); 
      ++itComplex_Chunks)
    {
      for(int a=0; a<static_cast<int>(itComplex_Chunks->size()); a++)
	std::cout<<itComplex_Chunks->at(a)<<" ";
      std::cout<<std::endl;
    }
}//Close Display_Complex_Chunks.

//PRIVATE.
std::vector<std::vector<int> > ComplexChunkGenerator::Find_Matching_BCs()
{
  int Maximum_Element_Value = Common_Basis_Alphas().at(0).Denominator();
  int Minimum_Element_Value = -1*Maximum_Element_Value;

  std::vector<std::vector<int> > Matching_BCs;
  Matching_BCs.push_back(Complex_Elements());

  for(int CBA=0; CBA<static_cast<int>(Common_Basis_Alphas().size()); CBA++)
    {
      std::vector<std::vector<int> > New_Matching_BCs;
      for(int MBC_Row=0; MBC_Row<static_cast<int>(Matching_BCs.size()); MBC_Row++)
	{
	  for(int E_Val=Minimum_Element_Value; E_Val<=Maximum_Element_Value; E_Val++)
	    {
	      std::vector<int> New_Matching_BCs_Loader;
	      for(int MBC_Col=0; MBC_Col<static_cast<int>(Matching_BCs.at(MBC_Row).size()); MBC_Col++)
		{
		  if(Common_Basis_Alphas().at(CBA).Numerator().at
		     (Matching_BCs.at(MBC_Row).at(MBC_Col)) == E_Val)
		    New_Matching_BCs_Loader.push_back
		      (Matching_BCs.at(MBC_Row).at(MBC_Col));
		}//Close for loop on Matching_BCs.at(MBC_Row).size().
	      if(New_Matching_BCs_Loader.size() != 0)
		New_Matching_BCs.push_back(New_Matching_BCs_Loader);
	      New_Matching_BCs_Loader.clear();
	    }//Close for loop on Maximum_Element_Value.
	}//Close for loop on Matching_BCs.size().
      Matching_BCs.swap(New_Matching_BCs);
      New_Matching_BCs.clear();
    }//Close for loop on Common_Basis_Alphas.
  return Matching_BCs;
}//Close Find_Matching_BCs.

void ComplexChunkGenerator::Generate_Complex_Chunks
(const std::vector<std::vector<int> >& Matching_BCs, int MBC_Row_Tracker, 
 int MBC_Col_Tracker, std::vector<int> Complex_Chunks_Loader)
{
  if(MBC_Row_Tracker<static_cast<int>(Matching_BCs.size()))
    {
      if(MBC_Col_Tracker<static_cast<int>(Matching_BCs.at(MBC_Row_Tracker).size()))
	{
	  int Full_BV_Adjust = Complex_Elements().at(0);
	  //Compensates for the size difference between the Complex_Chunks_Loader 
	  //vector and the Common_Basis_Alphas elements indexed b Matching_BCs.
	  int E_Val_Start = 0;
	  int E_Val_End = Order()-1;
	  if(First_Complex_Extension())
	    E_Val_End = Order()/2;//Rounds down to an integer by truncation.
	  if(MBC_Col_Tracker != 0)      
	    E_Val_Start = Complex_Chunks_Loader.at(Matching_BCs.at
						   (MBC_Row_Tracker).at
						   (MBC_Col_Tracker-1)-
						   Full_BV_Adjust);
	  for(int E_Val = E_Val_Start; E_Val<=E_Val_End; E_Val++)
	    {
	      Complex_Chunks_Loader.at(Matching_BCs.at(MBC_Row_Tracker).at
				       (MBC_Col_Tracker)-Full_BV_Adjust) = E_Val;
	      Complex_Chunks_Loader.at(Matching_BCs.at(MBC_Row_Tracker).at
				       (MBC_Col_Tracker+1)-Full_BV_Adjust) = E_Val;
	      Generate_Complex_Chunks(Matching_BCs, MBC_Row_Tracker, 
				      MBC_Col_Tracker+2, Complex_Chunks_Loader);
	    }//Close for loop on E_Val.
	}else//Close if/else on MBC_Col_Tracker.
	Generate_Complex_Chunks(Matching_BCs, MBC_Row_Tracker+1, 0, 
				Complex_Chunks_Loader);
    }else//Close if/else on MBC_Row_Tracker.
    Complex_Chunks_.insert(Complex_Chunks_Loader);
}//Close Generate_Complex_Chunks.
