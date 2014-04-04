#include "FermionModeMapBuilder.h"

//PUBLIC.
FermionModeMapBuilder::FermionModeMapBuilder()
{
  Large_ST_Dimensions_ = 0;
  Consistent_Pairings_ = true;
}//Close constructor.

FermionModeMapBuilder::FermionModeMapBuilder
(int Large_ST_Dimensions)
{
  Large_ST_Dimensions_ = Large_ST_Dimensions;
  Consistent_Pairings_ = true;
}//Close second constructor.

FermionModeMapBuilder::FermionModeMapBuilder
(const FermionModeMapBuilder& New_FermionModeMapBuilder)
{
  Large_ST_Dimensions_ = New_FermionModeMapBuilder.Large_ST_Dimensions();
  Fermion_Mode_Map_ = New_FermionModeMapBuilder.Fermion_Mode_Map();
  Consistent_Pairings_ = New_FermionModeMapBuilder.Consistent_Pairings();
}//Close copy constructor.

//INTERFACE.
void FermionModeMapBuilder::Build_Fermion_Mode_Map
(const std::vector<BasisAlpha>& Common_Basis_Alphas)
{
  if(Large_ST_Dimensions() == 0)
    Large_ST_Dimensions_ = Compute_Large_ST_Dimensions(Common_Basis_Alphas.at(0));

  Map_Complex_LM_Elements(Common_Basis_Alphas);
  Map_Complex_RM_Elements(Common_Basis_Alphas);
  std::vector<int> LR_Coordinates;
  
  //Add the LM to LR_Coordinates.
  int LM_Compact_Start = Large_ST_Dimensions();
  int LM_Compact_End = Common_Basis_Alphas.at(0).LM_Size();
  //Skip the complex ST Fermions -> Large_ST_Dimensions-2
  //Start at index 0 -> Large_ST_Dimensions-2-1
  //Start at the 'w' of the next triplet ->Large_ST_Dimensions -2 -1 +3 
  //->Large_ST_Dimensions.
  for(int a=LM_Compact_Start; a<LM_Compact_End; a+=3)
    {
      //x is a-2.
      LR_Coordinates.push_back(a-1);//y
      LR_Coordinates.push_back(a);//w
    }//Close for loop on Common_Basis_Alphas().at(0).LM_Size().

  //Add the Compact RMs to LR_Coordinates.
  int RM_Compact_Start = Common_Basis_Alphas.at(0).LM_Size()+16;
  //Skip the left movers and the Observable elements.
  int RM_Compact_End = RM_Compact_Start + 2*(10-Large_ST_Dimensions());
  for(int a=RM_Compact_Start; a<RM_Compact_End; a++)
    LR_Coordinates.push_back(a);

  Find_All_LR_Pairs(LR_Coordinates, Common_Basis_Alphas);
}//Close Build_Fermion_Mode_Map.

//DEBUG.
void FermionModeMapBuilder::Display_Fermion_Mode_Map() const
{
  std::cout<<"Fermion mode map:"<<std::endl;
  std::map<int, int>::const_iterator itFermion_Mode_Map = Fermion_Mode_Map_.begin();
  for(; itFermion_Mode_Map != Fermion_Mode_Map_.end(); 
      itFermion_Mode_Map++)
    std::cout<<itFermion_Mode_Map->first<<" "<<itFermion_Mode_Map->second<<std::endl;
  std::cout<<std::endl;
}//Close Display_Fermion_Mode_Map.

//PRIVATE.
//HELPERS.
int FermionModeMapBuilder::Compute_Large_ST_Dimensions(const BasisAlpha& 
							  Common_Basis_Alpha)
{
  return (14 - (Common_Basis_Alpha.LM_Size()/2));
}//Close Compute_Large_ST_Dimensions.

void FermionModeMapBuilder::Map_Complex_LM_Elements
(const std::vector<BasisAlpha>& Common_Basis_Alphas)
{
  //Complex ST Fermions.
  for(int a=1; a<(Large_ST_Dimensions()-2); a+=2)
      Fermion_Mode_Map_[a-1] = a;

  //x values.
  for(int a=(Large_ST_Dimensions()-2)+3; a<Common_Basis_Alphas.at(0).LM_Size(); 
      a+=6)
    Fermion_Mode_Map_[a-3] = a;
}//Close Map_Complex_LM_Elements.

void FermionModeMapBuilder::Map_Complex_RM_Elements
(const std::vector<BasisAlpha>& Common_Basis_Alphas)
{
  //First, the Observable elements.
  int RM_Observable_Start = Common_Basis_Alphas.at(0).LM_Size()+1;
  int RM_Observable_End = Common_Basis_Alphas.at(0).LM_Size()+16;
  for(int a=RM_Observable_Start; a<RM_Observable_End; a+=2)
    Fermion_Mode_Map_[a-1] = a;

  //Now the Hidden elements.
  int RM_Hidden_Start = Common_Basis_Alphas.at(0).LM_Size()+16+
    (2*(10-Large_ST_Dimensions()))+1;
  int RM_Hidden_End = Common_Basis_Alphas.at(0).Numerator().size();
  for(int a=RM_Hidden_Start; a<RM_Hidden_End; a+=2)
    Fermion_Mode_Map_[a-1] = a;
}//Close Map_Complex_RM_Elements.

void FermionModeMapBuilder::Find_All_LR_Pairs(std::vector<int> LR_Coordinates,
						 const std::vector<BasisAlpha>& 
						 Common_Basis_Alphas)
{
  //Firstly, pair the left movers.
  std::vector<int> Unpaired_LMs = Pair_LMs(LR_Coordinates, Common_Basis_Alphas);

  //Now pair the right movers.
  std::vector<int> Unpaired_RMs = Pair_RMs(LR_Coordinates, Common_Basis_Alphas);

  //Finally, pair the LMs and RMs that don't have matching boundary conditions
  //on their own side.
  std::vector<int> Unpaired_Mixed = Pair_Mixed(Unpaired_LMs, Unpaired_RMs, 
					       Common_Basis_Alphas);

  if(Unpaired_Mixed.size() != 0)
    Consistent_Pairings_ = false;
}//Close Find_All_LR_Pairs.

std::vector<int> FermionModeMapBuilder::Pair_LMs(std::vector<int> LR_Coordinates,
						    const std::vector<BasisAlpha>&
						    Common_Basis_Alphas)
{
  std::vector<std::vector<int> > Initial_Matching_BCs_LM;
  std::vector<int> Initial_Matching_BCs_LM_Loader;

  for(int a=0; a<(static_cast<int>(LR_Coordinates.size())/2); a++)
    Initial_Matching_BCs_LM_Loader.push_back(LR_Coordinates.at(a));
  Initial_Matching_BCs_LM.push_back(Initial_Matching_BCs_LM_Loader);
  Initial_Matching_BCs_LM_Loader.clear();

  std::vector<std::vector<int> > Final_Matching_BCs_LM = Find_Matching_BCs
    (Initial_Matching_BCs_LM, Common_Basis_Alphas);

  std::vector<int> Unpaired_LMs = Add_Pairs_To_Map(Final_Matching_BCs_LM);

  return Unpaired_LMs;
}//Close Pair_LMs.

std::vector<int> FermionModeMapBuilder::Pair_RMs(std::vector<int> LR_Coordinates,
						    const std::vector<BasisAlpha>& 
						    Common_Basis_Alphas)
{
 std::vector<std::vector<int> > Initial_Matching_BCs_RM;
  std::vector<int> Initial_Matching_BCs_RM_Loader;

  for(int a=(LR_Coordinates.size()/2); a<static_cast<int>(LR_Coordinates.size()); a++)
    Initial_Matching_BCs_RM_Loader.push_back(LR_Coordinates.at(a)); 
  Initial_Matching_BCs_RM.push_back(Initial_Matching_BCs_RM_Loader);
  Initial_Matching_BCs_RM_Loader.clear();

  std::vector<std::vector<int> > Final_Matching_BCs_RM = Find_Matching_BCs
    (Initial_Matching_BCs_RM, Common_Basis_Alphas);

  std::vector<int> Unpaired_RMs = Add_Pairs_To_Map(Final_Matching_BCs_RM);

  return Unpaired_RMs;
}//Close Pair_RMs.

std::vector<int> FermionModeMapBuilder::Pair_Mixed(std::vector<int> Unpaired_LMs,
						      std::vector<int> Unpaired_RMs,
						      const std::vector<BasisAlpha>&
						      Common_Basis_Alphas)
{
  std::vector<std::vector<int> > Initial_Matching_BCs_Mixed;
  std::vector<int> Initial_Matching_BCs_Mixed_Loader;

  for(int a=0; a<static_cast<int>(Unpaired_LMs.size()); a++)
    Initial_Matching_BCs_Mixed_Loader.push_back(Unpaired_LMs.at(a));
  for(int a=0; a<static_cast<int>(Unpaired_RMs.size()); a++)
    Initial_Matching_BCs_Mixed_Loader.push_back(Unpaired_RMs.at(a));

  Initial_Matching_BCs_Mixed.push_back(Initial_Matching_BCs_Mixed_Loader);
  Initial_Matching_BCs_Mixed_Loader.clear();

  std::vector<std::vector<int> > Final_Matching_BCs_Mixed = Find_Matching_BCs
    (Initial_Matching_BCs_Mixed, Common_Basis_Alphas);

  std::vector<int> Unpaired_Mixed = Add_Pairs_To_Map(Final_Matching_BCs_Mixed);

  return Unpaired_Mixed;
}//Close Pair_Mixed.

std::vector<std::vector<int> > FermionModeMapBuilder::Find_Matching_BCs
(std::vector<std::vector<int> > Matching_BCs, const std::vector<BasisAlpha>& 
 Common_Basis_Alphas)
{
  int Maximum_Element_Value = Common_Basis_Alphas.at(0).Denominator();
  int Minimum_Element_Value = -1*Maximum_Element_Value;
  for(int CBA=0; CBA<static_cast<int>(Common_Basis_Alphas.size()); CBA++)
    {
      std::vector<std::vector<int> > New_Matching_BCs;
      for(int MBC_Row=0; MBC_Row<static_cast<int>(Matching_BCs.size()); MBC_Row++)
	{
	  for(int E_Val=Minimum_Element_Value; E_Val<=Maximum_Element_Value; E_Val++)
	    {
	      std::vector<int> New_Matching_BCs_Loader;
	      for(int MBC_Col=0; MBC_Col<static_cast<int>(Matching_BCs.at(MBC_Row).size()); MBC_Col++)
		{
		  if(Common_Basis_Alphas.at(CBA).Numerator().at
		     (Matching_BCs.at(MBC_Row).at(MBC_Col)) == E_Val)
		    New_Matching_BCs_Loader.push_back
		      (Matching_BCs.at(MBC_Row).at(MBC_Col));
		}//Close for loop on Matching_BVs.at(b).size().
	      if(New_Matching_BCs_Loader.size()  != 0)
		New_Matching_BCs.push_back(New_Matching_BCs_Loader);
	      New_Matching_BCs_Loader.clear();
	    }//Close for loop on Maximum_Element_Value.
	}//Close for loop on Matching_BCs.size().
      Matching_BCs.swap(New_Matching_BCs);
      New_Matching_BCs.clear();
    }//Close for loop on Common_Basis_Alphas().size().
  return Matching_BCs;
}//Close Find_Pairs.

std::vector<int> FermionModeMapBuilder::Add_Pairs_To_Map
(std::vector<std::vector<int> > Final_Matching_BCs)
{
  std::vector<int> Unpaired_BCs;
  for(int a=0; a<static_cast<int>(Final_Matching_BCs.size()); a++)
    {
      for(int b=1; b<static_cast<int>(Final_Matching_BCs.at(a).size()); b+=2)
	Fermion_Mode_Map_[Final_Matching_BCs.at(a).at(b-1)] = 
	  Final_Matching_BCs.at(a).at(b);

      if(Final_Matching_BCs.at(a).size()%2 == 1)
	Unpaired_BCs.push_back(Final_Matching_BCs.at(a).back());
    }//Close for loop on Final_Matching_BCs.size().
  return Unpaired_BCs;
}//Close Add_Pairs_To_Map.


