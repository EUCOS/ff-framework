#include "ModelBuilder.h"

//PUBLIC.
ModelBuilder::ModelBuilder(int Large_ST_Dimensions)
{
  Model The_Model(Large_ST_Dimensions);
  FFHS_Model_ = The_Model;
  Consistent_Fermion_Mode_Pairs_ = true;
  Linearly_Independent_Alphas_ = true;
  Consistent_GSO_Matrix_ = true;
}//Close constructor.

ModelBuilder::ModelBuilder(const ModelBuilder& New_ModelBuilder)
{
  FFHS_Model_ = New_ModelBuilder.FFHS_Model();
  Consistent_Fermion_Mode_Pairs_ = 
		New_ModelBuilder.Consistent_Fermion_Mode_Pairs();
  Linearly_Independent_Alphas_ = 
		New_ModelBuilder.Linearly_Independent_Alphas();
  Consistent_GSO_Matrix_ = New_ModelBuilder.Consistent_GSO_Matrix();
  Basis_Alphas_ = New_ModelBuilder.Basis_Alphas();
  Common_Basis_Alphas_ = New_ModelBuilder.Common_Basis_Alphas();
  Alpha_Bosons_ = New_ModelBuilder.Alpha_Bosons();
  Alpha_Fermions_ = New_ModelBuilder.Alpha_Fermions();
  Alpha_SUSYs_ = New_ModelBuilder.Alpha_SUSYs();
  Fermion_Mode_Map_ = New_ModelBuilder.Fermion_Mode_Map();
  SUSY_States_ = New_ModelBuilder.SUSY_States();
}//Close copy constructor.

//INTERFACE.
void ModelBuilder::Load_S_Vector(int Large_ST_Dimensions)
{
  FFHS_Model_.Load_S_Vector(Large_ST_Dimensions);
}//Close Load_S_Vector.

void ModelBuilder::Load_Basis_Vector(const BasisVector& New_Basis_Vector)
{
  FFHS_Model_.Load_BV_Set(New_Basis_Vector);
}//Close Load_Basis_Vector.

void ModelBuilder::Load_k_ij_Row(const std::vector<int>& New_k_ij_Row)
{
  FFHS_Model_.Load_k_ij_Row(New_k_ij_Row);
}//Close Load_k_ij_Row.

void ModelBuilder::Load_Default_k_ij()
{
  std::vector<int> First_k_ij_Row(1, 1);
  Load_k_ij_Row(First_k_ij_Row);
  for(int a=1; a<static_cast<int>(FFHS_Model().BV_Set().size()); a++)
    {
      std::vector<int> k_ij_Row;
      for(int b=0; b<a; ++b)
	k_ij_Row.push_back(FFHS_Model().BV_Set().at(b).Order()/2);
      Load_k_ij_Row(k_ij_Row);
    }
}//Close Load_Default_k_ij.

void ModelBuilder::Load_NAHE_Set()
{
  NAHESetLoader NAHE_Loader;
  FFHS_Model_ = NAHE_Loader.Load_NAHE_Set(FFHS_Model(), 
					  (FFHS_Model().BV_Set().size() == 2));
}//Close Load_NAHE_Set.

void ModelBuilder::Load_NAHE_Variation()
{
  NAHEVariationLoader NAHE_Var_Loader;
  FFHS_Model_ = NAHE_Var_Loader.Load_NAHE_Variation(FFHS_Model(),
						    (FFHS_Model().BV_Set().size() 
						     == 2));
}//Close Load_NAHE_Variation.

bool ModelBuilder::Check_Modular_Invariance()
{
  if(Common_Basis_Alphas().size() == 0)
    Build_Basis_Alphas();
  ModularInvarianceChecker Checker;
  return Checker.Test_Modular_Invariance(Common_Basis_Alphas());
}//Close Check_Modular_Invariance.

bool ModelBuilder::Check_Linear_Independence()
{
  if(Basis_Alphas().size() == 0)
    { 
      Build_Basis_Alphas();
      Build_Alphas();
    }else if(Alpha_Bosons().size() == 0&&
	     Alpha_Fermions().size() == 0&&
	     Alpha_SUSYs().size() == 0)
    Build_Alphas();

  return Linearly_Independent_Alphas();
}//Close Check_Linear_Independence.

bool ModelBuilder::Check_k_ij_Consistency()
{
  if(Common_Basis_Alphas().size() == 0)
    Build_Basis_Alphas();
  if(FFHS_Model().k_ij().Numerators().front().size() !=
     FFHS_Model().k_ij().Numerators().size())
    Build_k_ij();
  return Consistent_GSO_Matrix();
}//Close Check_k_ij_Consistency.

bool ModelBuilder::Check_Model_Consistency()
{
  if(Check_Modular_Invariance())
    {
      if(Check_k_ij_Consistency())
	{   
	  if(Check_Linear_Independence())
	    return true;
	  else
	    return false;
	}else 
	return false;
    }else
    return false;
}//Close Check_Model_Consistency.

void ModelBuilder::Build_Gauge_Group_Model()
{
  if(Basis_Alphas().size() == 0)
    Build_Basis_Alphas();
  if(Alpha_SUSYs().size() == 0 && 
     Alpha_Fermions().size() == 0 &&
     Alpha_Bosons().size() == 0)
    Build_Alphas();
  if(Fermion_Mode_Map().size() == 0)
    Build_Fermion_Mode_Map();
  if(FFHS_Model().k_ij().Numerators().front().size() !=
     FFHS_Model().k_ij().Numerators().back().size())
    Build_k_ij();
  Build_Gauge_Groups();
  Compute_U1_Factors();
}//Close Build_Gauge_Group_Model.

void ModelBuilder::Build_Model()
{
  Build_Gauge_Group_Model();
  Build_SUSY_States();
  Build_MatterStates();
}//Close Build_Model.

//DEBUG.
void ModelBuilder::Display_Gauge_Group_Roots() const
{
  for(int a=0; a<static_cast<int>(FFHS_Model().Gauge_Groups().size()); a++)
    {
      std::list<State> Positive_Roots = 
	FFHS_Model().Gauge_Groups().at(a).Positive_Roots();
      FFHS_Model().Gauge_Groups().at(a).Display();
      std::cout<<"Positive roots for this gauge group: "<<Positive_Roots.size()<<
	std::endl;
      std::list<State>::iterator itPositive_Roots = Positive_Roots.begin();
      for(; itPositive_Roots != Positive_Roots.end(); 
	  ++itPositive_Roots)
	itPositive_Roots->Display();
    }
}//Close Display_Gauge_Group_Roots.

//PRIVATE.
//HELPERS.
void ModelBuilder::Build_Basis_Alphas()
{
  BasisAlphaBuilder BA_Builder(FFHS_Model());
  BA_Builder.Build_Basis_Alphas();
  BA_Builder.Build_Common_Basis_Alphas();

  Basis_Alphas_ = BA_Builder.Basis_Alphas();
  Common_Basis_Alphas_ = BA_Builder.Common_Basis_Alphas();
}//Close Build_Basis_Alphas.

void ModelBuilder::Build_Alphas()
{
  AlphaBuilder A_Builder(Common_Basis_Alphas(), Basis_Alphas());
  A_Builder.Build_Alphas();
  Alpha_Bosons_ = A_Builder.Alpha_Bosons();
  Alpha_Fermions_ = A_Builder.Alpha_Fermions();
  Alpha_SUSYs_ = A_Builder.Alpha_SUSYs();
  Linearly_Independent_Alphas_ = A_Builder.Linearly_Independent_Alphas();
}//Close Build_Alphas.

void ModelBuilder::Build_Fermion_Mode_Map()
{
  FermionModeMapBuilder Map_Builder;
  Map_Builder.Build_Fermion_Mode_Map(Common_Basis_Alphas());
  Fermion_Mode_Map_ = Map_Builder.Fermion_Mode_Map();
  Consistent_Fermion_Mode_Pairs_ = Map_Builder.Consistent_Pairings();
}//Close Build_Fermion_Mode_Map.

void ModelBuilder::Build_k_ij()
{
  GSOCoefficientMatrixBuilder k_ij_Builder(FFHS_Model(), Common_Basis_Alphas());
  k_ij_Builder.Build_Complete_GSO_Matrix();
  FFHS_Model_.Set_k_ij(k_ij_Builder.Complete_GSO_Matrix());
  Consistent_GSO_Matrix_ = k_ij_Builder.Consistent_GSO_Matrix();
}//Close Build_k_ij.

void ModelBuilder::Build_Gauge_Groups()
{
  //Consolidate the gauge states into a list.
  std::set<AlphaBoson>::iterator itBosons = Alpha_Bosons_.begin();
  std::list<State> Boson_States;
  for(; itBosons != Alpha_Bosons_.end(); ++itBosons)
    {
      StateBuilder Boson_Builder(*itBosons, Fermion_Mode_Map(), 
				  Common_Basis_Alphas(), FFHS_Model().k_ij());
      Boson_Builder.Build_States();

      std::list<State> Boson_State_Loader = Boson_Builder.States();
      std::list<State>::iterator itStates = Boson_State_Loader.begin();
      for(; itStates != Boson_State_Loader.end(); ++itStates)
	{
	  if(itStates->Is_Positive(Fermion_Mode_Map()))
	    {
	      State Boson_State = *itStates;
	      Boson_State.Calculate_Length_Squared(Fermion_Mode_Map());
	      Boson_States.push_back(Boson_State);
	    }//Close if statement on positivity of the state.
	}//Close for loop on states.
    }//Close for loop on the boson sectors.

  //Assemble the gauge states into groups of positive roots.
  while(Boson_States.size() != 0)
    {
      std::list<State>::iterator itBoson_States = Boson_States.begin();
      std::list<State> Positive_Roots;
      Positive_Roots.push_back(*itBoson_States);
      Boson_States.erase(itBoson_States);

      std::list<State>::iterator itPositive_Roots = Positive_Roots.begin();
      for(; itPositive_Roots != Positive_Roots.end(); 
	  ++itPositive_Roots)
			{
				for(itBoson_States = Boson_States.begin(); 
						itBoson_States != Boson_States.end(); )
				{
					if(Gauge_Dot(*itPositive_Roots, *itBoson_States) != 0)
					{
						Positive_Roots.push_back(*itBoson_States);
						itBoson_States = Boson_States.erase(itBoson_States);
					}else
						++itBoson_States;//Close if statement for dot product.
				}//Close for loop on Boson_States.
			}//Close for loop on Positive_Roots.
			//Now Identify the gauge group.
			GaugeGroupIdentifier ID(Positive_Roots, Fermion_Mode_Map());
			FFHS_Model_.Add_Gauge_Group(ID.Get_Group());
		}//Close while loop on Boson_States.
	FFHS_Model_.Sort_Gauge_Groups();
}//Close Build_Gauge_Groups.

void ModelBuilder::Build_SUSY_States()
{
  std::set<AlphaSUSY>::iterator itAlpha_SUSYs = Alpha_SUSYs_.begin();
  for(; itAlpha_SUSYs != Alpha_SUSYs_.end(); ++itAlpha_SUSYs)
    {
      StateBuilder SUSY_StateBuilder(*itAlpha_SUSYs, Fermion_Mode_Map(), 
				       Common_Basis_Alphas(), FFHS_Model().k_ij());
      SUSY_StateBuilder.Build_States();
      std::list<State> SUSY_States = SUSY_StateBuilder.States();
      std::list<State>::iterator itSUSY_States = SUSY_States.begin();
      //Add the built SUSY states to the model builder.
      for(; itSUSY_States != SUSY_States.end(); ++itSUSY_States)
	SUSY_States_.push_back(*itSUSY_States);
    }//Close for loop on SUSY sectors.
  //Copy the SUSY states to the model. Shouldn't take long since there can't be more
  //than four SUSY states.
  FFHS_Model_.Set_SUSY_States(SUSY_States());
}//Close Build_SUSY_States.

void ModelBuilder::Build_MatterStates()
{
  std::list<State> All_MatterStates;
  std::set<AlphaFermion>::iterator itAlpha_Fermions = 
		Alpha_Fermions_.begin();
	for(; itAlpha_Fermions != Alpha_Fermions_.end(); 
			++itAlpha_Fermions)
	{
		StateBuilder MatterStateBuilder(*itAlpha_Fermions, 
				Fermion_Mode_Map(),
				Common_Basis_Alphas(), FFHS_Model().k_ij());
		MatterStateBuilder.Build_States();

		std::list<State> MatterStates = MatterStateBuilder.States();
		std::list<State>::iterator itMatterStates = MatterStates.begin();
		for(; itMatterStates != MatterStates.end(); 
				++itMatterStates)
		{
			if(!Is_SUSY_Partner(*itMatterStates))
			{
				std::vector<GroupRepresentation> Representation_Dimensions;
				for(int a=0; a<static_cast<int>(FFHS_Model().Gauge_Groups().size()); a++)
				{
					GroupRepresentation 
						Rep_Dimension = rFFHS_Model().rGauge_Groups().at(a).
						Compute_Rep_Dimension(*itMatterStates, Fermion_Mode_Map());

					if(Rep_Dimension.Dimension() == 0)
						break;

					Representation_Dimensions.push_back(Rep_Dimension);
				}//Close for loop on gauge groups.
				if(Representation_Dimensions.size() == 
						FFHS_Model().Gauge_Groups().size())
					FFHS_Model_.Add_MatterState(MatterState
							(itMatterStates->Numerator(),
							 itMatterStates->Denominator(),
							 itMatterStates->LM_Size(),
							 Representation_Dimensions)); 
			}//Close if statement on SUSY partner.
		}//Close for loop on matter states from that sector.
	}//Close for loop on Fermion sectors.
	FFHS_Model_.Build_MatterRepresentations();
	if(FFHS_Model().MatterRepresentations().size() != 0)
	{
		FFHS_Model_.Set_Matter_Rep_Sign_Convention();
	}
}//Close Build_MatterStates.

void ModelBuilder::Compute_U1_Factors()
{
  int Total_Rank = (FFHS_Model().BV_Set().at(0).BV().size()/4) + 6;
  //This is the linear relationship between the size of the basis vector 
  //and the total rank.
  Total_Rank -= Find_Rank_Cuts();

  int Total_NA_Rank = 0;
  for(int a=0; a<static_cast<int>(FFHS_Model().Gauge_Groups().size()); a++)
    Total_NA_Rank += FFHS_Model().Gauge_Groups().at(a).Name().Rank();

  FFHS_Model_.Set_U1_Factors(Total_Rank - Total_NA_Rank);
}

int ModelBuilder::Gauge_Dot(const State& State1, const State& State2)
{
  int Dot = 0;
  std::map<int, int>::iterator itMap = Fermion_Mode_Map_.find(State1.LM_Size());
  for(; itMap != Fermion_Mode_Map_.end(); ++itMap)
    Dot += (State1.Numerator().at(itMap->first) * 
	    State2.Numerator().at(itMap->first));
  return Dot;
}//Close Gauge_Dot.

bool ModelBuilder::Is_SUSY_Partner(const State& New_MatterState)
{
	if(SUSY_States().size() == 0)
		return false;

	std::list<State>::iterator itSUSY_States = SUSY_States_.begin();
	for(; itSUSY_States != SUSY_States_.end(); ++itSUSY_States)
	{
		bool SUSY_State = true;
		for(int a=0; a<itSUSY_States->LM_Size(); a++)
		{
			if(New_MatterState.Numerator().at(a) != 
					(itSUSY_States->Numerator()).at(a))
			{
				SUSY_State = false;
				break;
			}
		}
		if(SUSY_State)
			return SUSY_State;
	}//Close for loop on the SUSY states.
	return false;
}//Close Is_SUSY_Partner.

int ModelBuilder::Find_Rank_Cuts()
{
  int Rank_Cuts = 0;
  int LM_Size = FFHS_Model().BV_Set().at(0).LM_Size();
  std::map<int, int>::iterator itFermion_Mode_Map = Fermion_Mode_Map_.begin();
  std::map<int, int>::iterator itMap_End = Fermion_Mode_Map_.end();
  for(; itFermion_Mode_Map != itMap_End; ++itFermion_Mode_Map)
    {
      if((itFermion_Mode_Map->first)<LM_Size && 
	 (itFermion_Mode_Map->second)>=LM_Size)
	Rank_Cuts++;
    }
  return Rank_Cuts/2;
}//Close Find_Rank_Cuts.
