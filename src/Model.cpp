#include "Model.h"

//PUBLIC.
Model::Model(int Large_ST_Dimensions)
{
  Load_One_Vector(Large_ST_Dimensions);
  NAHE_Loaded_ = false;
}//Close constructor.

Model::Model(const Model& New_Model)
{
  BV_Set_ = New_Model.BV_Set();
  k_ij_ = New_Model.k_ij();
  Gauge_Groups_ = New_Model.Gauge_Groups();
  MatterStates_ = New_Model.MatterStates();
  SUSY_States_ = New_Model.SUSY_States();
  NAHE_Loaded_ = New_Model.NAHE_Loaded();
}//Close copy constructor.

//INTERFACE.
void Model::Load_S_Vector(int Large_ST_Dimensions)
{
  std::vector<int> BV_S;
  //First load ST Fermions.
  for(int a=0; a<(Large_ST_Dimensions-2); a++)
    BV_S.push_back(1);
  //Now load compact triplets.
  for(int a=0; a<(10 - Large_ST_Dimensions); a++)
    {
      BV_S.push_back(1);
      BV_S.push_back(0);
      BV_S.push_back(0);
    }
  //Now load the right movers.
  for(int a=0; a<(52 - 2*Large_ST_Dimensions); a++)
    BV_S.push_back(0);

  Load_BV_Set(BasisVector(BV_S, 2, Large_ST_Dimensions));
}//Close Load_S_Vector.

void Model::Load_BV_Set(const BasisVector& New_BV)
{
  BV_Set_.push_back(New_BV);
  k_ij_.Load_GSOCoefficientMatrix_Order(New_BV.Order());
}//Close Load_BV_Set.

void Model::Load_k_ij_Row(const std::vector<int>& New_k_ij_Row)
{
  k_ij_.Load_GSOCoefficientMatrix_Row(New_k_ij_Row);
}//Close Load_k_ij_Matrix_Row.

void Model::Add_Gauge_Group(const GaugeGroup& New_Gauge_Group)
{
  Gauge_Groups_.push_back(New_Gauge_Group);
}//Close Add_Gauge_Group.

void Model::Sort_Gauge_Groups()
{
	std::sort(Gauge_Groups_.begin(), Gauge_Groups_.end());
}//Close Sort_Gauge_Groups.

void Model::Add_MatterState(const MatterState& New_MatterState)
{
  MatterStates_.push_back(New_MatterState);
}//Close Add_MatterState.

void Model::Build_MatterRepresentations()
{
	if(MatterStates().size() > 0)
	{
		MatterStates_.sort();
		std::reverse(MatterStates_.begin(), MatterStates_.end());
		std::list<MatterState>::iterator itMatterStates = 
			MatterStates_.begin();
		while(itMatterStates != MatterStates_.end())
		{
			int Duplicate_Representations = 
				std::count(itMatterStates, MatterStates_.end(), 
						*itMatterStates);
			MatterRepresentations_.push_back(MatterRepresentation
					(itMatterStates->Representations(), Duplicate_Representations));
			itMatterStates = std::find_if
				(itMatterStates, MatterStates_.end(), 
					bind2nd(std::not_equal_to<MatterState>(), *itMatterStates));
		}//Close while loop on finding equal representations.
	}//Close if statement on there being matter states.
}//Close Build_MatterRepresentations.

void Model::Set_Matter_Rep_Sign_Convention()
{
	for(int a=0;
    a<static_cast<int>(MatterRepresentations().front().Rep_Dimension().size());
    ++a)
	{
		std::vector<GroupRepresentation> Matter_Rep_Column;
		for(int b=0; b<static_cast<int>(MatterRepresentations().size()); ++b)
			Matter_Rep_Column.push_back(MatterRepresentations().at(b).
					Rep_Dimension().at(a));
		
		//Count the unique appearances of the absolute values of the rep 
		//dimensions.
		std::set<int> Complex_Rep_Dimensions = 
			Gauge_Groups().at(a).Complex_Rep_Dimensions();

		std::set<int>::iterator itComplex_Rep_Dimensions = 
			Complex_Rep_Dimensions.begin();
		//Only used for SO(8) Reps.
		bool v_Ordered = false;
		bool s_Ordered = false;
		bool t_Ordered = false;
		for(; itComplex_Rep_Dimensions != 
				Complex_Rep_Dimensions.end(); ++itComplex_Rep_Dimensions)
		{
			if(Gauge_Groups().at(a).Is_D4())
			{

				//Count vector reps.
				GroupRepresentation New_Group_Representation
					(*itComplex_Rep_Dimensions, 'v');
				int Distinct_Vector_Reps = 0;
				std::vector<GroupRepresentation>::iterator itMatter_Rep_Column = 
					Matter_Rep_Column.begin();
				while(itMatter_Rep_Column != Matter_Rep_Column.end())
				{
					itMatter_Rep_Column =std::find(itMatter_Rep_Column, 
							Matter_Rep_Column.end(), New_Group_Representation);
					int Rep_Index = int(itMatter_Rep_Column - 
							Matter_Rep_Column.begin());
					if(itMatter_Rep_Column != Matter_Rep_Column.end())
					{
						Distinct_Vector_Reps += 
							MatterRepresentations().at(Rep_Index).
							Duplicates();
						++itMatter_Rep_Column;
					}
				}
				//Count spinor reps.
				New_Group_Representation.Set_Triality(' ');
				itMatter_Rep_Column = Matter_Rep_Column.begin();
				int Distinct_Unbarred_Reps = 0;
				while(itMatter_Rep_Column != Matter_Rep_Column.end())
				{
					itMatter_Rep_Column =std::find(itMatter_Rep_Column, 
							Matter_Rep_Column.end(), New_Group_Representation);
					int Rep_Index = int(itMatter_Rep_Column - 
							Matter_Rep_Column.begin());
					if(itMatter_Rep_Column != Matter_Rep_Column.end())
					{
						Distinct_Unbarred_Reps += 
							MatterRepresentations().at(Rep_Index).
							Duplicates();
						++itMatter_Rep_Column;
					}
				}
				//COunt conjugate spinor reps.
				New_Group_Representation.Set_Dimension(-*itComplex_Rep_Dimensions);
				int Distinct_Barred_Reps = 0;
				itMatter_Rep_Column = Matter_Rep_Column.begin();
				while(itMatter_Rep_Column != Matter_Rep_Column.end())
				{
					itMatter_Rep_Column =std::find(itMatter_Rep_Column, 
							Matter_Rep_Column.end(), New_Group_Representation);
					int Rep_Index = int(itMatter_Rep_Column - 
							Matter_Rep_Column.begin());
					if(itMatter_Rep_Column != Matter_Rep_Column.end())
					{
						Distinct_Barred_Reps += MatterRepresentations().at(Rep_Index).
							Duplicates();
						++itMatter_Rep_Column;
					}
				}

				if(!v_Ordered && !s_Ordered && !t_Ordered)//Order all three.
				{
					if((Distinct_Vector_Reps != Distinct_Unbarred_Reps)&&
							(Distinct_Unbarred_Reps != Distinct_Barred_Reps))
					{
						v_Ordered = true;
						s_Ordered = true;
						t_Ordered = true;
						if((Distinct_Vector_Reps > Distinct_Unbarred_Reps)&&
								(Distinct_Vector_Reps > Distinct_Barred_Reps))
						{
							if(Distinct_Unbarred_Reps < Distinct_Barred_Reps)
								Invert_Trialities(a, 's', 't', Complex_Rep_Dimensions);

						}else if((Distinct_Vector_Reps < Distinct_Unbarred_Reps)&&
								(Distinct_Vector_Reps > Distinct_Barred_Reps))
						{
							Invert_Trialities(a, 'v', 's', Complex_Rep_Dimensions);	

						}else if((Distinct_Vector_Reps < Distinct_Unbarred_Reps)&&
								(Distinct_Vector_Reps < Distinct_Barred_Reps))
						{
							if(Distinct_Unbarred_Reps < Distinct_Barred_Reps)
								Invert_Trialities(a, 'v', 't', Complex_Rep_Dimensions);
							else if(Distinct_Unbarred_Reps < Distinct_Barred_Reps)
							{
								Invert_Trialities(a, 'v','t',Complex_Rep_Dimensions);
								Invert_Trialities(a, 'v', 's', Complex_Rep_Dimensions);
							}
						}//Close if/else on v != s != t/
					}else if((Distinct_Vector_Reps != Distinct_Unbarred_Reps)&&
							(Distinct_Unbarred_Reps == Distinct_Barred_Reps))
					{
						if(Distinct_Vector_Reps < Distinct_Unbarred_Reps)
						{
							Invert_Trialities(a, 'v', 's', Complex_Rep_Dimensions);
							Invert_Trialities(a, 's', 't', Complex_Rep_Dimensions);
							t_Ordered = true;
						}else if(Distinct_Vector_Reps > Distinct_Unbarred_Reps)
							v_Ordered = true;

					}else if((Distinct_Vector_Reps == Distinct_Unbarred_Reps)&&
							(Distinct_Vector_Reps != Distinct_Barred_Reps))
					{
						t_Ordered = true;
						if(Distinct_Vector_Reps < Distinct_Barred_Reps)
							Invert_Trialities(a, 'v', 't', Complex_Rep_Dimensions);

					}else if((Distinct_Vector_Reps == Distinct_Barred_Reps)&&
							(Distinct_Vector_Reps != Distinct_Unbarred_Reps))
					{
						if(Distinct_Vector_Reps < Distinct_Unbarred_Reps)
						{
							Invert_Trialities(a, 'v', 's', Complex_Rep_Dimensions);
							v_Ordered = true;
						}else if(Distinct_Vector_Reps > Distinct_Unbarred_Reps)
						{
							Invert_Trialities(a, 's', 't', Complex_Rep_Dimensions);
							t_Ordered = true;
						}
					}
				}else if(!v_Ordered && !s_Ordered)//Order s and v only.
				{
					if(Distinct_Vector_Reps < Distinct_Unbarred_Reps)
					{
						Invert_Trialities(a, 'v', 's', Complex_Rep_Dimensions);
						v_Ordered = true;
						s_Ordered = true;
					}else if(Distinct_Vector_Reps > Distinct_Unbarred_Reps)
					{
						v_Ordered = true;
						s_Ordered = true;
					}
				}else if(!s_Ordered && !t_Ordered)//Order s and t only.
				{
					if(Distinct_Vector_Reps < Distinct_Barred_Reps)
					{
						Invert_Trialities(a, 'v', 't', Complex_Rep_Dimensions);
						v_Ordered = true;
						t_Ordered = true;
					}else if(Distinct_Vector_Reps > Distinct_Barred_Reps)
					{
						v_Ordered = true;
						t_Ordered = true;
					}
				}//Close if/else on which elements are ordered.

				if(v_Ordered)
					Gauge_Groups_.at(a).Set_V_Ordered(true);
				if(s_Ordered || t_Ordered)
					Gauge_Groups_.at(a).Set_Ordered(true);

			}else//Group is not SO(8).
			{
				GroupRepresentation 
					New_Group_Representation(*itComplex_Rep_Dimensions, ' ');
				std::vector<GroupRepresentation>::iterator itMatter_Rep_Column = 
					Matter_Rep_Column.begin();
				int Distinct_Unbarred_Reps = 0;
				while(itMatter_Rep_Column != Matter_Rep_Column.end())
				{
					itMatter_Rep_Column =std::find(itMatter_Rep_Column, 
							Matter_Rep_Column.end(), New_Group_Representation);
					int Rep_Index = int(itMatter_Rep_Column - 
							Matter_Rep_Column.begin());
					if(itMatter_Rep_Column != Matter_Rep_Column.end())
					{
						Distinct_Unbarred_Reps += 
							MatterRepresentations().at(Rep_Index).
							Duplicates();
						++itMatter_Rep_Column;
					}
				}
				New_Group_Representation.Set_Dimension(-*itComplex_Rep_Dimensions);
				int Distinct_Barred_Reps = 0;
				itMatter_Rep_Column = Matter_Rep_Column.begin();
				while(itMatter_Rep_Column != Matter_Rep_Column.end())
				{
					itMatter_Rep_Column =std::find(itMatter_Rep_Column, 
							Matter_Rep_Column.end(), New_Group_Representation);
					int Rep_Index = int(itMatter_Rep_Column - 
							Matter_Rep_Column.begin());
					if(itMatter_Rep_Column != Matter_Rep_Column.end())
					{
						Distinct_Barred_Reps += MatterRepresentations().at(Rep_Index).
							Duplicates();
						++itMatter_Rep_Column;
					}
				}
				if(Distinct_Barred_Reps > Distinct_Unbarred_Reps)
				{
					Invert_Complex_Rep_Signs(a, Complex_Rep_Dimensions);
					Gauge_Groups_.at(a).Set_Ordered(true);
					break;
				}//CLose if statement on needing to 
				//invert the complex rep dimensions.
			}//Close if/else on SO(8).
		}//Close for loop over the unique rep dimensions.
	}//Close outer for loop on MatterRepresentations.
}//Close Set_Matter_Rep_Sign_Convention.

//DEBUG.
void Model::Display_BV_Set() const
{
  std::cout<<"BV Set: "<<BV_Set().size()<<std::endl;
  for(int a=0; a<static_cast<int>(BV_Set().size()); a++)
    BV_Set().at(a).Display();

  std::cout<<std::endl;
}//Close Display_BV_Set.

void Model::Display_k_ij() const
{
  k_ij_.Display();
}//Close Display_k_ij.

void Model::Display_Gauge_Groups() const
{
  for(int a=0; a<static_cast<int>(Gauge_Groups().size()); a++)
    Gauge_Groups().at(a).Display();
  std::cout<<std::endl;
}//Close Display_Gauge_Groups.

void Model::Display_MatterRepresentations() const
{
	for(int a=0; a<static_cast<int>(MatterRepresentations().size()); ++a)
		MatterRepresentations().at(a).Display();
    std::cout<<"Total matter representations: "
			<<MatterStates().size()<<std::endl;
}//Close Display_MatterRepresentations.

void Model::Display_Particle_Content() const
{
  Display_Gauge_Groups();
  Display_MatterRepresentations();
  std::cout<<"U(1)'s: "<<U1_Factors()<<std::endl;
  std::cout<<"ST SUSYs: "<<SUSY_States().size()<<std::endl;
}//Close Display_Particle_Content.

//PRIVATE.
//HELPERS.
void Model::Load_One_Vector(int Large_ST_Dimensions)
{
  std::vector<int> BV_One;
  for(int a=0; a<(80 - 4*Large_ST_Dimensions); a++)
    BV_One.push_back(1);

  Load_BV_Set(BasisVector(BV_One, 2, Large_ST_Dimensions));
}//Close Load_One_Vector.

void Model::Invert_Complex_Rep_Signs(int Gauge_Group_Index,
		const std::set<int>& Complex_Rep_Dimensions)
{
	std::set<int>::const_iterator itComplex_Rep_Dimensions = 
		Complex_Rep_Dimensions.begin();
	for(; itComplex_Rep_Dimensions != 
			Complex_Rep_Dimensions.end(); ++itComplex_Rep_Dimensions)
	{
		for(int a=0; a<static_cast<int>(MatterRepresentations().size()); ++a)
		{
			if(abs(MatterRepresentations().at(a).Rep_Dimension().
					at(Gauge_Group_Index).Dimension()) == *itComplex_Rep_Dimensions)
			{
				MatterRepresentations_.at(a).
					Switch_Dimension_Sign(Gauge_Group_Index);
			}//Close if statement on flipping the rep dimensions.
		}//Close for loop on MatterRepresentations().
	}//Close for loop on Complex_Rep_Dimensions.
}//Close Invert_Complex_Rep_Sign.

void Model::Invert_Trialities(int Gauge_Group_Index, char Triality1,
		char Triality2, const std::set<int>& Complex_Rep_Dimensions)
{
	std::set<int>::const_iterator itComplex_Rep_Dimensions = 
		Complex_Rep_Dimensions.begin();
	for(; itComplex_Rep_Dimensions !=
			Complex_Rep_Dimensions.end(); ++itComplex_Rep_Dimensions)
	{
		//Defaults to 's' state. 
		GroupRepresentation Group_Rep1(*itComplex_Rep_Dimensions,' ');
		if(Triality1 == 'v')
			Group_Rep1.Set_Triality('v');
		else if(Triality1 == 't')
			Group_Rep1.Set_Dimension(-*itComplex_Rep_Dimensions);

		GroupRepresentation Group_Rep2(*itComplex_Rep_Dimensions,' ');
		if(Triality2 == 'v')
			Group_Rep2.Set_Triality('v');
		else if(Triality2 == 't')
			Group_Rep2.Set_Dimension(-*itComplex_Rep_Dimensions);
		for(int a=0; a<static_cast<int>(MatterRepresentations().size()); ++a)
		{
			if(MatterRepresentations().at(a).
					Rep_Dimension().at(Gauge_Group_Index) == 
					Group_Rep1)
				MatterRepresentations_.at(a).
					Set_Rep_Dimension_at_Index(Group_Rep2, Gauge_Group_Index);
			else if(MatterRepresentations().at(a).
					Rep_Dimension().at(Gauge_Group_Index) ==
					Group_Rep2)
				MatterRepresentations_.at(a).Set_Rep_Dimension_at_Index
					(Group_Rep1, Gauge_Group_Index);
		}//Close for loop over MatterRepresentations.
	}//Close for loop over the different complex rep dimension states.
}//Close Invert_Trialities.
