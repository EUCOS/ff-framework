#include "GaugeGroup.h"

GaugeGroup::GaugeGroup(const std::list<State>& Positive_Roots, 
			 GaugeGroupName Name)
{
  Positive_Roots_ = Positive_Roots;
  Name_ = Name;

  Build_Weyl_Vector();
}//Close constructor.

GaugeGroup::GaugeGroup(const std::list<State>& Positive_Roots,
			 GaugeGroupName Name, 
			 const std::list<State>& Simple_Roots)
{
  Positive_Roots_ = Positive_Roots;
  Name_ = Name;
  Simple_Roots_ = Simple_Roots;

  Build_Weyl_Vector();
}//Close second constructor.

GaugeGroup::GaugeGroup(const GaugeGroup& New_GaugeGroup)
{
  Positive_Roots_ = New_GaugeGroup.Positive_Roots();
  Weyl_Vector_ = New_GaugeGroup.Weyl_Vector();
  Name_ = New_GaugeGroup.Name();
  Simple_Roots_ = New_GaugeGroup.Simple_Roots();
  Dynkin_Labels_ = New_GaugeGroup.Dynkin_Labels();
	Complex_Rep_Dimensions_ = New_GaugeGroup.Complex_Rep_Dimensions();
}//Close copy constructor.

//INTERFACE.
GroupRepresentation GaugeGroup::Compute_Rep_Dimension(
		const State& Weight, 
					     const std::map<int, int>& 
					     Fermion_Mode_Map)
{
	std::vector<int> Weight_Dynkin_Labels = Compute_Dynkin_Labels(Weight, 
			Fermion_Mode_Map);
	std::map<std::vector<int>, GroupRepresentation>
		::iterator itDynkin_Labels = Dynkin_Labels_.find(Weight_Dynkin_Labels);
	int Rep_Dimension = 0;

	if(Not_Highest_Weight(Weight_Dynkin_Labels))
		return GroupRepresentation(Rep_Dimension, ' ');

	//First, see if it's been mapped.
	if(itDynkin_Labels != Dynkin_Labels_.end())
		return (itDynkin_Labels->second);
	else
		Rep_Dimension = Apply_Weyl_Formula(Weight, Fermion_Mode_Map);

	//Now make sure it's a highest weight.
	if(Rep_Dimension == 0)
		return GroupRepresentation(Rep_Dimension, ' ');
	else if((Name().Class() == 'D') || 
			((Name().Class() == 'A')&&(Name().Rank()>1))||
			((Name().Class() == 'E')&&(Name().Rank()==6)))
		//Check if it's a barred rep, and add to Dynkin_Labels.
	{
		if(Is_Barred_Rep(Weight_Dynkin_Labels, Rep_Dimension))
			Rep_Dimension *= -1;
		char Triality = ' ';
		if(Is_D4())
			Triality = Compute_Triality(Weight_Dynkin_Labels);

		bool Is_Complex = Complex_Rep_Dimensions_.find(abs(Rep_Dimension)) !=
			Complex_Rep_Dimensions_.end();
		GroupRepresentation New_Group_Representation(Rep_Dimension, 
				Triality, Is_Complex);
		Dynkin_Labels_[Weight_Dynkin_Labels] = New_Group_Representation;
		return New_Group_Representation;
	}else
	{
		GroupRepresentation New_Group_Representation(Rep_Dimension, ' ');
		Dynkin_Labels_[Weight_Dynkin_Labels] = New_Group_Representation;
		return New_Group_Representation;
	}
}//Close Compute_Rep_Dimension.

//DEBUG.
void GaugeGroup::Display() const
{
  Name().Display();
}//Close Display.

void GaugeGroup::Display_Positive_Roots() const
{
  std::list<State>::const_iterator itPositive_Roots = 
		Positive_Roots_.begin();
  std::cout<<"Positive roots for this gauge group: "<<
		Positive_Roots().size()<<std::endl;
  for(; itPositive_Roots != Positive_Roots_.end(); 
      ++itPositive_Roots)
    {
      itPositive_Roots->Display();
    }
}//Close Display_Positive_Roots.

//PRIVATE.
//HELPERS.
void GaugeGroup::Build_Weyl_Vector()
{
	std::list<State>::iterator itPositive_Roots = Positive_Roots_.begin();
	std::vector<int> 
		Weyl_Vector_Num(itPositive_Roots->Numerator().size(), 0);
	for(; itPositive_Roots != Positive_Roots_.end(); 
			++itPositive_Roots)
	{
		for(int a=0; a<static_cast<int>(Weyl_Vector_Num.size()); a++)
			Weyl_Vector_Num.at(a) += (itPositive_Roots->Numerator().at(a));
	}
	Weyl_Vector_ = 
		State(Weyl_Vector_Num, 2*Positive_Roots().front().Denominator(),
				Positive_Roots().front().LM_Size());
}//Close Build_Weyl_Vector.

int GaugeGroup::Apply_Weyl_Formula(const State& Weight, 
				    const std::map<int, int>& Fermion_Mode_Map)
{
	std::list<State>::const_iterator itPositive_Roots = 
		Positive_Roots_.begin();
	double Full_Dimension = 1;
	long double Dimension_Numerator = 1;
	long double Dimension_Denominator = 1;
	for(; itPositive_Roots != Positive_Roots_.end(); 
			++itPositive_Roots)
	{
		int Weyl_Dot = Gauge_Dot(Weyl_Vector(), *itPositive_Roots, 
				Fermion_Mode_Map);
		int Weight_Dot = Gauge_Dot(Weight, *itPositive_Roots, 
				Fermion_Mode_Map);

		int Next_Dimension_Numerator = Weyl_Dot + (2*Weight_Dot);
		//2 is for the denominator of the Weyl vector.
		int Next_Dimension_Denominator = Weyl_Dot;

		Dimension_Numerator *= double(Next_Dimension_Numerator);
		Dimension_Denominator *= double(Next_Dimension_Denominator);

		Full_Dimension = Dimension_Numerator/Dimension_Denominator;
		//Check to see if it can be reduced.
		if(Full_Dimension - int(Full_Dimension) == 0)
		{
			Dimension_Numerator = Full_Dimension;
			Dimension_Denominator = 1;
		}

		//Break the loop if the dimension becomes zero.
		if(Full_Dimension == 0)
			return int(Full_Dimension);
	}//Close for loop on positive roots.

	if(Full_Dimension < 0)
		return 0;

	return int(Full_Dimension);
}//Close Apply_Weyl_Forumula.

int GaugeGroup::Gauge_Dot(const State& State1, const State& State2,
			   const std::map<int, int>& Fermion_Mode_Map) const
{
  
  int Dot = 0;
  std::map<int, int>::const_iterator itMap = 
		Fermion_Mode_Map.find(State1.LM_Size());
  for(; itMap != Fermion_Mode_Map.end(); ++itMap)
	{
    Dot += (State1.Numerator().at(itMap->first) *
	    State2.Numerator().at(itMap->first));
  
	}
  return Dot; 
}//Close Gauge_Dot.

std::vector<int> GaugeGroup::Compute_Dynkin_Labels(const State& Weight,
						    const std::map<int, int>& Fermion_Mode_Map)
{
	std::list<State>::iterator itSimple_Roots = Simple_Roots_.begin();
	std::list<State>::iterator itEnd_Simple_Roots = Simple_Roots_.end();
	std::vector<int> Weight_Dynkin_Labels;
	for(; itSimple_Roots != itEnd_Simple_Roots; 
			++itSimple_Roots)
	{
		int Weight_Dynkin_Label = Gauge_Dot(Weight, *itSimple_Roots, 
				Fermion_Mode_Map);
		if(Weight_Dynkin_Label<0)
			return Weight_Dynkin_Labels;
		else
			Weight_Dynkin_Labels.push_back(Weight_Dynkin_Label);
	}
	return Weight_Dynkin_Labels;
}//Close Compute_Dynkin_Labels.

bool GaugeGroup::Is_Barred_Rep(const std::vector<int>& 
		Weight_Dynkin_Labels, int Rep_Dimension)
{
	//The D and E6 are a little more difficult.
	if(Name().Class() == 'A')
	{
		std::vector<int> Reversed_Weight_Dynkin_Labels = 
			Weight_Dynkin_Labels;
		std::reverse(Reversed_Weight_Dynkin_Labels.begin(),
				Reversed_Weight_Dynkin_Labels.end());
		if(Weight_Dynkin_Labels != Reversed_Weight_Dynkin_Labels)
			Complex_Rep_Dimensions_.insert(Rep_Dimension);
		if(Weight_Dynkin_Labels < Reversed_Weight_Dynkin_Labels)
			return true;
		else 
			return false;
	}else if (Name().Class() == 'D')
	{
		std::vector<int> Reversed_Weight_Dynkin_Labels = 
			Weight_Dynkin_Labels;
		Reversed_Weight_Dynkin_Labels.at(1) = 
			Weight_Dynkin_Labels.at(2);
		Reversed_Weight_Dynkin_Labels.at(2) = 
			Weight_Dynkin_Labels.at(1);
		if(Weight_Dynkin_Labels != Reversed_Weight_Dynkin_Labels)
			Complex_Rep_Dimensions_.insert(Rep_Dimension);
		return (Weight_Dynkin_Labels < Reversed_Weight_Dynkin_Labels);
	}else//E6.
	{
		std::vector<int> Reversed_Weight_Dynkin_Labels = 
			Weight_Dynkin_Labels;
		Reversed_Weight_Dynkin_Labels.at(2) = 
			Weight_Dynkin_Labels.at(3);
		Reversed_Weight_Dynkin_Labels.at(3) = 
			Weight_Dynkin_Labels.at(2);
		Reversed_Weight_Dynkin_Labels.at(4) = 
			Weight_Dynkin_Labels.at(5);
		Reversed_Weight_Dynkin_Labels.at(5) = 
			Weight_Dynkin_Labels.at(4);
		if(Weight_Dynkin_Labels != Reversed_Weight_Dynkin_Labels)
			Complex_Rep_Dimensions_.insert(Rep_Dimension);
		return Weight_Dynkin_Labels < Reversed_Weight_Dynkin_Labels;
	}
}//Close Is_Barred_Rep.

char GaugeGroup::Compute_Triality(const std::vector<int>& 
		Weight_Dynkin_Labels)
{
	//ONLY WORKS FOR SO(8) GROUPS.
	//NO OTHER GROUP HAS TRIALITY.
	if(Weight_Dynkin_Labels.at(3) != 0)
		return 'v';
	else
		return ' ';
}//Close Compute_Triality.

bool GaugeGroup::Not_Highest_Weight
(const std::vector<int>& Weight_Dynkin_Labels)
{
  int Weight_DL_Size = Weight_Dynkin_Labels.size();
  if(Weight_DL_Size<Name().Rank())
   return true;

  for(int a=0; a<Weight_DL_Size; a++)
    {
      if(Weight_Dynkin_Labels.at(a)<0)
	return true;
    }
  return false;
}//Close Not_Highest_Weight.
