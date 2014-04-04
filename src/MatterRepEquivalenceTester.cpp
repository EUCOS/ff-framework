#include "MatterRepEquivalenceTester.h"

MatterRepEquivalenceTester::MatterRepEquivalenceTester(const std::
		vector<CondensedMatterRepresentation>& Matter_Rep_Class_A,
		const std::vector<CondensedMatterRepresentation>& 
		Matter_Rep_Class_B, const std::vector<GaugeGroupName>& Gauge_Groups)
{
	Matter_Rep_Class_A_ = Matter_Rep_Class_A;
	Matter_Rep_Class_B_ = Matter_Rep_Class_B;
	Gauge_Groups_ = Gauge_Groups;

	std::sort(Matter_Rep_Class_A_.begin(), Matter_Rep_Class_A_.end());
	std::sort(Matter_Rep_Class_B_.begin(), Matter_Rep_Class_B_.end());

	Equivalent_Matter_Classes_ = 
		(Matter_Rep_Class_A == Matter_Rep_Class_B);
	//The gauge groups must retain their absolute ordering, as the 
	//gauge group indices are stored in the matter rep classes.
}//Close constructor.

MatterRepEquivalenceTester::MatterRepEquivalenceTester(const 
		MatterRepEquivalenceTester& New_MatterRepEquivalenceTester)
{
	Matter_Rep_Class_A_ = 
		New_MatterRepEquivalenceTester.Matter_Rep_Class_A();
	Matter_Rep_Class_B_ = 
		New_MatterRepEquivalenceTester.Matter_Rep_Class_B();
	Gauge_Groups_ = New_MatterRepEquivalenceTester.Gauge_Groups();
	Equivalent_Matter_Classes_ = 
		New_MatterRepEquivalenceTester.Equivalent_Matter_Classes();
}//Close copy constructor.

void MatterRepEquivalenceTester::Test_Equivalence()
{
	Permute_Complex_Representations(0);
}//CLose Test_Equivalence.


//DEBUG.
void MatterRepEquivalenceTester::Display_Matter_Rep_Class_A() const
{
	std::cout<<"Matter_Rep_Class_A: "<<std::endl;
	for(int a=0; a<static_cast<int>(Matter_Rep_Class_A().size()); ++a)
		Matter_Rep_Class_A().at(a).Display();
}

//PRIVATE.
void MatterRepEquivalenceTester::Permute_Complex_Representations(int 
		Gauge_Group_Index)
{
	if((Gauge_Group_Index < static_cast<int>(Gauge_Groups().size()))&&
			!Equivalent_Matter_Classes())
	{
		if(Gauge_Groups().at(Gauge_Group_Index).Is_D4())
		{
			if(!Gauge_Groups().at(Gauge_Group_Index).V_Ordered() && 
					!Gauge_Groups().at(Gauge_Group_Index).Ordered())//v,s,t can all
				//be reordered.
			{
				Permute_Complex_Representations(Gauge_Group_Index+1);
				Invert_Trialities('s','t',Gauge_Group_Index);
				Permute_Complex_Representations(Gauge_Group_Index+1);
				Invert_Trialities('t','v',Gauge_Group_Index);
				Permute_Complex_Representations(Gauge_Group_Index+1);
				Invert_Trialities('s','v',Gauge_Group_Index);
				Permute_Complex_Representations(Gauge_Group_Index+1);
				Invert_Trialities('s','t',Gauge_Group_Index);
				Permute_Complex_Representations(Gauge_Group_Index+1);
				Invert_Trialities('v','t',Gauge_Group_Index);
				Permute_Complex_Representations(Gauge_Group_Index+1);
			}else if(!Gauge_Groups().at(Gauge_Group_Index).V_Ordered())
				//v and s can be reordered because s > t.
			{
				Permute_Complex_Representations(Gauge_Group_Index+1);
				Invert_Trialities('v','s',Gauge_Group_Index);
				Permute_Complex_Representations(Gauge_Group_Index+1);
			}else if(!Gauge_Groups().at(Gauge_Group_Index).Ordered())
				//s and t can be reordered.
			{
				Permute_Complex_Representations(Gauge_Group_Index+1);
				Invert_Trialities('s','t',Gauge_Group_Index);
				Permute_Complex_Representations(Gauge_Group_Index+1);
			}else
				Permute_Complex_Representations(Gauge_Group_Index+1);
		}else if(!Gauge_Groups().at(Gauge_Group_Index).Ordered())
		{
			Permute_Complex_Representations(Gauge_Group_Index+1);
			Invert_Complex_Reps(Gauge_Group_Index);
			Permute_Complex_Representations(Gauge_Group_Index+1);
		}else
			Permute_Complex_Representations(Gauge_Group_Index+1);
	}else if(!Equivalent_Matter_Classes())
	{
		std::sort(Matter_Rep_Class_A_.begin(), Matter_Rep_Class_A_.end());
		if(Matter_Rep_Class_A() == Matter_Rep_Class_B())
			Equivalent_Matter_Classes_ = true;
	}//Close if/else for ending the recursion.
}//Close Permute_Complex_Representations.

void MatterRepEquivalenceTester::Invert_Complex_Reps
(int Gauge_Group_Index)
{
	for(int a=0; a<static_cast<int>(Matter_Rep_Class_A().size()); ++a)
	{
		std::vector<int> Gauge_Group_Indices = 
			Matter_Rep_Class_A().at(a).Gauge_Group_Indices();
		//Check if the matter representation carried charge under the gauge
		//group and is a complex representation.
		if((std::find(Gauge_Group_Indices.begin(), Gauge_Group_Indices.end(),
					Gauge_Group_Index)!=Gauge_Group_Indices.end()))
		{
			int Rep_Index = int(std::find(Gauge_Group_Indices.begin(), 
						Gauge_Group_Indices.end(), Gauge_Group_Index) - 
					Gauge_Group_Indices.begin());
			if(Matter_Rep_Class_A().at(a).Rep_Dimension().at(Rep_Index).
					Is_Complex())
			{
				std::vector<GroupRepresentation> Rep_Dimension = 
					Matter_Rep_Class_A().at(a).Rep_Dimension();
				Rep_Dimension.at(Rep_Index) = 
					GroupRepresentation(-Rep_Dimension.at(Rep_Index).Dimension(),
							Rep_Dimension.at(Rep_Index).Triality(), 
							Rep_Dimension.at(Rep_Index).Is_Complex());
				Matter_Rep_Class_A_.at(a).Set_Rep_Dimension(Rep_Dimension);
			}
		}//Close if statement on the invertibility of the matter rep.
	}//CLose for loop over Matter_REp_Class_A.
}//Close Invert_Complex_Reps.

void MatterRepEquivalenceTester::Invert_Trialities(char Triality1,
		char Triality2, int Gauge_Group_Index)
{
	if(Triality1 < Triality2)
	{
		char New_Triality1 = Triality2;
		Triality2 = Triality1;
		Triality1 = New_Triality1;
	}
	for(int a=0; a<static_cast<int>(Matter_Rep_Class_A().size()); ++a)
	{
		std::vector<int> Gauge_Group_Indices = 
			Matter_Rep_Class_A().at(a).Gauge_Group_Indices();
		//Check if the matter representation carries charge under the 
		//gauge group and is a complex representation.
		if((std::find(Gauge_Group_Indices.begin(), Gauge_Group_Indices.end(),
						Gauge_Group_Index)!=Gauge_Group_Indices.end()))
		{

			int Rep_Index = int(std::find(Gauge_Group_Indices.begin(),
						Gauge_Group_Indices.end(), Gauge_Group_Index) - 
					Gauge_Group_Indices.begin());
			if(Matter_Rep_Class_A().at(a).Rep_Dimension().
					at(Rep_Index).Is_Complex())
			{
				std::vector<GroupRepresentation> Rep_Dimension = 
					Matter_Rep_Class_A().at(a).Rep_Dimension();
				//Adjust rep dimension based on which trialities are passed to 
				//the function.
				if(Triality1 == 'v')
				{
					if(Triality2 == 's')//Flips v and s.
					{
						if(Rep_Dimension.at(Rep_Index).Triality() == 'v')
							Rep_Dimension.at(Rep_Index).Set_Triality(' ');
						else if(Rep_Dimension.at(Rep_Index).Dimension()>0)
							Rep_Dimension.at(Rep_Index).Set_Triality('v');
					}else if(Triality2 == 't')//Flips v and t.
					{
						if(Rep_Dimension.at(Rep_Index).Triality() == 'v')
						{
							Rep_Dimension.at(Rep_Index).Set_Triality(' ');
							Rep_Dimension.at(Rep_Index).Set_Dimension
								(-Rep_Dimension.at(Rep_Index).Dimension());
						}else if(Rep_Dimension.at(Rep_Index).Dimension()<0)
						{
							Rep_Dimension.at(Rep_Index).Set_Dimension
								(-Rep_Dimension.at(Rep_Index).Dimension());
							Rep_Dimension.at(Rep_Index).Set_Triality('v');
						}
					}//Close if/else on Triality2.
				}else if(Triality1 == 't')//The ordering implies that the only other
					//triality is 's'.
				{
					//Switch 't' and 's'.
					if(Rep_Dimension.at(Rep_Index).Dimension()<0)
						Rep_Dimension.at(Rep_Index).Set_Dimension
							(-Rep_Dimension.at(Rep_Index).Dimension());
					else if(Rep_Dimension.at(Rep_Index).Triality()==' ')
						Rep_Dimension.at(Rep_Index).Set_Dimension
							(-Rep_Dimension.at(Rep_Index).Dimension());
				}//Close if/else on Triality1.
				//Replace existing Matter_Rep_Class element with inverted one.
				Matter_Rep_Class_A_.at(a).Set_Rep_Dimension(Rep_Dimension);
			}
		}//Close if/else on being a complex representation charged under the 
		//gauge group.
	}//Close for loop over Matter_Rep_Class_A.
}//Close Invert_Trialities.
