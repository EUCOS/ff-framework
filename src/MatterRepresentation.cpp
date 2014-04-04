#include "MatterRepresentation.h"

MatterRepresentation::MatterRepresentation
	(const std::vector<GroupRepresentation>& 
		Rep_Dimension, int Duplicates)
{
	Rep_Dimension_ = Rep_Dimension;
	Duplicates_ = Duplicates;
}//Close constructor.

MatterRepresentation::MatterRepresentation(const MatterRepresentation&
		New_MatterRepresentation)
{
	Rep_Dimension_ = New_MatterRepresentation.Rep_Dimension();
	Duplicates_ = New_MatterRepresentation.Duplicates();
}//Close copy constructor.

//INTERFACE.
void MatterRepresentation::Switch_Dimension_Sign(int Dimension_Index)
{
	Rep_Dimension_.at(Dimension_Index).Set_Dimension
		(-Rep_Dimension().at(Dimension_Index).Dimension());
}//Close Switch_Dimension_Sign.

void MatterRepresentation::Set_Rep_Dimension_at_Index
(const GroupRepresentation& New_Group_Representation, int Index)
{
	Rep_Dimension_.at(Index) = New_Group_Representation;
}//Close Set_Rep_Dimension_at_Index.

void MatterRepresentation::Display() const 
{
	std::cout<<Duplicates()<<": ";
	for(int a=0; a<static_cast<int>(Rep_Dimension().size()); ++a)
	{
		std::cout<<Rep_Dimension().at(a).Dimension();
		std::cout<<Rep_Dimension().at(a).Triality()<<" ";
	}//Close for loop over Rep_Dimension.
	std::cout<<std::endl;
}//Close Display.
