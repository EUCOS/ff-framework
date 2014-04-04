#include "CondensedMatterRepresentation.h"

CondensedMatterRepresentation::CondensedMatterRepresentation
(const std::vector<GroupRepresentation>& Rep_Dimension, int Duplicates,
 const std::vector<GaugeGroupName>& Gauge_Groups):
	MatterRepresentation(Rep_Dimension, Duplicates)
{
	//Reduce the matter representation and initialize Gauge_Group_Indices.
	Collapse_MatterRepresentation(Gauge_Groups);
}//Close constructor.

CondensedMatterRepresentation::CondensedMatterRepresentation
(const CondensedMatterRepresentation& New_CondensedMatterRepresentation)
	:MatterRepresentation(New_CondensedMatterRepresentation)
{
	Gauge_Groups_ = New_CondensedMatterRepresentation.Gauge_Groups();
	Gauge_Group_Indices_ = 
		New_CondensedMatterRepresentation.Gauge_Group_Indices();
}//Close copy constructor.
	
//DEBUG.
void CondensedMatterRepresentation::Display() const
{
	for(int a=0; a<static_cast<int>(Gauge_Groups().size()); ++a)
		Gauge_Groups().at(a).Display();
	std::cout<<std::endl;
	std::cout<<Duplicates()<<": ";
	for(int a=0; a<static_cast<int>(Rep_Dimension().size()); ++a)
		Rep_Dimension().at(a).Display();
	std::cout<<std::endl;
	for(int a=0; a<static_cast<int>(Gauge_Group_Indices().size()); ++a)
		std::cout<<Gauge_Group_Indices().at(a)<<" ";
	std::cout<<std::endl<<std::endl;
}//Close Display.

//HELPERS.
void CondensedMatterRepresentation::Collapse_MatterRepresentation(const
		std::vector<GaugeGroupName>& Gauge_Groups)
{
	std::vector<GroupRepresentation> New_Rep_Dimension;
	for(int a=0; a<static_cast<int>(Rep_Dimension().size()); ++a)
	{
		if(Rep_Dimension().at(a).Dimension() != 1)
		{
			New_Rep_Dimension.push_back(Rep_Dimension().at(a));
			Gauge_Group_Indices_.push_back(a);
			Gauge_Groups_.push_back(Gauge_Groups.at(a));
		}//Close if statment on collapsing the representation.
	}//Close for loop over the Rep_Dimension.
	

	Rep_Dimension_.swap(New_Rep_Dimension);
}//Close Collapse_MatterRepresentation.

