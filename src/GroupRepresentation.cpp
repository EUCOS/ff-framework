#include "GroupRepresentation.h"

GroupRepresentation::GroupRepresentation(int Dimension, char Triality)
{
	Dimension_ = Dimension;
	Triality_ = Triality;
	Is_Complex_ = false;
}//Close constructor.

GroupRepresentation::GroupRepresentation(int Dimension, char Triality, 
		bool Is_Complex)
{
	Dimension_ = Dimension;
	Triality_ = Triality;
	Is_Complex_ = Is_Complex;
}//Close second constructor.

GroupRepresentation::GroupRepresentation(const GroupRepresentation&
		New_GroupRepresentation)
{
	Dimension_ = New_GroupRepresentation.Dimension();
	Triality_ = New_GroupRepresentation.Triality();
	Is_Complex_ = New_GroupRepresentation.Is_Complex();
}//Close copy constructor.

//DEBUG.
void GroupRepresentation::Display() const
{
	std::cout<<Dimension()<<Triality()<<" ";
}//Close Display.
