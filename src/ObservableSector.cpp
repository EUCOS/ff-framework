#include "ObservableSector.h"

ObservableSector::ObservableSector(const ObservableSector& 
		New_ObservableSector)
{
	Generations_ = New_ObservableSector.Generations();

	Barred_Generations_ = New_ObservableSector.Barred_Generations();

	Anti_Generations_ = New_ObservableSector.Anti_Generations();

	Anti_Barred_Generations_ = New_ObservableSector.
		Anti_Barred_Generations();

}//Close copy constructor.
int ObservableSector::Count_Chiral_Generations() const
{
	//This only works when all of the vectors have the same size,
	//and that there are a maximum of two distinct reps in a generation.
	//There isn't error checking for this yet.
	//Put it on the to do list.

	int Generation_Count = 0;
	int Barred_Generation_Count = 0;

	//Count the generations.
	if(Generations().size() == 2)
		Generation_Count = std::min(Generations().at(0).Duplicates(),
				Generations().at(1).Duplicates());
	else
		Generation_Count = Generations().at(0).Duplicates();

	//Count the barred generations.
	if(Barred_Generations().size() == 2)
		Barred_Generation_Count = std::min(Barred_Generations().at(0).
				Duplicates(), Barred_Generations().at(1).Duplicates());
	else
		Barred_Generation_Count = Barred_Generations().at(0).Duplicates();

	return abs(Generation_Count - Barred_Generation_Count);
}//Close Count_Chiral_Generations.

int ObservableSector::Count_Chiral_Anti_Generations() const
{
	//This only works when all of the vectors have the same size,
	//and there is a maximum of one distinct rep in a generation.
	//There isn't error checking for this yet.
	//Put it on the to do list.

	if(Anti_Generations().size() == 0)
		return 0;

	int Anti_Generation_Count = Anti_Generations().at(0).Duplicates();
	int Anti_Barred_Generation_Count = Anti_Barred_Generations().at(0).
		Duplicates();

	return abs(Anti_Generation_Count - Anti_Barred_Generation_Count);
}//Close Count_Chiral_Anti_Generations.

void ObservableSector::Display() const
{
	std::cout<<"Generations: "<<std::endl;
	for(int a=0; a<static_cast<int>(Generations().size());++a)
	{
		Generations().at(a).Display();
	}//Close for loop on Generations.
	std::cout<<std::endl;

	std::cout<<"Barred Generations: "<<std::endl;
	for(int a=0; a<static_cast<int>(Barred_Generations().size()); ++a)
	{
		Barred_Generations().at(a).Display();
	}//Close for loop on Barred_Generations.
	std::cout<<std::endl;

	std::cout<<"Anti Generations: "<<std::endl;
	for(int a=0; a<static_cast<int>(Anti_Generations().size()); ++a)
	{
		Anti_Generations().at(a).Display();
	}//Close for loop on Anti_Generations.
	std::cout<<std::endl;

	std::cout<<"Anti Barred Generations: "<<std::endl;
	for(int a=0; a<static_cast<int>(Anti_Barred_Generations().size()); ++a)
	{
		Anti_Barred_Generations().at(a).Display();
	}//Close for loop on Anti_Barred_Generations.
	std::cout<<std::endl;
}//Close Display function.
