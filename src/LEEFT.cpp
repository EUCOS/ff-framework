#include "LEEFT.h"

LEEFT::LEEFT(std::vector<std::string> LEEFT_Data)
{
	Total_NA_Reps_ = 0;
	std::stringstream ss;
	ss<<LEEFT_Data.at(0);

	while(!ss.eof())
	{
		char Class = 'N';
		int Rank = 0;
		int KM_Level = 0;
		ss>>Class;
		ss>>Rank;
		ss>>KM_Level;

		if(KM_Level == 0)//Prevents overreading the whitespace.
			break;

		Gauge_Groups_.push_back(GaugeGroupName(Class, Rank, KM_Level));
	}//Close while loop.

	for(int a=1; a<static_cast<int>(LEEFT_Data.size())-2; a++)
	{
		ss.clear();
		ss<<LEEFT_Data.at(a);

		int Representation_QTY = 0;
		char junk; //for the colon.
		ss>>Representation_QTY;

		std::vector<GroupRepresentation> Representation_Loader;
		if(ss.peek() == ':')
			ss>>junk;
		while(!ss.eof())
		{
			std::string Representation;		
			ss>>Representation;

			if(Representation.size()==0)//Prevents overreading the whitespace.
				break;

			int Dimension=0;
			char Triality = ' ';
			if(Representation.at(Representation.size()-1)=='v')
			{
				Triality = Representation.at(Representation.size()-1);
				Representation = Representation.substr(0,Representation.size()-1);
				Dimension = atoi(Representation.c_str());
			}else
				Dimension = atoi(Representation.c_str());

			GroupRepresentation New_Group_Representation(Dimension, Triality);
			Representation_Loader.push_back(New_Group_Representation);
		}//Close while loop.
		MatterRepresentations_.push_back
			(MatterRepresentation(Representation_Loader, 
														 Representation_QTY));
		Total_NA_Reps_ += Representation_QTY;
		Representation_Loader.clear();
	}//Close for loop.
	std::string U1; // Another label in the file.
	ss.clear();
	ss<<LEEFT_Data.at(LEEFT_Data.size()-2);
	ss>>U1;
	ss>>U1_Factors_;

	std::string ST;
	std::string SUSY;
	ss.clear();
	ss<<LEEFT_Data.back();
	ss>>ST;
	ss>>SUSY;
	ss>>ST_SUSYs_;

	Build_Matter_Rep_Classes();
}//Close constructor.

LEEFT::LEEFT(const Model& FFHS_Model)
{
	Total_NA_Reps_ = 0;
  //Get the gauge group names.
  for(int a=0; a<static_cast<int>(FFHS_Model.Gauge_Groups().size()); ++a)
    Gauge_Groups_.push_back(FFHS_Model.Gauge_Groups().at(a).Name());


  //Get the matter representations.
	MatterRepresentations_ = FFHS_Model.MatterRepresentations();

  //Get the number of U(1) factors.
  U1_Factors_ = FFHS_Model.U1_Factors();

  //Get the number of ST SUSYs.
  ST_SUSYs_ = FFHS_Model.SUSY_States().size();

	//Calculate the total non-Abelian matter reps.
	for(int a=0; a<static_cast<int>(MatterRepresentations().size()); ++a)
		Total_NA_Reps_+=MatterRepresentations().at(a).Duplicates();

	//Organize the matter rep classes.
	Build_Matter_Rep_Classes();
}//Close second constructor.

LEEFT::LEEFT(const LEEFT& New_LEEFT)
{
	Gauge_Groups_ = New_LEEFT.Gauge_Groups();
	MatterRepresentations_ = New_LEEFT.MatterRepresentations();
	U1_Factors_ = New_LEEFT.U1_Factors();
	ST_SUSYs_ = New_LEEFT.ST_SUSYs();
	Total_NA_Reps_ = New_LEEFT.Total_NA_Reps();
	Matter_Rep_Classes_ = New_LEEFT.Matter_Rep_Classes();
	ObservableSectors_ = New_LEEFT.ObservableSectors();
}//Close copy constructor.

//INTERFACE.
	bool LEEFT::Equal_Matter_Rep_Classes
(const std::vector<CondensedMatterRepresentation>& Matter_Rep_Classes1,
 const std::vector<CondensedMatterRepresentation>& Matter_Rep_Classes2) 
const
{
	MatterRepEquivalenceTester New_Tester(Matter_Rep_Classes1,
			Matter_Rep_Classes2, Gauge_Groups());
	New_Tester.Test_Equivalence();
	return New_Tester.Equivalent_Matter_Classes();
}//Close Equal_Matter_Rep_Classes.

//DEBUG.
void LEEFT::Display() const
{
  std::cout<<"Gauge groups and matter representations: "<<std::endl;
  for(int a=0; a<static_cast<int>(Gauge_Groups().size()); a++)
    std::cout<<Gauge_Groups().at(a).Class()<<" "
	     <<Gauge_Groups().at(a).Rank()<<" "
	     <<Gauge_Groups().at(a).KM_Level()<<" ";
  std::cout<<std::endl;
	for(int a=0; a<static_cast<int>(MatterRepresentations().size()); a++)
		MatterRepresentations().at(a).Display();

	std::cout<<"ObservableSectors: "<<std::endl;
	for(int a=0; a<static_cast<int>(ObservableSectors().size()); ++a)
		ObservableSectors().at(a).Display();

  std::cout<<"ST SUSY: "<<ST_SUSYs()<<std::endl;
  std::cout<<"U(1)'s: "<<U1_Factors()<<std::endl;
  std::cout<<std::endl;
}//Close Display.

void LEEFT::Display_Matter_Rep_Classes() const
{
	std::cout<<"Matter_Rep_Classes: "<<Matter_Rep_Classes().size();
	std::cout<<std::endl;
	for(int a=0; a<static_cast<int>(Matter_Rep_Classes().size()); ++a)
		Matter_Rep_Classes().at(a).Display();
	std::cout<<std::endl;
}//Close DIsplay_MatterRepresentation_Classes.

void LEEFT::Build_Matter_Rep_Classes()
{
	for(int a=0; a<static_cast<int>(MatterRepresentations().size()); ++a)
		Matter_Rep_Classes_.push_back(CondensedMatterRepresentation
				(MatterRepresentations().at(a).Rep_Dimension(), 
				 MatterRepresentations().at(a).Duplicates(),
				 Gauge_Groups()));
	std::sort(Matter_Rep_Classes_.begin(), Matter_Rep_Classes_.end());

}//Close Build_Matter_Rep_Classes.
