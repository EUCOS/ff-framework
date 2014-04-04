#include "GaugeGroupName.h"

//PUBLIC.
GaugeGroupName::GaugeGroupName()
{
  Class_ = 'N';
  Rank_ = 0;
  KM_Level_ = 0;
	Ordered_ = false;
	V_Ordered_ = false;
}//Close constructor.

GaugeGroupName::GaugeGroupName(char Class, int Rank, int KM_Level)
{
  Class_ = Class;
  Rank_ = Rank;
  KM_Level_ = KM_Level;
	Ordered_ = false;
	V_Ordered_ = false;
}//Close second constructor.

GaugeGroupName::GaugeGroupName(const GaugeGroupName& New_GaugeGroupName)
{
  Class_ = New_GaugeGroupName.Class();
  Rank_ = New_GaugeGroupName.Rank();
  KM_Level_ = New_GaugeGroupName.KM_Level();
	Ordered_ = New_GaugeGroupName.Ordered();
	V_Ordered_ = New_GaugeGroupName.V_Ordered();
}//Close copy constructor.

void GaugeGroupName::Display() const
{
  std::cout<<Class()<<" "<<Rank()<<" "<<KM_Level()<<" | "<<std::flush;
  //std::cout<<std::endl;
}//Close Display.
