#include <cassert>
#include "GaugeGroupIdentifier.h"

//PUBLIC.
GaugeGroupIdentifier::GaugeGroupIdentifier(const std::list<State>& 
					       Positive_Roots, 
					       const std::map<int, int>& 
					       Fermion_Mode_Map)
{
  Positive_Roots_ = Positive_Roots;
  Fermion_Mode_Map_ = Fermion_Mode_Map;
  Long_Root_Length_Num_ = 0;
}//Close constructor.

GaugeGroupIdentifier::GaugeGroupIdentifier(const GaugeGroupIdentifier&
					       New_GaugeGroupIdentifier)
{
  Positive_Roots_ = New_GaugeGroupIdentifier.Positive_Roots();
  Simple_Roots_ = New_GaugeGroupIdentifier.Simple_Roots();
  Long_Root_Length_Num_ = New_GaugeGroupIdentifier.Long_Root_Length_Num();
}//Close copy constructor.

//INTERFACE.
GaugeGroup GaugeGroupIdentifier::Get_Group()
{
  
  Identify_KM_Level();
  char Class = Identify_Class();

  switch(Class)
    {
    case 'A':
      return Build_A_Class_Group();

      break;
    case 'B':
      return Build_B_Class_Group();

      break;
    case 'C':
      return Build_C_Class_Group();

      break;
    case 'D':
      return Build_D_Class_Group();

      break;
    case 'E':
      return Build_E_Class_Group();

      break;
    case 'F':
      return Build_F_Class_Group();

      break;
    case 'G':
      return Build_G_Class_Group();

      break;
    };
  std::cerr << "*** Error: gauge class unrecognized." << std::endl;
  assert(false);
}//Close Get_Group.

//DEBUG.
void GaugeGroupIdentifier::Display_Positive_Roots() const
{
  std::cout<<"Nonzero positive roots: "<<std::endl;
  std::list<State>::const_iterator itPositive_Roots = Positive_Roots_.begin();
  for(; itPositive_Roots != Positive_Roots_.end(); 
      ++itPositive_Roots)
    itPositive_Roots->Display();
  std::cout<<std::endl;
}//Close Display_Positive_Roots.

void GaugeGroupIdentifier::Display_Simple_Roots() const
{
  std::cout<<"Simple roots: "<<std::endl;
  std::list<State>::const_iterator itSimple_Roots = Simple_Roots_.begin();
  for(; itSimple_Roots != Simple_Roots_.end(); ++itSimple_Roots)
    itSimple_Roots->Display();
  std::cout<<std::endl;
}//Close Display_Simple_Roots.

//PRIVATE.
//HELPERS.

void GaugeGroupIdentifier::Identify_KM_Level()
{
	std::list<State>::iterator itPositive_Roots = Positive_Roots_.begin();
	int Max_Length_Squared_Num = 0;
	int Max_Length_Squared_Denom = 
		itPositive_Roots->Length_Squared_Denominator();
	for(; itPositive_Roots != Positive_Roots_.end(); 
			++itPositive_Roots)
	{
		int Next_Length_Squared_Num = 
			itPositive_Roots->Length_Squared_Numerator();
		if(Next_Length_Squared_Num == 2*Max_Length_Squared_Denom)
		{
			Long_Root_Length_Num_ = Next_Length_Squared_Num;
			break;
		}else if(Next_Length_Squared_Num > Max_Length_Squared_Num)
			Max_Length_Squared_Num = Next_Length_Squared_Num;
	}//Close for loop on itPositive_Roots.
	if(Long_Root_Length_Num() == 0)
		Long_Root_Length_Num_ = Max_Length_Squared_Num;
}//Close Identify_KM_Level.

char GaugeGroupIdentifier::Identify_Class()
{
	std::list<State>::iterator itPositive_Roots = Positive_Roots_.begin();
	for(; itPositive_Roots != Positive_Roots_.end(); 
			++itPositive_Roots)
	{
		if((itPositive_Roots->Length_Squared_Numerator()) != 
				Long_Root_Length_Num())
			return Identify_BCGF();
	}//Close for loop for checking ADE or BCGF.
	return Identify_ADE();
}//Close Identify_Class.

char GaugeGroupIdentifier::Identify_ADE()
{
  char Class = Identify_ADE_Degeneracies();
  if(Class == 'N')
    {
      double A_Rank = A_Class_Rank();

      if(A_Rank - int(A_Rank) == 0)
	return 'A';
      else
	return 'D';	
    }//Close if statement.
  else 
    return Class;
}//Close Identify_ADE.

char GaugeGroupIdentifier::Identify_BCGF()
{
  char Class = Identify_BCGF_Degeneracies();
  if(Class == 'N')
    {
      int Short_Root_Count = Count_Short_Roots();
      int Long_Root_Count = Positive_Roots().size() - Short_Root_Count;
      if(Short_Root_Count <= Long_Root_Count)//Only N=4. B2 = C2.
	return 'B';
      else
	return 'C';
    }//Close if statement.
  std::cerr << "*** Error: gauge class not recognized" << std::endl;
  assert(false);
}//Close Identify_BCGF.

char GaugeGroupIdentifier::Identify_ADE_Degeneracies()
{
  //First check the degeneracies.
  int NZPR = Positive_Roots().size();
  switch(NZPR)
    {
    case 36:
      Find_Simple_Roots();
      if(Simple_Roots().size() == 6)
	return 'E';
      else
	return 'A';
      break;

    case 63:
      return 'E';
      break;

    case 120:
      Find_Simple_Roots();
      if(Simple_Roots().size() == 8)
	return 'E';
      else
	return 'A';
      break;

    case 210:
      Find_Simple_Roots();
      if(Simple_Roots().size() == 15)
	return 'D';
      else
	return 'A';
      break;

    default:
      return 'N';
    };
}//Close Identify_ADE.

char GaugeGroupIdentifier::Identify_BCGF_Degeneracies()
{
  int NZPR = Positive_Roots().size();
  switch(NZPR)
    {
    case 6:
      return 'G';
      break;

    case 24:
      return 'F';
      break;

    default:
      return 'N';
      break;
    };
}//Close Identify_BCGF_Degeneracies.

void GaugeGroupIdentifier::Find_Simple_Roots()
{
	Positive_Roots_.sort();
	std::list<State>::iterator itPositive_Roots = Positive_Roots_.begin();
	for(; itPositive_Roots != Positive_Roots_.end();
			++itPositive_Roots)
	{
		//Build the initial list.
		std::list<State>::iterator itPositive_Roots2 = 
			Positive_Roots_.begin();

		std::list<State> New_Simple_Roots;
		New_Simple_Roots.push_back(*itPositive_Roots);
		for(; itPositive_Roots2 != Positive_Roots_.end();
				++itPositive_Roots2)
		{
			if(Gauge_Dot(*itPositive_Roots, *itPositive_Roots2) <= 0)
				New_Simple_Roots.push_back(*itPositive_Roots2);
		}//Close for loop on Positive_Roots.

		//Trim the list to form a set of mutually 
		//negative or zero dot products.
		std::list<State>::iterator itSimple_Roots = New_Simple_Roots.begin();
		++itSimple_Roots;
		for(; itSimple_Roots != New_Simple_Roots.end(); 
				++itSimple_Roots)
		{
			std::list<State>::iterator itSimple_Roots2 = 
				New_Simple_Roots.begin();
			for(; itSimple_Roots2 != New_Simple_Roots.end(); )
			{
				int Dot_Product = Gauge_Dot(*itSimple_Roots, *itSimple_Roots2);
				if(Dot_Product > 0 && (*itSimple_Roots != *itSimple_Roots2))
					itSimple_Roots2 = New_Simple_Roots.erase(itSimple_Roots2);
				else
					++itSimple_Roots2;

			}//Close inner for loop on Simple_Roots.
		}//Close for loop on Simple_Roots.

		if(Simple_Roots().size() < New_Simple_Roots.size() && 
				Connected_Simple_Roots(New_Simple_Roots))
			Simple_Roots_.swap(New_Simple_Roots);
	}//Close outer for loop on all positive roots.
}//Close Find_Simple_Roots.

void GaugeGroupIdentifier::Find_Simple_Roots(char Gauge_Group_Class, int Rank)
{
  //If the gauge group has alread been identified (in the case of non simply laced
  //groups) we need to make sure the proper number of short roots make it into 
  //the list of possible simple roots.
  int Short_Root_Count = 0;

  if((Gauge_Group_Class == 'B') || (Gauge_Group_Class == 'G'))
    Short_Root_Count = 1;
  else if(Gauge_Group_Class == 'F')
    Short_Root_Count = 2;
  else if(Gauge_Group_Class == 'C')
    Short_Root_Count = Rank - 1;

  Positive_Roots_.sort();

  std::list<State>::iterator itPositive_Roots = Positive_Roots_.begin();
  for(; itPositive_Roots != Positive_Roots_.end();
      ++itPositive_Roots)
    {
      //Build the initial list.
			std::list<State>::iterator itPositive_Roots2 = 
				Positive_Roots_.begin();
			std::list<State> New_Simple_Roots;
			New_Simple_Roots.push_back(*itPositive_Roots);
			int Short_Roots_Added = 0;
			if(itPositive_Roots->Length_Squared_Numerator() != 
					Long_Root_Length_Num())
				Short_Roots_Added++;
			for(; itPositive_Roots2 != Positive_Roots_.end();
					++itPositive_Roots2)
			{
				if((Short_Root_Count == 0) &&//Simply laced gauge group.
						(Gauge_Dot(*itPositive_Roots, *itPositive_Roots2) <= 0))
					New_Simple_Roots.push_back(*itPositive_Roots2);
				else
				{
					if(itPositive_Roots2->Length_Squared_Numerator() 
							== Long_Root_Length_Num())
					{
						if(Gauge_Dot(*itPositive_Roots, *itPositive_Roots2) <= 0)
							New_Simple_Roots.push_back(*itPositive_Roots2);
					}else if(Short_Roots_Added < Short_Root_Count)
					{
						if(Gauge_Dot(*itPositive_Roots, *itPositive_Roots2) <= 0)
						{
							New_Simple_Roots.push_back(*itPositive_Roots2);
							Short_Roots_Added++;
						}//Close if statement on gauge_Dot
					}//Close if/else on the corrector number of short roots.
				}//Close else on non simply laced gauge groups.
			}//Close for loop on Positive_Roots2.

      //Trim the list to form a set of mutually negative or zero dot products.
      std::list<State>::iterator itSimple_Roots = New_Simple_Roots.begin();
      ++itSimple_Roots;
      for(; itSimple_Roots != New_Simple_Roots.end(); ++itSimple_Roots)
	{
	  std::list<State>::iterator itSimple_Roots2 = New_Simple_Roots.begin();
	  for(; itSimple_Roots2 != New_Simple_Roots.end(); )
	    {
	      int Dot_Product = Gauge_Dot(*itSimple_Roots, *itSimple_Roots2);
	      if(Dot_Product > 0 && (*itSimple_Roots != *itSimple_Roots2))
		itSimple_Roots2 = New_Simple_Roots.erase(itSimple_Roots2);
	      else
		++itSimple_Roots2;

	    }//Close inner for loop on Simple_Roots.
	}//Close for loop on Simple_Roots.

      if(Simple_Roots().size() < New_Simple_Roots.size() && 
	 Connected_Simple_Roots(New_Simple_Roots) && 
	 (Short_Roots_Added == Short_Root_Count))
	Simple_Roots_.swap(New_Simple_Roots);
    }//Close outer for loop on all positive roots.
}//Close Find_Simple_Roots.

bool GaugeGroupIdentifier::Connected_Simple_Roots(const 
		std::list<State>& Simple_Roots)
{
	std::list<State> Unconnected_Simple_Roots = Simple_Roots;
	std::list<State> Connected_Simple_Roots;

	std::list<State>::iterator itUnconnected = 
		Unconnected_Simple_Roots.begin();
	Connected_Simple_Roots.push_back(*itUnconnected);
	Unconnected_Simple_Roots.erase(itUnconnected);

	std::list<State>::iterator itConnected = Connected_Simple_Roots.begin();
	for(; itConnected != Connected_Simple_Roots.end(); 
			++itConnected)
	{
		itUnconnected = Unconnected_Simple_Roots.begin();
		for(; itUnconnected != Unconnected_Simple_Roots.end();)
		{
			if(Gauge_Dot(*itConnected, *itUnconnected) != 0)
			{
				Connected_Simple_Roots.push_back(*itUnconnected);
				itUnconnected = Unconnected_Simple_Roots.erase(itUnconnected);
			}else
				++itUnconnected;
		}//Close for loop on unconnected simple roots.
	}//Close for loop on connected simple roots.
	if(Connected_Simple_Roots.size() == Simple_Roots.size())
		return true;
	else
		return false;
}//Close Connected_Simple_Roots.

void GaugeGroupIdentifier::Order_Simple_Roots(char Class)
{
	//Orders the simple roots of an A,D, or E6 gauge group.
	//This ordering convention is needed to distinguish the Dynkin
	//labels of the complex representations.

	//First, determine which roots to find.
	std::list<State>::iterator itSimple_Roots = Simple_Roots_.begin();
	State First_Simple_Root;
	std::vector<State> Special_Simple_Roots;//For the end roots of D class
	//and the "top" root of E 6. Not used for A class groups.
	if(Class == 'A')
	{
		for(; itSimple_Roots != Simple_Roots_.end(); 
				++itSimple_Roots)
		{
			int Simple_Connection_Count = 0;
			std::list<State>::iterator itSimple_Roots2 = Simple_Roots_.begin();
			for(; itSimple_Roots2 != Simple_Roots_.end();
					++itSimple_Roots2)
			{
				if((Gauge_Dot(*itSimple_Roots, *itSimple_Roots2) != 0) &&
						(*itSimple_Roots != *itSimple_Roots2))
					++Simple_Connection_Count;
			}//Close inner for loop on the simple roots.
			if(Simple_Connection_Count == 1)
			{
				First_Simple_Root = *itSimple_Roots;
				itSimple_Roots = Simple_Roots_.erase(itSimple_Roots);
				break;
			}//Close if statement on simple connections.
		}//Close outer for loop on the simple roots.
	}else 
	{
		for(; itSimple_Roots != Simple_Roots_.end();
				++itSimple_Roots)
		{
			int Simple_Connection_Count = 0;
			std::list<State>::iterator itSimple_Roots2 = Simple_Roots_.begin();
			for(; itSimple_Roots2 != Simple_Roots_.end();
					++itSimple_Roots2)
				{
					if((Gauge_Dot(*itSimple_Roots, *itSimple_Roots2) != 0) &&
							(*itSimple_Roots != *itSimple_Roots2))
					{
						Special_Simple_Roots.push_back(*itSimple_Roots2);
						++Simple_Connection_Count;
					}
				}//Close inner for loop on simple connections.
				if(Simple_Connection_Count == 3)//Found the "junction" root.
				{
					First_Simple_Root = *itSimple_Roots;
					itSimple_Roots = Simple_Roots_.erase(itSimple_Roots);
					break;
				}else//Close if statement on the simple connections.
					Special_Simple_Roots.clear();
		}//Close outer for loop on the simple roots.
		//Now order the roots which produce complex reps.
	}//Close if/else on the group class.

	//Now build a new set of simple roots using the "anchor" root.
	std::list<State> Ordered_Simple_Roots;
	Ordered_Simple_Roots.push_back(First_Simple_Root);
	//Add the special simple roots next.
	for(int a=0; a<static_cast<int>(Special_Simple_Roots.size()); ++a)
	{
		itSimple_Roots = Simple_Roots_.begin();
		bool Is_Connected = false;
		for(; itSimple_Roots != Simple_Roots_.end();
				++itSimple_Roots)
		{
			if((*itSimple_Roots != Special_Simple_Roots.at(a)) && 
					(Gauge_Dot(*itSimple_Roots, Special_Simple_Roots.at(a))!=0))
			{
				Is_Connected = true;
				break;
			}
		}//Close inner for loop over Simple_Roots.
		if(!Is_Connected)
		{
			itSimple_Roots = std::find(Simple_Roots_.begin(), 
					Simple_Roots_.end(), Special_Simple_Roots.at(a));
			Ordered_Simple_Roots.push_back(*itSimple_Roots);
			itSimple_Roots = Simple_Roots_.erase(itSimple_Roots);
		}//Close if statement on connectedness.
	}//Close for loop over Special_Simple_Roots.

	std::list<State>::iterator itOrdered_Simple_Roots = 
		Ordered_Simple_Roots.begin();
	for(; itOrdered_Simple_Roots != 
			Ordered_Simple_Roots.end(); ++itOrdered_Simple_Roots)
	{
		for(itSimple_Roots = Simple_Roots_.begin();
				itSimple_Roots != Simple_Roots_.end(); )
		{
			if(Gauge_Dot(*itOrdered_Simple_Roots, *itSimple_Roots) != 0)
			{
				Ordered_Simple_Roots.push_back(*itSimple_Roots);
				itSimple_Roots = Simple_Roots_.erase(itSimple_Roots);
			}else
				++itSimple_Roots;
		}//Close for loop on the inner simple roots.
	}//Close for loop over the ordered simple roots.
	Simple_Roots_.swap(Ordered_Simple_Roots);
}//Close Order_Simple_Roots.

int GaugeGroupIdentifier::Gauge_Dot(const State& State1, const State& State2)
{
  int Dot = 0;
  std::map<int, int>::iterator itMap = Fermion_Mode_Map_.find(State1.LM_Size());
  for(; itMap != Fermion_Mode_Map_.end(); ++itMap)
    Dot += (State1.Numerator().at(itMap->first) * 
	    State2.Numerator().at(itMap->first));
  return Dot;
}//Close Gauge_Dot.

double GaugeGroupIdentifier::A_Class_Rank()
{
  return (-1 + sqrt(1+8*Positive_Roots().size()))/2;
}//A_Class_Rank.

int GaugeGroupIdentifier::Count_Short_Roots()
{
  int Short_Root_Count = 0;
  std::list<State>::iterator itPositive_Roots = Positive_Roots_.begin();
  for(; itPositive_Roots != Positive_Roots_.end(); 
      ++itPositive_Roots)
    {
      if((itPositive_Roots->Length_Squared_Numerator()) != Long_Root_Length_Num())
	Short_Root_Count++;
    }
  return Short_Root_Count;
}//Count_Long_Roots.

GaugeGroup GaugeGroupIdentifier::Build_A_Class_Group()
{
  int KM_Level = (2*Positive_Roots().front().Denominator()*
		  Positive_Roots().front().Denominator())/
    Long_Root_Length_Num();

  if(Simple_Roots().size() != 0)//Simple roots have been found.
    {
      GaugeGroupName Name('A', Simple_Roots().size(), KM_Level);
			if(Name.Rank() > 1)
				Order_Simple_Roots('A');
      return GaugeGroup(Positive_Roots(), Name, Simple_Roots());
    }else//Simple roots have not been found.
    {
      int Rank = int(A_Class_Rank());
      GaugeGroupName Name('A', Rank, KM_Level);
      Find_Simple_Roots('A', Rank);
			if(Rank > 1)
				Order_Simple_Roots('A');
      return GaugeGroup(Positive_Roots(), Name, Simple_Roots());
    }
}//Close Build_A_Class_Group.

GaugeGroup GaugeGroupIdentifier::Build_B_Class_Group()
{
  int KM_Level = (2*Positive_Roots().front().Denominator()*
		  Positive_Roots().front().Denominator())/
    Long_Root_Length_Num();
  GaugeGroupName Name('B', int(sqrt(Positive_Roots().size())), KM_Level);
  if(Simple_Roots().size() == 0)
    Find_Simple_Roots(Name.Class(), Name.Rank());
  return GaugeGroup(Positive_Roots(), Name, Simple_Roots());
}//Close Build_B_Class_Group.

GaugeGroup GaugeGroupIdentifier::Build_C_Class_Group()
{
  int KM_Level = (2*Positive_Roots().front().Denominator()*
		  Positive_Roots().front().Denominator())/
    Long_Root_Length_Num();
  GaugeGroupName Name('C', int(sqrt(Positive_Roots().size())), KM_Level);
  if(Simple_Roots().size() == 0)
    Find_Simple_Roots(Name.Class(), Name.Rank());
  return GaugeGroup(Positive_Roots(), Name, Simple_Roots());
}//Close Build_C_Class_Group.

GaugeGroup GaugeGroupIdentifier::Build_D_Class_Group()
{
  int KM_Level = (2*Positive_Roots().front().Denominator()*
		  Positive_Roots().front().Denominator())/
    Long_Root_Length_Num();
  int Rank = int((1+sqrt(1 + 4*Positive_Roots().size()))/2);
  GaugeGroupName Name('D', Rank, KM_Level);
  if(Simple_Roots().size() == 0) 
    Find_Simple_Roots('D', Rank);
	Order_Simple_Roots('D');
  return GaugeGroup(Positive_Roots(), Name, Simple_Roots());
}//Close Build_D_Class_Group.

GaugeGroup GaugeGroupIdentifier::Build_E_Class_Group()
{
  int KM_Level = (2*Positive_Roots().front().Denominator()*
		  Positive_Roots().front().Denominator())/
    Long_Root_Length_Num();

	switch(Simple_Roots().size())
	{
		case 6:
			{
				GaugeGroupName Name('E', 6, KM_Level);
				Order_Simple_Roots('E');
				return GaugeGroup(Positive_Roots(), Name, Simple_Roots());
			}
			break;
		case 8:
			{
				GaugeGroupName Name('E', 8, KM_Level);
				return GaugeGroup(Positive_Roots(), Name, Simple_Roots());
			}

			break;
		default:
			{
				GaugeGroupName Name('E', 7, KM_Level);
				if(Simple_Roots().size() == 0)
					Find_Simple_Roots('E', 7);
				return GaugeGroup(Positive_Roots(), Name, Simple_Roots());
			}
	};
}//Close Build_E_Class_Group.

GaugeGroup GaugeGroupIdentifier::Build_F_Class_Group()
{
  int KM_Level = (2*Positive_Roots().front().Denominator())/
    Long_Root_Length_Num();
  GaugeGroupName Name('F', 4, KM_Level);
  if(Simple_Roots().size() == 0)
    Find_Simple_Roots(Name.Class(), Name.Rank());
  return GaugeGroup(Positive_Roots(), Name, Simple_Roots());
}//Close Build_F_Class_Group.

GaugeGroup GaugeGroupIdentifier::Build_G_Class_Group()
{
  int KM_Level = (2*Positive_Roots().front().Denominator()*
		  Positive_Roots().front().Denominator())/
    Long_Root_Length_Num();
  GaugeGroupName Name('G', 2, KM_Level);
  if(Simple_Roots().size() == 0)
    Find_Simple_Roots(Name.Class(), Name.Rank());
  return GaugeGroup(Positive_Roots(), Name, Simple_Roots());
}//Close Build_G_Class_Group.
