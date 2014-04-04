#include "Statistics.h"

Statistics::Statistics(const Statistics& New_Statistics)
{
  ;
}//Close copy constructor.

std::map<GaugeGroupName, int> Statistics::Count_Gauge_Groups
(const std::set<LEEFT>& Unique_LEEFTs)
{
  GaugeGroupName U1('U',1,1);
  std::map<GaugeGroupName, int> Gauge_Groups;
  std::set<LEEFT>::const_iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
  for(; itUnique_LEEFTs != Unique_LEEFTs.end(); ++itUnique_LEEFTs)
    {
      std::set<GaugeGroupName> Unique_GaugeGroupNames;
      //Unique gauge group factors for a given model.
      for(int a=0; a<static_cast<int>((itUnique_LEEFTs->Gauge_Groups()).size()); a++)
	Unique_GaugeGroupNames.insert((itUnique_LEEFTs->Gauge_Groups()).at(a));
      if(itUnique_LEEFTs->U1_Factors()>0)
	Unique_GaugeGroupNames.insert(U1);
      std::set<GaugeGroupName>::iterator itUnique_GaugeGroupNames = 
	Unique_GaugeGroupNames.begin();
      for(; itUnique_GaugeGroupNames != Unique_GaugeGroupNames.end(); 
	  ++itUnique_GaugeGroupNames)
	{
	  if(Gauge_Groups.find(*itUnique_GaugeGroupNames) != 
	     Gauge_Groups.end())
	    Gauge_Groups[*itUnique_GaugeGroupNames]++;
	  else
	    Gauge_Groups[*itUnique_GaugeGroupNames] = 1;

	}//Close for loop on Unique_GaugeGroupNames.
    }//Close for loop on Unique_LEEFTs.
  return Gauge_Groups;
}//Close Count_Gauge_Group_Factors.

std::map<int, int> Statistics::Count_ST_SUSYs(const std::set<LEEFT>& Unique_LEEFTs)
{
  std::map<int, int> ST_SUSYs;
  std::set<LEEFT>::const_iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
  for(; itUnique_LEEFTs != Unique_LEEFTs.end(); ++itUnique_LEEFTs)
    {
      if(ST_SUSYs.find(itUnique_LEEFTs->ST_SUSYs()) != ST_SUSYs.end())
	ST_SUSYs[itUnique_LEEFTs->ST_SUSYs()]++;
      else
	ST_SUSYs[itUnique_LEEFTs->ST_SUSYs()] = 1;
    }//Close for loop on Unique_LEEFTs.
  return ST_SUSYs;
}//Close Count_ST_SUSYs.

std::map<std::vector<GaugeGroupName>, int> 
Statistics::Count_Gauge_Group_Combinations(const std::set<LEEFT>& Unique_LEEFTs)
{
  std::set<LEEFT>::const_iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
  std::set<LEEFT>::const_iterator itLEEFT_End = Unique_LEEFTs.end();
  std::map<std::vector<GaugeGroupName>, int> Gauge_Group_Combinations;
  for(; itUnique_LEEFTs != itLEEFT_End; ++itUnique_LEEFTs)
    {
      //Build the gauge group combinations.
      std::set<std::vector<GaugeGroupName> > Unique_Gauge_Group_Combinations;
      std::vector<GaugeGroupName> Gauge_Groups = itUnique_LEEFTs->Gauge_Groups();
      if(itUnique_LEEFTs->U1_Factors()>0)
	Gauge_Groups.push_back(GaugeGroupName('U',1,1));
      int Gauge_Group_End = Gauge_Groups.size();
      for(int a=0; a<Gauge_Group_End; ++a)
	{
	  for(int b=(a+1); b<Gauge_Group_End; ++b)
	    {
	      std::vector<GaugeGroupName> Gauge_Group_Factors;
	      Gauge_Group_Factors.push_back
		(Gauge_Groups.at(a));
	      Gauge_Group_Factors.push_back
		(Gauge_Groups.at(b));
	      Unique_Gauge_Group_Combinations.insert(Gauge_Group_Factors);
	    }//Close inner for loop on gauge groups for a given model.
	}//Close outer for loop on gauge groups for a given model.
      std::set<std::vector<GaugeGroupName> >::iterator 
	itUnique_Gauge_Group_Combinations = Unique_Gauge_Group_Combinations.begin();
      std::set<std::vector<GaugeGroupName> >::iterator itUnique_Gauge_Group_End = 
	Unique_Gauge_Group_Combinations.end();
      for(; itUnique_Gauge_Group_Combinations != itUnique_Gauge_Group_End;
          ++itUnique_Gauge_Group_Combinations)
	{
	  if(Gauge_Group_Combinations.find(*itUnique_Gauge_Group_Combinations) != 
	     Gauge_Group_Combinations.end())
	    Gauge_Group_Combinations[*itUnique_Gauge_Group_Combinations]++;
	  else
	    Gauge_Group_Combinations[*itUnique_Gauge_Group_Combinations] = 1;
	}
    }//Close for loop on unique LEEFTs.
  return Gauge_Group_Combinations;
}//Close Count_Gauge_Group_Combinations.

std::map<int, int> Statistics::Count_U1_Factors(const std::set<LEEFT>& 
						Unique_LEEFTs)
{
  std::map<int, int> U1_Factors;
  std::set<LEEFT>::const_iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
  for(; itUnique_LEEFTs != Unique_LEEFTs.end(); ++itUnique_LEEFTs)
    {
      if(U1_Factors.find(itUnique_LEEFTs->U1_Factors()) != U1_Factors.end())
	U1_Factors[itUnique_LEEFTs->U1_Factors()]++;
      else
	U1_Factors[itUnique_LEEFTs->U1_Factors()] = 1;
    }//Close for loop on Unique_LEEFTs.
  return U1_Factors;
}//Close Count_U1_Factors.

std::map<int, int> Statistics::Count_Gauge_Group_Factors(const std::set<LEEFT>& 
							 Unique_LEEFTs)
{
  std::map<int, int> Gauge_Group_Factors;
  std::set<LEEFT>::const_iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
  std::set<LEEFT>::const_iterator itLEEFT_End = Unique_LEEFTs.end();
  for(; itUnique_LEEFTs != itLEEFT_End; ++itUnique_LEEFTs)
    {
      //Count the number of non-Abelian gauge groups.
      int Gauge_Group_Factor_Count = itUnique_LEEFTs->Gauge_Groups().size();
      //Add the U(1)'s.
      Gauge_Group_Factor_Count += itUnique_LEEFTs->U1_Factors();
      //Increment the proper spot on the map.
      if(Gauge_Group_Factors.find(Gauge_Group_Factor_Count) != 
	 Gauge_Group_Factors.end())
	Gauge_Group_Factors[Gauge_Group_Factor_Count]++;
      else
	Gauge_Group_Factors[Gauge_Group_Factor_Count] = 1;
    }
  return Gauge_Group_Factors;
}//Close Count_Gauge_Group_Factors.

std::map<int, int> Statistics::Count_NA_Singlets(const std::set<LEEFT>&
		Unique_LEEFTs)
{
	std::map<int, int> NA_Singlets;
	std::set<LEEFT>::const_iterator itUnique_LEEFTs = Unique_LEEFTs.begin();
	for(; itUnique_LEEFTs != Unique_LEEFTs.end();
			++itUnique_LEEFTs)
	{
		//First, count the number of Non-Abelian singlets in the LEEFT.
		int NA_Singlet_Count = 0;
		for(int a=0; a<static_cast<int>(itUnique_LEEFTs->MatterRepresentations().size());
				++a)
		{
			if(itUnique_LEEFTs->
					MatterRepresentations().at(a).Rep_Dimension() == 
					std::vector<GroupRepresentation> (itUnique_LEEFTs->
						MatterRepresentations().at(a).Rep_Dimension().size(), 
						GroupRepresentation(1,' ')))
			{
				NA_Singlet_Count+= 
					itUnique_LEEFTs->MatterRepresentations().at(a).Duplicates();
				break;
			}//Close if statement for finding the singlets.
		}//Close for loop on Representations.

		//Now increment the map.
		if(NA_Singlets.find(NA_Singlet_Count) != NA_Singlets.end())
			NA_Singlets[NA_Singlet_Count]++;
		else
			NA_Singlets[NA_Singlet_Count] = 1;
	}//Close for loop on the Unique_LEEFTs.
	return NA_Singlets;
}//Close Count_NA_Singlets.
