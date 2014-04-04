#include "InitialBuildParams.h"

//PUBLIC.
InitialBuildParams::InitialBuildParams(int Large_ST_Dimensions,
					   const std::vector<int>& Orders,
					   const std::vector
					   <std::vector<int> >&
					   Initial_BV_Set)
{
  Large_ST_Dimensions_ = Large_ST_Dimensions;
  Orders_.push_back(2);//For the '1' vector.
  for(int a=0; a<static_cast<int>(Orders.size()); a++)
    Orders_.push_back(Orders.at(a));
  Initial_BV_Set_ = Initial_BV_Set;
  k_ij_Extensions_.clear();
}//Close constructor.

InitialBuildParams::InitialBuildParams(const InitialBuildParams& 
					   New_InitialBuildParams)
{

  Large_ST_Dimensions_ = New_InitialBuildParams.Large_ST_Dimensions();
  Orders_ = New_InitialBuildParams.Orders();
  Initial_BV_Set_ = New_InitialBuildParams.Initial_BV_Set();
  k_ij_Extensions_ = New_InitialBuildParams.k_ij_Extensions();
}//Close copy constructor.

//DEBUG.
void InitialBuildParams::Display() const
{
  //Display the number of large ST dimensions.
  std::cout<<"Large ST Dimensions: "<<Large_ST_Dimensions()<<std::endl;
  //Display the orders.
  std::cout<<"Orders: ";
  for(int a=0; a<static_cast<int>(Orders().size()); a++)
    std::cout<<Orders().at(a)<<" ";
  std::cout<<std::endl;
  //Display the initial basis vectors.
  std::cout<<"Initial BV Set: "<<Initial_BV_Set().size()<<std::endl;
  for(int a=0; a<static_cast<int>(Initial_BV_Set().size()); a++)
    {
      for(int b=0; b<static_cast<int>(Initial_BV_Set().at(a).size()); b++)
	std::cout<<Initial_BV_Set().at(a).at(b)<<" ";
      std::cout<<std::endl;
    }//Close for loop on Initial_BV_Set.
  //Display the possible k_ij extensions.
  std::cout<<"k_ij matrix extensions: "<<k_ij_Extensions().size()<<std::endl;
  for(int a=0; a<static_cast<int>(k_ij_Extensions().size()); a++)
    {
      for(int b=0; b<static_cast<int>(k_ij_Extensions().at(a).size()); b++)
	std::cout<<k_ij_Extensions().at(a).at(b)<<" ";
      std::cout<<std::endl;
    }//Close for loop on k_ij_Extensions().size().
  std::cout<<std::endl;
}//Close Display.

//INTERFACE.
void InitialBuildParams::Build_k_ij_Extensions()
{
  k_ij_Extensions_.clear();
  std::vector<int> k_ij_Row(Initial_BV_Set().size()+1, 0);
  //The '+1' here is for the "1" vector, which is not included.
  Extend_k_ij(k_ij_Row, 0);
}//Close Build_k_ij_Extensions.

void InitialBuildParams::Load_Initial_BV_Set(const std::vector<int>& New_BV)
{
  Initial_BV_Set_.push_back(New_BV);
}//Close Load_Initial_BV_Set.

//PRIVATE.
//HELPERS.
void InitialBuildParams::Extend_k_ij(std::vector<int> k_ij_Row, 
				       int k_ij_index)
{
  if(k_ij_index != static_cast<int>(k_ij_Row.size()))
    {
      for(int a=0; a<=(Orders().at(k_ij_index)-1); ++a)
	{
	  k_ij_Row.at(k_ij_index) = a;
	  Extend_k_ij(k_ij_Row, k_ij_index+1);
	}//Close for loop on possible k_ij values. 
    }else
    k_ij_Extensions_.push_back(k_ij_Row);
}//Close Extend_k_ij.
