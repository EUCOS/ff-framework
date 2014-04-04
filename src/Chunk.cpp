#include "Chunk.h"

//PUBLIC.

Chunk::Chunk(const std::vector<int>& BV_Chunk, 
	     const std::vector<int>& MI_Dot_Products,
	     const std::vector<bool>& Simultaneous_Periodic_Modes)
{
  BV_Chunk_ = BV_Chunk;
  MI_Dot_Products_ = MI_Dot_Products;
  Simultaneous_Periodic_Modes_ = Simultaneous_Periodic_Modes;
}//Close constructor.

Chunk::Chunk(const std::vector<int>& BV_Chunk, 
	     const std::vector<int>& MI_Dot_Products)
{
  BV_Chunk_ = BV_Chunk;
  MI_Dot_Products_ = MI_Dot_Products;
}//Close second constructor.

Chunk::Chunk(const std::vector<int>& BV_Chunk)
{
  BV_Chunk_ = BV_Chunk;
  MI_Dot_Products_ = std::vector<int> (1,0);
}//Close third constructor.

Chunk::Chunk(const Chunk& New_Chunk)
{
  BV_Chunk_ = New_Chunk.BV_Chunk();
  MI_Dot_Products_ = New_Chunk.MI_Dot_Products();
  Simultaneous_Periodic_Modes_ = New_Chunk.Simultaneous_Periodic_Modes();
}//Close copy constructor.

//DEBUG.
void Chunk::Display_BV_Chunk() const 
{
  for(int a=0; a<static_cast<int>(BV_Chunk().size()); ++a)
    std::cout<<BV_Chunk().at(a)<<" ";
  std::cout<<std::endl;
}//Close Display_BV_Chunk.
void Chunk::Display_MI_Dot_Products() const
{
  for(int a=0; a<static_cast<int>(MI_Dot_Products().size()); ++a)
    std::cout<<MI_Dot_Products().at(a)<<" ";
  std::cout<<std::endl;
}//Close Display_MI_Dot_Products.

//PRIVATE.
//HELPERS.
