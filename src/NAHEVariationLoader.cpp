#include "NAHEVariationLoader.h"

NAHEVariationLoader::NAHEVariationLoader()
{
  b1_ = Build_b1();
  b2_ = Build_b2();
  b3_ = Build_b3();
}//Close constructor.

NAHEVariationLoader::NAHEVariationLoader
(const NAHEVariationLoader& Old_NAHEVariationLoader)
{
  b1_ = Old_NAHEVariationLoader.b1();
  b2_ = Old_NAHEVariationLoader.b2();
  b3_ = Old_NAHEVariationLoader.b3();
}//Close copy constructor.

Model NAHEVariationLoader::Load_NAHE_Variation(const Model& Old_Model, bool With_S)
{
  Model New_Model = Old_Model;

  //Load the basis vectors.
  New_Model.Load_BV_Set(b1());
  New_Model.Load_BV_Set(b2());
  New_Model.Load_BV_Set(b3());

  //Load the GSO coefficients.
  std::vector<int> k_ij_Loader;
  k_ij_Loader.push_back(0);
  New_Model.Load_k_ij_Row(k_ij_Loader);//1 vector.
  if(With_S)
    {
      New_Model.Load_k_ij_Row(k_ij_Loader);
      k_ij_Loader.at(0) = 1;
    }else
    k_ij_Loader.clear();
  for(int a=0; a<3; a++)
    {
      k_ij_Loader.push_back(1);
      New_Model.Load_k_ij_Row(k_ij_Loader);
    }//Close for loop for b1, b2, b3.

  New_Model.Set_NAHE_Loaded(true);

  return New_Model;
}//Close Load_NAHE_Variation.

//DEBUG.
void NAHEVariationLoader::Display_b1() const
{
  std::cout<<"b1: "<<std::endl;
  b1().Display();
}//Close Display_b1.

void NAHEVariationLoader::Display_b2() const
{
  std::cout<<"b2: "<<std::endl;
  b2().Display();
}//Close Display_b2.

void NAHEVariationLoader::Display_b3() const
{
  std::cout<<"b3: "<<std::endl;
  b3().Display();
}//Close Display_b3.

//PRIVATE.
BasisVector NAHEVariationLoader::Build_b1()
{
  int b1_Array [64] = {1,1,//ST1,2
		       1,0,0,//x1,y1,w1
		       1,0,0,//x2,y2,w2
		       0,1,0,//x3,y3,w3
		       0,1,0,//x4,y4,w4
		       0,1,0,//x5,y5,w5
		       0,1,0,//x6,y6,w6
		       1,1,1,1,1,1,1,1,1,1,//bar-psi1-10
		       1,1,0,0,0,0,//bar-eta1-6
		       0,0,1,1,1,1,//bar-y1-6
		       0,0,0,0,0,0,//bar-w1-6
		       0,0,0,0,0,0,//bar-etaprime1-6
		       0,0,0,0,0,0,0,0,0,0};//bar-psiprime1-10

  return BasisVector(Make_Vector(b1_Array), 2, 4);
}//Build_b1.

BasisVector NAHEVariationLoader::Build_b2()
{
  int b2_Array [64] = {1,1,//ST1,2
		       0,1,0,//x1,y1,w1
		       0,1,0,//x2,y2,w2
		       1,0,0,//x3,y3,w3
		       1,0,0,//x4,y4,w4
		       0,1,0,//x5,y5,w5
		       0,1,0,//x6,y6,w6
		       1,1,1,1,1,1,1,1,1,1,//bar-psi1-10
		       0,0,1,1,0,0,//bar-eta1-6
		       1,1,0,0,1,1,//bar-y1-6
		       0,0,0,0,0,0,//bar-w1-6
		       0,0,0,0,0,0,//bar-etaprime1-6
		       0,0,0,0,0,0,0,0,0,0};//bar-psiprime1-10

  return BasisVector(Make_Vector(b2_Array), 2, 4);
}//Build_b2.

BasisVector NAHEVariationLoader::Build_b3()
{
  int b3_Array [64] = {1,1,//ST1,2
		       0,1,0,//x1,y1,w1
		       0,1,0,//x2,y2,w2
		       0,1,0,//x3,y3,w3
		       0,1,0,//x4,y4,w4
		       1,0,0,//x5,y5,w5
		       1,0,0,//x6,y6,w6
		       1,1,1,1,1,1,1,1,1,1,//bar-psi1-10
		       0,0,0,0,1,1,//bar-eta1-6
		       1,1,1,1,0,0,//bar-y1-6
		       0,0,0,0,0,0,//bar-w1-6
		       0,0,0,0,0,0,//bar-etaprime1-6
		       0,0,0,0,0,0,0,0,0,0};//bar-psiprime1-10

  return BasisVector(Make_Vector(b3_Array), 2, 4);
}//Build_b3

std::vector<int> NAHEVariationLoader::Make_Vector(int Array[64])
{
  std::vector<int> New_Vector;
  for(int a=0; a<64; a++)
    New_Vector.push_back(Array[a]);

  return New_Vector;
}//Close Make_Vector.
