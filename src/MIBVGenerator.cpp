#include "MIBVGenerator.h"

MIBVGenerator::MIBVGenerator(int order, int large_st_dimensions)
{
  Order_ = order;
  Large_ST_Dimensions_ = large_st_dimensions;
  First_Complex_Extension_ = false;
  //ADD EXCEPTION HANDLING FOR ORDER = 1 AND LARGE_ST_DIMENSIONS BETWEEN 4
  //AND 10, INCLUSIVE
}//Close constructor.

MIBVGenerator::MIBVGenerator(const MIBVGenerator& New_MIBVGenerator)
{
  Order_=New_MIBVGenerator.Order();
  Large_ST_Dimensions_=New_MIBVGenerator.Large_ST_Dimensions();
  First_Complex_Extension_=New_MIBVGenerator.First_Complex_Extension();
  SP_LMs_=New_MIBVGenerator.SP_LMs();
  NSP_LMs_=New_MIBVGenerator.NSP_LMs();
  SP_RMs_Compact_=New_MIBVGenerator.SP_RMs_Compact();
  NSP_RMs_Compact_=New_MIBVGenerator.NSP_RMs_Compact();
  RMs_Observable_ = New_MIBVGenerator.RMs_Observable();
  RMs_Hidden_ = New_MIBVGenerator.RMs_Hidden();
  Observable_SU5_ = New_MIBVGenerator.Observable_SU5();
}//Close copy constructor.

//INTERFACE.
void MIBVGenerator::Build_Full_Chunks(const std::vector<BasisAlpha>& 
					Common_Basis_Alphas)
{
  First_Complex_Extension_ = Has_No_Complex_Elements(Common_Basis_Alphas);
  Build_LMs(Common_Basis_Alphas);
  if(Large_ST_Dimensions()<10)
    Build_RMs_Compact(Common_Basis_Alphas);
  Build_RMs_Observable(Common_Basis_Alphas);
  Build_RMs_Hidden(Common_Basis_Alphas);
}//Close Build_Chunks.

void MIBVGenerator::Build_Gauge_Chunks(const std::vector<BasisAlpha>& 
					 Common_Basis_Alphas)
{
  First_Complex_Extension_ = Has_No_Complex_Elements(Common_Basis_Alphas);
  //Add the LM.
  std::vector<int> Gauge_LM(Common_Basis_Alphas.at(0).LM_Size(), 0);
  SP_LMs_.push_back(Chunk(Gauge_LM, std::vector<int> 
			  (Common_Basis_Alphas.size()+1,0)));
  Build_RMs_Compact(Common_Basis_Alphas);
  Build_RMs_Observable(Common_Basis_Alphas);
  Build_RMs_Hidden(Common_Basis_Alphas);
}//Close Build_Gauge_Chunks.

void MIBVGenerator::Build_Observable_SU5_Chunks(const std::vector<BasisAlpha>&
						  Common_Basis_Alphas)
{
  if(!First_Complex_Extension())
    First_Complex_Extension_ = Has_No_Complex_Elements(Common_Basis_Alphas);
  int SU5_Limit = 0;
  int Eta_Limit = Order()-1;
  if(First_Complex_Extension())
    Eta_Limit = Order()/2;
  int Chunk_Start = 28-2*Large_ST_Dimensions();
  std::vector<int> Observable_SU5_Loader(16,0);
  if(Order()%2 == 0)
    SU5_Limit = Order()/2;
  else
    SU5_Limit = ((Order()+1)/2);
  for(int a=1; a<SU5_Limit; a++)
    {
      for(int b=0; b<10; b++)
	Observable_SU5_Loader.at(b) = a;

      for(int c=0; c<=Eta_Limit; c++)
	{
	  Observable_SU5_Loader.at(10) = c;
	  Observable_SU5_Loader.at(11) = c;
	  for(int d=0; d<=Eta_Limit; d++)
	    {
	      Observable_SU5_Loader.at(12) = d;
	      Observable_SU5_Loader.at(13) = d;
	      for(int e=0; e<=Eta_Limit; e++)
		{
		  Observable_SU5_Loader.at(14) = e;
		  Observable_SU5_Loader.at(15) = e;
		  Observable_SU5_.push_back(Chunk(Observable_SU5_Loader,
						  Compute_MI_Dot_Products
						  (Chunk_Start, 
						   Common_Basis_Alphas,
						   Observable_SU5_Loader)));
		  //No need for checking the simultaneous periodic modes
		  //for complex observables.
		}//Close for loop on bar-eta5,6
	    }//Close for loop on bar-eta3,4
	}//Close for loop on bar-eta1,2
    }//Close for loop on bar-psi1-10
}//Close Build_Observable_SU5_Chunks.



//DEBUG.
void MIBVGenerator::Display_Order() const
{
  std::cout<<"Order: "<<Order()<<std::endl<<std::endl;
}//Close Display_Order.

void MIBVGenerator::Display_Large_ST_Dimensions() const
{
  std::cout<<"Large_ST_Dimensions: "<<Large_ST_Dimensions()<<std::endl<<std::endl;
}//Close Display_Large_ST_Dimensions.

void MIBVGenerator::Display_LMs() const
{
  std::cout<<"Simply paired LMs: "<<SP_LMs().size()<<std::endl;
  std::list<Chunk>::const_iterator itSP_LMs = SP_LMs_.begin();
  for(; itSP_LMs != SP_LMs_.end(); ++itSP_LMs)
    {
      for(int b=0; b<static_cast<int>(itSP_LMs->BV_Chunk().size()); b++)
	std::cout<<itSP_LMs->BV_Chunk().at(b)<<" ";
      std::cout<<std::endl;
    }//Close for loop on SP_LMs.

  std::cout<<"Non simply paired LMs: "<<NSP_LMs().size()<<std::endl;
  std::list<Chunk>::const_iterator itNSP_LMs = NSP_LMs_.begin();
  for(; itNSP_LMs != NSP_LMs_.end(); ++itNSP_LMs)
    {
      for(int a=0; a<static_cast<int>(itNSP_LMs->BV_Chunk().size()); a++)
	std::cout<<itNSP_LMs->BV_Chunk().at(a)<<" ";
      std::cout<<std::endl;
    }//Close for loop on NSP_LMs.
  std::cout<<"Total LMs: "<<SP_LMs().size()+NSP_LMs().size()<<std::endl;
  std::cout<<std::endl;
}//Close Display_LMs.

void MIBVGenerator::Display_RMs_Compact() const
{
  std::cout<<"Simply paired compact RMs: "<<SP_RMs_Compact().size()<<std::endl;
  std::list<Chunk>::const_iterator itSP_RMs_Compact = 
    SP_RMs_Compact_.begin();
  for(; itSP_RMs_Compact != SP_RMs_Compact_.end(); 
      ++itSP_RMs_Compact)
    {
      for(int b=0; b<static_cast<int>(itSP_RMs_Compact->BV_Chunk().size()); b++)
	std::cout<<itSP_RMs_Compact->BV_Chunk().at(b)<<" ";
      std::cout<<std::endl;
    }//Close for loop on SP_RMs_Compact.
  
  std::cout<<"Non simply paired compact RMs: "<<NSP_RMs_Compact().size()<<std::endl;
  std::list<Chunk>::const_iterator itNSP_RMs_Compact = 
    NSP_RMs_Compact_.begin();
  for(; itNSP_RMs_Compact != NSP_RMs_Compact_.end(); 
      ++itNSP_RMs_Compact)
    {
      for(int a=0; a<static_cast<int>(itNSP_RMs_Compact->BV_Chunk().size()); a++)
	std::cout<<itNSP_RMs_Compact->BV_Chunk().at(a)<<" ";
      std::cout<<std::endl;
    }//Close for loop on NSP_RMs_Compact.

  std::cout<<"Total compact RMs: "<<SP_RMs_Compact().size()+NSP_RMs_Compact().size()
	   <<std::endl;
  std::cout<<std::endl;
}//Close Display_RMs_Compact.

void MIBVGenerator::Display_RMs_Observable() const
{
  std::cout<<"Observable RMs: "<<RMs_Observable().size()<<std::endl;
  std::list<Chunk>::const_iterator itRMs_Observable = 
    RMs_Observable_.begin();
  for(; itRMs_Observable != RMs_Observable_.end(); 
      ++itRMs_Observable)
    {
      for(int b=0; b<static_cast<int>(itRMs_Observable->BV_Chunk().size()); b++)
	std::cout<<itRMs_Observable->BV_Chunk().at(b)<<" ";
      std::cout<<std::endl;
    }//Close for loop on RMs_Observable.
  std::cout<<std::endl;
}//Close Display_RMs_Observable.

void MIBVGenerator::Display_RMs_Hidden() const
{
  std::cout<<"Hidden RMs: "<<RMs_Hidden().size()<<std::endl;
  std::list<Chunk>::const_iterator itRMs_Hidden = 
    RMs_Hidden_.begin();
  for(; itRMs_Hidden != RMs_Hidden_.end(); ++itRMs_Hidden)
    {
      for(int b=0; b<static_cast<int>(itRMs_Hidden->BV_Chunk().size()); b++)
	std::cout<<itRMs_Hidden->BV_Chunk().at(b)<<" ";
      std::cout<<std::endl;
    }//Close for loop on RMs_Hidden.
  std::cout<<std::endl;
}//Close Display_RMs_Hidden.

void MIBVGenerator::Display_Observable_SU5() const
{
  std::cout<<"Observable SU(5): "<<Observable_SU5().size()<<std::endl;
  std::list<Chunk>::const_iterator itObservable_SU5 = 
    Observable_SU5_.begin();
  for(; itObservable_SU5 != Observable_SU5_.end(); 
      ++itObservable_SU5)
    {
      for(int b=0; b<static_cast<int>(itObservable_SU5->BV_Chunk().size()); b++)
	std::cout<<itObservable_SU5->BV_Chunk().at(b)<<" ";
      std::cout<<std::endl;
    }//Close for loop on Observable_SU5.
  std::cout<<std::endl;
}//Close Display_Observable_SU5.

//PRIVATE.
void MIBVGenerator::Build_LMs(const std::vector<BasisAlpha>& Common_Basis_Alphas)
{
  LMGenerator The_LMGenerator(Large_ST_Dimensions(), Common_Basis_Alphas);
  The_LMGenerator.Build_LM_Chunks();
  std::list<std::vector<int> > SP_LMs_Set = The_LMGenerator.SP_LMs();
  std::list<std::vector<int> > NSP_LMs_Set = The_LMGenerator.NSP_LMs();

  std::list<std::vector<int> >::iterator itLMs_Set = SP_LMs_Set.begin();
  for(; itLMs_Set != SP_LMs_Set.end(); ++itLMs_Set)
    {
      SP_LMs_.push_back(Chunk(*itLMs_Set, Compute_MI_Dot_Products
			      (0, Common_Basis_Alphas, *itLMs_Set)));
      //We don't need to check the even simultaneous periodic modes since these
      //are simply paired.
    }//Close for loop on SP_LMs_Set.

  itLMs_Set = NSP_LMs_Set.begin();
  for(; itLMs_Set != NSP_LMs_Set.end(); ++itLMs_Set)
    {
      NSP_LMs_.push_back(Chunk(*itLMs_Set, Compute_MI_Dot_Products
			       (0, Common_Basis_Alphas, *itLMs_Set),
			       Compute_Simultaneous_Periodic_Modes
			       (0, Common_Basis_Alphas, *itLMs_Set)));
    }//Close for loop on NSP_LMs_Set.
}//Close Build_LMs.

void MIBVGenerator::Build_RMs_Compact(const std::vector<BasisAlpha>& 
					Common_Basis_Alphas)
{
  int Chunk_Start = 28-2*Large_ST_Dimensions()+16;
  RealChunkGenerator RC_Generator(Common_Basis_Alphas, Order(), 
				    First_Complex_Extension());
  RC_Generator.Build_Real_Chunks();

  std::set<std::vector<int> > SP_Real_Chunks_Set = 
    RC_Generator.SP_Real_Chunks();
  std::set<std::vector<int> >::iterator itSP_Real_Chunks_Set = 
    SP_Real_Chunks_Set.begin();
  for(; itSP_Real_Chunks_Set != SP_Real_Chunks_Set.end(); 
      ++itSP_Real_Chunks_Set)
    {
      SP_RMs_Compact_.push_back(Chunk(*itSP_Real_Chunks_Set, 
				      Compute_MI_Dot_Products
				      (Chunk_Start, Common_Basis_Alphas, 
				       *itSP_Real_Chunks_Set)));
      //There's no need to find simultaneous periodic modes for simply
      //paired chunks.
    }//Close for loop on itSP_Real_Chunks_Set.

  std::set<std::vector<int> > NSP_Real_Chunks_Set = 
    RC_Generator.NSP_Real_Chunks();
  std::set<std::vector<int> >::iterator itNSP_Real_Chunks_Set = 
    NSP_Real_Chunks_Set.begin();
  for(; itNSP_Real_Chunks_Set != NSP_Real_Chunks_Set.end();
      ++itNSP_Real_Chunks_Set)
    {
      NSP_RMs_Compact_.push_back(Chunk(*itNSP_Real_Chunks_Set,
				       Compute_MI_Dot_Products
				       (Chunk_Start,Common_Basis_Alphas, 
					*itNSP_Real_Chunks_Set), 
				       Compute_Simultaneous_Periodic_Modes
				       (Chunk_Start, Common_Basis_Alphas,
					*itNSP_Real_Chunks_Set)));
    }//Close for loop on itNSP_Real_Chunks_Set.
}//Close Build_RMs_Compact.

void MIBVGenerator::Build_RMs_Observable(const std::vector<BasisAlpha>& 
					   Common_Basis_Alphas)
{
  int Chunk_Start = 28-2*Large_ST_Dimensions();
  std::vector<int> Observable_BCs;
  for(int a=0; a<16; a++)
    Observable_BCs.push_back(a+Common_Basis_Alphas.at(0).LM_Size());

  ComplexChunkGenerator Obs_Generator(Common_Basis_Alphas, Order(),
					Observable_BCs, First_Complex_Extension());
  Obs_Generator.Build_Complex_Chunks();
  std::set<std::vector<int> > Observable_Chunks_Set = 
    Obs_Generator.Complex_Chunks();
  std::set<std::vector<int> >::iterator itObservable_Chunks_Set = 
    Observable_Chunks_Set.begin();
  for(; itObservable_Chunks_Set != Observable_Chunks_Set.end();
      ++itObservable_Chunks_Set)
    {
      RMs_Observable_.push_back(Chunk(*itObservable_Chunks_Set,
				      Compute_MI_Dot_Products
				      (Chunk_Start, Common_Basis_Alphas,
				       *itObservable_Chunks_Set)));
      //No need to compute the simultaneous periodic modes for 
      //complex obersvable chunk.
    }
}//Close Build_RMs_Observable.

void MIBVGenerator::Build_RMs_Hidden(const std::vector<BasisAlpha>&
				       Common_Basis_Alphas)
{
  int Chunk_Start = 28-2*Large_ST_Dimensions()+16+2*(10-Large_ST_Dimensions());
  std::vector<int> Hidden_BCs;
  for(int a=0; a<16; a++)
    Hidden_BCs.push_back(a + Common_Basis_Alphas.at(0).LM_Size() + 16 +
			 Common_Basis_Alphas.at(0).RM_Compact_Size());

  ComplexChunkGenerator Hid_Generator(Common_Basis_Alphas, Order(),
					Hidden_BCs, First_Complex_Extension());
  Hid_Generator.Build_Complex_Chunks();
  std::set<std::vector<int> > Hidden_Chunks_Set = 
    Hid_Generator.Complex_Chunks();
  std::set<std::vector<int> >::iterator itHidden_Chunks_Set = 
    Hidden_Chunks_Set.begin();
  for(; itHidden_Chunks_Set != Hidden_Chunks_Set.end();
      ++itHidden_Chunks_Set)
    {
      RMs_Hidden_.push_back(Chunk(*itHidden_Chunks_Set, Compute_MI_Dot_Products
				  (Chunk_Start, Common_Basis_Alphas, 
				   *itHidden_Chunks_Set)));
      //No need to compute the simultaneous periodic modes for 
      //complex hidden chunk.
    }//Close for loop on Hidden_Chunks_Set.
}//Close Build_RMs_Hidden.

std::vector<int> MIBVGenerator::Compute_MI_Dot_Products
(int Chunk_Start, const std::vector<BasisAlpha>& Common_Basis_Alphas,
 const std::vector<int>& BV_Chunk)
{
  std::vector<int> MI_Dot_Products;
  int Dot = 0;
  int N_i = Order();
  bool Is_LM = (Chunk_Start == 0);
  if(Is_LM && BV_Chunk.at(0)==1)//ST Fermion, LM Chunk
    N_i = FF::LCM(2, Order());
  else if(SP_LMs().front().BV_Chunk().at(0) == 1)//Also ST Fermion, RM Chunk. 
    N_i = FF::LCM(2, Order());
  int N_ij = FF::LCM(N_i, Common_Basis_Alphas.at(0).Denominator());
  std::vector<int> BA_Chunk = Make_Basis_Alpha(BV_Chunk, Is_LM);

  //First, do the self dot product.
  for(int a=0; a<static_cast<int>(BA_Chunk.size()); a++)
    Dot += (BA_Chunk.at(a)*BA_Chunk.at(a));
  MI_Dot_Products.push_back(N_i*Dot);

  //Now do the others.
  for(int a=0; a<static_cast<int>(Common_Basis_Alphas.size()); a++)
    {
      Dot = 0;
      for(int b=0; b<static_cast<int>(BA_Chunk.size()); b++)
	Dot+=(BA_Chunk.at(b)*Common_Basis_Alphas.at(a).Numerator().
	      at(b+Chunk_Start));
      MI_Dot_Products.push_back(N_ij*Dot);
    }//Close for loop on Common_Basis_Alphas.

  return MI_Dot_Products;
}//Close Compute_MI_Dot_Products.

std::vector<bool> MIBVGenerator::Compute_Simultaneous_Periodic_Modes
(int Chunk_Start, const std::vector<BasisAlpha>& Common_Basis_Alphas,
 const std::vector<int>& BV_Chunk)
{
  //This is to properly scale BV_Chunk, as it is in the integer coded form.
  int Scale_BV_Up = 2;
  //Left movers are order 2, and need to be scaled up to the RM order.
  //If the order is odd, then there must only be simply paired BVs,
  //so we don't need to check this condition in that case.
  if(Chunk_Start == 0 && Order()%2 == 0)
    Scale_BV_Up = FF::LCM(2, Order());
  std::vector<bool> Simultaneous_Periodic_Modes;
  for(int First_BA = 0; First_BA<static_cast<int>(Common_Basis_Alphas.size()); First_BA++)
    {
      for(int Second_BA = First_BA; Second_BA<static_cast<int>(Common_Basis_Alphas.size()); Second_BA++)
	{
	  int Simultaneous_Periodic_Mode_Counter = 0;
	  for(int a=0; a<static_cast<int>(BV_Chunk.size()); a++)
	    {
	      if((Common_Basis_Alphas.at(First_BA).Numerator().at(a+Chunk_Start) == 
		  Common_Basis_Alphas.at(0).Denominator())&&
		 (Common_Basis_Alphas.at(Second_BA).Numerator().at(a+Chunk_Start) == 
		  Common_Basis_Alphas.at(0).Denominator())&&
		 (Scale_BV_Up*BV_Chunk.at(a) == Order()))
		Simultaneous_Periodic_Mode_Counter++;
	    }//Close for loop on BV_Chunk.size().
	  Simultaneous_Periodic_Modes.push_back(Simultaneous_Periodic_Mode_Counter%2
						==0);
	}//Close for loop on Second_BA.
    }//Close for loop on First_BA.

  for(int First_BA = 0; First_BA<static_cast<int>(Common_Basis_Alphas.size()); ++First_BA)
    {
      int Simultaneous_Periodic_Mode_Counter = 0;
      for(int a=0; a<static_cast<int>(BV_Chunk.size()); ++a)
	{
	  if((Common_Basis_Alphas.at(First_BA).Numerator().at(a+Chunk_Start) == 
	      Common_Basis_Alphas.at(0).Denominator())&&
	     (Scale_BV_Up*BV_Chunk.at(a) == Order()))
	    Simultaneous_Periodic_Mode_Counter++;
	     }//Close for loop on BV_Chunk.size().
	  Simultaneous_Periodic_Modes.push_back(Simultaneous_Periodic_Mode_Counter%2
						==0);
	}//Close for loop on First_BA.

      int Simultaneous_Periodic_Mode_Counter = 0;
      for(int a=0; a<static_cast<int>(BV_Chunk.size()); ++a)
	{
	  if(Scale_BV_Up*BV_Chunk.at(a) == Order())
	    Simultaneous_Periodic_Mode_Counter++;
	} //Close for loop on BV_Chunk.size.
      Simultaneous_Periodic_Modes.push_back(Simultaneous_Periodic_Mode_Counter%2
					    ==0);

  return Simultaneous_Periodic_Modes;
}//Close Compute_Simultaneous_Periodic_Modes.

std::vector<int> MIBVGenerator::Make_Basis_Alpha
(const std::vector<int>& BV_Chunk, bool Is_LM)
{
  std::vector<int> BA_Chunk;
  int N = Order();
  if(Is_LM && BV_Chunk.at(0) == 1)//LM, ST Fermion
      N = FF::LCM(2, Order());
  else if(SP_LMs().front().BV_Chunk().at(0) == 1)//RM, ST Fermion.
    N = FF::LCM(2, Order());
  int Numerator_Conversion = (2*N)/Order();//Guaranteed to be an integer.
  //The proper form for a basis alpha is 2m/N, where N is the order and 
  //m is the basis chunk's integer code.
  if(Is_LM)//LM
    {
      for(int a=0; a<static_cast<int>(BV_Chunk.size()); a++)
	BA_Chunk.push_back(BV_Chunk.at(a)*N);
    }else//RM
    {
      for(int a=0; a<static_cast<int>(BV_Chunk.size()); a++)
	{
	  if(BV_Chunk.at(a)<=(double(Order())/double(2)))
	    BA_Chunk.push_back(Numerator_Conversion*BV_Chunk.at(a));
	  else
	    BA_Chunk.push_back(Numerator_Conversion*(BV_Chunk.at(a) - 
						     Order()));
	}//Close for loop on BV_Chunk.size().
    }//Close if/else on LM/RM.
  return BA_Chunk;
}//Close Make_Basis_Alpha.

bool MIBVGenerator::Has_No_Complex_Elements(const std::vector<BasisAlpha>& 
					    Common_Basis_Alphas)
{
  for(int a=0; a<static_cast<int>(Common_Basis_Alphas.size()); a++)
    {
      for(int b=Common_Basis_Alphas.at(a).LM_Size(); 
	  b<static_cast<int>(Common_Basis_Alphas.at(a).Numerator().size()); b++)
	{
	  if(Common_Basis_Alphas.at(a).Numerator().at(b) != 
	     Common_Basis_Alphas.at(a).Denominator() && 
	     Common_Basis_Alphas.at(a).Numerator().at(b) != 0)
	    return false;
	}//Close for loop on Common_Basis_Alphas.at(a).Numerator().size().
    }//Close for loop on Common_Basis_Alphas.size().
  return true;
}//Close Find_Complex_Elements.
