#include "OutputWriter.h"

OutputWriter::OutputWriter(const OutputWriter& New_OutputWriter)
{
  ;
}//Close copy constructor.

//MODEL CLASS INTERFACE.
void OutputWriter::Write_Model_Particle_Content(std::ofstream& Model_Out,
						const Model& FFHS_Model)
{
  //Write the gauge group names.
  for(int a=0; a<static_cast<int>(FFHS_Model.Gauge_Groups().size()); a++)
    Model_Out<<FFHS_Model.Gauge_Groups().at(a).Name().Class()<<" "
	     <<FFHS_Model.Gauge_Groups().at(a).Name().Rank()<<" "
	     <<FFHS_Model.Gauge_Groups().at(a).Name().KM_Level()<<" ";
  Model_Out<<std::endl;

	//Write the matter representations.
	for(int a=0; a<static_cast<int>(FFHS_Model.MatterRepresentations().size()); ++a)
	{
		Model_Out<<FFHS_Model.MatterRepresentations().at(a).Duplicates()<<": ";
		for(int b=0;
      b<static_cast<int>(FFHS_Model.MatterRepresentations().at(a).Rep_Dimension().size());
      ++b)
		{
			Model_Out<<FFHS_Model.MatterRepresentations().at(a).Rep_Dimension().
				at(b).Dimension()<<FFHS_Model.MatterRepresentations().at(a).
				Rep_Dimension().at(b).Triality()<<" ";
		}
		Model_Out<<std::endl;
	}//CLose for loop over the matter representations.

  //Now write the number of U(1) factors.
  Model_Out<<"U(1)'s: "<<FFHS_Model.U1_Factors()<<std::endl;
  //Now write the number of ST SUSYs.
  Model_Out<<"ST SUSY: "<<FFHS_Model.SUSY_States().size()<<std::endl;
}//Close Write_Model_Particle_Content.

void OutputWriter::Write_Model_BVs(std::ofstream& Model_Out, 
				    const Model& FFHS_Model)
{
  int Finish = FFHS_Model.BV_Set().size();
  for(int a=0; a<Finish; a++)
    {
      for(int b=0; b<static_cast<int>(FFHS_Model.BV_Set().at(a).BV().size()); b++)
	Model_Out<<FFHS_Model.BV_Set().at(a).BV().at(b)<<" ";
      Model_Out<<std::endl;
    }
}//Close Write_Model_BVs.

void OutputWriter::Write_Model_Extension_BVs(std::ofstream& Model_Out,
					      const Model& FFHS_Model, 
					      int Extensions)
{
  int Finish = FFHS_Model.BV_Set().size();
  int Start = Finish - Extensions;
  int BV_Size = FFHS_Model.BV_Set().at(0).BV().size();
  for(int a=Start; a<Finish; ++a)
    {
      for(int b=0; b<BV_Size; ++b)
	Model_Out<<FFHS_Model.BV_Set().at(a).BV().at(b)<<" ";
      Model_Out<<std::endl;
    }
}//Close Write_Model_Extension_BVs.


void OutputWriter::Write_Model_k_ij(std::ofstream& Model_Out, 
				     const Model& FFHS_Model)
{
  for(int a=0; a<static_cast<int>(FFHS_Model.k_ij().Numerators().size()); a++)
    {
      for(int b=0; b<static_cast<int>(FFHS_Model.k_ij().Numerators().at(a).size()); b++)
	Model_Out<<FFHS_Model.k_ij().Numerators().at(a).at(b)<<" ";
      Model_Out<<std::endl;
    }
}//Close Write_Model_k_ij.

void OutputWriter::Write_Model_Extension_k_ij(std::ofstream& Model_Out,
					       const std::vector<std::vector<int> >&
					       k_ij_Extensions)
{
  for(int a=0; a<static_cast<int>(k_ij_Extensions.size()); ++a)
    {
      for(int b=0; b<static_cast<int>(k_ij_Extensions.at(a).size()); ++b)
	Model_Out<<k_ij_Extensions.at(a).at(b)<<" ";
      Model_Out<<std::endl;
    }//Close for loop on k_ij_Extensions.
}//Close Write_Model_Extension_k_ij.

void OutputWriter::Write_Model_BV_Orders(std::ofstream& Model_Out,
					  const Model& FFHS_Model)
{
  for(int a=0; a<static_cast<int>(FFHS_Model.BV_Set().size()); a++)
    Model_Out<<FFHS_Model.BV_Set().at(a).Order()<<" ";
  Model_Out<<std::endl;

}//Close Write_Model_BV_Orders.

void OutputWriter::Write_BV_Search_Front_Info(std::ofstream& Model_Out,
						 const Model& FFHS_Model,
						 int Large_ST_Dimensions,
						 const std::vector<int>& 
						 Extension_Orders)
{
  Model_Out<<"Large ST Dimensions: "<<std::endl;
  Model_Out<<Large_ST_Dimensions<<std::endl;
  Model_Out<<"Initial Orders: "<<std::endl;
  Write_Model_BV_Orders(Model_Out, FFHS_Model);
  Model_Out<<"Extension Orders: "<<std::endl;
  for(int a=0; a<static_cast<int>(Extension_Orders.size()); a++)
    Model_Out<<Extension_Orders.at(a)<<" ";
  Model_Out<<std::endl;
  Model_Out<<"Extending: "<<std::endl;
  Write_Model_BVs(Model_Out, FFHS_Model);
  Model_Out<<std::endl;
  Model_Out<<"Initial k_ij: "
	   <<std::endl;
  Write_Model_k_ij(Model_Out, FFHS_Model);
  Model_Out<<std::endl;
  for(int a=0; a<static_cast<int>(FFHS_Model.BV_Set().at(0).BV().size()); a++)
    Model_Out<<"--";
  Model_Out<<std::endl;
}//Close Write_Syst_Search_Front_Info.

//LEEFT CLASS INTERFACE.
void OutputWriter::Write_LEEFT_Particle_Content(std::ofstream& LEEFT_Out,
						 const LEEFT& FFHS_LEEFT)
{
  for(int a=0; a<static_cast<int>(FFHS_LEEFT.Gauge_Groups().size()); a++)
    LEEFT_Out<<FFHS_LEEFT.Gauge_Groups().at(a).Class()<<" "
	     <<FFHS_LEEFT.Gauge_Groups().at(a).Rank()<<" "
	     <<FFHS_LEEFT.Gauge_Groups().at(a).KM_Level()<<" ";
  LEEFT_Out<<std::endl;
	
	//Write the matter representations.
	for(int a=0; a<static_cast<int>(FFHS_LEEFT.MatterRepresentations().size()); ++a)
	{
		LEEFT_Out<<FFHS_LEEFT.MatterRepresentations().at(a).
			Duplicates()<<": ";
		for(int b=0;
      b<static_cast<int>(FFHS_LEEFT.MatterRepresentations().at(a).Rep_Dimension().size());
      ++b)
			LEEFT_Out<<FFHS_LEEFT.MatterRepresentations().at(a).Rep_Dimension().
				at(b).Dimension()<<FFHS_LEEFT.MatterRepresentations().at(a).
				Rep_Dimension().at(b).Triality()<<" ";
		LEEFT_Out<<std::endl;
	}//CLose for loop over the matter representations.
  LEEFT_Out<<"U(1)'s: "<<FFHS_LEEFT.U1_Factors()<<std::endl;
  LEEFT_Out<<"ST SUSYs: "<<FFHS_LEEFT.ST_SUSYs()<<std::endl;


	LEEFT_Out<<std::endl;
	if(FFHS_LEEFT.ObservableSectors().size() > 0)
		LEEFT_Out<<"Observable sectors: "<<std::endl;
	for(int a=0; a<static_cast<int>(FFHS_LEEFT.ObservableSectors().size()); ++a)
	{
		LEEFT_Out<<std::endl;
		LEEFT_Out<<"Generations: "<<std::endl;
		for(int b=0;
      b<static_cast<int>(FFHS_LEEFT.ObservableSectors().at(a).Generations().size());
      ++b)
		{
			LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
				Generations().at(b).Duplicates()<<": ";
			for(int c=0;
        c<static_cast<int>(FFHS_LEEFT.ObservableSectors().at(a).Generations().at(b).Rep_Dimension().size());
        ++c)
			{
				LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
					Generations().at(b).Rep_Dimension().at(c).Dimension();
				LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
					Generations().at(b).Rep_Dimension().at(c).Triality()<<" ";
			}
			LEEFT_Out<<std::endl;
		}//Close for loop on Generations.

		//LEEFT_Out<<"Barred generations: "<<std::endl;
		for(int b=0; 
      b<static_cast<int>(FFHS_LEEFT.ObservableSectors().at(a).Barred_Generations().size());
      ++b)
		{
			LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
				Barred_Generations().at(b).Duplicates()<<": ";
			for(int c=0; 
				c<static_cast<int>(FFHS_LEEFT.ObservableSectors().at(a).Barred_Generations().at(b).Rep_Dimension().size());
        ++c)
			{
				LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
					Barred_Generations().at(b).Rep_Dimension().at(c).Dimension();
				LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
					Barred_Generations().at(b).Rep_Dimension().at(c).Triality()<<" ";
			}
			LEEFT_Out<<std::endl;
		}//Close for loop on Barred_Generations.

		//LEEFT_Out<<"Anti generations: "<<std::endl;
		for(int b=0;
      b<static_cast<int>(FFHS_LEEFT.ObservableSectors().at(a).Anti_Generations().size());
      ++b)
		{
			LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
				Anti_Generations().at(b).Duplicates()<<": ";
			for(int c=0;
        c<static_cast<int>(FFHS_LEEFT.ObservableSectors().at(a).Anti_Generations().at(b).Rep_Dimension().size());
        ++c)
			{
				LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
					Anti_Generations().at(b).Rep_Dimension().at(c).Dimension();
				LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
					Anti_Generations().at(b).Rep_Dimension().at(c).Triality()<<" ";
			}
			LEEFT_Out<<std::endl;
		}//Close for loop on Anti_Generations.

		//LEEFT_Out<<"Anti barred generations: "<<std::endl;
		for(int b=0;
      b<static_cast<int>(FFHS_LEEFT.ObservableSectors().at(a).Anti_Barred_Generations().size());
      ++b)
		{
			LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
				Anti_Barred_Generations().at(b).Duplicates()<<": ";
			for(int c=0;
        c<static_cast<int>(FFHS_LEEFT.ObservableSectors().at(a).Anti_Barred_Generations().at(b).Rep_Dimension().size());
        ++c)
			{
				LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
					Anti_Barred_Generations().at(b).Rep_Dimension().at(c).Dimension();
				LEEFT_Out<<FFHS_LEEFT.ObservableSectors().at(a).
					Anti_Barred_Generations().at(b).Rep_Dimension().at(c).Triality();
				LEEFT_Out<<" ";
			}
			LEEFT_Out<<std::endl;
		}//Close for loop on Anti_Barred_Generations.
	}//Close for loop on ObservableSectors.
  LEEFT_Out<<std::endl<<std::endl;
}//Close Write_LEEFT_Particle_Content.
