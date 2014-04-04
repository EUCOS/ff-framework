#include "LatexWriter.h"

//PUBLIC.
LatexWriter::LatexWriter(const LatexWriter& New_LatexWriter)
{
  ;
}//Close copy constructor.

//INTERFACE.
void LatexWriter::Begin_Latex_Document(std::ofstream& Latex_Out, 
					std::string Doc_Title)
{
  Latex_Out<<"\\documentclass [12pt] {article}"<<std::endl;
  Latex_Out<<"\\usepackage{tikz, pgfplots, supertabular,"; 	
	Latex_Out<<"longtable}"<<std::endl;
  Latex_Out<<"\\author{Timothy Renner}"<<std::endl;
  Latex_Out<<"\\title{"<<Doc_Title<<"}"<<std::endl;
  Latex_Out<<"\\begin{document}"<<std::endl;
  Latex_Out<<"\\maketitle"<<std::endl;
  Latex_Out<<"\\begin{flushleft}"<<std::endl;
}//Close Begin_Latex_Document.

void LatexWriter::Finish_Latex_Document(std::ofstream& Latex_Out)
{
  Latex_Out<<"\\end{flushleft}"<<std::endl;
  Latex_Out<<"\\end{document}"<<std::endl;
}//Close Finish_Latex_Document.

//MODEL INTERFACE.
void LatexWriter::Print_Model_BVs(std::ofstream& Latex_Out, const Model& FFHS_Model)
{
  std::vector<std::string> BV_Names = Get_BV_Names(FFHS_Model);

  Write_Left_Movers(BV_Names, Latex_Out, FFHS_Model);
  Latex_Out<<"\\newpage"<<std::endl;
  Write_Right_Movers(BV_Names, Latex_Out, FFHS_Model);
}//Close Print_Model_BVs.


void LatexWriter::Print_Model_k_ij(std::vector<std::string> BV_Names, 
						std::ofstream& Latex_Out, const Model& FFHS_Model)
{
  if(BV_Names.size() == 0)
    BV_Names = Get_BV_Names(FFHS_Model);

  Latex_Out<<"\\textbf{$k_{ij}\\ Matrix\\ \\times{"
		 <<FFHS_Model.k_ij().Denominators().at(0);
  Latex_Out<<"}:\\ $}"<<"\\\\"<<std::endl;
  Latex_Out<<"\\vspace{2 mm}"<<std::endl;
  Latex_Out<<"$\\left(\\begin{tabular}{c|";
  for(int a=0; a<static_cast<int>(FFHS_Model.k_ij().Numerators().size()); a++)
    Latex_Out<<"c";
  Latex_Out<<"}"<<std::endl;

  for(int a=0; a<static_cast<int>(BV_Names.size()); a++)
    Latex_Out<<"&"<<BV_Names.at(a);
  Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;

  for(int a=0; a<static_cast<int>(FFHS_Model.k_ij().Numerators().size()); a++)
    {
      Latex_Out<<BV_Names.at(a);
      for(int b=0; b<static_cast<int>(FFHS_Model.k_ij().Numerators().at(a).size()); b++)
	Latex_Out<<"&"<<FFHS_Model.k_ij().Numerators().at(a).at(b);
      Latex_Out<<"\\\\"<<std::endl;
    }
  Latex_Out<<"\\end{tabular}\\right)$"<<std::endl;
  Latex_Out<<"\\\\"<<std::endl;
  Latex_Out<<"\\vspace{2 mm}"<<std::endl;
}//Close Print_Model_k_ij.

void LatexWriter::Print_Model_Inputs(std::ofstream& Latex_Out, const Model&
							FFHS_Model)
{
  std::vector<std::string> BV_Names = Get_BV_Names(FFHS_Model);

  Write_Left_Movers(BV_Names, Latex_Out, FFHS_Model);
  Print_Model_k_ij(BV_Names,Latex_Out, FFHS_Model);
  Latex_Out<<"\\newpage"<<std::endl;
  Write_Right_Movers(BV_Names, Latex_Out, FFHS_Model);
    
}//Close Print_Model_Inputs.

void LatexWriter::Print_Model_Particle_Content(std::ofstream& Latex_Out, 
						const Model& FFHS_Model)
{
  std::vector<std::string> GaugeGroupNames =
		Get_GaugeGroupNames(FFHS_Model);

	Write_Matter_Reps(Latex_Out, GaugeGroupNames, 
			FFHS_Model.MatterRepresentations());
  Latex_Out<<"\\textbf{Total Matter Representations: "<<
		FFHS_Model.MatterStates().size()<<"}"
		 <<std::endl;
  Latex_Out<<"\\\\"<<std::endl<<"\\vspace{3 mm}"<<std::endl;
  Latex_Out<<"\\textbf{Number of $U(1)$'s: "<<FFHS_Model.U1_Factors()<<"}";
  Latex_Out<<"\\\\"<<std::endl;
  Latex_Out<<"\\textbf{Number of ST SUSYs: "
		<<FFHS_Model.SUSY_States().size()<<
    "}"<<std::endl;
}//Close Print_Model_Particle_Content.

void LatexWriter::Print_Full_Model(std::ofstream& Latex_Out, 
						const Model& FFHS_Model)
{
  Print_Model_Particle_Content(Latex_Out, FFHS_Model);
  Latex_Out<<"\\newpage"<<std::endl;
  Print_Model_Inputs(Latex_Out, FFHS_Model);
}//Close Print_Full_Model.

void LatexWriter::Print_Single_Model_Document(std::ofstream& Latex_Out, 
								 const Model& FFHS_Model)
{
  Begin_Latex_Document(Latex_Out, "Model Results");
  Print_Full_Model(Latex_Out, FFHS_Model);
  Finish_Latex_Document(Latex_Out);
}//Close Print_Single_Model_Document.

//LEEFT INTERFACE.
void LatexWriter::Print_LEEFT_Particle_Content(std::ofstream& Latex_Out,
						const LEEFT& FFHS_LEEFT)
{
  //First get the Gauge Group names.
 std::vector<std::string> GaugeGroupNames;
  for(int a=0; a<static_cast<int>(FFHS_LEEFT.Gauge_Groups().size()); a++)
		GaugeGroupNames.push_back(Cast_GaugeGroupName
				(FFHS_LEEFT.Gauge_Groups().at(a)));

	Write_Matter_Reps(Latex_Out, GaugeGroupNames, 
			FFHS_LEEFT.MatterRepresentations());
	
	int Total_Matter_Rep_Count = 0;
	for(int a=0; a<static_cast<int>(FFHS_LEEFT.MatterRepresentations().size()); ++a)
		Total_Matter_Rep_Count += FFHS_LEEFT.MatterRepresentations().at(a).
			Duplicates();

	Latex_Out<<"\\textbf{Total Matter Representations: "<<
		Total_Matter_Rep_Count<<"}\\\\"<<std::endl;
	Latex_Out<<"\\vspace{3 mm}"<<std::endl;
  Latex_Out<<"\\textbf{Number of $U(1)$'s: "<<FFHS_LEEFT.U1_Factors()<<"}";
  Latex_Out<<"\\\\"<<std::endl;
  Latex_Out<<"\\textbf{Number of ST SUSYs: "<<FFHS_LEEFT.ST_SUSYs()<<
    "}\\\\"<<std::endl;
  Latex_Out<<"\\vspace{3 mm}"<<std::endl;
}//Close Print_LEEFT_Particle_Content.

void LatexWriter::Print_Single_LEEFT_Document(std::ofstream& Latex_Out,
								 const LEEFT& FFHS_LEEFT)
{
  Begin_Latex_Document(Latex_Out, "Model Results");
  Print_LEEFT_Particle_Content(Latex_Out, FFHS_LEEFT);
  Finish_Latex_Document(Latex_Out);
}//Close Print_Single_LEEFT_Document.

//STATS INTERFACE.
void LatexWriter::Write_Gauge_Group_Count(std::ofstream& Latex_Out,
						 const std::map<GaugeGroupName, int>& 
						 Gauge_Group_Counts, 
						 int Total_Unique_Models)
{
  int GGC_Size = Gauge_Group_Counts.size();
  for(int a=0; a<GGC_Size; a++)
    {
      Latex_Out<<"\\begin{tabular}{||c|c|c||}"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;
      Latex_Out<<"Gauge Group & Number of Unique Models & \\% of Unique Models \\\\";
      Latex_Out<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;

      int Groups_in_Table = a;
      std::map<GaugeGroupName, int>::const_iterator itGGC = 
	Gauge_Group_Counts.begin();
      for(int b=0; b<a; b++)
	++itGGC;

      int Table_Loop_Limit = 0;
      if((GGC_Size - Groups_in_Table) < 20)
	Table_Loop_Limit = GGC_Size - Groups_in_Table;
      else
	Table_Loop_Limit = 20;

      for(int b=0; b<Table_Loop_Limit; b++)
	{
		Latex_Out<<"$";
		Latex_Out<<Cast_GaugeGroupName(itGGC->first);
		Latex_Out<<"$";
		Latex_Out<<"&";
		Latex_Out<<(itGGC->second)<<"&";
		double Percent = (double((itGGC->second))/double(Total_Unique_Models))*
			100;
		Latex_Out<<std::setprecision(4)<<Percent<<"\\%";
		Latex_Out<<"\\\\"<<std::endl;
		Latex_Out<<"\\hline"<<std::endl;

		//Increment.
		if(b != Table_Loop_Limit - 1)
			{
				a++;
				++itGGC;
			}
	}//Close for loop on Table_Loop_Limit.
      Latex_Out<<"\\end{tabular}"<<std::endl;
      Latex_Out<<"\\\\"<<std::endl;
    }//Close for loop on Gauge Group Factor Counts.
}//Close Write_Gauge_Group_Factor_Count.

void LatexWriter::Write_Gauge_Group_Product_Count(std::ofstream& Latex_Out,
							 const std::map<std::vector
							 <GaugeGroupName>, int>&
							 Gauge_Group_Products, 
							 int Total_Unique_Models)
{
 int GGP_Size = Gauge_Group_Products.size();
  for(int a=0; a<GGP_Size; a++)
    {
      Latex_Out<<"\\begin{tabular}{||c|c|c||}"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;
      Latex_Out<<"Gauge Group Product & ";
      Latex_Out<<"Number of Unique Models & \\% of Unique Models \\\\";
      Latex_Out<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;

      int Groups_in_Table = a;
      std::map<std::vector<GaugeGroupName>, int>::const_iterator itGGP = 
	Gauge_Group_Products.begin();
      for(int b=0; b<a; b++)
	++itGGP;

      int Table_Loop_Limit = 0;
      if((GGP_Size - Groups_in_Table) < 30)
	Table_Loop_Limit = GGP_Size - Groups_in_Table;
      else
	Table_Loop_Limit = 30;

      for(int b=0; b<Table_Loop_Limit; b++)
	{
		int GGP = (itGGP->first).size();
		Latex_Out<<"$";
		for(int c=0; c<GGP; ++c)
			{
				Latex_Out<<Cast_GaugeGroupName((itGGP->first).at(c));
				if(c != GGP-1)
		Latex_Out<<"\\otimes ";
			}
		Latex_Out<<"$";
		Latex_Out<<"&";
		Latex_Out<<(itGGP->second)<<"&";
		double Percent = (double((itGGP->second))/double(Total_Unique_Models))*
			100;
		Latex_Out<<std::setprecision(4)<<Percent<<"\\%";
		Latex_Out<<"\\\\"<<std::endl;
		Latex_Out<<"\\hline"<<std::endl;

		//Increment.
		if(b != Table_Loop_Limit - 1)
			{
				a++;
				++itGGP;
			}
	}//Close for loop on Table_Loop_Limit.
      Latex_Out<<"\\end{tabular}"<<std::endl;
      Latex_Out<<"\\\\"<<std::endl;
    }//Close for loop on Gauge Group Factor Counts.
}

void LatexWriter::Write_ST_SUSY_Count(std::ofstream& Latex_Out, 
							 const std::map<int, int>& ST_SUSY_Counts,
							 int Total_Unique_Models)
{
  Latex_Out<<"\\begin{tabular}{||c|c|c||}"<<std::endl;
  Latex_Out<<"\\hline"<<std::endl;
  Latex_Out<<"N&Number of Unique Models& \\% of Unique Models\\\\"<<std::endl;
  Latex_Out<<"\\hline"<<std::endl;
  std::map<int, int>::const_iterator itST_SUSY_Counts = ST_SUSY_Counts.begin();
  std::map<int, int>::const_iterator itSSCEnd = ST_SUSY_Counts.end();
  for(; itST_SUSY_Counts != itSSCEnd; ++itST_SUSY_Counts)
    {
      double Percent = (double(itST_SUSY_Counts->second)/
			double(Total_Unique_Models))*100;
      Latex_Out<<(itST_SUSY_Counts->first)<<"&";
      Latex_Out<<(itST_SUSY_Counts->second)<<"&";
      Latex_Out<<std::setprecision(4)<<Percent<<"\\%\\\\"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;
    }
  Latex_Out<<"\\end{tabular}"<<std::endl;
}//Close Write_ST_SUSY_Counts.

void LatexWriter::Plot_ST_SUSY_Count(std::ofstream& Latex_Out,
							const std::map<int, int>& ST_SUSY_Counts,
							int Large_ST_Dimensions)
{
  Latex_Out<<"\\begin{tikzpicture}"<<std::endl;
  Latex_Out<<"\\begin{axis}";
  Latex_Out<<" [ybar, ylabel = Number of Distinct Models,";
  Latex_Out<<" xlabel = Number of ST SUSYs]"<<std::endl;
  Latex_Out<<"\\addplot[draw=black, fill=black]";
  Latex_Out<<"coordinates{"<<std::endl;
  std::map<int, int>::const_iterator itST_SUSY_Counts = ST_SUSY_Counts.begin();
  int Loop_End = 0;
  switch(Large_ST_Dimensions)
    {
    case 10:
      Loop_End = 1;
      break;
    case 8:
      Loop_End = 1;
      break;
    case 6: 
      Loop_End = 2;
      break;
    case 4:
      Loop_End = 4;
      break;
    default:
      Loop_End = 4;
    };

  for(int a=0; a<=Loop_End; a++)
    {
      if(a != 3)
	{
		itST_SUSY_Counts = ST_SUSY_Counts.find(a);
		if(itST_SUSY_Counts != ST_SUSY_Counts.end())
			{
				Latex_Out<<"("<<a<<",";
				Latex_Out<<(itST_SUSY_Counts->second)<<") ";
			}else
			{
				Latex_Out<<"("<<a<<", ";
				Latex_Out<<0<<") ";
			}//Close if/else
	}//Close if(a != 3)
    }//Close for loop.
  Latex_Out<<"};"<<std::endl;
  Latex_Out<<"\\end{axis}"<<std::endl;
  Latex_Out<<"\\end{tikzpicture}"<<std::endl;
}//Close Plot_ST_SUSY_Count.

void LatexWriter::Write_U1_Factor_Count(std::ofstream& Latex_Out,
					 const std::map<int, int>& U1_Factor_Counts,
					 int Total_Unique_Models)
{
  Latex_Out<<"\\begin{tabular}{||c|c|c||}"<<std::endl;
  Latex_Out<<"\\hline"<<std::endl;
  Latex_Out<<"$U(1)$'s&Number of Unique Models& \\% of Unique Models\\\\"<<std::endl;
  Latex_Out<<"\\hline"<<std::endl;
  std::map<int, int>::const_iterator itU1_Factor_Counts = U1_Factor_Counts.begin();
  std::map<int, int>::const_iterator itU1_End = U1_Factor_Counts.end();
  for(; itU1_Factor_Counts != itU1_End; ++itU1_Factor_Counts)
    {
      double Percent = (double(itU1_Factor_Counts->second)/
			double(Total_Unique_Models))*100;
      Latex_Out<<(itU1_Factor_Counts->first)<<"&";
      Latex_Out<<(itU1_Factor_Counts->second)<<"&";
      Latex_Out<<std::setprecision(4)<<Percent<<"\\%\\\\"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;
    }
  Latex_Out<<"\\end{tabular}"<<std::endl;
}//Close Write_U1_Factors_Counts.

void LatexWriter::Plot_U1_Factor_Count(std::ofstream& Latex_Out,
					 const std::map<int, int>& U1_Factor_Counts)
{
  Latex_Out<<"\\begin{tikzpicture}"<<std::endl;
  Latex_Out<<"\\begin{axis}";
  Latex_Out<<" [ybar, ylabel = Number of Distinct Models,";
  Latex_Out<<" xlabel = Number of $U(1)$ Factors]"<<std::endl;
  Latex_Out<<"\\addplot[draw=black, fill=black]";
  Latex_Out<<"coordinates{"<<std::endl;
  std::map<int, int>::const_iterator itU1_Factor_Counts = U1_Factor_Counts.begin();
  std::map<int, int>::const_iterator itU1_End = U1_Factor_Counts.end();

  for(; itU1_Factor_Counts != itU1_End; ++itU1_Factor_Counts)
    {

      Latex_Out<<"("<<(itU1_Factor_Counts->first)<<",";
      Latex_Out<<(itU1_Factor_Counts->second)<<") ";
    }//Close for loop.
  Latex_Out<<"};"<<std::endl;
  Latex_Out<<"\\end{axis}"<<std::endl;
  Latex_Out<<"\\end{tikzpicture}"<<std::endl;
}//Close Plot_U1_Factor_Counts.

void LatexWriter::Write_Gauge_Group_Factor_Count(std::ofstream& Latex_Out,
							const std::map<int, int>& 
							Gauge_Group_Factors, 
							int Total_Unique_Models)
{
Latex_Out<<"\\begin{tabular}{||c|c|c||}"<<std::endl;
  Latex_Out<<"\\hline"<<std::endl;
  Latex_Out<<"$f$";
  Latex_Out<<"&Number of Unique Models& \\% of Unique Models\\\\"<<std::endl;
  Latex_Out<<"\\hline"<<std::endl;
  std::map<int, int>::const_iterator itGauge_Group_Factors = 
    Gauge_Group_Factors.begin();
  std::map<int, int>::const_iterator itGGF_End = Gauge_Group_Factors.end();
  for(; itGauge_Group_Factors != itGGF_End; 
      ++itGauge_Group_Factors)
    {
      double Percent = (double(itGauge_Group_Factors->second)/
			double(Total_Unique_Models))*100;
      Latex_Out<<(itGauge_Group_Factors->first)<<"&";
      Latex_Out<<(itGauge_Group_Factors->second)<<"&";
      Latex_Out<<std::setprecision(4)<<Percent<<"\\%\\\\"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;
    }
  Latex_Out<<"\\end{tabular}"<<std::endl;
}//Close Write_Gauge_Group_Factor_Count.

void LatexWriter::Plot_Gauge_Group_Factor_Count(std::ofstream& Latex_Out,
						 const std::map<int, int>&
						 Gauge_Group_Factors)
{
  Latex_Out<<"\\begin{tikzpicture}"<<std::endl;
  Latex_Out<<"\\begin{axis}";
  Latex_Out<<" [ybar, ylabel = Number of Distinct Models,";
  Latex_Out<<" xlabel = Number of Gauge Group Factors]"<<std::endl;
  Latex_Out<<"\\addplot[draw=black, fill=black]";
  Latex_Out<<"coordinates{"<<std::endl;
  std::map<int, int>::const_iterator itGauge_Group_Factors = 
    Gauge_Group_Factors.begin();
  std::map<int, int>::const_iterator itGGF_End = Gauge_Group_Factors.end();

  for(; itGauge_Group_Factors != itGGF_End; 
      ++itGauge_Group_Factors)
    {

      Latex_Out<<"("<<(itGauge_Group_Factors->first)<<",";
      Latex_Out<<(itGauge_Group_Factors->second)<<") ";
    }//Close for loop.
  Latex_Out<<"};"<<std::endl;
  Latex_Out<<"\\end{axis}"<<std::endl;
  Latex_Out<<"\\end{tikzpicture}"<<std::endl;
}//Close Plot_Gauge_Group_Factor_Count.

void LatexWriter::Write_Matter_Generations(std::ofstream& Latex_Out,
		const std::map<int, int>& Generations)
{
	if(Generations.size() != 0)
	{
		Latex_Out<<"\\tablehead{\\hline Number of Chiral Matter Generations";
		Latex_Out<<" & Number of Unique Observable Sectors\\\\}"<<std::endl;
		Latex_Out<<"\\tabletail{\\hline}"<<std::endl;
		Latex_Out<<"\\tablelasttail{\\hline}"<<std::endl;
		Latex_Out<<"\\begin{supertabular}{||c|c||}"<<std::endl;
		Latex_Out<<"\\hline"<<std::endl;

		std::map<int, int>::const_iterator itGenerations = Generations.begin();
		for(; itGenerations != Generations.end(); 
				++itGenerations)
		{
			Latex_Out<<itGenerations->first<<"&";
			Latex_Out<<itGenerations->second<<"\\\\"<<std::endl;
			Latex_Out<<"\\hline"<<std::endl;
		}//Close for loop on Generations.
		Latex_Out<<"\\end{supertabular}"<<std::endl;
	}
		else
			Latex_Out<<"No models with chiral matter reps."<<std::endl;
}//Close Write_Matter_Generations.

void LatexWriter::Plot_Matter_Generations(std::ofstream& Latex_Out,
		const std::map<int, int>& Generations)
{
	if(Generations.size() > 0)
	{
		Latex_Out<<"\\begin{tikzpicture}"<<std::endl;
		Latex_Out<<"\\begin{axis}";
		Latex_Out<<" [ybar, ylabel = Number of Distinct Observable Sectors, ";
		Latex_Out<<" xlabel = Number of Chiral Matter Generations]"<<std::endl;
		Latex_Out<<"\\addplot[draw=black, fill=black]";
		Latex_Out<<"coordinates{"<<std::endl;
		std::map<int, int>::const_iterator itGenerations = Generations.begin();
		for(; itGenerations != Generations.end(); ++itGenerations)
		{
			Latex_Out<<"("<<(itGenerations->first)<<",";
			Latex_Out<<(itGenerations->second)<<")";
		}//Close for loop on Generations.
		Latex_Out<<"};"<<std::endl;
		Latex_Out<<"\\end{axis}"<<std::endl;
		Latex_Out<<"\\end{tikzpicture}"<<std::endl;
	}else
		Latex_Out<<"No models with chiral matter reps."<<std::endl;
}//Close Plot_matter_Generations.

void LatexWriter::Write_NA_Singlets(std::ofstream& Latex_Out,
		const std::map<int, int>& NA_Singlets)
{
	if(NA_Singlets.size() != 0)
	{
		Latex_Out<<"\\tablehead{\\hline Number of NA Singlets";
		Latex_Out<<"& Number of Unique Models\\\\}"<<std::endl;
		Latex_Out<<"\\tabletail{\\hline}"<<std::endl;
		Latex_Out<<"\\tablelasttail{\\hline}"<<std::endl;
		Latex_Out<<"\\begin{supertabular}{||c|c||}"<<std::endl;
		Latex_Out<<"\\hline"<<std::endl;

		std::map<int, int>::const_iterator itNA_Singlets = 
			NA_Singlets.begin();
		for(; itNA_Singlets != NA_Singlets.end();
				++itNA_Singlets)
		{
			Latex_Out<<itNA_Singlets->first<<"&";
			Latex_Out<<itNA_Singlets->second<<"\\\\"<<std::endl;
			Latex_Out<<"\\hline"<<std::endl;
		}//Close for loop on NA_Singlets.
		Latex_Out<<"\\end{supertabular}"<<std::endl;
	}else
		Latex_Out<<"No models with NA singlets."<<std::endl;
}//Close Write_NA_Singlets.

void LatexWriter::Plot_NA_Singlets(std::ofstream& Latex_Out,
		const std::map<int, int>& NA_Singlets)
{
	if(NA_Singlets.size() != 0)
	{
		Latex_Out<<"\\begin{tikzpicture}"<<std::endl;
		Latex_Out<<"\\begin{axis}";
		Latex_Out<<" [ybar, ylabel = Number of Distinct Models, ";
		Latex_Out<<" xlabel = Number of NA Singlets]"<<std::endl;
		Latex_Out<<"\\addplot[draw=black, fill=black]";
		Latex_Out<<"coordinates{"<<std::endl;
		std::map<int, int>::const_iterator itNA_Singlets = NA_Singlets.begin();
		for(; itNA_Singlets != NA_Singlets.end(); 
				++itNA_Singlets)
		{
			Latex_Out<<"("<<(itNA_Singlets->first)<<",";
			Latex_Out<<(itNA_Singlets->second)<<")";
		}//Close for loop on NA_Singlets.
		Latex_Out<<"};"<<std::endl;
		Latex_Out<<"\\end{axis}"<<std::endl;
		Latex_Out<<"\\end{tikzpicture}"<<std::endl;
	}else
		Latex_Out<<"No models with NA singlets."<<std::endl;
}//Close Plot_NA_Singlets.

void LatexWriter::Write_OS_Charged_Exotics(std::ofstream& Latex_Out,
		const std::map<int, int>& OS_Charged_Exotics)
{
	if(OS_Charged_Exotics.size() != 0)
	{
		Latex_Out<<"\\begin{longtable}{||c|c||}"<<std::endl;
		Latex_Out<<"\\hline"<<std::endl;
		Latex_Out<<"Number of charged exotics";
		Latex_Out<<"& Number of unique models\\\\"<<std::endl;
		Latex_Out<<"\\hline"<<std::endl;
		Latex_Out<<"\\endhead"<<std::endl;
		Latex_Out<<"\\hline"<<"\\hline"<<std::endl;
		Latex_Out<<"\\endfoot"<<std::endl;

		std::map<int, int>::const_iterator itOS_Charged_Exotics = 
			OS_Charged_Exotics.begin();
		for(; itOS_Charged_Exotics !=
				OS_Charged_Exotics.end(); ++itOS_Charged_Exotics)
		{
			Latex_Out<<itOS_Charged_Exotics->first<<"&";
			Latex_Out<<itOS_Charged_Exotics->second<<"\\\\";
			Latex_Out<<std::endl;
		}//Close for loop on OS_Charged_Exotics.
		Latex_Out<<"\\end{longtable}"<<std::endl;
	}
		else
			Latex_Out<<"No charged exotics."<<std::endl;
}//Close Write_OS_Charged_Exotics.

void LatexWriter::Plot_OS_Charged_Exotics(std::ofstream& Latex_Out,
		const std::map<int, int>& OS_Charged_Exotics)
{

	if(OS_Charged_Exotics.size() != 0)
	{
		Latex_Out<<"\\begin{tikzpicture}"<<std::endl;
		Latex_Out<<"\\begin{axis}";
		Latex_Out<<" [ybar, ylabel = Number of Distinct Models, ";
		Latex_Out<<" xlabel = Number of Charged Exotics]"<<std::endl;
		Latex_Out<<"\\addplot[draw=black, fill=black]";
		Latex_Out<<"coordinates{"<<std::endl;
		std::map<int, int>::const_iterator itOS_Charged_Exotics = 
			OS_Charged_Exotics.begin();
		for(; itOS_Charged_Exotics != 
				OS_Charged_Exotics.end(); ++itOS_Charged_Exotics)
		{
			Latex_Out<<"("<<(itOS_Charged_Exotics->first)<<",";
			Latex_Out<<(itOS_Charged_Exotics->second)<<")";
		}//Close for loop on OS_Charged_Exotics.
		Latex_Out<<"};"<<std::endl;
		Latex_Out<<"\\end{axis}"<<std::endl;
		Latex_Out<<"\\end{tikzpicture}"<<std::endl;
	}else
		Latex_Out<<"No models with charged exotics."<<std::endl;
}//Close Plot_OS_Charged_Exotics.

//PRIVATE.
//HELPERS.
std::vector<std::string> LatexWriter::Get_BV_Names(const Model& FFHS_Model)
{
  std::vector<std::string> Names;
  int Extension_Start = 0;
  if(FFHS_Model.NAHE_Loaded())
    {
      Names.push_back("1");
      if(Has_S_Vector(FFHS_Model))//This class assumes the second vector is the
	//S vector, so load it accordingly when building the model.
	{
		Names.push_back("S");
		Extension_Start = 5;
	}else
	Extension_Start = 4;
      Names.push_back("$b_1$");
      Names.push_back("$b_2$");
      Names.push_back("$b_3$");
      
      for(int a=Extension_Start; a<static_cast<int>(FFHS_Model.BV_Set().size()); a++)
	{
		std::ostringstream ss;
		ss<<"$\\alpha_{"<<a-4<<"}$";
		std::string BV_Name = ss.str();
		Names.push_back(BV_Name);
	}//Close for loop on the rest of the BVs.
    }else
    {
      Names.push_back("1");
      if(Has_S_Vector(FFHS_Model))//This class assumes the second vector is the 
	//S vector, so load it accordingly when building the model.
	{
		Names.push_back("S");
		Extension_Start = 2;
	}else
	Extension_Start = 1;
      for(int a=Extension_Start; a<static_cast<int>(FFHS_Model.BV_Set().size()); a++)
	{
		std::ostringstream ss;
		ss<<"\\alpha_{"<<a-1<<"}$";
		std::string BV_Name = ss.str();
		Names.push_back(BV_Name);
	}//Close for loop on the rest of the BVs.
    }//Close if/else on NAHE_Loaded.
  return Names;
}//Close Get_BV_Names.

std::vector<std::string> LatexWriter::Get_GaugeGroupNames(const Model& FFHS_Model)
{
  std::vector<std::string> GaugeGroupNames;
  for(int a=0; a<static_cast<int>(FFHS_Model.Gauge_Groups().size()); a++)
		GaugeGroupNames.push_back(Cast_GaugeGroupName
				(FFHS_Model.Gauge_Groups().at(a).Name()));

  return GaugeGroupNames;
}//Close Get_GaugeGroupNames.

void LatexWriter::Write_Left_Movers(std::vector<std::string> BV_Names,
						 std::ofstream& Latex_Out, 
						 const Model& FFHS_Model)
{
  int ST_LM = 12 - (FFHS_Model.BV_Set().at(0).LM_Size()/2);
  int Real_LM = (3*FFHS_Model.BV_Set().at(0).LM_Size()/2) - 12;

  if(Real_LM <= 12)
    {
      //ST left movers first.
      Latex_Out<<"\\begin{tabular}{||c|c||";
      for(int a=0; a<ST_LM; a++)
		Latex_Out<<"c|";
      Latex_Out<<"|}"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;
      Latex_Out<<"\\textbf{N}&\\textbf{O}";
      for(int a=0; a<ST_LM; a+=2)
	{
		Latex_Out<<"&$\\psi^"<<((a/2)+1)<<"$";
		Latex_Out<<"&$\\psi^{"<<((a/2)+1)<<"*}$";
	}
      Latex_Out<<"\\\\"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;
      for(int a=0; a<static_cast<int>(FFHS_Model.BV_Set().size()); a++)
	{
		Latex_Out<<BV_Names.at(a)<<"&"<<FFHS_Model.BV_Set().at(a).Order();
		for(int b=0; b<ST_LM; b++)
			Latex_Out<<"&"<<FFHS_Model.BV_Set().at(a).BV().at(b);
		Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
	}
      Latex_Out<<"\\hline"<<std::endl<<"\\end{tabular}"<<std::endl;
      Latex_Out<<"\\\\"<<std::endl<<"\\vspace{3 mm}"<<std::endl;
      //Now for the compactified LMs.
      if(Real_LM != 0)
	{
		Latex_Out<<"\\begin{tabular}{||c|c||";
		for(int a=0; a<Real_LM; a+=3)
			Latex_Out<<"c|c|c||";
		Latex_Out<<"}"<<std::endl;
		Latex_Out<<"\\hline"<<std::endl<<"\\hline"<<std::endl;
		Latex_Out<<"\\textbf{N}&\\textbf{O}";
		for(int a=0; a<Real_LM; a+=3)
			{
				Latex_Out<<"&$x^"<<((a/3)+1)<<"$";
				Latex_Out<<"&$y^"<<((a/3)+1)<<"$";
				Latex_Out<<"&$w^"<<((a/3)+1)<<"$";
			}
		Latex_Out<<"\\\\"<<std::endl;
		Latex_Out<<"\\hline"<<std::endl;
		for(int a=0; a<static_cast<int>(FFHS_Model.BV_Set().size()); a++)
			{
				Latex_Out<<BV_Names.at(a)<<"&"<<FFHS_Model.BV_Set().at(a).Order();
				for(int b=ST_LM; b<Real_LM+ST_LM; b++)
		Latex_Out<<"&"<<FFHS_Model.BV_Set().at(a).BV().at(b);
				Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
			}
		Latex_Out<<"\\hline"<<std::endl<<"\\end{tabular}"<<std::endl;
		Latex_Out<<"\\\\"<<std::endl<<"\\vspace{3 mm}"<<std::endl;
	}//Close if statement on Real_LM != 0.
    }else
    {
      Latex_Out<<"\\begin{tabular}{||c|c||";
      for(int a=0; a<ST_LM; a++)//ST
	Latex_Out<<"c|";
      Latex_Out<<"|";
      for(int a=0; a<6; a+=3)//First two compactified triplets.
	Latex_Out<<"c|c|c||";
      Latex_Out<<"}"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl<<"\\hline"<<std::endl;
      Latex_Out<<"\\textbf{N}&\\textbf{O}";
      for(int a=0; a<ST_LM; a+=2)//ST
	{
		Latex_Out<<"&$\\psi^"<<(a/2)+1<<"$";
		Latex_Out<<"&$\\psi^{"<<(a/2)+1<<"*}$";
	}
      for(int a=0; a<6; a+=3)//x,y,w 1,2
	{
		Latex_Out<<"&$x^"<<(a/3)+1<<"$"<<std::endl;
		Latex_Out<<"&$y^"<<(a/3)+1<<"$"<<std::endl;
		Latex_Out<<"&$w^"<<(a/3)+1<<"$"<<std::endl;
	}
      Latex_Out<<"\\\\"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;
      for(int a=0; a<static_cast<int>(FFHS_Model.BV_Set().size()); a++)
	{
		Latex_Out<<BV_Names.at(a)<<"&"<<FFHS_Model.BV_Set().at(a).Order();
		for(int b=0; b<ST_LM+6; b++)
			Latex_Out<<"&"<<FFHS_Model.BV_Set().at(a).BV().at(b);
		Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
	}
      Latex_Out<<"\\hline"<<std::endl<<"\\end{tabular}"<<std::endl;
      Latex_Out<<"\\\\"<<std::endl<<"\\vspace{3 mm}"<<std::endl;
      //Next do the rest of the compactified LM.
      Latex_Out<<"\\begin{tabular}{||c|c||";
      for(int a=6; a<Real_LM; a+=3)
	{
		Latex_Out<<"c|c|c||";
	}
      Latex_Out<<"}"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl<<"\\hline"<<std::endl;
      Latex_Out<<"\\textbf{N}&\\textbf{O}";
      for(int a=6; a<Real_LM; a+=3)
	{
		Latex_Out<<"&$x^"<<((a/3)+1)<<"$";
		Latex_Out<<"&$y^"<<((a/3)+1)<<"$";
		Latex_Out<<"&$w^"<<((a/3)+1)<<"$";
	}
      Latex_Out<<"\\\\"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl;
      for(int a=0; a<static_cast<int>(FFHS_Model.BV_Set().size()); a++)
	{
		Latex_Out<<BV_Names.at(a)<<"&"<<FFHS_Model.BV_Set().at(a).Order();
		for(int b=ST_LM+6; b<Real_LM+ST_LM; b++)
			Latex_Out<<"&"<<FFHS_Model.BV_Set().at(a).BV().at(b);
		Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
	}
      Latex_Out<<"\\hline"<<std::endl<<"\\end{tabular}"<<std::endl;
      Latex_Out<<"\\\\"<<std::endl<<"\\vspace{3 mm}"<<std::endl;
    }//Close if/else on splitting left movers.
}//Close Write_Left_Movers.

void LatexWriter::Write_Right_Movers(std::vector<std::string> BV_Names, 
							std::ofstream& Latex_Out, 
							const Model& FFHS_Model)
{
  int LM_Size = FFHS_Model.BV_Set().at(0).LM_Size();
  int Compact_Size = LM_Size - 8;
  //Do the observable elements first.
  Latex_Out<<"\\begin{tabular}{||c|c||";
  for(int a=0; a<10; a++)
    Latex_Out<<"c|";
  Latex_Out<<"|";
  for(int a=0; a<6; a++)
    Latex_Out<<"c|";
  Latex_Out<<"|}"<<std::endl;
  Latex_Out<<"\\hline"<<std::endl<<"\\hline"<<std::endl;
  Latex_Out<<"\\textbf{N}&\\textbf{O}";
  for(int a=0; a<10; a+=2)
    {
      Latex_Out<<"&$\\overline{\\psi}^{"<<(a/2)+1<<"}$";
      Latex_Out<<"&$\\overline{\\psi}^{"<<(a/2)+1<<"*}$";
    }
  for(int a=0; a<6; a+=2)
    {
      Latex_Out<<"&$\\overline{\\eta}^{"<<(a/2)+1<<"}$";
      Latex_Out<<"&$\\overline{\\eta}^{"<<(a/2)+1<<"*}$";
    }
  Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
  for(int a=0; a<static_cast<int>(FFHS_Model.BV_Set().size()); a++)
    {
      Latex_Out<<BV_Names.at(a)<<"&"<<FFHS_Model.BV_Set().at(a).Order();
      for(int b=LM_Size; b<LM_Size+16; b++)
	Latex_Out<<"&"<<FFHS_Model.BV_Set().at(a).BV().at(b);
      Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
    }
  Latex_Out<<"\\end{tabular}"<<std::endl<<"\\\\"<<std::endl;
  Latex_Out<<"\\vspace{3 mm}"<<std::endl;

  //Now for the compact parts.
  if(Compact_Size != 0)
    {
      Latex_Out<<"\\begin{tabular}{||c|c||";
      for(int a=0; a<Compact_Size/2; a++)
	Latex_Out<<"c|";
      Latex_Out<<"|";
      for(int a=0; a<Compact_Size/2; a++)
	Latex_Out<<"c|";
      Latex_Out<<"|}"<<std::endl;
      Latex_Out<<"\\hline"<<std::endl<<"\\hline"<<std::endl;
      Latex_Out<<"\\textbf{N}&\\textbf{O}";
      for(int a=0; a<Compact_Size/2; a++)
	Latex_Out<<"&$\\overline{y}^"<<a+1<<"$";
      for(int a=0; a<Compact_Size/2; a++)
	Latex_Out<<"&$\\overline{w}^"<<a+1<<"$";
      Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
      for(int a=0; a<static_cast<int>(FFHS_Model.BV_Set().size());a++)
	{
		Latex_Out<<BV_Names.at(a)<<"&"<<FFHS_Model.BV_Set().at(a).Order();
		for(int b=LM_Size+16; b<LM_Size+16+Compact_Size; b++)
			Latex_Out<<"&"<<FFHS_Model.BV_Set().at(a).BV().at(b);
		Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
	}
      Latex_Out<<"\\end{tabular}"<<std::endl<<"\\\\"<<std::endl;
      Latex_Out<<"\\vspace{3 mm}"<<std::endl;
    }//Close if statement on Compact_Size.

  //Finally, the hidden elements.
  Latex_Out<<"\\begin{tabular}{||c|c||";
  for(int a=0; a<16; a++)
    Latex_Out<<"c|";
  Latex_Out<<"|}"<<std::endl;
  Latex_Out<<"\\hline"<<std::endl<<"\\hline"<<std::endl;
  Latex_Out<<"\\textbf{N}&\\textbf{O}";
  for(int a=0; a<16; a+=2)
    {
      Latex_Out<<"&$\\overline{\\phi}^{"<<(a/2)+1<<"}$";
      Latex_Out<<"&$\\overline{\\phi}^{"<<(a/2)+1<<"*}$";
    }
  Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
  for(int a=0; a<static_cast<int>(FFHS_Model.BV_Set().size()); a++)
    {
      Latex_Out<<BV_Names.at(a)<<"&"<<FFHS_Model.BV_Set().at(a).Order();
      for(int b=LM_Size+16+Compact_Size; 
		b<static_cast<int>(FFHS_Model.BV_Set().at(a).BV().size()); b++)
	Latex_Out<<"&"<<FFHS_Model.BV_Set().at(a).BV().at(b);
      Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
    }
  Latex_Out<<"\\end{tabular}"<<std::endl<<"\\\\"<<std::endl;
  Latex_Out<<"\\vspace{3 mm}"<<std::endl;
}//Close Write_Right_Movers.

bool LatexWriter::Has_S_Vector(const Model& FFHS_Model)
{
  //This function assumes the S vector is placed in the spot after the
  //1 vector, so load accordingly when building the model.
  std::vector<int> S_Vector;
  int Large_ST_Dimensions = 14 - (FFHS_Model.BV_Set().at(0).LM_Size()/2);
  for(int a=0; a<Large_ST_Dimensions-2; a++)
    S_Vector.push_back(1);
  for(int a=0; a<(10-Large_ST_Dimensions); a++)
    {
      S_Vector.push_back(1);
      S_Vector.push_back(0);
      S_Vector.push_back(0);
    }//Close for loop on building S vector.
  for(int a=FFHS_Model.BV_Set().at(0).LM_Size(); 
      a<static_cast<int>(FFHS_Model.BV_Set().at(0).BV().size()); a++)
    S_Vector.push_back(0);

  if(FFHS_Model.BV_Set().size() > 1)
    return (FFHS_Model.BV_Set().at(1).BV() == S_Vector);
  else
    return false;
}//Close Has_S_Vector.

std::string LatexWriter::Cast_GaugeGroupName(const GaugeGroupName& Name)
{
  std::stringstream ss;
  switch(Name.Class())
    {
    case 'A':
      ss<<"SU("<<Name.Rank()+1<<")";
      break;
    case 'B':
      ss<<"SO("<<(2*Name.Rank())+1<<")";
      break;
    case 'C':
      ss<<"Sp("<<2*Name.Rank()<<")";
      break;
    case 'D':
      ss<<"SO("<<2*Name.Rank()<<")";
      break;
    case 'E':
      ss<<"E_"<<Name.Rank();
      break;
    case 'F':
      ss<<"F_"<<Name.Rank();
      break;
    case 'G':
      ss<<"G_"<<Name.Rank();
      break;
    case 'U':
      ss<<"U(1)";
    };
  if(Name.KM_Level()>1)
    ss<<"^{("<<Name.KM_Level()<<")}";
  return ss.str();
}//Close Cast_GaugeGroupName.

void LatexWriter::Write_Matter_Reps(std::ofstream& Latex_Out, 
		const std::vector<std::string>& GaugeGroupNames, 
		const std::vector<MatterRepresentation>& MatterRepresentations)
{
	if(MatterRepresentations.size() == 0)
	{
		Latex_Out<<"\\textbf{Gauge Groups: ";
		for(int a=0; a<static_cast<int>(GaugeGroupNames.size()); a++)
			Latex_Out<<"$"<<GaugeGroupNames.at(a)<<"$"<<" ";
		Latex_Out<<"}\\\\"<<std::endl;
		Latex_Out<<"\\textbf{No matter representations.}\\\\"<<std::endl;
	}else
	{
		Latex_Out<<"\\begin{longtable}{||c|";
		for(int a=0; a<static_cast<int>(GaugeGroupNames.size()); a++)
			Latex_Out<<"|c";
		Latex_Out<<"||}"<<std::endl;
		Latex_Out<<"\\hline"<<std::endl<<"\\hline"<<std::endl;
		Latex_Out<<"\\textbf{QTY}";
		for(int a=0; a<static_cast<int>(GaugeGroupNames.size()); a++)
			Latex_Out<<"&"<<"$"<<GaugeGroupNames.at(a)<<"$";
		Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
		Latex_Out<<"\\hline"<<std::endl<<"\\endhead"<<std::endl;
		Latex_Out<<"\\hline"<<std::endl<<"\\endfoot"<<std::endl;
		for(int a=0; a<static_cast<int>(MatterRepresentations.size()); a++)
		{
			Latex_Out<<MatterRepresentations.at(a).Duplicates();
			for(int b=0; b<static_cast<int>(MatterRepresentations.at(a).Rep_Dimension().size()); ++b)
			{
				int MatterRepresentation_Dimension = 
					MatterRepresentations.at(a).Rep_Dimension().at(b).Dimension();
				char MatterRepresentation_Triality = 
					MatterRepresentations.at(a).Rep_Dimension().at(b).Triality();
				if(MatterRepresentation_Dimension<0)
				{
					Latex_Out<<"&$\\overline{"<<abs(MatterRepresentation_Dimension);
					Latex_Out<<"}$";
				}else
				{
					Latex_Out<<"&"<<"$"<<MatterRepresentation_Dimension;
					Latex_Out<<"_{"<<MatterRepresentation_Triality<<"}$";
				}
			}
			Latex_Out<<"\\\\"<<std::endl<<"\\hline"<<std::endl;
		}//Close for loop on a.
		Latex_Out<<"\\end{longtable}"<<std::endl;
	}//Close if/else on MatterRepresentations.size().
	Latex_Out<<"\\hspace{\\fill}\\\\"<<std::endl;
}//Close Write_Matter_Reps.
