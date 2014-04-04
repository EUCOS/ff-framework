/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/LatexWriter.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the LatexWrite class which is designed
 * to facilitate writing models out to LaTeX source files.
 */
#pragma once

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include <LEEFT.h>
#include <MatterRepresentation.h>
#include <Model.h>

/*!
 * This class writes a model in LaTeX markup language for easy examination.
 * There are interfaces for both the Model and LEEFT classes.
 *
 * Bash scripts are available with the FF Framework for easy command line
 * rendering of pdf documents.
 */
class LatexWriter {
  public:
    LatexWriter() {}
    LatexWriter(const LatexWriter& New_LatexWriter);

    ~LatexWriter() {}

    void Begin_Latex_Document(std::ofstream& Latex_Out, std::string Doc_Title);
    void Finish_Latex_Document(std::ofstream& Latex_Out);

    void Print_Model_BVs(std::ofstream& Latex_Out, const Model& FFHS_Model);
    void Print_Model_k_ij(std::vector<std::string> BV_Names, std::ofstream& Latex_Out, const Model& FFHS_Model);
    void Print_Model_Inputs(std::ofstream& Latex_Out, const Model& FFHS_Model);
    void Print_Model_Particle_Content(std::ofstream& Latex_Out, const Model& FFHS_Model);
    void Print_Full_Model(std::ofstream& Latex_Out, const Model& FFHS_Model);
    void Print_Single_Model_Document(std::ofstream& Latex_Out, const Model& FFHS_Model);

    void Print_LEEFT_Particle_Content(std::ofstream& LEEFT_Out, const LEEFT& FFHS_LEEFT);
    void Print_Single_LEEFT_Document(std::ofstream& LEEFT_Out, const LEEFT& FFHS_LEEFT);

    void Write_Gauge_Group_Count(std::ofstream& Latex_Out, const std::map<GaugeGroupName, int>& Gauge_Group_Counts, int Total_Unique_Models);
    void Write_Gauge_Group_Product_Count(std::ofstream& Latex_Out, const std::map<std::vector<GaugeGroupName>, int>& Gauge_Group_Products, int Total_Unique_Models);
    void Write_ST_SUSY_Count(std::ofstream& Latex_Out, const std::map<int, int>& ST_SUSY_Counts, int Total_Unique_Models);
    void Plot_ST_SUSY_Count(std::ofstream& Latex_Out, const std::map<int, int>& ST_SUSY_Counts, int Large_ST_Dimensions);
    void Write_U1_Factor_Count(std::ofstream& Latex_Out, const std::map<int, int>& U1_Factor_Counts, int Total_Unique_Models);
    void Plot_U1_Factor_Count(std::ofstream& Latex_Out, const std::map<int, int>& U1_Factor_Counts);
    void Write_Gauge_Group_Factor_Count(std::ofstream& Latex_Out, const std::map<int, int>& Gauge_Group_Factors, int Total_Unique_Models);
    void Plot_Gauge_Group_Factor_Count(std::ofstream& Latex_Out, const std::map<int, int>& Gauge_Group_Factors);
    void Write_Matter_Generations(std::ofstream& Latex_Out, const std::map<int, int>& Generations);
    void Plot_Matter_Generations(std::ofstream& Latex_Out, const std::map<int, int>& Generations);
    void Write_NA_Singlets(std::ofstream& Latex_Out, const std::map<int, int>& NA_Singlets);
    void Plot_NA_Singlets(std::ofstream& Latex_Out, const std::map<int, int>& NA_Singlets);
    void Write_OS_Charged_Exotics(std::ofstream& Latex_Out, const std::map<int, int>& OS_Charged_Exotics);
    void Plot_OS_Charged_Exotics(std::ofstream& Latex_Out, const std::map<int, int>& OS_Charged_Exotics);

  private:
    std::vector<std::string> Get_BV_Names(const Model& FFHS_Model);
    std::vector<std::string> Get_GaugeGroupNames(const Model& FFHS_Model);
    void Write_Left_Movers(std::vector<std::string> BV_Names, std::ofstream& Latex_Out, const Model& FFHS_Model);
    void Write_Right_Movers(std::vector<std::string> BV_Names, std::ofstream& Latex_Out, const Model& FFHS_Model);
    bool Has_S_Vector(const Model& FFHS_Model);
    std::string Cast_GaugeGroupName(const GaugeGroupName& Name);
    void Write_Matter_Reps(std::ofstream& Latex_Out, const std::vector<std::string>& GaugeGroupNames, const std::vector<MatterRepresentation>& MatterRepresentations);
};
