/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/OutputWriter.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file declares a class for efficient output handling.
 */
#pragma once

#include <fstream>
#include <sstream>
#include <vector>

#include <GaugeGroup.h>
#include <GaugeGroupName.h>
#include <GroupRepresentation.h>
#include <LEEFT.h>
#include <Model.h>

class OutputWriter {
  public:
    OutputWriter() {}
    OutputWriter(const OutputWriter& New_OutputWriter);

    ~OutputWriter() {}

    void Write_Model_Particle_Content(std::ofstream& Model_Out, const Model& FFHS_Model);
    void Write_Model_BVs(std::ofstream& Model_Out, const Model& FFHS_Model);
    void Write_Model_Extension_BVs(std::ofstream& Model_Out, const Model& FFHS_Model, int Extensions);
    void Write_Model_k_ij(std::ofstream& Model_Out, const Model& FFHS_Model);
    void Write_Model_Extension_k_ij(std::ofstream& Model_Out, const std::vector<std::vector<int> >& k_ij_Extensions);
    void Write_Model_BV_Orders(std::ofstream& Model_Out, const Model& FFHS_Model);
    void Write_BV_Search_Front_Info(std::ofstream& Model_Out, const Model& FFHS_Model, int Large_ST_Dimensions, const std::vector<int>& Extension_Orders);
    void Write_LEEFT_Particle_Content(std::ofstream& LEEFT_Out, const LEEFT& FFHS_LEEFT);
};
