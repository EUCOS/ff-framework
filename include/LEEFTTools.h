/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/LEEFTTools.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * A namespace of useful LEEFT related functions is provided.
 */
#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

/*!
 * This namespace has miscellaenous tools for reading the LEEFT information
 * from standanrd FF Framework output files.
 */
namespace LEEFTTools {
  void Read_Front_Info(std::ifstream& LEEFT_In);
  std::vector<std::string> Read_LEEFT_Data(std::ifstream& LEEFT_In);
  std::vector<std::string> Read_Model_Inputs(std::ifstream& LEEFT_In, int Layers);
}
