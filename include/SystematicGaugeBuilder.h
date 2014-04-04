/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/SystematicGaugeBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Here we define a derived class of SystematicBuilder specializing in gauge
 * models.
 */
#pragma once

#include <SystematicBuilder.h>

class SystematicGaugeBuilder : public SystematicBuilder {
  public:
    SystematicGaugeBuilder(const ModelBuilder& Initial_Build, int Large_ST_Dimensions, const std::vector<int>& Extension_Orders, const std::string& Output_File_Name);
    SystematicGaugeBuilder(const SystematicGaugeBuilder& New_SystematicGaugeBuilder);

    ~SystematicGaugeBuilder() {}

    void Perform_Search();

  protected:
    void Build_Extensions(int Layer, std::vector<BasisVector> Basis_Vectors);
    void Build_Basis_Vectors(int Layer, std::vector<BasisVector> Basis_Vectors);
    void Build_Models(std::vector<BasisVector> Basis_Vectors);
};
