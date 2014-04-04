/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/SystematickijBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 */
#pragma once

#include <SystematicBuilder.h>

class SystematickijBuilder : public SystematicBuilder {
  public:
    SystematickijBuilder(const ModelBuilder& Initial_Build, int Large_ST_Dimensions, const std::vector<int>& Extension_Orders, const std::string& Output_File_Name, const std::vector<BasisVector>& Basis_Vectors);
    SystematickijBuilder(const SystematickijBuilder& New_SystematickijBuilder);

    ~SystematickijBuilder() {}

    void Perform_Search();

    const std::vector<BasisVector>& Basis_Vectors() const {
      return Basis_Vectors_;
    }

    void Set_Basis_Vectors(const std::vector<BasisVector>& New_Basis_Vectors) {
      Basis_Vectors_ = New_Basis_Vectors;
    }

  private:
    std::vector<BasisVector> Basis_Vectors_;
};
