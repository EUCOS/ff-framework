/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/SystematicD10Builder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Here we define a derived class of SystematicBuilder specializing in D=10
 * models.
 */
#include <SystematicBuilder.h>

class SystematicD10Builder : public SystematicBuilder {
  public:
    SystematicD10Builder(const ModelBuilder& Initial_Build, int Large_ST_Dimensions, const std::vector<int>& Extension_Orders, const std::string& Output_File_Name);
    SystematicD10Builder(const SystematicD10Builder& New_SystematicD10Builder);

    ~SystematicD10Builder() {}

    void Perform_Search();

  private:
    void Build_Extensions(int Layer, std::vector<BasisVector> Basis_Vectors);
    void Build_Basis_Vectors(int Layer, std::vector<BasisVector> Basis_Vectors);
};
