/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/MidStartSystematicBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 */
#pragma once

#include <SystematicBuilder.h>

class MidStartSystematicBuilder : public SystematicBuilder {
  public:
    MidStartSystematicBuilder(const ModelBuilder& Initial_Build, int Large_ST_Dimensions, const std::vector<int>& Extension_Orders, const std::string& Output_File_Name, const std::vector<BasisVector>& Starting_Basis_Vectors);
    MidStartSystematicBuilder(const ModelBuilder& Initial_Build, int Large_ST_Dimensions, const std::vector<int>& Extension_Orders, const std::string& Output_File_Name, const std::vector<bool>& Simply_Paired_Layers, const std::vector<BasisVector>& Starting_Basis_Vectors);
    MidStartSystematicBuilder(const MidStartSystematicBuilder& New_MidStartSystematicBuilder);

    ~MidStartSystematicBuilder() {}

    void Perform_Search();

    const std::vector<BasisVector>& Starting_Basis_Vectors() const {
      return Starting_Basis_Vectors_;
    }

    const std::vector<bool>& Layer_Built() const {
      return Layer_Built_;
    }

    void Set_Starting_Basis_Vectors(const std::vector<BasisVector>& New_Starting_Basis_Vectors) {
      Starting_Basis_Vectors_ = New_Starting_Basis_Vectors;
    }

    void Set_Layer_Build(const std::vector<bool>& New_Layer_Built) {
      Layer_Built_ = New_Layer_Built;
    }

  private:
    std::vector<BasisVector> Starting_Basis_Vectors_;
    std::vector<bool> Layer_Built_;

    void Build_Extensions(int Layer, std::vector<BasisVector> Basis_Vectors);
    void Build_Basis_Vectors(int Layer, std::vector<BasisVector> Basis_Vectors);
    std::list<Chunk>::const_iterator Set_LM_Start_Point (const std::list<Chunk>& LM, const BasisVector& Basis_Vector);
    std::list<Chunk>::const_iterator Set_Comp_Start_Point (const std::list<Chunk>& Comp, const BasisVector& Basis_Vector);
    std::list<Chunk>::const_iterator Set_Obs_Start_Point (const std::list<Chunk>& Obs, const BasisVector& Basis_Vector);
    std::list<Chunk>::const_iterator Set_Hid_Start_Point (const std::list<Chunk>& Hid, const BasisVector& Basis_Vector);
};
