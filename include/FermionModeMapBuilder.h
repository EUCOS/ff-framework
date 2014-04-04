/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/FermionModeMapBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the FermionModeMap class, a working
 * class designed to pair complex fermions.
 */
#pragma once

#include <iostream>
#include <map>
#include <vector>

#include <BasisAlpha.h>
#include <Model.h>

/*!
 * Finds the sets of simultaneously matching boundary conditions for a set of
 * basis alphas with a common denominator. Those boundary conditions are then
 * placed into complex pairs, and a map between the first and second fermion
 * modes in the pair is created.
 */
class FermionModeMapBuilder {
  public:
    FermionModeMapBuilder();
    FermionModeMapBuilder(int Large_ST_Dimensions);
    FermionModeMapBuilder(const FermionModeMapBuilder& New_FermionModeMapBuilder);

    ~FermionModeMapBuilder() {}

    void Build_Fermion_Mode_Map(const std::vector<BasisAlpha>& Common_Basis_Alphas);

    void Display_Fermion_Mode_Map() const;

    int Large_ST_Dimensions() const {
      return Large_ST_Dimensions_;
    }

    const std::map<int, int>& Fermion_Mode_Map() const {
      return Fermion_Mode_Map_;
    }

    bool Consistent_Pairings() const {
      return Consistent_Pairings_;
    }

    void Set_Large_ST_Dimensions(int New_Large_ST_Dimensions) {
      Large_ST_Dimensions_ = New_Large_ST_Dimensions;
    }

    void Set_Fermion_Mode_Map(const std::map<int, int>& New_Fermion_Mode_Map) {
      Fermion_Mode_Map_ = New_Fermion_Mode_Map;
    }

    void Set_Consistent_Pairings(bool New_Consistent_Pairings) {
      Consistent_Pairings_ = New_Consistent_Pairings;
    }

  private:
    int Large_ST_Dimensions_;

    std::map<int, int> Fermion_Mode_Map_;
    bool Consistent_Pairings_;

    int Compute_Large_ST_Dimensions(const BasisAlpha& Common_Basis_Alpha);
    void Map_Complex_LM_Elements(const std::vector<BasisAlpha>& Common_Basis_Alphas);
    void Map_Complex_RM_Elements(const std::vector<BasisAlpha>& Common_Basis_Alphas);
    void Find_All_LR_Pairs(std::vector<int> LR_Coordinates, const std::vector<BasisAlpha>& Common_Basis_Alphas);
    std::vector<int> Pair_LMs(std::vector<int> LR_Coordinates, const std::vector<BasisAlpha>& Common_Basis_Alphas);
    std::vector<int> Pair_RMs(std::vector<int> LR_Coordinates, const std::vector<BasisAlpha>& Common_Basis_Alphas);
    std::vector<int> Pair_Mixed(std::vector<int> Unpaired_LMs, std::vector<int> Unpaired_RMs, const std::vector<BasisAlpha>& Common_Basis_Alphas);
    std::vector<std::vector<int> > Find_Matching_BCs (std::vector<std::vector<int> > Matching_BCs, const std::vector<BasisAlpha>& Common_Basis_Alphas);
    std::vector<int> Add_Pairs_To_Map(std::vector<std::vector<int> > Final_Pairs);
};
