/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/GSOProjector.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Herein we declare the GSOProjector class to facilitate the application of
 * the GSO projection on the particle states.
 */
#pragma once

#include <iostream>
#include <list>
#include <map>
#include <vector>

#include <Alpha.h>
#include <BasisAlpha.h>
#include <GSOCoefficientMatrix.h>
#include <State.h>

/*!
 * This class performs the GSO projections onto the states from a specified
 * sector.  To perform the calculation, all of the basis alphas with common
 * denominators are needed, as well as the GSO coefficient matrix and the
 * fermion mode map.
 */
class GSOProjector {
  public:
    GSOProjector() {}
    GSOProjector(const std::vector<BasisAlpha>& Common_Basis_Alphas, char Alpha_Type, const GSOCoefficientMatrix& k_ij, const std::vector<int>& Coefficients, const std::map<int, int>& Fermion_Mode_Map);
    GSOProjector(const std::vector<BasisAlpha>& Common_Basis_Alphas, const Alpha& The_Alpha, const GSOCoefficientMatrix& k_ij, const std::map<int, int>& Fermion_Mode_Map);
    GSOProjector(const GSOProjector& New_GSOProjector);

    ~GSOProjector() {}

    bool GSOP(const State& The_State) const;

    const std::vector<BasisAlpha>& Common_Basis_Alphas() const {
      return Common_Basis_Alphas_;
    }

    char Alpha_Type() const {
      return Alpha_Type_;
    }

    const GSOCoefficientMatrix& k_ij() const {
      return k_ij_;
    }

    const std::vector<int>& Coefficients() const {
      return Coefficients_;
    }

    const std::map<int, int>& Fermion_Mode_Map() const {
      return Fermion_Mode_Map_;
    }

    void Set_Common_Basis_Alphas(std::vector<BasisAlpha> New_Common_Basis_Alphas) {
      Common_Basis_Alphas_ = New_Common_Basis_Alphas;
    }

    void Set_Alpha_Type(char New_Alpha_Type) {
      Alpha_Type_ = New_Alpha_Type;
    }

    void Set_k_ij(GSOCoefficientMatrix New_k_ij) {
      k_ij_ = New_k_ij;
    }

    void Set_Coefficients(std::vector<int> New_Coefficients) {
      Coefficients_ = New_Coefficients;
    }

    void Set_Fermion_Mode_Map(std::map<int, int> New_Fermion_Mode_Map) {
      Fermion_Mode_Map_ = New_Fermion_Mode_Map;
    }

  private:
    std::vector<BasisAlpha> Common_Basis_Alphas_;
    char Alpha_Type_;
    GSOCoefficientMatrix k_ij_;
    std::vector<int> Coefficients_;
    std::map<int, int> Fermion_Mode_Map_;

    bool GSOP_Boson(const State& The_State) const;
    bool GSOP_Fermion(const State& The_State) const;
};
