/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/StateBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Here we define a worker class to construct model states.
 */
#pragma once

#include <iostream>
#include <list>
#include <map>
#include <vector>

#include <Alpha.h>
#include <AlphaBoson.h>
#include <AlphaFermion.h>
#include <AlphaSUSY.h>
#include <GSOCoefficientMatrix.h>
#include <GSOProjector.h>
#include <State.h>
#include <StateLMBuilder.h>
#include <StateRMBuilder.h>

/*!
 * This class builds the physically consistent massless states for a model.
 * First, all of the massless states are generated, then the GSO projections
 * are performed.
 */
class StateBuilder {
  public:
    StateBuilder(const Alpha& The_Alpha, const std::map<int, int>& Fermion_Mode_Map, const std::vector<BasisAlpha>& Common_Basis_Alphas, const GSOCoefficientMatrix& k_ij);
    StateBuilder(const StateBuilder& New_StateBuilder);

    ~StateBuilder() {}

    void Build_States();

    void Display_The_Alpha() const;
    void Display_States() const;

    const Alpha& The_Alpha() const {
      return The_Alpha_;
    }

    char Alpha_Type() const {
      return Alpha_Type_;
    }

    const std::map<int, int>& Fermion_Mode_Map() const {
      return Fermion_Mode_Map_;
    }

    const GSOProjector& GSO() const {
      return GSO_;
    }

    const std::list<State>& States() const {
      return States_;
    }

    void Set_The_Alpha(const Alpha& New_Alpha) {
      The_Alpha_ = New_Alpha;
    }

    void Set_Alpha_Type(char New_Alpha_Type) {
      Alpha_Type_ = New_Alpha_Type;
    }

    void Set_Fermion_Mode_Map(const std::map<int, int>& New_Fermion_Mode_Map) {
      Fermion_Mode_Map_ = New_Fermion_Mode_Map;
    }

    void Set_GSO(const GSOProjector& New_GSO) {
      GSO_ = New_GSO;
    }

    void Set_States(const std::list<State>& New_States) {
      States_ = New_States;
    }

  private:
    Alpha The_Alpha_;
    char Alpha_Type_;
    std::map<int, int> Fermion_Mode_Map_;
    GSOProjector GSO_;
    std::list<State> States_;

    void Build_Fermion_States();
    void Build_Boson_States();
    void Build_SUSY_States();
    int Calculate_Large_ST_Dimensions(int LM_Size);
};
