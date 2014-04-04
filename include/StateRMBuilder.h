/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/StateRMBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Here we declare the StateLMBuilder to build massless right-movers for our
 * states.
 */
#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
#include <vector>

#include <Alpha.h>

/*!
 * This class builds the massless right movers for states corresponding to a
 * boson or fermion sector.
 */
class StateRMBuilder {
  public:
    StateRMBuilder(const std::vector<int>& Alpha_RM_Numerator, int Alpha_RM_Denominator, int Large_ST_Dimensions, const std::map<int, int>& Fermion_Mode_Map);
    StateRMBuilder(const Alpha& The_Alpha, int Large_ST_Dimensions, const std::map<int, int>& Fermion_Mode_Map);
    StateRMBuilder(const StateRMBuilder& New_StateRMBuilder);

    ~StateRMBuilder() {}

    void Build_Massless_State_RMs();

    void Display_Massless_State_RMs() const;

    int Large_ST_Dimensions() const {
      return Large_ST_Dimensions_;
    }

    const std::vector<int>& Alpha_RM_Numerator() const {
      return Alpha_RM_Numerator_;
    }

    int Alpha_RM_Denominator() const {
      return Alpha_RM_Denominator_;
    }

    const std::map<int, int>& Fermion_Mode_Map() const {
      return Fermion_Mode_Map_;
    }

    int LM_Size() const {
      return LM_Size_;
    }

    int Mass_Limit() const {
      return Mass_Limit_;
    }

    const std::map<int, int>& Reverse_Fermion_Mode_Map() const {
      return Reverse_Fermion_Mode_Map_;
    }

    const std::list<std::vector<int> >& Massless_State_RMs() const {
      return Massless_State_RMs_;
    }

    void Set_Large_ST_Dimensions(int New_Large_ST_Dimensions) {
      Large_ST_Dimensions_ = New_Large_ST_Dimensions;
    }

    void Set_Alpha_RM_Numerator(std::vector<int> New_Alpha_RM_Numerator) {
      Alpha_RM_Numerator_ = New_Alpha_RM_Numerator;
    }

    void Set_Alpha_RM_Denominator(int New_Alpha_RM_Denominator) {
      Alpha_RM_Denominator_ = New_Alpha_RM_Denominator;
    }

    void Set_Fermion_Mode_Map(std::map<int, int> New_Fermion_Mode_Map) {
      Fermion_Mode_Map_ = New_Fermion_Mode_Map;
    }

    void Set_LM_Size(int New_LM_Size) {
      LM_Size_ = New_LM_Size;
    }

    void Set_Mass_Limit(int New_Mass_Limit) {
      Mass_Limit_ = New_Mass_Limit;
    }

    void Set_Reverse_Fermion_Mode_Map(std::map<int, int> New_Reverse_Fermion_Mode_Map) {
      Reverse_Fermion_Mode_Map_ = New_Reverse_Fermion_Mode_Map;
    }

    void Set_Massless_State_RMs(std::list<std::vector<int> > New_Massless_State_RMs) {
      Massless_State_RMs_ = New_Massless_State_RMs;
    }

  private:
    int Large_ST_Dimensions_;
    std::vector<int> Alpha_RM_Numerator_;
    int Alpha_RM_Denominator_;
    std::map<int, int> Fermion_Mode_Map_;
    int LM_Size_;
    int Mass_Limit_;
    std::map<int, int> Reverse_Fermion_Mode_Map_;
    std::list<std::vector<int> > Massless_State_RMs_;

    void Build_Reverse_Fermion_Mode_Map();
    bool In_Map(int element);
    bool Real_RM_Mode(int element);
    void Select_F_Operator(int element, std::vector<int>& RM, int Mass);
    int Compute_Mass(std::vector<int>& RM);
    int Mass_Increase_Raise(int element);
    int Mass_Increase_Lower(int element);
    void Apply_Real_F_Operator(int element, std::vector<int> RM, int Mass);
    void Apply_Complex_F_Operator(int element, std::vector<int> RM, int F_Lower, int F_Raise, int Mass);
};
