/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/GaugeGroupIdentifier.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the GaugeGroupIdentifier class which
 * is designed to identify gauge group instances.
 */
#pragma once

#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <vector>

#include <GaugeGroup.h>
#include <GaugeGroupName.h>
#include <State.h>

/*!
 * This class takes the nonzero positive roots for a gauge group and identifies
 * the group's class, rank, and Kac-Moody level. It also identifies the simple
 * roots of the gauge group, which are needed to determine the Dynkin labels of
 * matter representations.
 */
class GaugeGroupIdentifier {
  public:
    GaugeGroupIdentifier(const std::list<State>& Positive_Roots, const std::map<int, int>& Fermion_Mode_Map);
    GaugeGroupIdentifier(const GaugeGroupIdentifier& New_GaugeGroupIdentifier);

    ~GaugeGroupIdentifier() {}

    GaugeGroup Get_Group();

    void Display_Positive_Roots() const;
    void Display_Simple_Roots() const;

    const std::list<State>& Positive_Roots() const {
      return Positive_Roots_;
    }

    const std::map<int, int>& Fermion_Mode_Map() const {
      return Fermion_Mode_Map_;
    }

    const std::list<State>& Simple_Roots() const {
      return Simple_Roots_;
    }

    int Long_Root_Length_Num() const {
      return Long_Root_Length_Num_;
    }

    void Set_Positive_Roots(std::list<State> New_Positive_Roots) {
      Positive_Roots_ = New_Positive_Roots;
    }

    void Set_Fermion_Mode_Map(std::map<int, int> New_Fermion_Mode_Map) {
      Fermion_Mode_Map_ = New_Fermion_Mode_Map;
    }

    void Set_Simple_Roots(std::list<State> New_Simple_Roots) {
      Simple_Roots_ = New_Simple_Roots;
    }

    void Set_Long_Root_Length_Num(int New_Long_Root_Length_Num) {
      Long_Root_Length_Num_ = New_Long_Root_Length_Num;
    }

  private:
    std::list<State> Positive_Roots_;
    std::map<int, int> Fermion_Mode_Map_;

    std::list<State> Simple_Roots_;
    int Long_Root_Length_Num_;

    void Identify_KM_Level();
    char Identify_Class();
    char Identify_ADE();
    char Identify_BCGF();
    char Identify_ADE_Degeneracies();
    char Identify_BCGF_Degeneracies();
    void Find_Simple_Roots();
    void Find_Simple_Roots(char Gauge_Group_Class, int Rank);
    bool Connected_Simple_Roots(const std::list<State>& Simple_Roots);
    void Order_Simple_Roots(char Class);
    int Gauge_Dot(const State& State1, const State& State2);
    double A_Class_Rank();
    int Count_Short_Roots();

    GaugeGroup Build_A_Class_Group();
    GaugeGroup Build_B_Class_Group();
    GaugeGroup Build_C_Class_Group();
    GaugeGroup Build_D_Class_Group();
    GaugeGroup Build_E_Class_Group();
    GaugeGroup Build_F_Class_Group();
    GaugeGroup Build_G_Class_Group();
};
