/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/GaugeGroup.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the GaugeGroup class to represent a
 * gauge group instance.
 */
#pragma once

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
#include <set>

#include <State.h>
#include <GaugeGroupName.h>
#include <GroupRepresentation.h>

/*!
 * This class holds the information needed to completely specify the behavior of a 
 * gauge group in the model. It also has a method for computing the dimension 
 * of a state's representation for the group.
 */
class GaugeGroup {
  public:
    GaugeGroup(const std::list<State>& Positive_Roots, GaugeGroupName Name);
    GaugeGroup(const std::list<State>& Positive_Roots, GaugeGroupName Name, const std::list<State>& Simple_Roots);
    GaugeGroup(const GaugeGroup& New_GaugeGroup);

    ~GaugeGroup() {}

    GroupRepresentation Compute_Rep_Dimension(const State& Weight, const std::map<int, int>& Fermion_Mode_Map);

    bool Is_D4() const {
      return Name().Is_D4();
    }

    void Set_Ordered(bool New_Ordered) {
      Name_.Set_Ordered(New_Ordered);
    }

    void Set_V_Ordered(bool New_V_Ordered) {
      Name_.Set_V_Ordered(New_V_Ordered);
    }

    void Display() const;
    void Display_Positive_Roots() const;

    friend bool operator<(const GaugeGroup& GaugeGroup1, const GaugeGroup& GaugeGroup2) {
      return GaugeGroup1.Name()<GaugeGroup2.Name();
    }

    friend bool operator==(const GaugeGroup& GaugeGroup1, const GaugeGroup& GaugeGroup2) {
      return GaugeGroup1.Name() == GaugeGroup2.Name();
    }

    const std::list<State>& Positive_Roots() const {
      return Positive_Roots_;
    }

    const State& Weyl_Vector() const {
      return Weyl_Vector_;
    }

    const GaugeGroupName& Name() const {
      return Name_;
    }

    const std::list<State>& Simple_Roots() const {
      return Simple_Roots_;
    }

    const std::map<std::vector<int>, GroupRepresentation>& Dynkin_Labels() const {
      return Dynkin_Labels_;
    }

    const std::set<int>& Complex_Rep_Dimensions() const {
      return Complex_Rep_Dimensions_;
    }

    void Set_Positive_Roots(std::list<State> New_Positive_Roots) {
      Positive_Roots_ = New_Positive_Roots;
    }

    void Set_Weyl_Vector(State New_Weyl_Vector) {
      Weyl_Vector_ = New_Weyl_Vector;
    }

    void Set_Name(GaugeGroupName New_Name) {
      Name_ = New_Name;
    }

    void Set_Simple_Roots(std::list<State> New_Simple_Roots) {
      Simple_Roots_ = New_Simple_Roots;
    }

    void Set_Dynkin_Labels(const std::map<std::vector<int>, GroupRepresentation>& New_Dynkin_Labels) {
      Dynkin_Labels_ = New_Dynkin_Labels;
    }

    void Set_Complex_Rep_Dimensions(const std::set<int>& New_Complex_Rep_Dimensions) {
      Complex_Rep_Dimensions_ = New_Complex_Rep_Dimensions;
    }

  private:
    std::list<State> Positive_Roots_;
    State Weyl_Vector_;

    GaugeGroupName Name_;
    std::list<State> Simple_Roots_;
    std::map<std::vector<int>, GroupRepresentation> Dynkin_Labels_;
    std::set<int> Complex_Rep_Dimensions_;

    void Build_Weyl_Vector();
    int Apply_Weyl_Formula(const State& Weight, const std::map<int, int>& Fermion_Mode_Map);
    int Gauge_Dot(const State& State1, const State& State2, const std::map<int, int>& Fermion_Mode_Map) const;
    std::vector<int> Compute_Dynkin_Labels(const State& Weight, const std::map<int, int>& Fermion_Mode_Map);
    bool Is_Barred_Rep(const std::vector<int>& Weight_Dynkin_Labels, int Rep_Dimension);
    char Compute_Triality(const std::vector<int>& Weight_Dynkin_Labels);
    bool Not_Highest_Weight(const std::vector<int>& Weight_Dynkin_Labels);
};
