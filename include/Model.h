/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/Alpha.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the Alpha class; a superclass for the
 * various types of BasisAlpha instances.
 */
#pragma once

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

#include <BasisVector.h>
#include <GSOCoefficientMatrix.h>
#include <GaugeGroup.h>
#include <Math.h>
#include <MatterRepresentation.h> 
#include <MatterState.h>
#include <State.h>

/*!
 * This class holds the information needed to specify a model, as well as methods
 * for loading an incomplete model.
 */
class Model {
  public:
    Model() {}
    Model(int Large_ST_Dimensions);
    Model(const Model& New_Model);

    ~Model(){}

    void Load_S_Vector(int Large_ST_Dimensions);
    void Load_BV_Set(const BasisVector& New_BV);
    void Load_k_ij_Row(const std::vector<int>& New_k_ij_Row);
    void Add_Gauge_Group(const GaugeGroup& New_Gauge_Group);
    void Sort_Gauge_Groups();
    void Add_MatterState(const MatterState& New_MatterState);
    void Build_MatterRepresentations();
    void Set_Matter_Rep_Sign_Convention();

    void Display_BV_Set() const;
    void Display_k_ij() const;
    void Display_Gauge_Groups() const;
    void Display_MatterRepresentations() const;
    void Display_Particle_Content() const;

    const std::vector<BasisVector>& BV_Set() const {
      return BV_Set_;
    }

    const GSOCoefficientMatrix& k_ij() const {
      return k_ij_;
    }

    const std::vector<GaugeGroup>& Gauge_Groups() const {
      return Gauge_Groups_;
    }

    std::vector<GaugeGroup>& rGauge_Groups() {
      return Gauge_Groups_;
    }

    const std::list<MatterState>& MatterStates() const {
      return MatterStates_;
    }

    const std::vector<MatterRepresentation>& MatterRepresentations() const {
      return MatterRepresentations_;
    }

    const std::list<State>& SUSY_States() const {
      return SUSY_States_;
    }

    int U1_Factors() const {
      return U1_Factors_;
    }

    bool NAHE_Loaded() const {
      return NAHE_Loaded_;
    }

    void Set_BV_Set(std::vector<BasisVector> New_BV_Set) {
      BV_Set_ = New_BV_Set;
    }

    void Set_k_ij(GSOCoefficientMatrix New_k_ij) {
      k_ij_ = New_k_ij;
    }

    void Set_Gauge_Groups(std::vector<GaugeGroup> New_Gauge_Groups) {
      Gauge_Groups_ = New_Gauge_Groups;
    }

    void Set_MatterStates(std::list<MatterState> New_MatterStates) {
      MatterStates_ = New_MatterStates;
    }

    void Set_MatterRepresentations(const std::vector<MatterRepresentation>& New_MatterRepresentations) {
      MatterRepresentations_ = New_MatterRepresentations;
    }

    void Set_SUSY_States(std::list<State> New_SUSY_States) {
      SUSY_States_ = New_SUSY_States;
    }

    void Set_U1_Factors(int New_U1_Factors) {
      U1_Factors_ = New_U1_Factors;
    }

    void Set_NAHE_Loaded(bool New_NAHE_Loaded) {
      NAHE_Loaded_ = New_NAHE_Loaded;
    }

  private:
    std::vector<BasisVector> BV_Set_;
    GSOCoefficientMatrix k_ij_;
    std::vector<GaugeGroup> Gauge_Groups_;
    std::list<MatterState> MatterStates_;
    std::vector<MatterRepresentation> MatterRepresentations_;
    std::list<State> SUSY_States_;
    int U1_Factors_;
    bool NAHE_Loaded_;

    void Load_One_Vector(int Large_ST_Dimensions);
    void Invert_Complex_Rep_Signs(int Gauge_Group_Index, const std::set<int>& Complex_Rep_Dimensions);
    void Invert_Trialities(int Gauge_Group_Index, char Triality1, char Triality2, const std::set<int>& Complex_Rep_Dimensions);
};
