/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/ModelBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Here we define the primary working class, ModelBuilder, which is used to
 * construct models.
 */
#pragma once

#include <iostream>
#include <list>
#include <map>
#include <vector>

#include <Alpha.h>
#include <AlphaBoson.h>
#include <AlphaBuilder.h>
#include <AlphaFermion.h>
#include <AlphaSUSY.h>
#include <BasisAlpha.h>
#include <BasisAlphaBuilder.h>
#include <BasisVector.h>
#include <FermionModeMapBuilder.h>
#include <GSOCoefficientMatrix.h>
#include <GSOCoefficientMatrixBuilder.h>
#include <GaugeGroup.h>
#include <GaugeGroupIdentifier.h>
#include <MatterState.h>
#include <Model.h>
#include <ModularInvarianceChecker.h>
#include <NAHESetLoader.h>
#include <NAHEVariationLoader.h>
#include <State.h>
#include <StateBuilder.h>

/*!
 * This class uses the other elements of the FF Framework to construct a free
 * fermionic heterotic string model. Has an interface for setting the inputs,
 * checking for modular invariance, constructing and altering the model.
 */
class ModelBuilder {
  public:
    ModelBuilder(int Large_ST_Dimensions);
    ModelBuilder(const ModelBuilder& New_ModelBuilder);

    ~ModelBuilder() {}

    void Load_S_Vector(int Large_ST_Dimensions);
    void Load_Basis_Vector(const BasisVector& New_Basis_Vector);
    void Load_k_ij_Row(const std::vector<int>& New_k_ij_Row);
    void Load_Default_k_ij();
    void Load_NAHE_Set();
    void Load_NAHE_Variation();
    bool Check_Modular_Invariance();
    bool Check_Linear_Independence();
    bool Check_k_ij_Consistency();
    bool Check_Model_Consistency();
    void Build_Gauge_Group_Model();
    void Build_Model();

    void Display_Gauge_Group_Roots() const;

    const Model& FFHS_Model() const {
      return FFHS_Model_;
    }

    Model& rFFHS_Model() {
      return FFHS_Model_;
    }

    bool Consistent_Fermion_Mode_Pairs() const {
      return Consistent_Fermion_Mode_Pairs_;
    }

    bool Linearly_Independent_Alphas() const {
      return Linearly_Independent_Alphas_;
    }

    bool Consistent_GSO_Matrix() const {
      return Consistent_GSO_Matrix_;
    }

    const std::vector<BasisAlpha>& Basis_Alphas() const {
      return Basis_Alphas_;
    }

    const std::vector<BasisAlpha>& Common_Basis_Alphas() const {
      return Common_Basis_Alphas_;
    }

    const std::set<AlphaBoson>& Alpha_Bosons() const {
      return Alpha_Bosons_;
    }

    const std::set<AlphaFermion>& Alpha_Fermions() const {
      return Alpha_Fermions_;
    }

    const std::set<AlphaSUSY>& Alpha_SUSYs() const {
      return Alpha_SUSYs_;
    }

    const std::map<int, int>& Fermion_Mode_Map() const {
      return Fermion_Mode_Map_;
    }

    const std::list<State>& SUSY_States() const {
      return SUSY_States_;
    }

    void Set_FFHS_Model(Model New_FFHS_Model) {
      FFHS_Model_ = New_FFHS_Model;
    }

    void Set_Model_k_ij(GSOCoefficientMatrix New_k_ij) {
      FFHS_Model_.Set_k_ij(New_k_ij);
    }

    void Set_Consistent_Fermion_Mode_Pairs(bool New_Consistent_Fermion_Mode_Pairs) {
      Consistent_Fermion_Mode_Pairs_ = New_Consistent_Fermion_Mode_Pairs;
    }

    void Set_Linearly_Independent_Alphas(bool New_Linearly_Independent_Alphas) {
      Linearly_Independent_Alphas_ = New_Linearly_Independent_Alphas;
    }

    void Set_Consistent_GSO_Matrix(bool New_Consistent_GSO_Matrix) {
      Consistent_GSO_Matrix_ = New_Consistent_GSO_Matrix;
    }

    void Set_Basis_Alphas(std::vector<BasisAlpha> New_Basis_Alphas) {
      Basis_Alphas_ = New_Basis_Alphas;
    }

    void Set_Common_Basis_Alphas(std::vector<BasisAlpha> New_Common_Basis_Alphas) {
      Common_Basis_Alphas_ = New_Common_Basis_Alphas;
    }

    void Set_Alpha_Bosons(std::set<AlphaBoson> New_Alpha_Bosons) {
      Alpha_Bosons_ = New_Alpha_Bosons;
    }

    void Set_Alpha_Fermion(std::set<AlphaFermion> New_Alpha_Fermions) {
      Alpha_Fermions_ = New_Alpha_Fermions;
    }

    void Set_Alpha_SUSY(std::set<AlphaSUSY> New_Alpha_SUSYs) {
      Alpha_SUSYs_ = New_Alpha_SUSYs;
    }

    void Set_Fermion_Mode_Map(std::map<int, int> New_Fermion_Mode_Map) {
      Fermion_Mode_Map_ = New_Fermion_Mode_Map;
    }

    void Set_SUSY_States(std::list<State> New_SUSY_States) {
      SUSY_States_ = New_SUSY_States;
    }

  private:
    Model FFHS_Model_;
    bool Consistent_Fermion_Mode_Pairs_;
    bool Linearly_Independent_Alphas_;
    bool Consistent_GSO_Matrix_;
    std::vector<BasisAlpha> Basis_Alphas_;
    std::vector<BasisAlpha> Common_Basis_Alphas_;
    std::set<AlphaBoson> Alpha_Bosons_;
    std::set<AlphaFermion> Alpha_Fermions_;
    std::set<AlphaSUSY> Alpha_SUSYs_;
    std::map<int, int> Fermion_Mode_Map_;
    std::list<State> SUSY_States_;

    void Build_Basis_Alphas();
    void Build_Alphas();
    void Build_Fermion_Mode_Map();
    void Build_k_ij();
    void Build_Gauge_Groups();
    void Build_SUSY_States();
    void Build_MatterStates();
    void Compute_U1_Factors();

    int Gauge_Dot(const State& State1, const State& State2);
    bool Is_SUSY_Partner(const State& New_MatterState);
    int Find_Rank_Cuts();
};
