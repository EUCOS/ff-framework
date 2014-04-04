/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/AlphaBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the AlphaBuilder class; a class
 * designed to build and handle Alphas of various types.
 */
#pragma once

#include <iostream>
#include <set>
#include <vector>

#include <Alpha.h>
#include <AlphaBoson.h>
#include <AlphaFermion.h>
#include <AlphaSUSY.h>
#include <BasisAlpha.h>

/*!
 * This class builds the linear combinations of alphas from a set of basis
 * alphas with matching denominators, checks whether or not they can produce
 * massless states, then classifies them as Alpha_Boson, Alpha_Fermion, and
 * Alpha SUSY.
 */
class AlphaBuilder {
  public:
    AlphaBuilder(const std::vector<BasisAlpha>& Common_BasisAlphas, const std::vector<BasisAlpha>& BasisAlphas);
    AlphaBuilder(const AlphaBuilder& New_BasisAlphaBuilder);

    ~AlphaBuilder() {}

    //INTERFACE.
    void Build_Alphas();

    //DEBUG.
    void Display_Common_BasisAlphas() const;
    void Display_Coefficient_Limits() const;
    void Display_Alpha_Bosons() const;
    void Display_Alpha_Fermions() const;
    void Display_Alpha_SUSYs() const;
    void Display_All_Alphas() const;

    //ACCESSORS.
    const std::vector<BasisAlpha>& Common_BasisAlphas() const {
      return Common_BasisAlphas_;
    }

    const std::vector<int>& Coefficient_Limits() const {
      return Coefficient_Limits_;
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

    bool Linearly_Independent_Alphas() const {
      return Linearly_Independent_Alphas_;
    }

    //SETTERS.
    void Set_Common_BasisAlphas(const std::vector<BasisAlpha>& New_Common_BasisAlphas) {
      Common_BasisAlphas_ = New_Common_BasisAlphas;
    }

    void Set_Coefficient_Limits(const std::vector<int>& New_Coefficient_Limits) {
      Coefficient_Limits_ = New_Coefficient_Limits;
    }

    void Set_Alpha_Bosons(const std::set<AlphaBoson>& New_Alpha_Bosons) {
      Alpha_Bosons_ = New_Alpha_Bosons;
    }

    void Set_Alpha_Fermions(const std::set<AlphaFermion>& New_Alpha_Fermions) {
      Alpha_Fermions_ = New_Alpha_Fermions;
    }

    void Set_Alpha_SUSYs(const std::set<AlphaSUSY>& New_Alpha_SUSYs) {
      Alpha_SUSYs_ = New_Alpha_SUSYs;
    }

    void Set_Linearly_Independent_Alphas(bool New_Linearly_Independent_Alphas) {
      Linearly_Independent_Alphas_ = New_Linearly_Independent_Alphas;
    }

  private:
    //FUNDAMENTAL MEMBERS.
    std::vector<BasisAlpha> Common_BasisAlphas_;
    std::vector<int> Coefficient_Limits_;

    //DERIVED MEMBERS.
    std::set<AlphaBoson> Alpha_Bosons_;
    std::set<AlphaFermion> Alpha_Fermions_;
    std::set<AlphaSUSY> Alpha_SUSYs_;

    bool Linearly_Independent_Alphas_;

    //HELPERS.
    std::vector<int> Get_Coefficient_Limits(const std::vector<BasisAlpha>& BasisAlphas);
    void Add_Alphas(int Layer, Alpha Last_Alpha);
    std::vector<int> Adjust_Alpha_Range(std::vector<int> Alpha_Numerator);
    bool Linearly_Independent_Alpha(const Alpha& Last_Alpha);
};
