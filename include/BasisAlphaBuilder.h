/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/BasisAlphaBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the BasisAlphaBuilder class, a class
 * for constructing BasisAlphas from BasisVectors.
 */
#pragma once

#include <vector>

#include <BasisAlpha.h>
#include <BasisVector.h>
#include <Math.h>
#include <Model.h>

/*!
 * This class takes the basis vectors (BasisVector's) for a model and converts
 * them to actual phase values (BasisAlpha's). It builds them with and without
 * a common denominator. Both are needed to construct an FFHS model.
 */
class BasisAlphaBuilder {
  public:
    BasisAlphaBuilder(const Model& FFHS_Model);
    BasisAlphaBuilder(const std::vector<BasisVector>& Basis_Vectors);
    BasisAlphaBuilder(const BasisAlphaBuilder& New_BasisAlphaBuilder);

    ~BasisAlphaBuilder() {}

    void Build_Basis_Alphas();
    void Build_Common_Basis_Alphas();

    void Display_Basis_Alphas() const;
    void Display_Common_Basis_Alphas() const;

    const std::vector<BasisVector>& Basis_Vectors() const {
      return Basis_Vectors_;
    }

    const std::vector<BasisAlpha>& Basis_Alphas() const {
      return Basis_Alphas_;
    }

    const std::vector<BasisAlpha>& Common_Basis_Alphas() const {
      return Common_Basis_Alphas_;
    }

    void Set_Basis_Vectors(const std::vector<BasisVector>& New_Basis_Vectors) {
      Basis_Vectors_ = New_Basis_Vectors;
    }

    void Set_Basis_Alphas(const std::vector<BasisAlpha>& New_Basis_Alphas) {
      Basis_Alphas_ = New_Basis_Alphas;
    }

    void Set_Common_Basis_Alphas(const std::vector<BasisAlpha>& New_Common_Basis_Alphas) {
      Common_Basis_Alphas_ = New_Common_Basis_Alphas;
    }

  private:
    std::vector<BasisVector> Basis_Vectors_;

    std::vector<BasisAlpha> Basis_Alphas_;
    std::vector<BasisAlpha> Common_Basis_Alphas_;
};
