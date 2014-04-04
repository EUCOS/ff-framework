/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/ComplexChunkGenerator.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the ComplexChunkGenerator which generates
 * the complex observable or hidden right moving chunks.
 */
#pragma once

#include <iostream>
#include <set>
#include <vector>

#include <BasisAlpha.h>

/*!
 * This class can generate the complex observable or hidden right moving chunks
 * of a basis vector for an FFHS model (bar-psis, bar-etas or bar-phis). It
 * uses the matching boundary conditions of the previous layers of basis
 * vectors in a model to determine the largest nonredundant combinations of
 * values for these elements.
 */
class ComplexChunkGenerator {
  public:
    ComplexChunkGenerator(const std::vector<BasisAlpha>& Common_Basis_Alphas, int Order, const std::vector<int>& Complex_Elements);
    ComplexChunkGenerator(const std::vector<BasisAlpha>& Common_Basis_Alphas, int Order, const std::vector<int>& Complex_Elements, bool First_Complex_Extension);
    ComplexChunkGenerator(const ComplexChunkGenerator& New_ComplexChunkGenerator);

    ~ComplexChunkGenerator() {}

    void Build_Complex_Chunks();

    void Display_Complex_Chunks() const;

    const std::vector<BasisAlpha>& Common_Basis_Alphas() const {
      return Common_Basis_Alphas_;
    }

    int Order() const {
      return Order_;
    }

    bool First_Complex_Extension() const {
      return First_Complex_Extension_;
    }

    const std::vector<int>& Complex_Elements() const {
      return Complex_Elements_;
    }

    const std::set<std::vector<int> >& Complex_Chunks() const {
      return Complex_Chunks_;
    }

    inline void Set_Common_Basis_Alphas(std::vector<BasisAlpha> New_Common_Basis_Alphas) {
      Common_Basis_Alphas_ = New_Common_Basis_Alphas;
    }

    inline void Set_Order(int New_Order) {
      Order_ = New_Order;
    }

    inline void Set_First_Complex_Extension(bool New_First_Complex_Extension) {
      First_Complex_Extension_ = New_First_Complex_Extension;
    }

    inline void Set_Complex_Elements(std::vector<int> New_Complex_Elements) {
      Complex_Elements_ = New_Complex_Elements;
    }

    inline void Set_Complex_Chunks(std::set<std::vector<int> > New_Complex_Chunks) {
      Complex_Chunks_ = New_Complex_Chunks;
    }

  private:
    std::vector<BasisAlpha> Common_Basis_Alphas_;
    int Order_;
    std::vector<int> Complex_Elements_;
    bool First_Complex_Extension_;

    std::set<std::vector<int> > Complex_Chunks_;

    std::vector<std::vector<int> > Find_Matching_BCs();
    void Generate_Complex_Chunks(const std::vector<std::vector<int> >& Matching_BCs, int MBC_Row_Tracker, int MBC_Col_Tracker, std::vector<int> Real_Chunks_Loader);
};
