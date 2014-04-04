/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/RealChunkGenerator.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 */
#pragma once

#include <iostream>
#include <list>
#include <set>
#include <vector>

#include <BasisAlpha.h>

/*!
 * This class generates the real chunk of a basis vector for an FFHS model
 * (bar-ys and bar-ws) without redundancies.
 *
 * It requires the common basis alphas for a model prior to the layer being
 * built by this class. It finds the matching boundary conditions for that set
 * of basis alphas and then within those boundaries, builds the possible
 * nonredundant boundary conditions available.
 */
class RealChunkGenerator {
  public:
    RealChunkGenerator(std::vector<BasisAlpha> Common_Basis_Alphas, int Order);
    RealChunkGenerator(const std::vector<BasisAlpha>& Common_Basis_Alphas, int Order, bool First_Complex_Extension);
    RealChunkGenerator(const RealChunkGenerator& New_RealChunkGenerator);

    ~RealChunkGenerator() {}

    void Build_Real_Chunks();

    void Display_Real_Chunks() const;

    const std::vector<BasisAlpha>& Common_Basis_Alphas() const {
      return Common_Basis_Alphas_;
    }

    int Order() const {
      return Order_;
    }

    bool First_Complex_Extension() const {
      return First_Complex_Extension_;
    }

    const std::set<std::vector<int> >& SP_Real_Chunks() const {
      return SP_Real_Chunks_;
    }

    const std::set<std::vector<int> >& NSP_Real_Chunks() const {
      return NSP_Real_Chunks_;
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

    inline void Set_SP_Real_Chunks(std::set<std::vector<int> > New_SP_Real_Chunks) {
      SP_Real_Chunks_ = New_SP_Real_Chunks;
    }

    inline void Set_NSP_Real_Chunks(std::set<std::vector<int> > New_NSP_Real_Chunks) {
      NSP_Real_Chunks_ = New_NSP_Real_Chunks;
    }


  private:
    std::vector<BasisAlpha> Common_Basis_Alphas_;
    int Order_;
    bool First_Complex_Extension_;
    std::set< std::vector<int> > SP_Real_Chunks_;
    std::set< std::vector<int> > NSP_Real_Chunks_;

    std::vector<std::vector<int> > Find_Matching_BCs();
    void Generate_Real_Chunks(const std::vector<std::vector<int> >& Matching_BCs, int MBC_Row_Tracker, int MBC_Col_Tracker, std::vector<int> Real_Chunks_Loader);
    void Add_Real_Chunk(const std::vector<std::vector<int> >& Matching_BCs, const std::vector<int>& Real_Chunks_Loader);
    bool Is_Simply_Paired(const std::vector<std::vector<int> >& Matching_BCs, const std::vector<int>& Real_Chunks_Loader);
};
