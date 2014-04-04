/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/LMGenerator.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Here we declare the LMGenerator class.
 */
#pragma once

#include <iostream>
#include <list>
#include <vector>

#include <BasisAlpha.h>

/*!
 * This class construct valid left-movers satisfying all of the extra
 * constraints imposed on them.
 */
class LMGenerator {
  public:
    LMGenerator(int Large_ST_Dimensions, const std::vector<BasisAlpha>& Common_Basis_Alphas);
    LMGenerator(const LMGenerator& New_LMGenerator);

    ~LMGenerator() {}

    void Build_LM_Chunks();

    void Display_Large_ST_Dimensions() const;
    void Display_Compact_LM_Triplets() const;
    void Display_Matching_BCs() const;
    void Display_LMs() const;

    int Large_ST_Dimensions() const {
      return Large_ST_Dimensions_;
    }

    const std::vector< std::vector<int> >& Compact_LM_Triplets() const {
      return Compact_LM_Triplets_;
    }

    const std::vector< std::vector<int> >& Matching_BCs() const {
      return Matching_BCs_;
    }

    const std::list< std::vector<int> >& SP_LMs() const {
      return SP_LMs_;
    }

    const std::list< std::vector<int> >& NSP_LMs() const {
      return NSP_LMs_;
    }

    inline void Set_Large_ST_Dimensions(int New_Large_ST_Dimensions) {
      Large_ST_Dimensions_ = New_Large_ST_Dimensions;
    }

    inline void Set_Compact_LM_Triplets(std::vector< std::vector<int> > New_Compact_LM_Triplets) {
      Compact_LM_Triplets_ = New_Compact_LM_Triplets;
    }

    inline void Set_Matching_BCs(std::vector< std::vector<int> > New_Matching_BCs) {
      Matching_BCs_ = New_Matching_BCs;
    }

    inline void Set_SP_LMs(std::list< std::vector<int> > New_SP_LMs) {
      SP_LMs_ = New_SP_LMs;
    }

    inline void Set_NSP_LMs(std::list< std::vector<int> > New_NSP_LMs) {
      NSP_LMs_ = New_NSP_LMs;
    }

  private:
    int Large_ST_Dimensions_;
    std::vector< std::vector<int> > Compact_LM_Triplets_;
    std::vector< std::vector<int> > Matching_BCs_;

    std::list< std::vector<int> > SP_LMs_;
    std::list< std::vector<int> > NSP_LMs_;

    void Build_Compact_LM_Triplets();
    void Build_LMs(int Tracker, std::vector<int> LMs_Loader);
    void Find_Matching_BCs(const std::vector<BasisAlpha>& Common_Basis_Alphas);
    void Add_LM_Chunk(const std::vector<int>& LMs_Loader);
    bool Is_Simply_Paired(const std::vector<int>& LMs_Loader);
};
