/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/ChunkConsistencyChecker.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the ChunkConsistencyChecker class.
 */
#pragma once

#include <vector>

#include <Chunk.h>

class ChunkConsistencyChecker {
  public:
    ChunkConsistencyChecker(int BV_Denom, int CBA_Denom);
    ChunkConsistencyChecker(const ChunkConsistencyChecker& New_ChunkConsistencyChecker);

    ~ChunkConsistencyChecker() {}

    bool Check_Modular_Invariance(const Chunk& LM, const Chunk& Obs, const Chunk& Comp, const Chunk& Hid);
    bool Check_D10_Modular_Invariance(const Chunk& LM, const Chunk& Obs, const Chunk& Hid);
    bool Check_Simultaneous_Periodic_Modes(const Chunk& LM, const Chunk& Comp);
    bool Check_BV_Consistency(const Chunk& LM, const Chunk& Obs, const Chunk& Comp, const Chunk& Hid);
    bool Check_D10_BV_Consistency(const Chunk& LM, const Chunk& Obs, const Chunk& Hid);

    int BV_Denom() const {return BV_Denom_;}
    int CBA_Denom() const {return CBA_Denom_;}

    inline void Set_BV_Denom(int New_BV_Denom) {
      BV_Denom_ = New_BV_Denom;
    }

    inline void Set_CBA_Denom(int New_CBA_Denom) {
      CBA_Denom_ = New_CBA_Denom;
    }

  private:
    int BV_Denom_;
    int CBA_Denom_;
};
