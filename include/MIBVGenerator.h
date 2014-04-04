/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/Math.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Here we define a class to generate modularly invariant basis vectors.
 */
#pragma once

#include <iostream>
#include <list>
#include <vector>

#include <BasisAlpha.h>
#include <Chunk.h>
#include <ComplexChunkGenerator.h>
#include <LMGenerator.h>
#include <Math.h>
#include <RealChunkGenerator.h>

class MIBVGenerator {
  public:
    MIBVGenerator(int order, int large_st_dimensions);
    MIBVGenerator(const MIBVGenerator& New_MIBVGenerator);

    ~MIBVGenerator() {}

    void Build_Full_Chunks(const std::vector<BasisAlpha>& Common_Basis_Alphas);
    void Build_Gauge_Chunks(const std::vector<BasisAlpha>& Common_Basis_Alphas);
    void Build_Observable_SU5_Chunks(const std::vector<BasisAlpha>& Common_Basis_Alphas);

    void Display_Order() const;
    void Display_Large_ST_Dimensions() const;
    void Display_LMs() const;
    void Display_RMs_Compact() const;
    void Display_RMs_Observable() const;
    void Display_RMs_Hidden() const;
    void Display_Observable_SU5() const;

    int Large_ST_Dimensions() const {
      return Large_ST_Dimensions_;
    }

    int Order() const {
      return Order_;
    }

    bool First_Complex_Extension() const {
      return First_Complex_Extension_;
    }

    const std::list<Chunk>& SP_LMs() const {
      return SP_LMs_;
    }

    const std::list<Chunk>& NSP_LMs() const {
      return NSP_LMs_;
    }

    const std::list<Chunk>& SP_RMs_Compact() const {
      return SP_RMs_Compact_;
    }

    const std::list<Chunk>& NSP_RMs_Compact() const {
      return NSP_RMs_Compact_;
    }

    const std::list<Chunk>& RMs_Observable() const {
      return RMs_Observable_;
    }

    const std::list<Chunk>& RMs_Hidden() const {
      return RMs_Hidden_;
    }

    const std::list<Chunk>& Observable_SU5() const {
      return Observable_SU5_;
    }

    inline void Set_Large_ST_Dmensions(int New_Large_ST_Dimensions) {
      Large_ST_Dimensions_ = New_Large_ST_Dimensions;
    }

    inline void  Set_Order(int New_Order) {
      Order_ = New_Order;
    }

    inline void Set_First_Complex_Extension(bool New_First_Complex_Extension) {
      First_Complex_Extension_ = New_First_Complex_Extension;
    }

    inline void Set_SP_LMs(std::list<Chunk> New_SP_LMs) {
      SP_LMs_ = New_SP_LMs;
    }

    inline void Set_NSP_LMs(std::list<Chunk> New_NSP_LMs) {
      NSP_LMs_ = New_NSP_LMs;
    }

    inline void Set_SP_RMs_Compact(std::list<Chunk> New_SP_RMs_Compact) {
      SP_RMs_Compact_ = New_SP_RMs_Compact;
    }

    inline void Set_NSP_RMs_Compact(std::list<Chunk> New_NSP_RMs_Compact) {
      NSP_RMs_Compact_ = New_NSP_RMs_Compact;
    }

    inline void Set_RMs_Observable(std::list<Chunk> New_RMs_Observable) {
      RMs_Observable_ = New_RMs_Observable;
    }

    inline void Set_RMs_Hidden(std::list<Chunk> New_RMs_Hidden) {
      RMs_Hidden_ = New_RMs_Hidden;
    }

    inline void Set_Observable_SU5(std::list<Chunk> New_Observable_SU5) {
      Observable_SU5_ = New_Observable_SU5;
    }

  private:
    int Large_ST_Dimensions_;
    int Order_;
    bool First_Complex_Extension_;
    std::list<Chunk> SP_LMs_;
    std::list<Chunk> NSP_LMs_;
    std::list<Chunk> SP_RMs_Compact_;
    std::list<Chunk> NSP_RMs_Compact_;
    std::list<Chunk> RMs_Observable_;
    std::list<Chunk> RMs_Hidden_;
    std::list<Chunk> Observable_SU5_;

    void Build_LMs(const std::vector<BasisAlpha>& Common_Basis_Alphas);
    void Build_RMs_Compact(const std::vector<BasisAlpha>& Common_Basis_Alphas);
    void Build_RMs_Observable(const std::vector<BasisAlpha>& Common_Basis_Alphas);
    void Build_RMs_Hidden(const std::vector<BasisAlpha>& Common_Basis_Alphas);

    std::vector<int> Compute_MI_Dot_Products(int Chunk_Start,  const std::vector<BasisAlpha>& Common_Basis_Alphas, const std::vector<int>& BV_Chunk);
    std::vector<bool> Compute_Simultaneous_Periodic_Modes (int Chunk_Start, const std::vector<BasisAlpha>& Common_Basis_Alphas, const std::vector<int>& BV_Chunk);
    std::vector<int> Make_Basis_Alpha(const std::vector<int>& BV_Chunk, bool Is_LM);
    bool Has_No_Complex_Elements(const std::vector<BasisAlpha>& Common_Basis_Alphas);
};
