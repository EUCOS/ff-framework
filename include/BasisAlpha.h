/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/BasisAlpha.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the BasisAlpha class representing the
 * phases fermions pick up under parallel transport.
 */
#pragma once

#include <iostream>
#include <vector>

#include <BasisVector.h>

/*!
 * This class holds the actual phases of the fermions modes, serving as a basis
 * set for the excitation space.
 */
class BasisAlpha {
  public:
    BasisAlpha() {}
    BasisAlpha(const std::vector<int>& Numerator, int Denominator);
    BasisAlpha(const std::vector<int>& Numerator, int Denominator, const BasisVector& BV);
    BasisAlpha(const std::vector<int>& Numerator, int Denominator, int Large_ST_Dimensions);
    BasisAlpha(const BasisAlpha& New_BasisAlpha);

    ~BasisAlpha() {}

    int Lorentz_Dot(const BasisAlpha& BasisAlpha_2) const;

    void Display() const;

    const std::vector<int>& Numerator() const {
      return Numerator_;
    }

    int Denominator() const {
      return Denominator_;
    }

    int LM_Size() const {
      return LM_Size_;
    }

    int RM_Compact_Size() const {
      return RM_Compact_Size_;
    }

    void Set_Numerator(std::vector<int> New_Numerator) {
      Numerator_ = New_Numerator;
    }

    void Set_Denominator(int New_Denominator) {
      Denominator_ = New_Denominator;
    }

    void Set_LM_Size(int New_LM_Size) {
      LM_Size_ = New_LM_Size;
    }

    void Set_RM_Compact_Size(int New_RM_Compact_Size) {
      RM_Compact_Size_ = New_RM_Compact_Size;
    }

  private:
    //Fundamental members.
    std::vector<int> Numerator_;
    int Denominator_;
    int LM_Size_;
    int RM_Compact_Size_;

    //Helpers.
    void Calculate_LM_Size();
    void Calculate_RM_Compact_Size();
};
