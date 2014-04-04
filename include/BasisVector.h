/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/BasisVector.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the BasisVector class representing
 * the integer-encoded basis vector.
 */
#pragma once

#include <iostream>
#include <vector>

/*!
 * This class holds the integer coded information for the basis vector, its
 * order, and the number of left moving modes, as well as the number of compact
 * right moving modes.
 */
class BasisVector {
  public:
    BasisVector() {}
    BasisVector(const std::vector<int>& BV, int Order);
    BasisVector(const std::vector<int>& BV, int Order, int Large_ST_Dimensions);
    BasisVector(const BasisVector& New_BasisVector);

    ~BasisVector() {}

    void Display() const;

    const std::vector<int>& BV() const {
      return BV_;
    }

    int Order() const {
      return Order_;
    }

    int LM_Size() const {
      return LM_Size_;
    }

    int RM_Compact_Size() const {
      return RM_Compact_Size_;
    }

    void Set_BV(std::vector<int> New_BV) {
      BV_ = New_BV;
    }

    void Set_Order(int New_Order) {
      Order_ = New_Order;
    }

    void Set_LM_Size(int New_LM_Size) {
      LM_Size_ = New_LM_Size;
    }

    void Set_RM_Compact_Size(int New_RM_Compact_Size) {
      RM_Compact_Size_ = New_RM_Compact_Size;
    }

  private:
    //FUNDAMENTAL MEMBERS.
    std::vector<int> BV_;
    int Order_;
    int LM_Size_;
    int RM_Compact_Size_;

    //HELPERS.
    void Calculate_LM_Size();
    void Calculate_RM_Compact_Size();
};
