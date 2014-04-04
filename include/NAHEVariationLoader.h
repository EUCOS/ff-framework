/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/NAHEVariationLoader.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Within we provide a class to load the NAHE Variation.
 */
#pragma once

#include <vector>
#include <iostream>

#include <BasisVector.h>
#include <Model.h>

class NAHEVariationLoader {
  public:
    NAHEVariationLoader();
    NAHEVariationLoader(const NAHEVariationLoader& Old_NAHEVariationLoader);

    ~NAHEVariationLoader() {}

    Model Load_NAHE_Variation(const Model& Old_Model, bool With_S);

    void Display_b1() const;
    void Display_b2() const;
    void Display_b3() const;

    const BasisVector& b1() const {
      return b1_;
    }

    const BasisVector& b2() const {
      return b2_;
    }

    const BasisVector& b3() const {
      return b3_;
    }

    inline void Set_b1(BasisVector New_b1) {
      b1_ = New_b1;
    }

    inline void Set_b2(BasisVector New_b2) {
      b2_ = New_b2;
    }

    inline void Set_b3(BasisVector New_b3) {
      b3_ = New_b3;
    }

  private:
    BasisVector b1_;
    BasisVector b2_;
    BasisVector b3_;

    BasisVector Build_b1();
    BasisVector Build_b2();
    BasisVector Build_b3();

    std::vector<int> Make_Vector(int Array[64]);
};
