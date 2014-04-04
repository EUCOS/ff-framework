/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/Alpha.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the Alpha class; a superclass for the
 * various types of BasisAlpha instances.
 */
#pragma once

#include <iostream>
#include <vector>

#include <BasisAlpha.h>

/*!
 * This class serves as a bridge between the basis alphas, which are derived
 * from the basis vectors in a one to one fashion, and the boson, fermion, and
 * SUSY sectors.
 */
class Alpha : public BasisAlpha {
  public:
    Alpha() {}
    Alpha(const std::vector<int>& Numerator, int Denominator, const std::vector<int>& Coefficients);
    Alpha(const Alpha& New_Alpha);

    ~Alpha() {}

    int Mass_Left();
    int Mass_Right();
    virtual char Type() const;

    bool operator<(const Alpha& Other_Alpha) const {
      return Numerator() < Other_Alpha.Numerator();
    }

    void Display_Coefficients() const;

    const std::vector<int>& Coefficients() const {
      return Coefficients_;
    }

    void Set_Coefficients(const std::vector<int>& New_Coefficients) {
      Coefficients_ = New_Coefficients;
    }

  protected:
    //FUNDAMENTAL MEMBERS.
    std::vector<int> Coefficients_;
};
