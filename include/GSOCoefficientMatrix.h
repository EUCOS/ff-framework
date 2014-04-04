/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/GSOCoefficientMatrix.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the GSOCoefficientMatrix class for
 * representing a GSO projection matrix.
 */
#pragma once

#include <iostream>
#include <vector>

/*!
 * This class holds the GSO coefficient matrix for the model.
 */
class GSOCoefficientMatrix {
  public:
    GSOCoefficientMatrix() {}
    GSOCoefficientMatrix(const std::vector<std::vector<int> >& Numerators, const std::vector<int>& Denominators);
    GSOCoefficientMatrix(const GSOCoefficientMatrix& New_GSOCoefficientMatrix);

    ~GSOCoefficientMatrix() {}

    void Load_GSOCoefficientMatrix_Row(const std::vector<int>& New_Row);
    void Load_GSOCoefficientMatrix_Order(int New_Order);

    void Display() const;

    const std::vector<std::vector<int> >& Numerators() const {
      return Numerators_;
    }

    const std::vector<int>& Denominators() const {
      return Denominators_;
    }

    void Set_Numerators(std::vector<std::vector<int> > New_Numerators) {
      Numerators_ = New_Numerators;
    }

    void Set_Denominators(std::vector<int> New_Denominators) {
      Denominators_ = New_Denominators;
    }

  private:
    std::vector<std::vector<int> > Numerators_;
    std::vector<int> Denominators_;
};
