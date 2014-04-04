/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/GSOCoefficientMatrixGenerator.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file declares the GSOCoefficientMatrixGenerator class.
 */
#pragma once

#include <iostream>
#include <list>
#include <vector>

/*!
 * This class generates GSO matrix coefficients.
 */
class GSOCoefficientMatrixGenerator {
  public:
    GSOCoefficientMatrixGenerator(std::vector<int> Orders, int Layer);
    GSOCoefficientMatrixGenerator(const GSOCoefficientMatrixGenerator& New_GSOCoefficientMatrixGenerator);

    ~GSOCoefficientMatrixGenerator() {}

    void Build_GSO_Coefficient_Extensions();

    void Display_GSO_Coefficient_Matrix_Extensions() const;

    std::vector<int> Orders() const {
      return Orders_;
    }

    int Layer() const {
      return Layer_;
    }

    const std::list<std::vector<int> >& GSO_Coefficient_Matrix_Extensions() const {
      return GSO_Coefficient_Matrix_Extensions_;
    }


    void Set_Orders(std::vector<int> New_Orders) {
      Orders_ = New_Orders;
    }

    void Set_Layer(int New_Layer) {
      Layer_ = New_Layer;
    }

    void Set_GSO_Coefficient_Matrix_Extensions(const std::list<std::vector<int> >& New_GSO_Coefficient_Matrix_Extensions) {
      GSO_Coefficient_Matrix_Extensions_ = New_GSO_Coefficient_Matrix_Extensions;
    }

  private:
    std::vector<int> Orders_;
    int Layer_;
    std::list<std::vector<int> > GSO_Coefficient_Matrix_Extensions_;

    void Extend_GSO_Coefficient_Matrix(std::vector<int> Row, int Index);
};
