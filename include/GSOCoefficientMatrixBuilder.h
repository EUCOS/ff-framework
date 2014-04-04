/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/GSOCoefficientMatrixBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 */
#pragma once

#include <iostream>
#include <vector>

#include <BasisAlpha.h>
#include <GSOCoefficientMatrix.h>
#include <Math.h>
#include <Model.h>

/*!
 * This class takes the half of the GSO coefficient matrix in integer coded
 * form which is specified by the user, and builds the other half based on the
 * modular invariance constraints. It also checks the matrix to ensure the
 * proper values have been built for the model.
 */
class GSOCoefficientMatrixBuilder {
  public:
    GSOCoefficientMatrixBuilder(const GSOCoefficientMatrix& Half_GSO_Matrix, const std::vector<BasisAlpha>& Common_Basis_Alphas);
    GSOCoefficientMatrixBuilder(const Model& FFHS_Model, const std::vector<BasisAlpha>& Common_Basis_Alphas);
    GSOCoefficientMatrixBuilder(const std::vector<std::vector<int> >& Half_Numerators, const std::vector<int>& Half_Denominators, const std::vector<BasisAlpha>& Common_Basis_Alphas);
    GSOCoefficientMatrixBuilder(const GSOCoefficientMatrixBuilder& New_GSOCoefficientMatrixBuilder);

    ~GSOCoefficientMatrixBuilder() {}

    void Build_Complete_GSO_Matrix();

    void Display_Half_GSO_Matrix() const;
    void Display_Common_Basis_Alphas() const;
    void Display_Complete_GSO_Matrix() const;

    const GSOCoefficientMatrix& Half_GSO_Matrix() const {
      return Half_GSO_Matrix_;
    }

    const std::vector<BasisAlpha>& Common_Basis_Alphas() const {
      return Common_Basis_Alphas_;
    }

    const std::vector<int>& GSO_Coefficient_Orders() const {
      return GSO_Coefficient_Orders_;
    }

    const GSOCoefficientMatrix& Complete_GSO_Matrix() const {
      return Complete_GSO_Matrix_;
    }

    bool Consistent_GSO_Matrix() const {
      return Consistent_GSO_Matrix_;
    }

    void Set_Half_GSO_Matrix(GSOCoefficientMatrix New_Half_GSO_Matrix) {
      Half_GSO_Matrix_ = New_Half_GSO_Matrix;
    }

    void Set_Common_Basis_Alphas(std::vector<BasisAlpha> New_Common_Basis_Alphas) {
      Common_Basis_Alphas_ = New_Common_Basis_Alphas;
    }

    void Set_GSO_Coefficient_Orders(const std::vector<int>& New_GSO_Coefficient_Orders) {
      GSO_Coefficient_Orders_ = New_GSO_Coefficient_Orders;
    }

    void Set_Complete_GSO_Matrix(GSOCoefficientMatrix New_Complete_GSO_Matrix) {
      Complete_GSO_Matrix_ = New_Complete_GSO_Matrix;
    }

    void Set_Consistent_GSO_Matrix(bool New_Consistent_GSO_Matrix) {
      Consistent_GSO_Matrix_ = New_Consistent_GSO_Matrix;
    }

  private:
    GSOCoefficientMatrix Half_GSO_Matrix_;
    std::vector<BasisAlpha> Common_Basis_Alphas_;
    std::vector<int> GSO_Coefficient_Orders_;

    GSOCoefficientMatrix Complete_GSO_Matrix_;
    bool Consistent_GSO_Matrix_;

    void Convert_Half_GSO_Matrix();
    void Complete_Half_GSO_Matrix();
    int Compute_Diagonal_GSO_Element(const BasisAlpha& The_Alpha, int First_GSO_Coefficient_Row);
    int Compute_Off_Diagonal_GSO_Element(const BasisAlpha& Alpha_a, const BasisAlpha& Alpha_b, int Half_GSO_Coefficient);
    int Reset_GSO_Coefficient_Range(int Converted_GSO_Coefficient, int Converted_GSO_Denominator);
    bool Check_GSO_Element_Consistency(int GSO_Element, int column);
};
