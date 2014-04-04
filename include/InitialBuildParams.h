/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/InitialBuildParams.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 */
#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

class InitialBuildParams {
  public:
    InitialBuildParams() {}
    InitialBuildParams(int Large_ST_Dimensions, const std::vector<int>& Orders, const std::vector<std::vector<int> >& Initial_BV_Set);
    InitialBuildParams(const InitialBuildParams& New_InitialBuildParams);

    ~InitialBuildParams() {}

    void Display() const;

    void Build_k_ij_Extensions();
    void Load_Initial_BV_Set(const std::vector<int>& New_BV);

    int Large_ST_Dimensions() const {
      return Large_ST_Dimensions_;
    }

    const std::vector<int>& Orders() const {
      return Orders_;
    }

    const std::vector<std::vector<int> >& Initial_BV_Set() const {
      return Initial_BV_Set_;
    }

    const std::vector<std::vector<int> >& k_ij_Extensions() const {
      return k_ij_Extensions_;
    }


    void Set_Large_ST_Dimensions(int New_Large_ST_Dimensions) {
      Large_ST_Dimensions_ = New_Large_ST_Dimensions;
    }

    void Set_Orders(const std::vector<int>& New_Orders) {
      Orders_ = New_Orders;
    }

    void Set_Initial_BV_Set(const std::vector<std::vector<int> >& New_Initial_BV_Set) {
      Initial_BV_Set_ = New_Initial_BV_Set;
    }

    void Set_k_ij_Extensions(const std::vector<std::vector<int> >& New_k_ij_Extensions) {
      k_ij_Extensions_ = New_k_ij_Extensions;
    }

  private:
    int Large_ST_Dimensions_;
    std::vector<int> Orders_;
    std::vector<std::vector<int> > Initial_BV_Set_;
    std::vector<std::vector<int> > k_ij_Extensions_;

    void Extend_k_ij(std::vector<int> k_ij_Row, int k_ij_index);
};
