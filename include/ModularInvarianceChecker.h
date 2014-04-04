/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/ModularInvarianceChecker.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * We define a working class for checking modular invariance.
 */
#pragma once

#include <iostream>
#include <vector>

#include <BasisAlpha.h>
#include <Math.h>

/*!
 * This class holds the functions for checking modular invariance for a given
 * free fermionic heterotic string model.
 */
class ModularInvarianceChecker {
  public:
    ModularInvarianceChecker() {}
    ModularInvarianceChecker(const ModularInvarianceChecker& New_ModularInvarianceChecker);

    ~ModularInvarianceChecker() {}

    bool Test_Modular_Invariance(const std::vector<BasisAlpha>& Basis_Alpha_Set);

  private:
    bool Check_Dot_Products(const std::vector<BasisAlpha>& Basis_Alpha_Set);
    bool Check_Simultaneous_Periodic_Modes(const std::vector<BasisAlpha>& Basis_Alpha_Set);
};
