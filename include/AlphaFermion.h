/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/AlphaFermion.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the AlphaFermion class; a derived class
 * of the Alpha class.
 */
#pragma once

#include <iostream>
#include <vector>

#include <Alpha.h>

/*!
 * This class holds the data for a fermion sector, as well as a character to
 * indicate which subclass of Alpha is being passed to StateBuilder.
 */
class AlphaFermion : public Alpha {
  public:
    AlphaFermion(const std::vector<int>& Numerator, int Denominator, const std::vector<int>& Coefficients);
    AlphaFermion(const AlphaFermion& New_AlphaFermion);

    ~AlphaFermion() {}

    char Type() const;
};
