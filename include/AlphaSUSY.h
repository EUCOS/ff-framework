/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/AlphaSUSY.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the AlphaSUSY class; a derived class
 * of the Alpha class.
 */
#pragma once

#include <iostream>
#include <vector>

#include <Alpha.h>

/*!
 * This class holds the data for a SUSY generating sector, as well as a character
 * to indicate which subclass of Alpha is being passed to StateBuilder.
 */
class AlphaSUSY : public Alpha {
  public:
    AlphaSUSY(const std::vector<int>& Numerator, int Denominator, const std::vector<int>& Coefficients);
    AlphaSUSY(const AlphaSUSY& New_AlphaSUSY);

    ~AlphaSUSY() {}

    char Type() const;
};
