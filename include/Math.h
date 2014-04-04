/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/Math.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Provides several relavent mathematical functions.
 */
#pragma once

#include <cmath>
#include <cstdlib>

/*!
 * We include several functions for computing some mathematical quantities.
 */
namespace FF {
  int GCD(int a, int b);
  int LCM(int a, int b);
  bool abseq(int a, int b);
}
