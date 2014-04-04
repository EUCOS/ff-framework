/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/Statistics.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * We define the Statistics class for calculating landscape statistics.
 */
#pragma once

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <vector>

#include <LEEFT.h>
#include <GaugeGroupName.h>

class Statistics {
  public:
    Statistics() {}
    Statistics(const Statistics& New_Statistics);

    ~Statistics() {}

    std::map<GaugeGroupName, int> Count_Gauge_Groups (const std::set<LEEFT>& Unique_LEEFTs);
    std::map<int, int> Count_ST_SUSYs(const std::set<LEEFT>& Unique_LEEFTs);
    std::map<std::vector<GaugeGroupName>, int> Count_Gauge_Group_Combinations(const std::set<LEEFT>&  Unique_LEEFTs);
    std::map<int, int> Count_U1_Factors(const std::set<LEEFT>& Unique_LEEFTs);
    std::map<int, int> Count_Gauge_Group_Factors (const std::set<LEEFT>& Unique_LEEFTs);
    std::map<int, int> Count_NA_Singlets(const std::set<LEEFT>& Unique_LEEFTs);
};
