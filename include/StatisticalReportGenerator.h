/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/StatisticalReportGenerator.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Herein we define a class for generating statistical reports.
 */
#pragma once

#include <fstream>
#include <map>
#include <set>
#include <vector>

#include <GaugeGroupName.h>
#include <LEEFT.h>
#include <LatexWriter.h>
#include <Statistics.h>

/*!
 * This program generates the statistics for a given set of LEEFTs and writes
 * them to a Latex document.
 */
class StatisticalReportGenerator {
  public:
    StatisticalReportGenerator() {}
    StatisticalReportGenerator(const std::set<LEEFT>& Unique_LEEFTs);
    StatisticalReportGenerator(const StatisticalReportGenerator& New_StatisticalReportGenerator);

    ~StatisticalReportGenerator() {}

    void Generate_Statistical_Report(std::ofstream& Latex_Out, int Large_ST_Dimensions);

    const std::set<LEEFT>& Unique_LEEFTs() const {
      return Unique_LEEFTs_;
    }

    void Set_Unique_LEEFTs(const std::set<LEEFT>& New_Unique_LEEFTs) {
      Unique_LEEFTs_ = New_Unique_LEEFTs;
    }

  private:
    std::set<LEEFT> Unique_LEEFTs_;
    Statistics Stats_;
    LatexWriter Writer_;
};
