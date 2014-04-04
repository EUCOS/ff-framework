/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/CondensedMatterRepresentation.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 */
#pragma once

#include <GaugeGroupName.h>
#include <MatterRepresentation.h>

class CondensedMatterRepresentation : public MatterRepresentation {
  public:
    CondensedMatterRepresentation(const std::vector<GroupRepresentation>& Rep_Dimension, int Duplicates, const std::vector<GaugeGroupName>& Gauge_Groups);
    CondensedMatterRepresentation(const CondensedMatterRepresentation& New_CondensedMatterRepresentation);

    ~CondensedMatterRepresentation() {}

    friend bool operator<(const CondensedMatterRepresentation& Rep1, const CondensedMatterRepresentation& Rep2) {
      if(Rep1.Gauge_Groups() != Rep2.Gauge_Groups())
        return Rep1.Gauge_Groups()<Rep2.Gauge_Groups();
      if(Rep1.Duplicates() != Rep2.Duplicates())
        return Rep1.Duplicates()<Rep2.Duplicates();
      else
        return Rep1.Rep_Dimension()<Rep2.Rep_Dimension();
    }

    friend bool operator==(const CondensedMatterRepresentation& Rep1, const CondensedMatterRepresentation& Rep2) {
      return ((Rep1.Duplicates() == Rep2.Duplicates()) &&
          (Rep1.Rep_Dimension() == Rep2.Rep_Dimension()) &&
          (Rep1.Gauge_Groups() == Rep2.Gauge_Groups()));
    }

    void Display() const;

    const std::vector<GaugeGroupName>& Gauge_Groups() const {
      return Gauge_Groups_;
    }

    const std::vector<int> Gauge_Group_Indices() const {
      return Gauge_Group_Indices_;
    }

    void Set_Gauge_Groups(const std::vector<GaugeGroupName>& New_Gauge_Groups) {
      Gauge_Groups_ = New_Gauge_Groups;
    }

    void Set_Gauge_Group_Indices(const std::vector<int>& New_Gauge_Group_Indices) {
      Gauge_Group_Indices_ = New_Gauge_Group_Indices;
    }

  protected:
    std::vector<GaugeGroupName> Gauge_Groups_;
    std::vector<int> Gauge_Group_Indices_;

    void Collapse_MatterRepresentation(const std::vector<GaugeGroupName>& Gauge_Groups);
};
