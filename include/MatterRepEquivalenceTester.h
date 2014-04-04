/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/MatterRepEquivalenceTester.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 */
#pragma once

#include <algorithm>
#include <vector>

#include <CondensedMatterRepresentation.h>
#include <GaugeGroupName.h>

class MatterRepEquivalenceTester {
  public:
    MatterRepEquivalenceTester(const std:: vector<CondensedMatterRepresentation>& Matter_Rep_Class_A, const std::vector<CondensedMatterRepresentation>& Matter_Rep_Class_B, const std::vector<GaugeGroupName>& Gauge_Groups);
    MatterRepEquivalenceTester(const MatterRepEquivalenceTester& New_MatterRepEquivalenceTester);

    ~MatterRepEquivalenceTester() {}

    void Test_Equivalence();

    void Display_Matter_Rep_Class_A() const;

    const std::vector<CondensedMatterRepresentation>& Matter_Rep_Class_A() const {
      return Matter_Rep_Class_A_;
    }

    const std::vector<CondensedMatterRepresentation>& Matter_Rep_Class_B() const {
      return Matter_Rep_Class_B_;
    }

    const std::vector<GaugeGroupName>& Gauge_Groups() const {
      return Gauge_Groups_;
    }

    bool Equivalent_Matter_Classes() const {
      return Equivalent_Matter_Classes_;
    }

    void Set_Matter_Rep_Class_A(const std:: vector<CondensedMatterRepresentation>& New_Matter_Rep_Class_A) {
      Matter_Rep_Class_A_ = New_Matter_Rep_Class_A;
    }

    void Set_Matter_Rep_Class_B(const std:: vector<CondensedMatterRepresentation>& New_Matter_Rep_Class_B) {
      Matter_Rep_Class_B_ = New_Matter_Rep_Class_B;
    }

    void Set_Gauge_Groups(const std::vector<GaugeGroupName>& New_Gauge_Groups) {
      Gauge_Groups_ = New_Gauge_Groups;
    }

    void Set_Equivalent_Matter_Classes(bool New_Equivalent_Matter_Classes) {
      Equivalent_Matter_Classes_ = New_Equivalent_Matter_Classes;
    }

  private:
    std::vector<CondensedMatterRepresentation> Matter_Rep_Class_A_;
    std::vector<CondensedMatterRepresentation> Matter_Rep_Class_B_;
    std::vector<GaugeGroupName> Gauge_Groups_;

    bool Equivalent_Matter_Classes_;

    void Permute_Complex_Representations(int Gauge_Group_Index);
    void Invert_Complex_Reps(int Gauge_Group_Index);
    void Invert_Trialities(char Triality1, char Triality2, int Gauge_Group_Index);
};
