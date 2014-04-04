/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/LEEFT.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the LEEFT class for representing the
 * Low Energy Effective Field Theory.
 */
#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <CondensedMatterRepresentation.h>
#include <GaugeGroupName.h>
#include <MatterRepEquivalenceTester.h>
#include <Model.h>
#include <ObservableSector.h>

/*!
 * The Low Energy Effective Field Theory (LEEFT) includes the gauge group,
 * matter representations, number of U(1) group factors, number of spacetime
 * supersymmetries and observable sectors.
 */
class LEEFT {
  public:
    LEEFT() {}
    LEEFT(std::vector<std::string> LEEFT_Data);
    LEEFT(const Model& FFHS_Model);
    LEEFT(const LEEFT& New_LEEFT);

    ~LEEFT() {}

    bool Equal_Matter_Rep_Classes(const std:: vector<CondensedMatterRepresentation>& Matter_Rep_Classes1, const std::vector<CondensedMatterRepresentation>& Matter_Rep_Classes2) const;

    inline bool operator<(const LEEFT& LEEFT2) const {
      if(Gauge_Groups() != LEEFT2.Gauge_Groups()) {
        return Gauge_Groups()<LEEFT2.Gauge_Groups();
      } else if(ST_SUSYs() != LEEFT2.ST_SUSYs()) {
        return ST_SUSYs()<LEEFT2.ST_SUSYs();
      } else if(U1_Factors() != LEEFT2.U1_Factors()) {
        return U1_Factors()<LEEFT2.U1_Factors();
      } else if(Total_NA_Reps()!=LEEFT2.Total_NA_Reps()) {
        return Total_NA_Reps()<LEEFT2.Total_NA_Reps();
      } else if(Matter_Rep_Classes() != LEEFT2.Matter_Rep_Classes()) {
        if(!Equal_Matter_Rep_Classes(Matter_Rep_Classes(), LEEFT2.Matter_Rep_Classes())) {
          return Matter_Rep_Classes()<LEEFT2.Matter_Rep_Classes();
        } else {
          return false;
        }
      } else {
        return false;
      }
    }

    inline bool operator==(const LEEFT& LEEFT2) const {
      return ((Gauge_Groups() == LEEFT2.Gauge_Groups()) && (ST_SUSYs() == LEEFT2.ST_SUSYs()) && (U1_Factors() == LEEFT2.U1_Factors()) && (Total_NA_Reps() == LEEFT2.Total_NA_Reps())&& (Equal_Matter_Rep_Classes(Matter_Rep_Classes(), LEEFT2.Matter_Rep_Classes())));
    }

    void Display() const;
    void Display_Matter_Rep_Classes() const;

    const std::vector<GaugeGroupName>& Gauge_Groups() const {
      return Gauge_Groups_;
    }

    const std::vector<MatterRepresentation>& MatterRepresentations() const {
      return MatterRepresentations_;
    }

    int U1_Factors() const {
      return U1_Factors_;
    }

    int ST_SUSYs() const {
      return ST_SUSYs_;
    }

    int Total_NA_Reps() const {
      return Total_NA_Reps_;
    }

    const std::vector<CondensedMatterRepresentation>& Matter_Rep_Classes() const {
      return Matter_Rep_Classes_;
    }

    const std::vector<ObservableSector>& ObservableSectors() const {
      return ObservableSectors_;
    }

    void Set_Gauge_Groups(std::vector<GaugeGroupName> New_Gauge_Groups) {
      Gauge_Groups_ = New_Gauge_Groups;
    }

    void Set_MatterRepresentations(const std::vector<MatterRepresentation> New_MatterRepresentations) {
      MatterRepresentations_ = New_MatterRepresentations;
    }

    void Set_U1_Factors(int New_U1_Factors) {
      U1_Factors_ = New_U1_Factors;
    }

    void Set_ST_SUSYs(int New_ST_SUSYs) {
      ST_SUSYs_ = New_ST_SUSYs;
    }

    void Set_Total_NA_Reps(int New_Total_NA_Reps) {
      Total_NA_Reps_ = New_Total_NA_Reps;
    }

    void Set_Matter_Rep_Classes(const std:: vector<CondensedMatterRepresentation> New_Matter_Rep_Classes) {
      Matter_Rep_Classes_ = New_Matter_Rep_Classes;
    }

    void Set_ObservableSectors(const std::vector<ObservableSector>& New_ObservableSectors) {
      ObservableSectors_ = New_ObservableSectors;
    }

  private:
    std::vector<GaugeGroupName> Gauge_Groups_;
    std::vector<MatterRepresentation> MatterRepresentations_;
    int U1_Factors_;
    int ST_SUSYs_;
    int Total_NA_Reps_;
    std::vector<CondensedMatterRepresentation> Matter_Rep_Classes_;
    std::vector<ObservableSector> ObservableSectors_;

    void Build_Matter_Rep_Classes();
};
