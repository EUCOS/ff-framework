/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/ObservableSector.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Here we define a container class for observable matter generations.
 */
#pragma once

#include <cstdlib>
#include <iostream>
#include <vector>

#include <MatterRepresentation.h>

/*!
 * This class holds the information for all matter generations for a choice of
 * observable sector.
 */
class ObservableSector {
  public:
    ObservableSector() {}
    ObservableSector(const ObservableSector& New_ObservableSector);

    ~ObservableSector() {}

    int Count_Chiral_Generations() const;
    int Count_Chiral_Anti_Generations() const;

    void Display() const;

    const std::vector<MatterRepresentation>& Generations() const {
      return Generations_;
    }

    const std::vector<MatterRepresentation>& Barred_Generations() const {
      return Barred_Generations_;
    }

    const std::vector<MatterRepresentation>& Anti_Generations() const {
      return Anti_Generations_;
    }

    const std::vector<MatterRepresentation>& Anti_Barred_Generations() const {
      return Anti_Barred_Generations_;
    }

    void Set_Generations(const std::vector<MatterRepresentation>& New_Generations) {
      Generations_ = New_Generations;
    }

    void Set_Barred_Generations(const std::vector<MatterRepresentation>& New_Barred_Generations) {
      Barred_Generations_ = New_Barred_Generations;
    }

    void Set_Anti_Generations(const std::vector<MatterRepresentation>& New_Anti_Generations) {
      Anti_Generations_ = New_Anti_Generations;
    }

    void Set_Anti_Barred_Generations(const std::vector<MatterRepresentation>& New_Anti_Barred_Generations) {
      Anti_Barred_Generations_ = New_Anti_Barred_Generations;
    }

  private:
    std::vector<MatterRepresentation> Generations_;
    std::vector<MatterRepresentation> Barred_Generations_;
    std::vector<MatterRepresentation> Anti_Generations_;
    std::vector<MatterRepresentation> Anti_Barred_Generations_;
};
