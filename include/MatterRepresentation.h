/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/MatterRepresentation.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Herein we define the MatterRepresentation class.
 */
#pragma once

#include <iostream>
#include <vector>

#include <GroupRepresentation.h>

/*!
 * This class holds the information related to the matter representations.
 * It has members for matter representations and quantities.
 */
class MatterRepresentation {
  public:
    MatterRepresentation(const std::vector<GroupRepresentation>& Rep_Dimension, int Duplicates);
    MatterRepresentation(const MatterRepresentation& New_MatterRepresentation);

    ~MatterRepresentation() {}

    friend bool operator<(const MatterRepresentation& MatterRepresentation1, const MatterRepresentation& MatterRepresentation2) {
      if(MatterRepresentation1.Duplicates() !=  MatterRepresentation2.Duplicates()) {
        return MatterRepresentation1.Duplicates()< MatterRepresentation2.Duplicates();
      } else {
        return (MatterRepresentation1.Rep_Dimension()< MatterRepresentation2.Rep_Dimension());
      }
    }

    friend bool operator==(const MatterRepresentation& MatterRepresentation1, const MatterRepresentation& MatterRepresentation2) {
      return ((MatterRepresentation1.Rep_Dimension() == MatterRepresentation2.Rep_Dimension()) && (MatterRepresentation1.Duplicates() == MatterRepresentation2.Duplicates()));
    }

    void Switch_Dimension_Sign(int Dimension_Index);
    void Set_Rep_Dimension_at_Index(const GroupRepresentation& New_Group_Representation, int Index);

    virtual void Display() const;

    const std::vector<GroupRepresentation> Rep_Dimension() const {
      return Rep_Dimension_;
    }

    int Duplicates() const {
      return Duplicates_;
    }

    void Set_Rep_Dimension(const std::vector<GroupRepresentation>& New_Rep_Dimension) {
      Rep_Dimension_ = New_Rep_Dimension;
    }

    void Set_Duplicates(int New_Duplicates) {
      Duplicates_ = New_Duplicates;
    }

  protected:
    std::vector<GroupRepresentation> Rep_Dimension_;
    int Duplicates_;
};
