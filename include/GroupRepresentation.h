/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/GroupRepresentation.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the GroupRepresentation class
 * providing information regarding the group representation under which the
 * state transforms.
 */
#pragma once

#include <iostream>

/*!
 * This class provides an object representing a non-abelian group
 * representation including dimension, reality and triality.
 */
class GroupRepresentation {
  public:
    GroupRepresentation() {}
    GroupRepresentation(int Dimension, char Triality);
    GroupRepresentation(int Dimension, char Triality, bool Is_Complex);
    GroupRepresentation(const GroupRepresentation& New_GroupRepresentation);

    ~GroupRepresentation() {}

    void Display() const;

    friend bool operator<(const GroupRepresentation& GroupRepresentation1, const GroupRepresentation& GroupRepresentation2) {
      if(GroupRepresentation1.Dimension() == GroupRepresentation2.Dimension()) {
        if((GroupRepresentation2.Triality() == 'v')&& (GroupRepresentation1.Triality() != 'v')) {
          return true;
        } else {
          return false;
        }
      } else {
        return GroupRepresentation1.Dimension()< GroupRepresentation2.Dimension();
      }
    }

    friend bool operator==(const GroupRepresentation& GroupRepresentation1, const GroupRepresentation& GroupRepresentation2) {
      return ((GroupRepresentation1.Dimension() == GroupRepresentation2.Dimension())&& (GroupRepresentation1.Triality() == GroupRepresentation2.Triality()));
    }

    int Dimension() const {
      return Dimension_;
    }

    char Triality() const {
      return Triality_;
    }

    bool Is_Complex() const {
      return Is_Complex_;
    }

    void Set_Dimension(int New_Dimension) {
      Dimension_ = New_Dimension;
    }

    void Set_Triality(char New_Triality) {
      Triality_ = New_Triality;
    }

    void Set_Is_Complex(bool New_Is_Complex) {
      Is_Complex_ = New_Is_Complex;
    }

  private:
    int Dimension_;
    char Triality_;
    bool Is_Complex_;
};
