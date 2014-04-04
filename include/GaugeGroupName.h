/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/GaugeGroupName.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the GaugeGroupName class which
 * represents the necessary information to specify a gauge group.
 */
#pragma once

#include <iostream>

/*!
 * This class holds the data needed to specify a gauge group's name - the class
 * (in Cartan notation), the rank, and the Kac-Moody level.
 */
class GaugeGroupName {
  public:
    GaugeGroupName();
    GaugeGroupName(char Class, int Rank, int KM_Level);
    GaugeGroupName(const GaugeGroupName& New_GaugeGroupName);

    ~GaugeGroupName() {}

    bool Is_D4() const {
      return ((Class()=='D')&&(Rank()==4));
    }

    void Display() const;

    bool operator<(const GaugeGroupName& GaugeGroupName2) const {
      if(Class() != GaugeGroupName2.Class())
        return (Class()<GaugeGroupName2.Class());
      else if(Rank() != GaugeGroupName2.Rank())
        return (Rank()<GaugeGroupName2.Rank());
      else
        return (KM_Level()<GaugeGroupName2.KM_Level());
    }

    bool operator==(const GaugeGroupName& GaugeGroupName2) const {
      return (Class() == GaugeGroupName2.Class())&&
        (Rank() == GaugeGroupName2.Rank())&&
        (KM_Level() == GaugeGroupName2.KM_Level());
    }

    bool operator!=(const GaugeGroupName& GaugeGroupName2) const {
      return !(*this==GaugeGroupName2);
    }

    char Class() const {
      return Class_;
    }

    int Rank() const {
      return Rank_;
    }

    int KM_Level() const {
      return KM_Level_;
    }

    bool Ordered() const {
      return Ordered_;
    }

    bool V_Ordered() const {
      return V_Ordered_;
    }

    void Set_Class(char New_Class) {
      Class_ = New_Class;
    }

    void Set_Rank(int New_Rank) {
      Rank_ = New_Rank;
    }

    void Set_KM_Level(int New_KM_Level) {
      KM_Level_ = New_KM_Level;
    }

    void Set_Ordered(bool New_Ordered) {
      Ordered_ = New_Ordered;
    }

    void Set_V_Ordered(bool New_V_Ordered) {
      V_Ordered_ = New_V_Ordered;
    }

  private:
    char Class_;
    int Rank_;
    int KM_Level_;
    bool Ordered_;
    bool V_Ordered_;
};
