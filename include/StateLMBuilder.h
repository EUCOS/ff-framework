/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/StateLMBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Here we declare the StateLMBuilder to build massless left-movers for our
 * states.
 */
#pragma once

#include <iostream>
#include <list>
#include <map>
#include <vector>

#include <Alpha.h>

/*!
 * This class builds the massless left movers for states corresponding to a
 * fermion or SUSY sector.
 */
class StateLMBuilder {
  public:
    StateLMBuilder(const std::vector<int>& Alpha_LM_Numerator, int Alpha_LM_Denominator, int Large_ST_Dimensions, const std::map<int, int>& Fermion_Mode_Map);
    StateLMBuilder(const Alpha& The_Alpha, int Large_ST_Dimensions, const std::map<int, int>& Fermion_Mode_Map);
    StateLMBuilder(const StateLMBuilder& New_StateLMBuilder);

    ~StateLMBuilder() {}

    void Build_Massless_State_LMs();

    void Display_Massless_State_LMs() const;
    void Display_Fermion_Mode_Map() const;

    int Large_ST_Dimensions() const {
      return Large_ST_Dimensions_;
    }

    const std::vector<int>& Alpha_LM_Numerator() const {
      return Alpha_LM_Numerator_;
    }

    int Alpha_LM_Denominator() const {
      return Alpha_LM_Denominator_;
    }

    const std::map<int, int>& Fermion_Mode_Map() const {
      return Fermion_Mode_Map_;
    }

    const std::list<std::vector<int> >& ST_LM_Modes() const {
      return ST_LM_Modes_;
    }

    const std::list<std::vector<int> >& Compact_LM_Modes() const {
      return Compact_LM_Modes_;
    }

    const std::list<std::vector<int> >& Massless_State_LMs() const {
      return Massless_State_LMs_;
    }

    void Set_Large_ST_Dimensions(int New_Large_ST_Dimensions) {
      Large_ST_Dimensions_ = New_Large_ST_Dimensions;
    }

    void Set_Alpha_LM_Numerator(std::vector<int> New_Alpha_LM_Numerator) {
      Alpha_LM_Numerator_ = New_Alpha_LM_Numerator;
    }

    void Set_Alpha_LM_Denominator(int New_Alpha_LM_Denominator) {
      Alpha_LM_Denominator_ = New_Alpha_LM_Denominator;
    }

    void Set_ST_LM_Modes(std::list<std::vector<int> > New_ST_LM_Modes) {
      ST_LM_Modes_ = New_ST_LM_Modes;
    }

    void Set_Compact_LM_Modes(std::list<std::vector<int> > New_Compact_LM_Modes) {
      Compact_LM_Modes_ = New_Compact_LM_Modes;
    }

    void Set_Fermion_Mode_Map(std::map<int, int> New_Fermion_Mode_Map) {
      Fermion_Mode_Map_ = New_Fermion_Mode_Map;
    }

    void Set_Massles_State_LMs(std::list<std::vector<int> > New_Massless_State_LMs) {
      Massless_State_LMs_ = New_Massless_State_LMs;
    }

  private:
    int Large_ST_Dimensions_;
    std::vector<int> Alpha_LM_Numerator_;
    int Alpha_LM_Denominator_;
    std::map<int, int> Fermion_Mode_Map_;
    std::list<std::vector<int> > ST_LM_Modes_;
    std::list<std::vector<int> > Compact_LM_Modes_;
    std::list<std::vector<int> > Massless_State_LMs_;

    void Build_ST_LM_Modes();
    void Build_Compact_LM_Modes(int element, std::vector<int>& LM_Compact);
    bool Real_LM_Mode(int element);
    bool In_Map(int element);
};
