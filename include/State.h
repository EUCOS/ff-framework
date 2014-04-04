/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/State.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * We declare a container class for states.
 */
#pragma once

#include <iostream>
#include <map>
#include <vector>

/*!
 * This class holds the data needed to specify a state in the model.
 */
class State {
  public:
    State() {}
    State(const std::vector<int>& Numerator, int Denominator, int LM_Size);
    State(const State& New_State);

    ~State() {}

    void Calculate_Length_Squared(const std::map<int, int>& Fermion_Mode_Map);
    bool Is_Positive(const std::map<int, int>& Fermion_Mode_Map);

    void Display() const;

    virtual  bool operator==(const State& State2) const {
      return Numerator() == State2.Numerator();
    }

    virtual  bool operator!=(const State& State2) const {
      return !(*this == State2);
    }

    virtual  bool operator<(const State& State2) const {
      return Numerator()<State2.Numerator();
    }

    const std::vector<int>& Numerator() const {
      return Numerator_;
    }

    int Denominator() const {
      return Denominator_;
    }

    int LM_Size() const {
      return LM_Size_;
    }

    int Length_Squared_Numerator() const {
      return Length_Squared_Numerator_;
    }

    int Length_Squared_Denominator() const {
      return Length_Squared_Denominator_;
    }

    void Set_Numerator(std::vector<int> New_Numerator) {
      Numerator_ = New_Numerator;
    }

    void Set_Denominator(int New_Denominator) {
      Denominator_ = New_Denominator;
    }

    void Set_LM_Size(int New_LM_Size) {
      LM_Size_ = New_LM_Size;
    }

    void Set_Length_Squared_Numerator(int New_Length_Squared_Numerator) {
      Length_Squared_Numerator_ = New_Length_Squared_Numerator;
    }

    void Set_Length_Squared_Denominator(int New_Length_Squared_Denominator) {
      Length_Squared_Denominator_ = New_Length_Squared_Denominator;
    }


  protected:
    std::vector<int> Numerator_;
    int Denominator_;
    int LM_Size_;

    int Length_Squared_Numerator_;
    int Length_Squared_Denominator_;
};
