/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/MatterState.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Here is declared a class to represent matter states.
 */
#pragma once

#include <GroupRepresentation.h>
#include <State.h>

/*!
 * MatterState is a derived class of State used to represent matter states in a
 * model.
 */
class MatterState : public State {
  public:
    MatterState() {}
    MatterState(const std::vector<int>& Numerator, int Denominator, int LM_Size, const std::vector<GroupRepresentation>& Representations);
    MatterState(const MatterState& New_MatterState);

    ~MatterState() {}

    void Display_Representations() const;

    bool operator<(const MatterState& MatterState2) const {
      return Representations()<MatterState2.Representations();
    }

    bool operator==(const MatterState& MatterState2) const {
      return Representations() == MatterState2.Representations();
    }

    bool operator!=(const MatterState& MatterState2) const {
      return !(*this==MatterState2);
    }

    const std::vector<GroupRepresentation>& Representations() const {
      return Representations_;
    }

    void Set_Representations(std::vector<GroupRepresentation> New_Representations) {
      Representations_ = New_Representations;
    }

  private:
    std::vector<GroupRepresentation> Representations_;
};
