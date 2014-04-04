/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/Chunk.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * This file provides the declaration of the Chunk class which is designed to
 * systematically construct BasisVectors.
 */
#pragma once

#include <iostream>
#include <vector>

class Chunk {
  public:
    Chunk(const std::vector<int>& BV_Chunk, const std::vector<int>& MI_Dot_Products, const std::vector<bool>& Simultaneous_Periodic_Modes);
    Chunk(const std::vector<int>& BV_Chunk, const std::vector<int>& MI_Dot_Products);
    Chunk(const std::vector<int>& BV_Chunk);
    Chunk(const Chunk& New_Chunk);

    ~Chunk() {}

    void Display_BV_Chunk() const;
    void Display_MI_Dot_Products() const;

    inline bool operator<(const Chunk& Chunk2) const {
      return BV_Chunk() < Chunk2.BV_Chunk();
    }

    const std::vector<int>& BV_Chunk() const {
      return BV_Chunk_;
    }

    const std::vector<int>& MI_Dot_Products() const {
      return MI_Dot_Products_;
    }

    const std::vector<bool>& Simultaneous_Periodic_Modes() const {
      return Simultaneous_Periodic_Modes_;
    }

    inline void Set_BV_Chunk(std::vector<int> New_BV_Chunk) {
      BV_Chunk_ = New_BV_Chunk;
    }

    inline void Set_MI_Dot_Products(std::vector<int> New_MI_Dot_Products) {
      MI_Dot_Products_ = New_MI_Dot_Products;
    }

    inline void Set_Simultaneous_Periodic_Modes(std::vector<bool> New_Simultaneous_Periodic_Modes) {
      Simultaneous_Periodic_Modes_ = New_Simultaneous_Periodic_Modes;
    }

  private:
    std::vector<int> BV_Chunk_;
    std::vector<int> MI_Dot_Products_;
    std::vector<bool> Simultaneous_Periodic_Modes_;
};
