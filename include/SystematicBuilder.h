/* Copyright 2009-2014 Baylor University, EUCOS Group.
 * All rights reserved. Use of this source code is governed
 * by the MIT license that can be found in the LICENSE file.
 */

/*!
 * @file include/SystematicBuilder.h
 * @author Timothy Renner <timothy_renner@baylor.edu>
 *
 * Here we define the systematic builder class for performing systematic
 * surveys.
 */
#pragma once

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <BasisAlpha.h>
#include <BasisVector.h>
#include <ChunkConsistencyChecker.h>
#include <GSOCoefficientMatrixGenerator.h>
#include <LEEFT.h>
#include <MIBVGenerator.h>
#include <Model.h>
#include <ModelBuilder.h>
#include <OutputWriter.h>

class SystematicBuilder {
  public:
    SystematicBuilder(const ModelBuilder& Initial_Build, int Large_ST_Dimensions, const std::vector<int>& Extension_Orders, const std::string& Output_File_Name);
    SystematicBuilder(const ModelBuilder& Initial_Build, int Large_ST_Dimensions, const std::vector<int>& Extension_Orders, const std::string& Output_File_Name, std::vector<bool> Simply_Paired_Layers);
    SystematicBuilder(const SystematicBuilder& New_SystematicBuilder);

    ~SystematicBuilder() {}

    virtual void Perform_Search();

    void Display_k_ij_Extensions() const;

    const Model& Initial_Model() const {
      return Initial_Model_;
    }

    const std::vector<int>& Extension_Orders() const {
      return Extension_Orders_;
    }

    const std::vector<BasisAlpha>& Initial_Common_Basis_Alphas() const {
      return Initial_Common_Basis_Alphas_;
    }

    int Large_ST_Dimensions() const {
      return Large_ST_Dimensions_;
    }

    std::vector<bool> Simply_Paired_Layers() const {
      return Simply_Paired_Layers_;
    }

    int Consistent_Models() const {
      return Consistent_Models_;
    }

    int Consistent_BVs() const {
      return Consistent_BVs_;
    }

    int Total_BVs() const {
      return Total_BVs_;
    }

    const std::set<LEEFT>& Unique_LEEFTs() const {
      return Unique_LEEFTs_;
    }

    const std::list<std::vector<std::vector<int> > >& k_ij_Extensions() const {
      return k_ij_Extensions_;
    }

    void Set_Initial_Model(const Model& New_Initial_Model) {
      Initial_Model_ = New_Initial_Model;
    }

    void Set_Extension_Orders(const std::vector<int>& New_Extension_Orders) {
      Extension_Orders_ = New_Extension_Orders;
    }

    void Set_Initial_Common_Basis_Alphas(const std::vector<BasisAlpha>& New_Initial_Common_Basis_Alphas) {
      Initial_Common_Basis_Alphas_ = New_Initial_Common_Basis_Alphas;
    }

    void Set_Large_ST_Dimensions(int New_Large_ST_Dimensions) {
      Large_ST_Dimensions_ = New_Large_ST_Dimensions;
    }

    void Set_Simply_Paired_Layers(std::vector<bool> New_Simply_Paired_Layers) {
      Simply_Paired_Layers_ = New_Simply_Paired_Layers;
    }

    void Set_Consistent_Models(int New_Consistent_Models) {
      Consistent_Models_ = New_Consistent_Models;
    }

    void Set_Consistent_BVs(int New_Consistent_BVs) {
      Consistent_BVs_ = New_Consistent_BVs;
    }

    void Set_Total_BVs(int New_Total_BVs) {
      Total_BVs_ = New_Total_BVs;
    }

    void Set_Unique_LEEFTs(const std::set<LEEFT>& New_Unique_LEEFTs) {
      Unique_LEEFTs_ = New_Unique_LEEFTs;
    }

    void Set_k_ij_Extensions(const std::list<std::vector<std::vector<int> > >& New_k_ij_Extensions) {
      k_ij_Extensions_ = New_k_ij_Extensions;
    }

  protected:
    Model Initial_Model_;
    std::vector<int> Extension_Orders_;
    std::vector<BasisAlpha> Initial_Common_Basis_Alphas_;
    int Large_ST_Dimensions_;
    std::vector<bool> Simply_Paired_Layers_;

    //Members with no accessors. These members are not copied either.
    //TODO: Figure out a way to copy the file streams, and/or initialize them
    //in the copy constructor.
    OutputWriter Data_Writer_;
    std::ofstream Model_Out_;
    std::ofstream MD_Out_;

    unsigned Consistent_Models_;
    unsigned Consistent_BVs_;
    unsigned Total_BVs_;

    std::set<LEEFT> Unique_LEEFTs_;
    std::list<std::vector<std::vector<int> > > k_ij_Extensions_;

    void Build_k_ij_Extensions();
    void Form_k_ij_Layer_Extensions(const std::vector<std::list<std::vector<int> > >& k_ij_Layer_Extensions, int Layer, std::vector<std::vector<int> > k_ij_Extension);
    virtual void Build_Extensions(int Layer, std::vector<BasisVector> Basis_Vectors);
    virtual void Build_Basis_Vectors(int Layer, std::vector<BasisVector> Basis_Vectors);
    virtual void Build_Models(std::vector<BasisVector> Basis_Vectors);
    void Write_End_MD(double Total_Time);
};
