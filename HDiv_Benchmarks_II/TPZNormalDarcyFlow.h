//
//  TPZNormalDarcyFlow.hpp
//  HDiv
//
//  Created by Omar Durán on 4/8/19.
//

#ifndef TPZNormalDarcyFlow_h
#define TPZNormalDarcyFlow_h

#include <stdio.h>
#include "TPZMaterial.h"
#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include "pzvec.h"


class TPZNormalDarcyFlow : public TPZDiscontinuousGalerkin {
    
protected:
    
    /// Material identifier
    int m_mat_id;
    
    /// Material dimension
    int m_dim;
    
    /// Normal permeability
    REAL m_kappa_normal;
    
public:
    
    /// Constuctor based on material identifier and dimension
    TPZNormalDarcyFlow(int mat_id, int dim);
    
    /// @brief Default destructor
    ~TPZNormalDarcyFlow();
    
    /// Copy constructor
    TPZNormalDarcyFlow(const TPZNormalDarcyFlow & other);
    
    /// Copy constructor
    TPZNormalDarcyFlow & operator=(const TPZNormalDarcyFlow & other);
    
    /// Returns problem dimension
    virtual int Dimension() const override { return m_dim; }
    
    /// Sets problem dimension
    virtual void SetDimension(int dim) { m_dim = dim; }
    
    /// Returns number of state variables
    virtual int NStateVariables() const override { return 1; }
    
    /// Creates another material of the same type
    virtual TPZMaterial *NewMaterial() override
    {
        return new TPZNormalDarcyFlow(*this);
    }
    
    /// Sets normal permeability
    void SetKappaNormal(REAL kappa_normal) { m_kappa_normal = kappa_normal; }
    
    /// Gets normal permeability
    REAL GetKappaNormal() { return m_kappa_normal; }
    
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override ;
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override { DebugStop();}
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override { DebugStop();}
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override { DebugStop();}
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)  override{DebugStop();}
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {DebugStop();}
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)  override {DebugStop();}
    
};

#endif /* TPZNormalDarcyFlow_h */
