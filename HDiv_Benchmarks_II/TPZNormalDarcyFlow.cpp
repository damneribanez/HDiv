//
//  TPZNormalDarcyFlow.cpp
//  HDiv
//
//  Created by Omar Durán on 4/8/19.
//

#include "TPZNormalDarcyFlow.h"


TPZNormalDarcyFlow::TPZNormalDarcyFlow(int mat_id, int dim) : TPZDiscontinuousGalerkin(mat_id){
    
    m_mat_id = mat_id;
    m_dim = dim;
    m_kappa_normal = 0.0;
    
}

TPZNormalDarcyFlow::~TPZNormalDarcyFlow(){
    
}

TPZNormalDarcyFlow::TPZNormalDarcyFlow(const TPZNormalDarcyFlow & other) : TPZDiscontinuousGalerkin(other){
    m_mat_id        = other.m_mat_id;
    m_dim           = other.m_dim;
    m_kappa_normal  = other.m_kappa_normal;
}

TPZNormalDarcyFlow & TPZNormalDarcyFlow::operator=(const TPZNormalDarcyFlow & other){
 
    if (this != & other) // prevent self-assignment
    {
        TPZDiscontinuousGalerkin::operator=(other);
        m_mat_id        = other.m_mat_id;
        m_dim           = other.m_dim;
        m_kappa_normal  = other.m_kappa_normal;
    }
    return *this;
}

void TPZNormalDarcyFlow::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    TPZFNMatrix<10,REAL> phi = data.phi;
    
    int n_phi = phi.Rows();
    
    if (n_phi != ek.Rows() && n_phi != ek.Cols()) {
        std::cout << "TPZNormalDarcyFlow:: ek and phi has incompatible number of positions." << std::endl;
    }
    
    for (int i = 0; i < n_phi; i++) {
        for (int j = 0; j < n_phi; j++) {
            ek(i,j) += weight * (1.0/m_kappa_normal) * phi(i,0) * phi(j,0);
        }
    }
}
