//
//  TPZCompElDiscScaled.cpp
//  MixedElasticityGirkmann
//
//  Created by Philippe Devloo on 03/05/18.
//

#include "TPZCompelDiscScaled.h"
#include <math.h>


void TPZCompElDiscScaled::ComputeScale() {
    TPZGeoEl * ref = this->Reference();
    fScale = ref->CharacteristicSize();
}

/** @brief Compute shape functions multiplied by scale value*/
void TPZCompElDiscScaled::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                          TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                          REAL &detjac, TPZFMatrix<REAL> &jacinv,
                          TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphix){

    TPZCompElDisc::ComputeShape(intpoint,X,jacobian,axes,detjac,jacinv,phi,dphi,dphix);
    phi*=fScale;
    dphi*=fScale;
    dphix*=fScale;
}

TPZCompElDiscScaled::~TPZCompElDiscScaled()
{
}



