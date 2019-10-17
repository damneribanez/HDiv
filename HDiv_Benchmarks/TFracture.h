//
//  TFracture.h
//  HDiv
//
//  Created by Omar Durán and José Villegas on 3/25/19.
//

#ifndef TFracture_h
#define TFracture_h

#include <stdio.h>
#include "pzreal.h"
#include <set>

class TFracture {
    
public:
    
    std::set<int>  m_id;
    int     m_dim;
    REAL    m_kappa_normal;
    REAL    m_kappa_tangential;
    REAL    m_d_opening;
    REAL    m_porosity;
    
    /// Default constructor
    TFracture();
    
    /// Copy constructor
    TFracture(const TFracture &other);
    
    /// Assignment constructor
    TFracture &operator=(const TFracture &other);
    
};

#endif /* TFracture_h */
