//
//  TFracture.cpp
//  HDiv
//
//  Created by Omar Durán and José Villegas on 3/25/19.
//

#include "TFracture.h"


TFracture::TFracture() : m_id(), m_dim(0), m_kappa_normal(0), m_kappa_tangential(0), m_d_opening(0), m_porosity()
{
    
}

TFracture::TFracture(const TFracture &other) : m_id(other.m_id), m_dim(other.m_dim), m_kappa_normal(other.m_kappa_normal), m_kappa_tangential(other.m_kappa_tangential), m_d_opening(other.m_d_opening), m_porosity(other.m_porosity)

{
    
}

TFracture & TFracture::operator=(const TFracture &other)
{
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_id = other.m_id;
    m_dim = other.m_dim;
    m_kappa_normal = other.m_kappa_normal;
    m_kappa_tangential = other.m_kappa_tangential;
    m_d_opening = other.m_d_opening;
    m_porosity  = other.m_porosity;
    
    return *this;
    }
