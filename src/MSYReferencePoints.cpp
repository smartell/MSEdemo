/**
 * 
 * @file MSYReferencePoints.cpp
 * @author Steven Martell <stevem@iphc.int>
 * 
 * 
 */

#include <admodel.h>
#include "MSYReferencePoints.h"

/**
 * \brief Calculate MSY-based reference points for LRGS model
 * \author Steve Martell
 * \remarks 
 */
void msy_reference_points::calcReferencePoints()
{
	m_bmsy = (m_bo*(-1.+sqrt(m_k))/(m_k-1.)); 
	m_msy  = (m_bmsy * (-1. + m_s) * (-1. + m_k)*(-m_bo + m_bmsy)) 
			 / (m_bo + (-1. + m_k)*m_bmsy);
	
	m_fmsy = m_msy / m_bmsy;	

	// if( m_msy<m_bmsy )
	// {
	// 	m_fmsy = -log((-m_msy+m_bmsy)/m_bmsy);
	// }
	// else
	// {
	// 	m_fmsy = m_msy/m_bmsy;
	// }

	// cout<<"Bmsy = "<<m_bmsy<<endl;
	// cout<<"Bmsy/Bo = "<<m_bmsy/m_bo<<endl;
	// cout<<m_msy<<endl;
	// cout<<m_fmsy<<endl;	
}