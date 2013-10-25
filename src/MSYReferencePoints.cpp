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
	m_bmsy = (m_bo*(sqrt(m_k)-1.0)/(m_k-1.0)); 

	m_msy  = m_bmsy * (m_s-1.0) * (m_k-1.0) * (m_bmsy - m_bo);
	m_msy /= m_bo + (m_k-1.0) * m_bmsy;
	
	m_fmsy = m_msy / m_bmsy;	

	// verified calculations on Oct 24.
}