#include <admodel.h>
#undef MAXITER
#define MAXITER  17

#ifndef MSY_REFERENCE_POINTS_H
#define MSY_REFERENCE_POINTS_H


/**
 * This class calculates the MSY-based reference points for the
 * Lagged Recruit Growth Survival Model from the Ecological detective.
 * 
 * @file MSYReferencePoints.h
 * @author Steven Martell <stevem@iphc.int>
 * 
 * 
 */

 /**
  * \brief Calculate MSY-based reference points
  * \author Steve Martell
 **/
class msy_reference_points
{
private:
	double m_k;
	double m_s;
	double m_bo;

	double m_fmsy;
	double m_msy;
	double m_bmsy;

public:
	~msy_reference_points(){}

	// default constructor
	msy_reference_points(const double &k, const double &s, const double &bo)
		:m_k(k), m_s(s), m_bo(bo)
		{
			calcReferencePoints();
		}

	// getters
	double get_fmsy()  { return m_fmsy; }
	double get_bmsy()  { return m_bmsy; }
	double get_msy()   { return m_msy;  }
	// member functions
	void calcReferencePoints();


};

#endif

