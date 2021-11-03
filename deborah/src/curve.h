/***************************************************************************
 *   Copyright (C) 2008 by Luca Bisti   *
 *   luca.bisti@iet.unipi.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef DEBORAHCURVE_H
#define DEBORAHCURVE_H

#include <vector>
#include "linearsegment.h"

namespace deborah {

// Some defines for function return values
#define CURVE_ERROR	(-2147483647-1)

#define CURVE_EPSILON		0.000001
#define CURVE_POS_INFINITY	DBL_MAX
#define CURVE_NEG_INFINITY	-DBL_MAX

#define IS_POS_INFINITY(x)	((x)==CURVE_POS_INFINITY)
#define IS_NEG_INFINITY(x)	((x)==CURVE_NEG_INFINITY)
// Macro for testing is x ~= 0.0
#define IS_ZERO(x)	(fabs(x)<CURVE_EPSILON)
// Macro for testing if x ~= y
#define IS_EQUAL(x,y)	(fabs((x)-(y))<CURVE_EPSILON)
// Macro for testing if x > y
#define IS_GREATER(x,y)	((x)>(y)+CURVE_EPSILON)
#define IS_GREATER_OR_EQUAL(x,y) ((x)+CURVE_EPSILON>=(y))
// Macro for testing if x < y
#define IS_SMALLER(x,y)	((x)+CURVE_EPSILON<(y))
#define IS_SMALLER_OR_EQUAL(x,y) ((x)<=(y)+CURVE_EPSILON)
	
	
	
/**
This class implements a piecewise linear curve object and some common operations
that can be performed on it, such as convolution, summation, etc.
The curve is assumed to be monotone non-decreasing.

	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/

class Curve{
public:

    	Curve();
	Curve(int nsegs);
    	~Curve();
	
	int num_segments() const;
	void Zero();
	bool add(Curve &c2); //TODO: rewrite
	bool add_segment(LinearSegment &s);
	bool add_segment(double x, double y, double rate);
	bool create_latency_rate(double latency, double rate);
	bool create_line_burst(double x1, double x2, double rate, double burst, uint flow_id);
	bool create_burst_line(double x1, double x2, double rate, double burst, uint flow_id);
	bool create_token_bucket(double burst, double rate);
		
	double f(double x1);
	double f(double x1, bool right_limit);
	double f_inv(double y1, bool rightmost);
	double f_inv2(double y1, bool rightmost);
	int getSegmentDefining(double x);
		
	double getFirstBit(uint flow);
	double getLastBit(uint flow);
	bool ConvolveWithLatencyRate(Curve &beta);
	
	double getLatency();
	bool removeLatency();
	bool shiftRight(double L);
	LinearSegment &getSeg(uint s);
	double getSlope(double x_start, bool right);
	double getTotalTraffic();
	int has_BPX_at(double x);
	int find_gaps_at(double x, int &num_gaps);
	int has_BPY_at(double y);
	int add_BPX(double x);
	int add_BPY(double y);
	
	bool recompute_slopes();
	bool remove_colinear_segments();
	bool remove_infinitesimal_segments();
	bool move_burst_id(uint flow_id, bool forward);
	
	Curve& operator=(const Curve &c2);
	
	void Print();
	bool SanityCheck();
	
private:
	std::vector<LinearSegment> segments;
	int getSegmentFirstAtValue(double y1);
	bool convolution(Curve &c2);
};

}

#endif
