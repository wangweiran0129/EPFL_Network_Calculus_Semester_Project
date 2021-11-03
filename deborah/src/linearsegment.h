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
#ifndef DEBORAHLINEARSEGMENT_H
#define DEBORAHLINEARSEGMENT_H

#include <float.h>

namespace deborah {

/**
This class represent a linear segment of a curve
For instance, it can represent a latency-rate service curve

	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/


#define LINEARSEGMENT_POSITIVE_INFINITE DBL_MAX
#define LINEARSEGMENT_NEGATIVE_INFINITE -DBL_MAX
#define LINEARSEGMENT_EPSILON 0.000000001

class LinearSegment{
public:
    LinearSegment();
	LinearSegment(double x0, double y0, double r);
	LinearSegment(double x0, double r);
	LinearSegment(const LinearSegment &ls);
    ~LinearSegment();
	
	LinearSegment& operator=(const LinearSegment &s2);
	
	void setLatencyRate(double l, double r);
	void setAffine(double b, double r);
	void add(LinearSegment &s2, double xpos, bool leftopen);
	void sub(LinearSegment &s2, double xpos, bool leftopen);
	void Zero();
	
	//Return the value of the function at the given x-coordinate
	double value(double x1);
	
	//Return the x-coordinate of the intersection between this segment and the specified one
	double get_X_intersection(LinearSegment &ls);
	
	//Modifies the current segment by convolving it with the specified one
	bool convolve(LinearSegment &ls);
	bool convolve(double latency, double rate);
	
	double x;
	double y;
	double rate;
	bool leftopen;
	int id;
	double gap;
	
	void Print();

};

}

#endif
