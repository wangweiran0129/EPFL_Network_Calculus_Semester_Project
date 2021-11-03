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
#include "linearsegment.h"
#include <iostream>
#include <math.h>
#include <stdio.h>

namespace deborah {

LinearSegment::LinearSegment()
{
	x = 0.0;
	y = 0.0;
	rate = LINEARSEGMENT_POSITIVE_INFINITE;
	leftopen = false;
	id = -1;
	gap = 0.0;
}


LinearSegment::LinearSegment(double x0, double r)
{
	x = x0;
	y = 0.0;
	rate = r;
	leftopen = false;
	id = -1;
	gap = 0.0;
}


LinearSegment::LinearSegment(double x0, double y0, double r)
{
	x = x0;
	y = y0;
	rate = r;
	leftopen = false;
	id = -1;
	gap = 0.0;
}


LinearSegment& LinearSegment::operator=(const LinearSegment &s2)
{
	if(this != &s2)
	{
		x = s2.x;
		y = s2.y;
		rate = s2.rate;
		leftopen = s2.leftopen;
		id = s2.id;
		gap = s2.gap;
	}
	return *this;
}


LinearSegment::LinearSegment(const LinearSegment &ls)
{
	x = ls.x;
	y = ls.y;
	rate = ls.rate;
	leftopen = ls.leftopen;
	id = ls.id;
	gap = ls.gap;
}


LinearSegment::~LinearSegment()
{
}


void LinearSegment::setLatencyRate(double l, double r)
{
	x = l;
	y = 0.0;
	rate = r;
}


void LinearSegment::setAffine(double b, double r)
{
	x = 0.0;
	y = b;
	rate = r;
}


double LinearSegment::value(double x1)
{
	if(fabs(x1 - x) < LINEARSEGMENT_EPSILON) return y;
	return (x1 - x) * rate + y;
}


double LinearSegment::get_X_intersection(LinearSegment &other)
{
	if(other.rate == LINEARSEGMENT_POSITIVE_INFINITE)
		return LINEARSEGMENT_POSITIVE_INFINITE;
	if(fabs(rate - other.rate) < LINEARSEGMENT_EPSILON)
		return LINEARSEGMENT_POSITIVE_INFINITE;
	double y1 = y - x * rate;
	double y2 = other.y - other.x * other.rate;
	return (y2 - y1) / (rate - other.rate);
}


bool LinearSegment::convolve(double x1, double rate1)
{
	if(y != 0.0) return false;
	x += x1;
	if((rate1 < rate) || (rate == LINEARSEGMENT_POSITIVE_INFINITE)) rate = rate1;	
}


bool LinearSegment::convolve(LinearSegment &other)
{
	//supports only latency-rate segments
	if((y != 0.0) || (other.y != 0.0)) return false;
	return convolve(other.x, other.rate);
}

void LinearSegment::add(LinearSegment &s2, double xpos, bool lopen)
{
	x = xpos;
	y = value(xpos) + s2.value(xpos);
	rate = rate + s2.rate;
	leftopen = lopen;
}

void LinearSegment::sub(LinearSegment &s2, double xpos, bool lopen)
{
	x = xpos;
	y = value(xpos) - s2.value(xpos);
	rate = rate - s2.rate;
	leftopen = lopen;
}


void LinearSegment::Zero()
{
	x = y = rate = 0.0;
	leftopen = false;
	id = -1;
	gap = 0.0;
}

void LinearSegment::Print()
{
	printf("(%.6lf, %.6lf)  R = ",x,y);
	if(rate == LINEARSEGMENT_POSITIVE_INFINITE) printf("+INF");
	else if(rate == LINEARSEGMENT_NEGATIVE_INFINITE) printf("-INF");
	else printf("%.6lf",rate);
	//printf("  leftopen=%s  flow_burst_id=%i  gap=%.3lf",(leftopen?"true":"false"),id,gap);
	if(id>=0) printf("   bursting_flow = %02i   burst_size = %.3lf", id, gap);
}

}
