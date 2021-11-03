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

#include "types.h"
#include "curve.h"
#include "linearsegment.h"
#include <iostream>
#include <math.h>
#include <stdio.h>

using namespace std;

namespace deborah {

//#define CURVE_DEBUG
//#define CURVE_DEBUG_CONV
//#define CURVE_DEBUG_CHECKS
//#define CURVE_DEBUG_ADD


Curve::Curve()
{
	segments.clear();
}


Curve::Curve(int nsegs)
{
	LinearSegment s;
	s.Zero();
	segments.clear();
	for(int i=0; i<nsegs; i++)
		segments.push_back(s);
}


void Curve::Zero()
{
	segments.clear();
}


int Curve::num_segments() const
{
	return segments.size();
}


bool Curve::add_segment(LinearSegment &s)
{
	LinearSegment ns(s.x, s.y, s.rate);
	ns.leftopen = s.leftopen;
	ns.id = s.id;
	ns.gap = s.gap;
	segments.push_back(ns);
	//if(segments.size()>1000) printf("CURVE IS BIG\n");
	return true;
}


bool Curve::add_segment(double x, double y, double rate)
{
	LinearSegment s(x,y,rate);
	return add_segment(s);
}


bool Curve::create_latency_rate(double l, double r)
{
	LinearSegment s;
	s.Zero();
	Zero();
	if(l == 0.0)
	{
		s.x = 0.0;
		s.y = 0.0;
		s.rate = r;
		return add_segment(s);
	}
	s.x = s.y = s.rate = 0.0;
	add_segment(s);
	s.x = l;
	s.y = 0.0;
	s.rate = r;
	return add_segment(s);
}


bool Curve::create_line_burst(double x1, double x2, double r, double burst, uint flow_id)
{
	LinearSegment s;
	Zero();
	// zero flat segment
	s.Zero();
	add_segment(s);
	// linear segment
	s.x = x1;
	s.y = 0.0;
	s.rate = r;
	add_segment(s);
	// burst segment
	if(burst > 0.0)
	{
		s.y = s.value(x2);
		s.x = x2;
		s.rate = CURVE_POS_INFINITY;
		s.id = flow_id;
		s.gap = burst;
		s.leftopen = true;
		add_segment(s);
	}
	// horizontal segment
	s.y = s.value(x2) + burst;
	s.x = x2;
	s.rate = 0.0;
	s.leftopen = true; // check this???
	s.id = -1;
	s.gap = 0.0;
	return add_segment(s);
}


bool Curve::create_burst_line(double x1, double x2, double r, double burst, uint flow_id)
{
	LinearSegment s;
	segments.clear();
	if(x1 > 0.0)
	{
		// zero flat segment
		s.Zero();
		add_segment(s);
	}
	// burst
	if(burst > 0.0)
	{
		s.id = flow_id;
		s.gap = burst;
		s.rate = CURVE_POS_INFINITY;
		s.x = x1;
		s.y = 0.0;
		s.leftopen = true;
		add_segment(s);
	}
	s.id = -1;
	s.gap = 0.0;
	s.leftopen = false;
	if(IS_GREATER(x2, x1))
	{
		// linear segment at vertical offset
		s.x = x1;
		s.y = burst;
		s.rate = r;
		s.leftopen = true;
		add_segment(s);
		s.y = s.value(x2);
	}
	else s.y = burst;
	// horizontal segment
	s.x = x2;
	s.rate = 0.0;
	s.leftopen = false; // check this???
	if((x1==0.0) && (x2==0.0)) s.leftopen=true;
	s.id = -1;
	return add_segment(s);
}


bool Curve::create_token_bucket(double b, double r)
{
	LinearSegment s;
	segments.clear();
	s.Zero();
	add_segment(s);
	s.x = 0.0;
	s.y = b;
	s.rate = r;
	s.leftopen = true; // left-continuos
	return add_segment(s);
}


Curve& Curve::operator=(const Curve &c2)
{
	if(&c2 != this)
	{
		int ns2 = c2.num_segments();
		segments.clear();
		for(int i=0; i<ns2; i++)
		{
			segments.push_back(c2.segments[i]);
		}
	}
	return *this;
}


Curve::~Curve()
{
	segments.clear();
	//printf("Curve() destructor\n");
}


void Curve::Print()
{
	double dy;
	int size = segments.size();
	printf("Curve contains %i segment(s):\n",size);
	for(int i=0; i<size; i++)
	{
		printf("   ");
		segments[i].Print();
		dy = segments[i].y - ((i==0)?0.0:segments[i-1].value(segments[i].x));
		if(dy > CURVE_EPSILON) printf("  dy=%lf",dy);
		printf("\n");
	}
}


/* Returns the index of the LinearSegment to which x-coordinate x1 belongs
 *
 */
int Curve::getSegmentDefining(double x1)
{
	int size = segments.size();

	for(int i=size-1; i>=0; i--)
	{
		if(segments[i].leftopen)
		{
			if(segments[i].x < x1) return i;
		}
		else {
			if(segments[i].x <= x1) return i;
		}
	}
	return -1;
}


/*	Find the first segment which defines the y-coordinate y1
	Return the segment index, or a negative value if the curve is not defined for y1
*/
int Curve::getSegmentFirstAtValue(double y1)
{
	int size = segments.size();

	if(size == 0) return -1;
	// check if the curve starts already above y1 (and quit, since it's monotonic)
	if(segments[0].y > y1 + LINEARSEGMENT_EPSILON) return -1;
	for (int i = 0; i < size-1; i++)
	{
		if (segments[i+1].y >= y1+LINEARSEGMENT_EPSILON)
			return i;
	}
	// handle the last segment separately
	if(segments[size-1].y >= y1-LINEARSEGMENT_EPSILON) return size-1;
	// finally check whether the curve diverges (but we expect it to never happen!)
	if(segments[size-1].rate > LINEARSEGMENT_EPSILON) return size-1;
	return -1;
}


/* Returns the value of the curve at the x-coordinate x1
 */
double Curve::f(double x1)
{
	int i = getSegmentDefining(x1);
	if (i < 0) return CURVE_NEG_INFINITY;
	return (x1 - segments[i].x) * segments[i].rate + segments[i].y;
}


double Curve::f(double x1, bool right_limit)
{
	int size = num_segments();
	double x;
	if(right_limit)
	{
		int s = size-1;
		while(s>=0)
		{
			x = segments[s].x;
			if(IS_EQUAL(x,x1)) return segments[s].y;
			if(IS_SMALLER(x,x1))
			{
				if(IS_ZERO(segments[s].rate)) return segments[s].y;
				else return (x1 - x) * segments[s].rate + segments[s].y;
			}
			else s--;
		}
		return CURVE_NEG_INFINITY; // undefined in x1!!!
	}
	else {
		int s = 0;
		x = segments[0].x;
		if(IS_GREATER(x,x1)) return CURVE_NEG_INFINITY; // undefined in x1
		while(s<size)
		{
			x = segments[s].x;
			if(IS_EQUAL(x,x1)) return segments[s].y;
			if(IS_GREATER(x,x1))
			{
				if(IS_ZERO(segments[s-1].rate)) return segments[s-1].y;
				else return (x1 - segments[s-1].x) * segments[s-1].rate + segments[s-1].y;
			}
			else s++;
		}
		s = size-1;
		if(IS_ZERO(segments[s].rate)) return segments[s].y;
		else return (x1 - segments[s].x) * segments[s].rate + segments[s].y;
	}
}


/* Returns the x value at which the function value is equal to y1.
 * If rightmost is true then it returns the rightmost x-coordinate,
 * otherwise the leftmost coordinate.
 */
double Curve::f_inv2(double y1, bool right_limit)
{
	int size = num_segments();
	if(size < 1) return CURVE_NEG_INFINITY;
	int s = 0;
	while(s < size)
	{
		if(IS_GREATER(y1,segments[s].y)) s++;
		else break;
	}
	if(s == 0) return CURVE_NEG_INFINITY; // undefined for y=y1
	if(right_limit) while((s < size) && IS_EQUAL(y1,segments[s].y)) s++;
	if(s == size)
	{
		if(IS_ZERO(segments[s-1].rate)) return CURVE_POS_INFINITY; // undefined for y=y1
		else return segments[s-1].x + (y1-segments[s-1].y) / segments[s-1].rate;
	}
	if(right_limit) if(IS_EQUAL(y1,segments[s-1].y)) return segments[s-1].x;
	else if(IS_EQUAL(y1,segments[s].y)) return segments[s].x;
	return segments[s-1].x + (y1-segments[s-1].y) / segments[s-1].rate;
}

double Curve::f_inv(double y1, bool rightmost)
{
	int size = segments.size();
	int i = getSegmentFirstAtValue(y1);

	if (i < 0) {
		//error: segment never assumes value y1
		return CURVE_NEG_INFINITY;
	}

	if (rightmost) {
		while(i < size && segments[i].rate == 0.0) {
			i++;
		}
		if (i >= size) {
			return CURVE_POS_INFINITY;
		}
	}
	if (segments[i].rate != 0.0) {
		return segments[i].x + (y1 - segments[i].y) / segments[i].rate;
	}
	else return segments[i].x;
}



double Curve::getLatency()
{
	int nsegs = segments.size();
	if(nsegs<=0)
	{
		printf("Curve::getLatency() called on empty curve!!!\n");
		return 0.0;
	}

	if (segments[0].y > 0.0) {
		return 0.0;
	}
	for (int i = 0; i < nsegs; i++) {
		double y0 = segments[i].y;
		if (y0 < 0.0 && y0 > -CURVE_EPSILON) {
			y0 = 0.0;
		}
		if (y0 > 0.0 || ((y0 >= 0.0) && (segments[i].rate > 0.0)))
		{
			return segments[i].x;
		}
		if ((y0 < 0.0) || (segments[i].rate < 0.0))
		{
			printf("Curve::getLatency(): NEGATIVE rate or y-coordinate somewhere!!!\n");
			return CURVE_NEG_INFINITY;
		}
	}
	return CURVE_POS_INFINITY;
}


bool Curve::removeLatency()
{
	Curve result = *this;
	int nsegs;

	// printf("removeLatency(): input "); result.Print();
	nsegs = result.num_segments();
	while(nsegs > 0)
	{
		if ((result.segments[0].y > CURVE_EPSILON) || (result.segments[0].rate > CURVE_EPSILON))
		{
			break;
		}
		if ((result.segments[0].y < 0.0) || (result.segments[0].rate < 0.0)) {
			printf("Curve::removeLatency(): NEGATIVE rate or y-coordinate somewhere!!!\n");
			return false;
		}
		//removeSegment(0);
		result.segments.erase(result.segments.begin());
		nsegs--;
	}

	// in case that we've removed everything, the curve had infinite
	// latency, so return the NULL curve.
	nsegs = result.num_segments();
	if (nsegs <= 0) {
		segments.clear();
		printf("Curve::removeLatency(): removed all segments, bug somewhere!\n");
		return false;
	}
	// shift remaining segments left by latency
	double L = result.segments[0].x;
	result.segments[0].x = 0.0;
	for (int i = 1; i < nsegs; i++)
	{
		result.segments[i].x -= L;
	}
	*this = result;
	return true;
}


bool Curve::ConvolveWithLatencyRate(Curve &beta)
{
	Curve tmp1, tmp2;
	double L1 = getLatency();
	double L2 = beta.getLatency();
	bool result;
	//double K = 0.0;
	if (IS_POS_INFINITY(L1))
	{
		*this = beta;
		return true;
	}
	if (IS_POS_INFINITY(L2))
	{
		return true;
	}

	tmp1= *this;
	tmp2 = beta;
	//K = tmp1.f(0); // initial offset
	//if(K > 0.0) tmp1.shiftUp(-K);
	tmp1.removeLatency(); // the cumulative arrival curve
	tmp2.removeLatency(); // the node's service curve

#ifdef VERBOSE
	printf("Convolution: A = "); tmp1.Print();
	printf("Convolution: Beta = "); tmp2.Print();
#endif

	// perform the convolution
	result = tmp1.convolution(tmp2);

#ifdef VERBOSE
	printf("Convolution: RESULT = "); tmp1.Print();
#endif

	tmp1.shiftRight(L1+L2);
	//if(K > 0.0) tmp1.shiftUp(K);
	//tmp1.beautify();
	*this = tmp1;

	remove_infinitesimal_segments();
	remove_colinear_segments();

	return result;
}


/*	Shift the curve to the right of the specified amount
*/
bool Curve::shiftRight(double L)
{
	LinearSegment ls;
	int i;
	int nsegs = segments.size();

	if(nsegs <= 0) return false;

	i = 0;
	while(i<nsegs)
	{
		segments[i].x += L;
		if(segments[i].x >= 0.0)
		{
			i++;
			continue;
		}
		segments.erase(segments.begin() + i);
		nsegs--;
	}
	if(segments[0].x == 0.0) return true;
	ls.Zero();
	segments.insert(segments.begin(),ls);
	return true;
}



/*	Check if the curve has a break-point with an X-coordinate very close to the specified one
	Returns the segment index (>=0) if there is such a point, -1 otherwise
*/
int Curve::has_BPX_at(double x)
{
	int ns = segments.size();
	for(int i=0; i<ns; i++)
	{
		if(IS_EQUAL(segments[i].x, x)) return i;
	}
	return -1;
}


/*	Check if the curve has a break-point with an Y-coordinate very close to the specified one
	Returns the segment index (>=0) if there is such a point, -1 otherwise
*/
int Curve::has_BPY_at(double y)
{
	int ns = segments.size();
	for(int i=0; i<ns; i++)
	{
		if(IS_EQUAL(segments[i].y, y)) return i;
	}
	return -1;
}


/*	This function inserts a breakpoint at the specified y-coordinate.
	If the curve already has a breakpoint at the specified point, no action is taken.
	Input: y = the y-coordinate of the point
	Output: an integer containing the index of the segment added, or a negative value
*/
int Curve::add_BPY(double y)
{
	LinearSegment ns;
	int i, bpy;
	double x;

	bpy = has_BPY_at(y);
	if(bpy >= 0) return -1;
	x = f_inv(y, true);
	i = getSegmentDefining(x);
	if(i < 0) return CURVE_ERROR;
	ns.x = x;
	ns.y = y;
	ns.rate = segments[i].rate;
	ns.gap = 0.0;
	ns.id = -1;
	ns.leftopen = false;
	segments.insert(segments.begin()+i+1, ns);
	return i+1;
}


/*	This function inserts a breakpoint at the specified x-coordinate.
	If the curve already has a breakpoint at the specified point, no action is taken.
	Input: x = the x-coordinate of the point
	Output: an integer containing the index of the segment added, or a negative value
*/
int Curve::add_BPX(double x)
{
	LinearSegment ns;
	int i,bpx;
	double y;

	bpx = has_BPX_at(x);
	if(bpx >= 0) return -1;
	i = getSegmentDefining(x);
	if(i < 0) return CURVE_ERROR;
	ns.x = x;
	ns.y = f(ns.x, false); // the curve is surely continuous here
	ns.rate = segments[i].rate;
	ns.gap = 0.0;
	ns.id = -1;
	ns.leftopen = false;
	segments.insert(segments.begin()+i+1, ns);
	return i+1;
}


/*	This method recomputes the slopes of the segments in order to enforce the linear
	continuity of the curve; it leaves x and y breakpoint coordinates intact.
	No approximations here!
*/
bool Curve::recompute_slopes()
{
	int nsegs = num_segments();
	if(nsegs < 2) return false;

	for(int i=0; i<nsegs-1; i++)
	{
		double delta_y = segments[i+1].y - segments[i].y;
		double delta_x = segments[i+1].x - segments[i].x;
		if(delta_x == 0.0) segments[i].rate = CURVE_POS_INFINITY;
		else segments[i].rate = delta_y / delta_x;
	}
	segments[nsegs-1].rate = 0.0;
	return true;
}


/*	This method checks for discontinuities at the specified x0.
	If there is at least one, it returns the segment index of the discontinuity and
	it stores the number of consecutive gaps in "num_gaps"
*/
int Curve::find_gaps_at(double x0, int &num_gaps)
{
	int size = num_segments();
	int s0 = has_BPX_at(x0);
	if((s0<0) || (s0>=size))
	{
		num_gaps = 0;
		return -1;
	}
	int s = s0;
	while(s < size)
	{
		if(IS_POS_INFINITY(segments[s].rate)) s++;
		else break;
	}
	num_gaps = s-s0;
	return (num_gaps > 0) ? s0 : -1;
}


/*	This methods sums the current curve to the specified one.
	In DEBORAH it is mostly used to compute the cumulative CAF at the node's ingress
*/
bool Curve::add(Curve &c2)
{
	Curve result;
	int nsegs1 = num_segments();
	int nsegs2 = c2.num_segments();
	int ngaps1, ngaps2, ngaps;
	int sg1, sg2;
	double x1, x2, x, y;
	LinearSegment ns;
	int size;

	result.Zero();
	std::vector<double> result_bp;

	// first some checks for trivial situations: empty curves.
	if(nsegs2 == 0) return true;
	if(nsegs1 == 0)
	{
		*this = c2;
		return true;
	}

#ifdef CURVE_DEBUG_ADD
	printf("----------- add: entering\n");
	printf("add(): c1 = "); Print();
	printf("add(): c2 = "); c2.Print();
#endif

	// compute the union of the unique x-coordinates of the breakpoints
	x = -1.0;
	for(int i=0; i<nsegs2; i++)
	{
		x2 = c2.segments[i].x;
		if(IS_EQUAL(x2,x)) continue;
		result_bp.push_back(x2);
		x = x2;
	}
	size = result_bp.size();
	for(int i=0; i<nsegs1; i++)
	{
		x1 = segments[i].x;
		// check if x1 is already in result_bp
		int j = 0;
		while(j < size)
		{
			x = result_bp[j];
			if(IS_EQUAL(x1,x))
			{
				// snap breakpoint and then ignore it
				segments[i].x = x1;
				break;
			}
			j++;
		}
		if(j >= size)
		{
			// ordered insert into breakpoints list
			int k=0;
			while(k<size) if(result_bp[k]>x1) break;
			else k++;
			if(k>=size) result_bp.push_back(x1);
			else result_bp.insert(result_bp.begin()+k,x1);
			size++;
		}
	}
	// the common breakpoints between c1 and c2 have been snapped.
	// now the two curves either have totally distinct breakpoints, or perfectly coincident.
	size = result_bp.size();
#ifdef CURVE_DEBUG_ADD
	printf("add(): unique breakpoint x-coordinates are:  [");
	for(int i=0; i<size; i++) printf(" %lf ",result_bp[i]);
	printf("]\n");
#endif

	// now compute the sum
	for(int i=0; i<size; i++)
	{
		ns.x = result_bp[i];
		//printf("x=%lf   f1(left)=%lf   f2(left)=%lf   f1(right)=%lf   f2(right)=%lf\n",ns.x,f(ns.x,false),c2.f(ns.x,false),f(ns.x,true),c2.f(ns.x,true));
		ns.y = f(ns.x, false) + c2.f(ns.x,false);
		sg1 = find_gaps_at(ns.x, ngaps1);
		sg2 = c2.find_gaps_at(ns.x, ngaps2);
		if((ngaps1 == 0) && (ngaps2 == 0))
		{
			// both curves are continuous in x0
			ns.gap = 0.0;
			ns.id = -1;
			result.add_segment(ns);
			continue;
		}
		// handle discontinuities
		while(ngaps1 > 0)
		{
			ns.id = segments[sg1].id;
			ns.gap = segments[sg1].gap;
			result.add_segment(ns);
			ns.y += ns.gap;
			ngaps1--;
			sg1++;
		}
		while(ngaps2 > 0)
		{
			ns.id = c2.segments[sg2].id;
			ns.gap = c2.segments[sg2].gap;
			result.add_segment(ns);
			ns.y += ns.gap;
			ngaps2--;
			sg2++;
		}
		// add the post-gap segment
		ns.gap = 0.0;
		ns.id = -1;
		ns.y = f(ns.x,true) + c2.f(ns.x, true);
		result.add_segment(ns);
	}

	// recompute slopes
	result.recompute_slopes();

#ifdef CURVE_DEBUG_ADD
	printf("add(): result "); result.Print();
	printf("----------- add: leaving\n");
#endif
	*this = result;
	return true;
}


double Curve::getFirstBit(uint flow)
{
	if(segments.size()==0) return 0.0;
	return getLatency();
}


double Curve::getLastBit(uint flow)
{
	//TODO CHECK if it's correct
	int nsegs = segments.size();
	if(nsegs==0) return 0.0;
	return segments[nsegs-1].x;
}


/*	This method performs the convolution between the current curve A (non-convex, non-concave) and a rate-latency curve
	If we are at a point where beta's slope is greater than A's, then we must copy A until we reach a break-point
	where A's slope becomes greater than beta's (this includes a case of discontinuity, where A's slope becomes +INF)
	If we are at a point where beta's slope is less than A's, then we must copy beta (centered on the start point) until
	we reach an intersection between A and beta (where beta's slope has necessarily become greater than A's)
*/
bool Curve::convolution(Curve &beta)
{
	Curve cdf;
	LinearSegment scdf, b;
	double x_cross, x_next,x_cross_invalid;
	bool beta_below;
	int cseg;
	int ngaps;

	int nsegs = num_segments();

	if(beta.num_segments() != 1)
	{
		printf("convolution(): beta must be sigle-segment!\n");
		return false;
	}
	b = beta.segments[0];
	if(nsegs <= 0)
	{
		printf("convolution(): CAF is empty!\n");
		return false;
	}

	// smooth CAF (remove colinear segments and recompute rates)
	remove_colinear_segments();
	// adjust segments which are colinear with beta (preserve Y's)
	for(int i=0; i<nsegs; i++)
	{
		if(IS_EQUAL(segments[i].rate, b.rate))
			segments[i].rate = b.rate;
		if(IS_ZERO(segments[i].rate))
			continue;
		if(i >= nsegs-1) break;
		if(IS_POS_INFINITY(segments[i].rate)) continue;
		segments[i+1].x = segments[i].x + (segments[i+1].y - segments[i].y) / segments[i].rate;
		//segments[i+1].y = segments[i].y + segments[i].rate * (segments[i+1].x - segments[i].x);
	}

	b.x = segments[0].x; b.y = segments[0].y;
	beta_below = (b.rate < segments[0].rate);
	cseg = 0;
	ngaps = 0;

#ifdef CURVE_DEBUG_CONV
	printf("convolution(): beautified CAF = ");
	Print();
	printf("Initial status is 'beta %s CAF'\n",(beta_below?"below":"above"));
#endif

	while(cseg < nsegs)
	{
		if(beta_below)
		{
			// count discontinuities
			if(IS_POS_INFINITY(segments[cseg].rate))
			{
				ngaps++;
				cseg++;
				continue;
			}

			// find intersection between beta and CAF
			double y1 = segments[cseg].y - segments[cseg].x * segments[cseg].rate;
			double y2 = b.y - b.x * b.rate;
			// discontinuity in CAF?
			if(IS_POS_INFINITY(segments[cseg].rate))
			{
				cseg++;
				continue;
			}
			// parallel segments?
			if(segments[cseg].rate == b.rate)
			{
				if(b.y == segments[cseg].y) x_cross = segments[cseg].x;
				else {
					cseg++;
					continue;
				}
			}
			else x_cross = (y2 - y1) / (segments[cseg].rate - b.rate);
			x_next = (cseg < nsegs-1) ? segments[cseg+1].x : CURVE_POS_INFINITY;
#ifdef CURVE_DEBUG_CONV
			printf("test intersection between seg#%i = ",cseg);segments[cseg].Print();
			printf(" and beta = "); b.Print(); printf(" --> x_cross = %lf (%s)\n",x_cross,(x_cross==segments[cseg].x)?"==":"!=");
#endif
			// ignore intersections at the breakpoint-x == beta.x
			if(IS_EQUAL(x_cross,segments[cseg].x) && (b.x==segments[cseg].x)) //if(IS_EQUAL(x_cross,segments[cseg].x))
				x_cross = CURVE_POS_INFINITY;
			// snap to next breakpoint-x
			if(cseg < nsegs-1)
			{
				if(IS_EQUAL(x_cross,x_next))
					x_cross = x_next;
			}
			// check x_cross within segment's range
			if((x_cross < segments[cseg].x - CURVE_EPSILON) || (x_cross > x_next + CURVE_EPSILON))
			{
				cseg++;
				continue;
			}
#ifdef CURVE_DEBUG_CONV
			printf("convolution(): found intersection at x_cross=%.6lf\n",x_cross);
#endif
			// we have a valid intersection at this segment's x_cross coordinate
			// write beta segment
			cdf.add_segment(b);
#ifdef CURVE_DEBUG_CONV
			printf("write beta: "); b.Print(); printf("\n");
#endif
			b.id = -1;
			b.gap = 0.0;
			beta_below = false;
			// center beta at intersection point
			b.y = b.y + (x_cross - b.x) * b.rate;
			b.x = x_cross;
		}
		// x_cross = x-coordinate of intersection
		// b is now centered in intersection point
		// beta is above -- scan CAF for slope higher than beta's
		if(segments[cseg].rate <= b.rate) //if(segments[cseg].rate < b.rate)
		{
			if(!IS_POS_INFINITY(x_cross))
			{
				scdf.x = x_cross;
				x_cross = CURVE_POS_INFINITY;
				scdf.y = b.y;
			}
			else {
				scdf.x = segments[cseg].x;
				scdf.y = segments[cseg].y;
			}
			scdf.rate = segments[cseg].rate;
			cdf.add_segment(scdf);
#ifdef CURVE_DEBUG_CONV
			printf("copy CAF segment: "); scdf.Print(); printf("\n");
#endif
			cseg++;
			continue;
		}
		// here the CAF slope is grater than beta's
		b.x = segments[cseg].x;
		b.y = segments[cseg].y;
		beta_below = true;
#ifdef CURVE_DEBUG_CONV
		printf("convolution(): found CAF slope change at x=%.6lf (CAF slope=%.3lf, segment %i)\n",b.x,segments[cseg].rate,cseg);
#endif
	}
	// we should never exit with beta_below = true, otherwise CAF was unlimited!
#ifndef CURVE_DEBUG_CHECKS
	if(beta_below) printf("convolution(): warning: exiting with beta below the CAF (unlimited CAF or bug?).\n");
#endif

	// post-process the CDF, inserting breakpoints to keep track of the CAF discontinuities
	/*
	if(ngaps > 0)
	{
		int k;
		for(int i=0; i<nsegs; i++)
		{
			if(!IS_POS_INFINITY(segments[i].rate)) continue;
			// CAF has a discontinuity here
			k = cdf.add_BPY(segments[i].y + segments[i].gap);
			k = cdf.add_BPY(segments[i].y);
			if(k == CURVE_ERROR) continue;
			if(k < 0) k = -k;
			cdf.segments[k].gap = segments[i].gap;
			cdf.segments[k].id = segments[i].id;
#ifdef CURVE_DEBUG_CONV
			printf("+ add_breakpoint(y=%lf) --> segment %i:  ",segments[i].y,k); cdf.segments[k].Print(); printf("\n");
#endif
		}
	}
	*/

#ifdef CURVE_DEBUG_CONV
	printf("convolution(): CDF = "); cdf.Print();
#endif
	*this = cdf;
	return true;
}


LinearSegment& Curve::getSeg(uint s)
{
	int ns = num_segments();
	if(s>=ns) return segments[ns-1];
	return segments[s];
}


// Find the slope of the curve at x_start, right or left limit
double Curve::getSlope(double x_start, bool right)
{
	int ns = num_segments();
	if(ns == 0)
	{
		printf("Curve::getSlope(): ERROR, this curve is empty!\n");
		return 0.0;
	}
#ifdef CURVE_DEBUG
	printf("Curve::getSlope() at x=%.3lf (%s-limit)\n",x_start,(right?"right":"left"));
#endif
	if(right)
	{
		while(ns >= 1)
		{
			if(IS_POS_INFINITY(segments[ns].rate))
			{
				ns--;
				continue;
			}
			if((segments[ns].x > x_start) && (segments[ns-1].x > x_start))
			{
				ns--;
				continue;
			}
			return segments[ns-1].rate;
		}
		return segments[0].rate;
	}
	else
		for(int i=ns-1; i>=0; i--)
		{
			if(IS_POS_INFINITY(segments[i].rate)) continue;
			if(segments[i].x > x_start) continue;
			return segments[i].rate;
		}
#ifdef CURVE_DEBUG
	printf("Curve::getSlope(): ERROR, x-coordinate %.3lf NOT FOUND\n",x_start);
#endif
	return 0.0;
}


bool Curve::remove_colinear_segments()
{
	int size = num_segments();
	if(size < 2) return false;
	for(int i=0; i<size-1; i++)
	{
		if(IS_POS_INFINITY(segments[i+1].rate)) continue;
		if(!IS_EQUAL(segments[i].rate, segments[i+1].rate)) continue;
		segments.erase(segments.begin()+i+1);
		size--;
		i--;
	}
	return recompute_slopes();
}


bool Curve::remove_infinitesimal_segments()
{
	int size = num_segments();
	if(size < 2) return false;
	for(int i=0; i<size-1; i++)
	{
		if(IS_ZERO(segments[i+1].x-segments[i].x))
		{
			if(IS_ZERO(segments[i+1].y-segments[i].y))
			{
				//printf("remove infinitesimal segment: ");segments[i+1].Print();printf("\n");
				segments.erase(segments.begin()+i+1);
				size--;
				i--;
			}
		}
	}
	recompute_slopes();
	size = num_segments();
	for(int i=0; i<size; i++)
	{
		if(IS_POS_INFINITY(segments[i].rate) && (segments[i].id < 0))
		{
			segments.erase(segments.begin()+i);
			size--;
			i--;
		}
	}
	return true;
}


/*	This function returns the total amount of traffic transmitted by the flow
*/
double Curve::getTotalTraffic()
{
	int ns = segments.size();
	if(ns <= 0) return 0.0;
	if(fabs(segments[ns-1].rate) < CURVE_EPSILON)
		return segments[ns-1].y;
	return CURVE_POS_INFINITY;
}


/*	This method performs some basic tests to check if the curve contains nonsense data
*/
bool Curve::SanityCheck()
{
	int ns = segments.size();

	if(ns == 0) return true;

	// check for final rate = 0.0
	if(fabs(segments[ns-1].rate) <= CURVE_EPSILON)
		segments[ns-1].rate = 0.0;
	else {
#ifdef CURVE_DEBUG_CHECKS
		printf("warning: curve ends with non-null rate\n");
#endif
	}

	// check for vertical continuity
	int i = 0;
	double dx, dy;
	while(i < ns-1)
	{
		dy = segments[i+1].y - segments[i].y;
		dx = segments[i+1].x - segments[i].x;
		if(IS_POS_INFINITY(segments[i].rate))
		{
			if(dy == segments[i].gap) {
				i++;
				continue;
			}
			break;
		}
		if(IS_ZERO(dx))
		{
#ifdef CURVE_DEBUG_CHECKS
			printf("warning: X-distance between segment %i and %i is %s\n",i,i+1,((dx==0.0)?"zero":"epsilon"));
#endif
		}
		if(IS_ZERO(dy))
		{
#ifdef CURVE_DEBUG_CHECKS
			printf("warning: Y-distance between segment %i and %i is %s\n",i,i+1,((dy==0.0)?"zero":"epsilon"));
#endif
		}

		if(segments[i].rate * dx != dy)
		{
			if(fabs(segments[i].rate * dx - dy) < CURVE_EPSILON)
				segments[i+1].y = segments[i].y + dx * segments[i].rate;
			else break;
		}
		i++;
	}
	if(i == ns-1) return true;
#ifdef CURVE_DEBUG_CHECKS
	printf("Curve :: sanity check fails at segment %i: dx=%.6lf  *  rate=%.6lf  !=  dy=%.6lf\n",i,dx,segments[i].rate,dy);
#endif
	return false;
}


bool Curve::move_burst_id(uint flow_id, bool forward)
{
	int size = num_segments();
	if(size <= 0) return false;
	int i = 0;
	// find the flow's burst transmission
	while(i < size-1)
	{
		if(IS_POS_INFINITY(segments[i].rate) && (segments[i].id == flow_id))
			break;
		i++;
	}
	if(i >= size-1) return false;
	int j = i;
	while(j < size-1)
	{
		if(IS_POS_INFINITY(segments[j+1].rate)) j++;
		else break;
	}
	if(j == i) return false;
	LinearSegment ts = segments[j];
	segments[j] = segments[i];
	segments[i] = ts;
	return true;
}


}
