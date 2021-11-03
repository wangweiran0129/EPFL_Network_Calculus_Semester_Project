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
#ifndef DEBORAHNODE_H
#define DEBORAHNODE_H

#include "types.h"
#include "curve.h"

namespace deborah {

/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/

class Node{
public:
    Node();

    ~Node();
	
	Node(uint src, uint exit, double theta, double R);
		
	double latency;
	double rate;
	
	void clear_curves(uint n_flows);
	
	Curve CDF;
	Curve CAF;
	Curve *CDF_i;
	Curve *CAF_i;
	Curve beta;
	
	bool set_output_rates(uint num_flows, double *out_rates, bool rate_red);
	double get_output_rate(uint flow);
	bool set_output_rate(uint flow, double rate);
	bool is_rate_reducing();
	
private:
	double *max_output_rates;
	bool rate_reduction;	
	
};

}

#endif
