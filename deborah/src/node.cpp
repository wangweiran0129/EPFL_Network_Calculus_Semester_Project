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
#include "node.h"

namespace deborah {

Node::Node()
{
	latency = 0.0;
	rate  = 0.0;
	CDF_i = NULL;
	CAF_i = NULL;
	max_output_rates = NULL;
}

Node::Node(uint src, uint exit, double theta, double R)
{
	latency = theta;
	rate = R;
	max_output_rates = NULL;
}


Node::~Node()
{
	if(CDF_i != NULL) delete [] CDF_i;
	CDF_i = NULL;
	if(CAF_i != NULL) delete [] CAF_i;
	CAF_i = NULL;
	if(max_output_rates != NULL) delete[] max_output_rates;
	max_output_rates = NULL;
}


void Node::clear_curves(uint n_flows)
{
	CDF.Zero();
	CAF.Zero();
	beta.create_latency_rate(latency, rate);
	if(CDF_i != NULL)
		for(int i=0;i<n_flows;i++) CDF_i[i].Zero();
	else CDF_i = new Curve[n_flows];
	if(CAF_i != NULL)
		for(int i=0;i<n_flows;i++) CAF_i[i].Zero();
	else CAF_i = new Curve[n_flows];
}


bool Node::set_output_rates(uint n_flows, double *out_rates, bool rate_red)
{
	if(max_output_rates == NULL)
		max_output_rates = new double [n_flows];
	if(max_output_rates == NULL) return false;
	if(out_rates == NULL)
	{
		for(int i=0; i<n_flows; i++)
			max_output_rates[i] = -1.0;
	}
	else
	{
		for(int i=0; i<n_flows; i++)
			max_output_rates[i] = out_rates[i];
	}
	rate_reduction = rate_red;
	return true;
}


/*	Return the max output rate of the specified flow at this node
	A negative value means that the flow does not interfere with the tagged flow at this node
*/
double Node::get_output_rate(uint flow)
{
	if(max_output_rates == NULL) return 0.0;
	return max_output_rates[flow];
}


bool Node::set_output_rate(uint flow, double rate)
{
	if(max_output_rates == NULL) return false;
	max_output_rates[flow] = rate;
}


bool Node::is_rate_reducing()
{
	return rate_reduction;
}

}
