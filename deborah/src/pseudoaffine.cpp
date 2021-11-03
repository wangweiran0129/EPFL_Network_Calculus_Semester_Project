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
#include "pseudoaffine.h"
#include <iostream>
#include <stdio.h>

namespace deborah {

PseudoAffine::PseudoAffine(uint num_variables): delay(num_variables)
{
	poly_variables = num_variables;
	num_stages = 0;
	//stages.clear();
	stages2.clear();
}


PseudoAffine::~PseudoAffine()
{
	//stages.clear();
	stages2.clear();
}


PseudoAffine PseudoAffine::operator=(PseudoAffine &p)
{
	if(&p != this)
	{
		Zero();
		delay = p.delay;
		num_stages = p.num_stages;
		for(int i=0; i<num_stages; i++)
		{
			//stages.push_back(p.stages[i]);
			stages2.push_back(p.stages2[i]);
		}
	}
	return *this;
}


bool PseudoAffine::MakeLatencyRate(LinearSegment &ls)
{
	//LeakyBucket lb;
	LeakyBucket2 lb2(poly_variables);
	
	//stages.clear();
	stages2.clear();
	delay.Zero();
	delay.SetB(ls.x);
	//lb.sigma = new Poly(poly_variables);
	//lb.rho = ls.rate;
	lb2.rho = ls.rate;
	//stages.push_back(lb);
	stages2.push_back(lb2);
	num_stages = stages2.size();
	return true;
}


/* Convolve the current curve with the passed one:
 * sum delays and add leaky bucket stages
 */
bool PseudoAffine::Convolve(PseudoAffine &p)
{
	LeakyBucket2 newstage(p.delay.getNumVariables());
	delay = delay + p.delay;

	for(int i=0; i<p.num_stages; i++)
	{
		newstage.sigma = p.stages2[i].sigma;
		newstage.rho = p.stages2[i].rho;
		//stages.push_back(p.stages[i]);
		stages2.push_back(newstage);
	}
	num_stages = stages2.size();
	return true;
}


void PseudoAffine::Zero()
{
	//stages.clear();
	stages2.clear();
	num_stages = stages2.size();
	delay.Zero();
}


bool PseudoAffine::AddStage(LeakyBucket2 &lb)
{
	stages2.push_back(lb);
	num_stages = stages2.size();
}


void PseudoAffine::Print()
{
	num_stages = stages2.size();
	printf("%i leaky-bucket stage(s):\n",num_stages);
	printf("Delay = "); delay.Print(); printf("\n");
	for(int i=0;i<num_stages;i++)
	{
		printf("Stage %02i: sigma = ",i);
		stages2[i].sigma.Print();
		printf("\n            rho = %.2lf\n",stages2[i].rho);
	}
}

}
