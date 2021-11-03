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
#ifndef DEBORAHPSEUDOAFFINE_H
#define DEBORAHPSEUDOAFFINE_H

#include "types.h"
#include "poly.h"
#include "linearsegment.h"
#include <vector>
#include "leakybucket2.h"

using namespace deborah;


namespace deborah {
	
	/*
	typedef struct LeakyBucket {
		Poly *sigma;
		double rho;
	} LeakyBucket;
	*/
	
/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/
class PseudoAffine{
	public:
    PseudoAffine(uint variables);
	bool MakeLatencyRate(LinearSegment &ls);
	bool Convolve(PseudoAffine &p);
	void Zero();
	bool AddStage(LeakyBucket2 &lb);

	PseudoAffine operator=(PseudoAffine &p);
	
    ~PseudoAffine();
	void Print();
	
	Poly delay;
	uint num_stages;
	//std::vector<LeakyBucket> stages; // leaky bucket stages
	std::vector<LeakyBucket2> stages2; // leaky bucket stages

	private:
	uint poly_variables;

};

}

#endif
