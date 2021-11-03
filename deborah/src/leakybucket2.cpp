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
#include "leakybucket2.h"
#include <iostream>

namespace deborah {

LeakyBucket2::LeakyBucket2(uint num_variables): sigma(num_variables)
{
	sigma.Zero();
	rho = 0.0;
	nvars=num_variables;
}


LeakyBucket2::LeakyBucket2(const LeakyBucket2 &lb2): sigma(lb2.nvars)
{
	sigma = lb2.sigma;
	rho = lb2.rho;
	nvars=sigma.getNumVariables();
}


LeakyBucket2::~LeakyBucket2()
{
}


LeakyBucket2& LeakyBucket2::operator=(const LeakyBucket2 &lb2)
{
	if(&lb2 != this)
	{
		sigma = lb2.sigma;
		rho = lb2.rho;
	}
	return *this;
}


}
