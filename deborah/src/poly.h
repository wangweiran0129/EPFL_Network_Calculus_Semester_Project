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
#ifndef DEBORAHPOLY_H
#define DEBORAHPOLY_H

#include "types.h"
#include "linearsegment.h"

namespace deborah {

/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/
class Poly{
public:
    Poly(uint num_variables);
	Poly(uint num_variables, double *init_coeff);
	Poly(const Poly &p);
	//Poly(uint num_variables, LinearSegment &s);
    ~Poly();
	
	void Zero();
	bool SetB(double b);
	bool SetCoeff(uint idx, double a);
	bool SumCoeff(uint idx, double a);
	double GetCoeff(uint idx);
	double GetB();
	uint getNumVariables();
	bool isConstant();
	void Negate();
	
	Poly operator+(const Poly &p);
	Poly operator+(double b);
	Poly operator-(const Poly &p);
	Poly operator-(double b);
	Poly operator*(const double f);
	Poly operator/(const double f);
	Poly &operator=(const Poly &p);
	
	void Print();
	void PrintCoefficients();
	
	private:
	double *coeff;
	uint num_coeff;

};

}

#endif
