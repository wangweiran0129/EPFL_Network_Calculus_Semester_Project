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
#include "poly.h"
#include <iostream>
#include <stdio.h>

namespace deborah {

void Poly::Zero()
{
	for(int i=0;i<num_coeff; i++) coeff[i] = 0.0;
}
	
	
Poly::Poly(uint num_variables)
{
	num_coeff = num_variables + 1;
	coeff = new double [num_coeff];
	Zero();
}


Poly::Poly(uint num_variables, double *init_coeff)
{
	num_coeff = num_variables + 1;
	coeff = new double [num_coeff];
	for(int i=0;i<=num_coeff; i++) coeff[i] = init_coeff[i];
}


Poly::Poly(const Poly &p)
{
	num_coeff = p.num_coeff;
	coeff = new double [num_coeff];
	for(int i=0;i<=num_coeff; i++) coeff[i] = p.coeff[i];
}


uint Poly::getNumVariables()
{
	return num_coeff - 1;
}

/*
Poly::Poly(uint num_variables, LinearSegment &ls)
{
	num_coeff = num_variables + 1;
	coeff = new double [num_coeff];
	Zero();
	if(ls.y == 0.0)
		SetB(ls.x);
	else SetB(ls.y);
	
}
*/

bool Poly::SetCoeff(uint idx, double a)
{
	if(idx >= num_coeff) return false;
	coeff[idx] = a;
	return true;
}


bool Poly::SumCoeff(uint idx, double a)
{
	if(idx >= num_coeff) return false;
	coeff[idx] += a;
	return true;
}


double Poly::GetCoeff(uint idx)
{
	if(idx >= num_coeff) return false;
	return coeff[idx];
}


bool Poly::SetB(double b)
{
	return SetCoeff(num_coeff-1, b);
}


double Poly::GetB()
{
	return coeff[num_coeff-1];
}


Poly::~Poly()
{
	if(coeff != NULL) delete [] coeff;
	coeff = NULL;
	num_coeff = 0;
}


bool Poly::isConstant()
{
	bool result = true;
	for(int i=0; i<num_coeff-1; i++)
		if(coeff[i] != 0.0) return false;
	return true;
}


Poly &Poly::operator=(const Poly &p)
{	
	if(&p != this)
	{
		if(coeff != NULL) delete [] coeff;
		num_coeff = p.num_coeff;
		coeff = new double [num_coeff];
		for(int i=0; i<num_coeff; i++)
			coeff[i] = p.coeff[i];
	}	
	return *this;
}


Poly Poly::operator+(const Poly &p)
{
	Poly result(num_coeff-1);
	
	for(int i=0; i<result.num_coeff; i++)
		result.coeff[i] = this->coeff[i] + p.coeff[i];
	
	return result;
}


Poly Poly::operator+(double b)
{
	Poly result(num_coeff-1);
	for(int i=0; i<result.num_coeff; i++)
		result.coeff[i] = coeff[i];
	result.coeff[num_coeff-1] += b;
	return result;
}


Poly Poly::operator-(const Poly &p)
{
	Poly result(num_coeff-1);
	result = *this;
	
	for(int i=0; i<result.num_coeff; i++)
		result.coeff[i] = this->coeff[i] - p.coeff[i];
	
	return result;
}


Poly Poly::operator-(double b)
{
	Poly result(num_coeff-1);
	for(int i=0; i<result.num_coeff; i++)
		result.coeff[i] = coeff[i];
	result.coeff[num_coeff-1] -= b;
	return result;
}


Poly Poly::operator*(const double f)
{
	Poly result(num_coeff-1);
	
	for(int i=0; i<result.num_coeff; i++)
		result.coeff[i] = this->coeff[i] * f;
		
	return result;
}


Poly Poly::operator/(const double f)
{
	Poly result(num_coeff-1);
	
	if(f == 0.0)
	{
		for(int i=0; i<result.num_coeff; i++)
			result.coeff[i] = 999999999999.9999999999999;		
		return result;
	}
	
	for(int i=0; i<result.num_coeff; i++)
		result.coeff[i] = this->coeff[i] / f;
		
	return result;
}


void Poly::Negate()
{
	for(int i=0; i<num_coeff; i++)
		coeff[i] = -coeff[i];
}


void Poly::Print()
{
	if(num_coeff == 0) return;
	for(int i=0;i<num_coeff-1;i++)
	{
		printf("%.2lf*X%i + ",coeff[i],i);
	}
	printf("%.2lf",coeff[num_coeff-1]);
}

void Poly::PrintCoefficients()
{
	if(num_coeff == 0) return;
	for(int i=0;i<num_coeff-1;i++)
	{
		printf("%+.02lf ",coeff[i]);
	}
	printf("| %+.02lf",coeff[num_coeff-1]);
}

}
