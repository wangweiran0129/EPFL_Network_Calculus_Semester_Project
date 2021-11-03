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
#include "cpoly.h"
#include <iostream>
#include <stdio.h>

namespace deborah {

CPoly::CPoly(uint num_variables): inequality(num_variables)
{
	greater_than_zero = true;
	nvars = num_variables;
	tnode_id = CPOLY_TNODEID_NULL;
}

CPoly::CPoly(uint num_variables, uint tnode): inequality(num_variables)
{
	greater_than_zero = true;
	nvars = num_variables;
	tnode_id = tnode;
}

CPoly::CPoly(const CPoly &cp): inequality(cp.nvars)
{
	inequality = cp.inequality;
	greater_than_zero = cp.greater_than_zero;
	nvars = inequality.getNumVariables();
	tnode_id = cp.tnode_id;
}


CPoly::~CPoly()
{
}


CPoly& CPoly::operator=(const CPoly &cp)
{
	if(&cp != this)
	{
		inequality = cp.inequality;
		greater_than_zero = cp.greater_than_zero;
		tnode_id = cp.tnode_id;
	}
	return *this;
}


CPoly& CPoly::operator=(const Poly &p)
{
	if(&p != &inequality)
	{
		inequality = p;
		greater_than_zero = true;
	}
	return *this;
}


double CPoly::GetB()
{
	return inequality.GetB();
}


double CPoly::GetCoeff(uint idx)
{
	return inequality.GetCoeff(idx);
}


bool CPoly::SetB(double b)
{
	return inequality.SetB(b);
}


bool CPoly::SetCoeff(uint idx, double a)
{
	return inequality.SetCoeff(idx,a);
}


bool CPoly::isTautology()
{
	// only one type of check: const >= 0   or  const <= 0
	if(!inequality.isConstant()) return false;
	if (greater_than_zero) return (inequality.GetB() >= 0);
	else return (inequality.GetB() <= 0);
}


bool CPoly::isOxymoron()
{
	// only one type of check: const >= 0   or  const <= 0
	if(!inequality.isConstant()) return false;
	if (greater_than_zero) return (inequality.GetB() < 0);
	else return (inequality.GetB() > 0);
}


// This method flips the direction of the inequality by negating the coefficients
void CPoly::InvertDirection()
{
	inequality.Negate();
	greater_than_zero = !greater_than_zero;
}


/*	Checks whether the constraint is respected by the given point
*/
bool CPoly::isSatisfied(Poly &p)
{
	if(nvars != p.getNumVariables())
	{
		printf("CPoly::isSatisfied(): num_variables mismatch\n");
		return false;
	}
	double val = inequality.GetB();
	for(int i=0; i<nvars; i++)
	{
		val += inequality.GetCoeff(i) * p.GetCoeff(i);
	}
	if(greater_than_zero) return (val >= 0.0);
	else return (val <= 0.0);
}


void CPoly::Print()
{
	inequality.Print();
	if(greater_than_zero) printf(" >= 0");
	else printf(" <= 0");
	printf("  (id=%i)",tnode_id);
}

}
