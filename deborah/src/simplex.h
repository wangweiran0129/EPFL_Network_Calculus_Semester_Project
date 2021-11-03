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

namespace Simplex {
	
	
#define SIMPLEX_RESULT_OK			 0
#define SIMPLEX_RESULT_UNFEASIBLE	-1
#define SIMPLEX_RESULT_UNBOUNDED	-2
#define SIMPLEX_RESULT_UNK_ERROR	-999
	
#define SIMPLEX_ERROR_OK			 0
#define SIMPLEX_ERROR_MEMORY		-1
	
	
	//MMAX = Max number of constraints
	#define  MMAX  320
	
	//NMAX = Max number of variables
	#define  NMAX  63

	extern int test_simplex();
	extern int minimum(Poly &EcFun, ConstraintSet &c_set, const uint dead_var, Poly &Optimum);
	extern void test_es4();
	extern int reset();
}
