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

#ifndef DEBORAH_H
#define DEBORAH_H

/*	Revision history:
 *
 * 		31 Jul 2009 (v0.90):	Added support for LUDB analysis on non-nested tandems
 * 								Minor optimizations and bugfixes in LUDB-related code
 *
 * 		23 Oct 2008 (v0.85):	First public release of Deborah
 */


#define DEBORAH_VERSION 	"0.91.11"

// Global application error codes
#define DEBORAH_ERROR_CONFIG	-1
#define DEBORAH_ERROR_PROV		-2
#define DEBORAH_ERROR_CLI		-3
#define DEBORAH_ERROR_SAFECHECK	-4
#define DEBORAH_ERROR_ALLOC		-5
#define DEBORAH_ERROR_SIZE		-6

#define DEBORAH_ERROR(x)		(x<0)
#define DEBORAH_OK(x)			(x>=0)

#endif
