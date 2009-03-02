/*BHEADER**********************************************************************

  Copyright (c) 1995-2009, Lawrence Livermore National Security,
  LLC. Produced at the Lawrence Livermore National Laboratory. Written
  by the Parflow Team (see the CONTRIBUTORS file)
  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.

  This file is part of Parflow. For details, see
  http://www.llnl.gov/casc/parflow

  Please read the COPYRIGHT file or Our Notice and the LICENSE file
  for the GNU Lesser General Public License.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License (as published
  by the Free Software Foundation) version 2.1 dated February 1999.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
  and conditions of the GNU General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
**********************************************************************EHEADER*/
/******************************************************************************
 *
 * Member functions for the pf_module class.
 *
 *****************************************************************************/

#include "parflow.h"
#include "pf_module.h"


/*--------------------------------------------------------------------------
 * NewPFModule
 *--------------------------------------------------------------------------*/

PFModule  *NewPFModule(call,
		       init_instance_xtra, free_instance_xtra,
		       new_public_xtra, free_public_xtra,
		       sizeof_temp_data,
		       instance_xtra, public_xtra)
void      *call;
void      *init_instance_xtra;
void      *free_instance_xtra;
void      *new_public_xtra;
void      *free_public_xtra;
void      *sizeof_temp_data;
void      *instance_xtra;
void      *public_xtra;
{
    PFModule         *new;

    new = talloc(PFModule, 1);

    (new -> call)               = (void (*)())call;
    (new -> init_instance_xtra) = (void (*)())init_instance_xtra;
    (new -> free_instance_xtra) = (void (*)())free_instance_xtra;
    (new -> new_public_xtra)    = (void (*)())new_public_xtra;
    (new -> free_public_xtra)   = (void (*)())free_public_xtra;
    (new -> sizeof_temp_data)   = (int  (*)())sizeof_temp_data;

    PFModuleInstanceXtra(new)     = instance_xtra;
    PFModulePublicXtra(new)       = public_xtra;

    return new;
}


/*--------------------------------------------------------------------------
 * DupPFModule
 *--------------------------------------------------------------------------*/

PFModule  *DupPFModule(pf_module)
PFModule  *pf_module;
{
    return  NewPFModule((void *)(pf_module -> call),
			(void *)(pf_module -> init_instance_xtra),
			(void *)(pf_module -> free_instance_xtra),
			(void *)(pf_module -> new_public_xtra),
			(void *)(pf_module -> free_public_xtra),
			(void *)(pf_module -> sizeof_temp_data),
			PFModuleInstanceXtra(pf_module),
			PFModulePublicXtra(pf_module));
}


/*--------------------------------------------------------------------------
 * FreePFModule
 *--------------------------------------------------------------------------*/

void            FreePFModule(pf_module)
PFModule       *pf_module;
{
   tfree(pf_module);
}

