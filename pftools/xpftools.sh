#!/bin/sh
# the next line restarts using wish \
exec wish "$0" "$@"

#BHEADER**********************************************************************
#
#  Copyright (c) 1995-2009, Lawrence Livermore National Security,
#  LLC. Produced at the Lawrence Livermore National Laboratory. Written
#  by the Parflow Team (see the CONTRIBUTORS file)
#  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.
#
#  This file is part of Parflow. For details, see
#  http://www.llnl.gov/casc/parflow
#
#  Please read the COPYRIGHT file or Our Notice and the LICENSE file
#  for the GNU Lesser General Public License.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License (as published
#  by the Free Software Foundation) version 2.1 dated February 1999.
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
#  and conditions of the GNU General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
#  USA
#**********************************************************************EHEADER

#
# Load in the required parflow packages 
#
lappend auto_path $env(PARFLOW_DIR)/bin/

    
if [catch { package require parflow } ] {
    puts "Error: Could not find parflow TCL library"
    exit
}

if [catch { package require xparflow } ] {
    puts "Error: Could not find xparflow TK library"
    exit
}


if [catch { package require vtktcl } ] {
    set XParflow::haveVTK 0
} {
    set XParflow::haveVTK 1
}

#
# To avoid typing import the namespaces from the packages
#
namespace import Parflow::*
namespace import XParflow::*

#
# Start up the main script
#
XParflow::MainWin_Init

