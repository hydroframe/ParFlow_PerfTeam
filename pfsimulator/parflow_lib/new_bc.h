#ifndef _NEW_BC_H
#define _NEW_BC_H

#include "bc_util.h"

#define Do_BCContrib(i, j, k, ival, bc_struct, ipatch, is, bc_patch_values, body) \
  ForBCStructNumPatches(bc_struct, ipatch)                                     \
  {                                                                     \
    bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);       \
    {                                                                   \
      switch(BCStructBCType(bc_struct, ipatch))                         \
      {                                                                 \
        body;                                                           \
      }                                                                 \
    }                                                                   \
  }



#endif // _NEW_BC_H
