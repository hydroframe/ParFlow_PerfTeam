#ifndef _BC_RICHARDS_H
#define _BC_RICHARDS_H

#define Do_Richards_SymmContrib(i, j, k, ival, bc_struct, ipatch, is, body) \
  ForBCStructNumPatches(ipatch, bc_struct)                              \
  {                                                                     \
    body;                                                               \
  }

#define Do_Richards_BuildJCMatrix(i, j, k, ival, bc_struct, ipatch, is, body) \
  ForBCStructNumPatches(ipatch, bc_struct)                              \
  {                                                                     \
    body;                                                               \
  }

#endif // _BC_RICHARDS_H
