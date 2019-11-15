#include "bc_util.h"
#include "def_util.h"
#include "problem_bc.h"
#include "grgeometry.h"


BcPatchLoop(bc_struct, ipatch, is, bc_patch_values,
            (i, j, k, ival, bc_struct, ipatch, is),
            ApplyBCPatch(DirichletBC,
                         PROLOGUE({
                             ip = 1;
                           }),
                         EPILOGUE({
                             ip = 2;
                           }),
                         FACE(Left, {
                             op = wp;
                           })),
            ApplyBCPatch(FluxBC,
                         PROLOGUE({
                             ip = 3;
                           }),
                         EPILOGUE({
                             ip = 4;
                           }),
                         FACE(Right, {
                             op = ep;
                           }))
            );
