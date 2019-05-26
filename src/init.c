#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(cinc)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(crrf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(crrfit)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(crrfsv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(crrsr)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(crrvv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(crstm)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(tpoi)(void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"cinc",   (DL_FUNC) &F77_NAME(cinc),    7},
    {"crrf",   (DL_FUNC) &F77_NAME(crrf),   16},
    {"crrfit", (DL_FUNC) &F77_NAME(crrfit), 16},
    {"crrfsv", (DL_FUNC) &F77_NAME(crrfsv), 20},
    {"crrsr",  (DL_FUNC) &F77_NAME(crrsr),  17},
    {"crrvv",  (DL_FUNC) &F77_NAME(crrvv),  25},
    {"crstm",  (DL_FUNC) &F77_NAME(crstm),  18},
    {"tpoi",   (DL_FUNC) &F77_NAME(tpoi),    5},
    {NULL, NULL, 0}
};

void R_init_cmprsk(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
