#ifndef __C_CONNECTOR_H
#define __C_CONNECTOR_H

#pragma once
#ifdef __cplusplus
extern "C"
{
#endif

int MMTbleile(int ncount, int ecount, int *elist, int *ncolors,
                 COLORset **colorclasses, int L, int T, int time_limit, int lb);

#ifdef __cplusplus
}
#endif

#endif
