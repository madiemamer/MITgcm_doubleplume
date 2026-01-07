
C     *==========================================================*
C     | DOUBLEPLUME_OPTIONS.h
C     | o CPP options file for ICEFRONT package.
C     *==========================================================*
C     | Use this file for selecting options within the ICEFRONT
C     | package.
C     *==========================================================*

#ifndef DOUBLEPLUME_OPTIONS_H
#define DOUBLEPLUME_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_DOUBLEPLUME
C     Package-specific Options & Macros go here

C     ALLOW FOR DETACHED PLUME - requires different array sizes
#undef DOUBLEPLUME_ALLOW_DETACHED_PLUME

#endif /* ALLOW_DOUBLEPLUME */
#endif /* DOUBLEPLUME_OPTIONS_H */
