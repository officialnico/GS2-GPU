! Suppress compiler warnings about unused dummy variables.
!
! This is useful for procedures are only functional for certain build
! options, or when a certain signature is required for compatibility
! with a library or other routines.
!
! Tested with most modern compilers. Compiler explorer shows that this
! is entirely elided at any level of optimisation .

#ifndef UNUSED_DUMMY
#define UNUSED_DUMMY(x) if (.false.) print*, shape(x)
#endif
