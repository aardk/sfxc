!
! COMMON FOR PROGRAMS SKIPPING DBCAT, AND THEREFORE NEEDING ALTERNATE
! WAY TO LOCK DATA BASE CATALOG.  COPIED FROM CTGCM, CATLG'S COMMON.
! KDB 6/27/88
!
      INTEGER*2 ICATLU
      INTEGER*4 CATSTART, CATSIZE
!
      COMMON /CMDSK/ ICATLU, CATSTART, CATSIZE
