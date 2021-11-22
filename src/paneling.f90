! === PANELING ===
! Get paneling parameters for parallel computation.
! Style 1: Continuous equally divided remainder.
! Style 2: Block cyclic
!
! INVARIENTS:
! If this thread participates, nFirst>0. Else nFirst==0
! If I have columns left: nP > 0
!
! Loop prototype:
!
! nThd=omp_get_max_threads()
! !$OMP PARALLEL PRIVATE(iThd,j,nP) SHARED(style,n,jCmpt,nPanel,nThd,Panels,ldPanels)
! iThd=omp_get_thread_num()
! call getPanels(style,n,jCmpt,nPanel,iThd,nThd,Panels,ldPanels)
! call firstPanel(Panels,ldPanels, iThd, j, nP)
! do while (nP>0)
!    ! PROCESS j:j+nP-1
!    call nextPanel(Panels,ldPanels, iThd, j, nP)
! enddo
! !$OMP END PARALLEL
! --- GETPANELS ---
subroutine getPanels(style,mPanel,n,jCmpt,nPanel,nUsed,iThd,nThd,Panels,ldPanels)
   implicit none

   ! .. Scalar Arguments ..
   integer, intent(in) :: style, mPanel, n, jCmpt, nPanel, nUsed, iThd, nThd, ldPanels
   ! style: 1 - Equally divided continuous blocking
   !        2 - Block cyclic
   ! mPanel: Rows referenced in panel
   ! n: Total columns in matrix
   ! jCmpt: Last complete / processed column
   ! nPanel: Column limit for one panel
   ! nUsed: Fixed number of doubles used.
   ! iThd: This thread index, counting from 0 to nThd-1
   ! nThd: Total number of threads

   ! .. Array Arguments ..
   integer, intent(out) :: Panels(ldPanels,nThd)

   ! .. Local Scalars ..
   integer :: ind, jFirst, nFirst, nSkip, jLast, nThdUsed, nP, nC
   ! jFirst: First column this thread will process.
   ! nFirst: Number of columns this thread will process in next block
   ! nSkip: Columns skipped when moving to next block
   ! jLast: Last column this thread could process. It does not have to be inside block that will be processed.
   
   ! .. Parameters ..
   integer, parameter :: nGroup=12, cacheMB=30,nPanelMin=24 !Megabytes of highest level of cache shared among all processors in group.
   real*8, parameter :: cacheUse=0.6 !Fraction of cache to fill.
   logical, parameter :: checkError=.false.

   ! ===========================
   ! .. Executable Statements ..

   if (style==1) then ! Continuous equally divided remainder

      nP = nPanel !Use nPanel suggestion.

      ! Remove processors if they don't have enough work
      nThdUsed=min(1+(n-jCmpt-1)/nP,nThd)

      jFirst=jCmpt+1+(iThd*(n-jCmpt))/nThdUsed
      jLast=jCmpt+((iThd+1)*(n-jCmpt))/nThdUsed
      if (jFirst>n) then !This thread does not participate in paneling
         jFirst=n+1
         jLast=n
         nFirst=0
      elseif (jLast>n) then !This shouldn't happen
         !write(0,'(a)') 'Error: jFirst<=n, jLast>n'
      else
         nFirst=min(nP,jLast-jFirst+1)
      endif
      nSkip=0

   elseif (style==2) then

      nP = nPanel !Use nPanel suggestion.

      jFirst=1+((jCmpt/(nThd*nP))*nThd+iThd)*nP ! First column in this superblock.
      nFirst=nP
      nSkip=nP*(nThd-1)
      jLast=n

      if (jFirst+nP-1<=jCmpt) then ! This entire panel in this superblock is already complete
         jFirst=jFirst+nFirst+nSkip ! Move to next superblock
      elseif (jFirst<=jCmpt) then ! Complete columns must be removed from this panel.
         nFirst=jFirst+nP-1-jCmpt
         ! jFirst+nP-1 > jCmpt   -> nFirst > 0
         ! jFirst <= jCmpt           -> nFirst <= nP-1
         jFirst=jCmpt+1
      endif

      if (jFirst>jLast) then ! This thread does not participate
         jFirst=n+1
         nFirst=0
      else
         nFirst=min(nFirst,jLast-jFirst+1) !Truncate first panel to fit remaining columns
      endif

   elseif (style==3) then

      !    fraction utilized        dbls/MB  divided among prc in group
      nC = cacheUse*cacheMB*131072/min(nThd,nGroup)!Cache available.
      ! Get maximum panel size to fill cache.
      nP = (max(nPanelMin,min(nPanel,(nC-nUsed)/mPanel))/8)*8 ! Don't let panels get too small or large. Watch sdWS.
      ! Remove processors if they don't have enough work
      nThdUsed=min(1+(n-jCmpt-1)/nP,nThd)

      if (iThd==0.and.checkError) then
         !write (*,'(a,i5,a1,i3)')'Panel=',mPanel,'*',nP
      endif

      jFirst=jCmpt+1+(iThd*(n-jCmpt))/nThdUsed
      jLast=jCmpt+((iThd+1)*(n-jCmpt))/nThdUsed
      if (jFirst>n) then !This thread does not participate in paneling
         jFirst=n+1
         jLast=n
         nFirst=0
      elseif (jLast>n) then !This shouldn't happen
         !write(0,'(a)') 'Error: jFirst<=n, jLast>n'
      else
         nFirst=min(nP,jLast-jFirst+1)
      endif
      nSkip=0


   endif
   
   ! INVARIENTS:
   ! If this thread participates, nFirst>0. Else nFirst==0

   ind=iThd+1 ! Processor column of panelling data.
   Panels(1,ind)=jFirst !First column to process
   Panels(2,ind)=nFirst !Number of columns in first panel. Allows for partial panel to be complete.
   Panels(3,ind)=nP !Number of columns in subsequent full panels.
   Panels(4,ind)=nSkip  !Number of columns to skip between panels - Used for block cyclic.
   Panels(5,ind)=jLast  !Last column to process for this processor.
   Panels(6,ind)=n      !Total columns in original matrix.

end subroutine getPanels



! --- FIRSTPANEL ---
! INVARIENTS:
! If I have columns left: nP > 0
subroutine firstPanel(Panels,ldPanels,iThd,j,nP)
   implicit none

   ! .. Scalar Arguements ..
   integer, intent(in) :: ldPanels,iThd
   integer, intent(out) :: j, nP

   ! .. Array Arguements ..
   integer, intent(in) :: Panels(ldPanels,*)

   ! .. Local Scalars ..
   integer :: ind

   ind=iThd+1
   j=Panels(1,ind) !jFirst
   nP=Panels(2,ind) !nFirst

end subroutine firstPanel



! --- NEXTPANEL ---
! INVARIENTS:
! If I have columns left: nP > 0
subroutine nextPanel(Panels,ldPanels,iThd,j,nP)
   implicit none

   ! .. Scalar Arguements ..
   integer, intent(in) :: ldPanels,iThd
   integer, intent(inout) :: j, nP

   ! .. Array Arguements ..
   integer, intent(in) :: Panels(ldPanels,*)

   ! .. Local Scalars ..
   integer :: ind

   if (nP==0) return !This thread is finished. Quick exit.

   ind=iThd+1

   if (nP<0) then
      !write(0,'(a)') 'Error: nP<0'
      !write(0,'(a,i,a,i,a,i,a,i,a,i,a,i,a,i)') 'iThd=', iThd,' jFirst=', Panels(1,ind),' nFirst=', Panels(2,ind),' nPanel=', Panels(3,ind),' nSkip=', Panels(4,ind),' jLast=', Panels(5,ind),' n=', Panels(6,ind)
      !write(0,'(a,i,a,i,a,i)') 'nP=', nP,' j=', j
      stop
   endif

   j=j+nP+Panels(4,ind) !nSkip=Panels(4,ind)
   if (j>Panels(5,ind)) then !jLast=Panels(5,ind)
      j=Panels(6,ind)+1 !n=Panels(6,ind);
      nP=0
      return
   endif

   nP=min(Panels(5,ind)-j+1, Panels(3,ind)) !jLast=Panels(5,ind), nPanel=Panels(3,ind);

end subroutine nextPanel



! --- REMOVEPROCESSED ---
! INVARIENTS:
! If I have columns left: nP > 0
subroutine removeProcessed(jPrc,iThd,Panels,ldPanels)
   implicit none

   ! .. Scalar Arguements ..
   integer, intent(in) :: jPrc, iThd, ldPanels

   ! .. Array Arguements ..
   integer, intent(inout) :: Panels(ldPanels,*)

   ! .. Local Scalars ..
   integer :: ind, jFirst, nFirst, nPanel, nSkip, jLast, n

   ind=iThd+1
   jFirst=Panels(1,ind)
   nFirst=Panels(2,ind)
   nPanel=Panels(3,ind)
   nSkip=Panels(4,ind)
   jLast=Panels(5,ind)
   n=Panels(6,ind)

   if (jFirst>jPrc) return !This thread paneling schedule doesn't change.

   if (jFirst<jPrc) then
      !write(0,'(a)') 'Error: jFirst<jPrc'
      !write(0,'(a,i,a,i,a,i,a,i,a,i,a,i,a,i)') 'iThd=', iThd,' jFirst=', Panels(1,ind),' nFirst=', Panels(2,ind),' nPanel=', Panels(3,ind),' nSkip=', Panels(4,ind),' jLast=', Panels(5,ind),' n=', Panels(6,ind)
      !write(0,'(a,i)') 'jPrc=', jPrc
      stop
   endif

   ! jFirst==jPrc. This thread owns processed column.

   if (nFirst<=0) then
      !write(0,'(a)') 'Error: nFirst<=0'
      !write(0,'(a,i,a,i,a,i,a,i,a,i,a,i,a,i)') 'iThd=', iThd,' jFirst=', Panels(1,ind),' nFirst=', Panels(2,ind),' nPanel=', Panels(3,ind),' nSkip=', Panels(4,ind),' jLast=', Panels(5,ind),' n=', Panels(6,ind)
      !write(0,'(a,i)') 'jPrc=', jPrc
      stop
   endif

   ! nFirst>0

   jFirst=jFirst+1
   nFirst=nFirst-1

   if (nFirst==0) then
      jFirst=jFirst+nSkip
      if (jFirst>jLast) then
         jFirst=n+1
         nFirst=0
      else
         nFirst=min(jLast-jFirst+1, nPanel)
      endif
   endif

   Panels(1,ind)=jFirst
   Panels(2,ind)=nFirst

end subroutine removeProcessed
