! --- UPDATEA ---
subroutine updateA(style,m,n,jCmpt,nBlock,nPanel,nThd,A,ldA,T,ldT,W,ldW,WS,ldWS,sdWS,Panels,ldPanels)
   implicit none
   ! A := (I - Y*T'*Y')*A
   ! T is formed from Y and tau

   ! .. Scalar Arguments ..
   integer, intent(in) :: style, m, n, jCmpt, nBlock, nPanel, nThd, ldA, ldT, ldW, ldWS, sdWS, ldPanels
   ! style: paneling style.
   ! m: rows A
   ! n: cols A
   ! nCmpt: cols complete
   ! nBlock: number of reflectors to apply
   ! nPanel: max cols in a panel

   ! .. Array Arguments ..
   real*8, intent(in) :: T(ldT,nBlock)
   ! tau: reflector coefficients
   real*8, intent(inout) :: A(ldA,n)
   ! A: matrix to be updated
   real*8, intent(out) :: W(ldW,nBlock), WS(ldWS,*)
   ! T: Reflector connection matrix workspace
   ! W: (T'Y'A)' for pass back to updateB
   ! WS: work array nBlock by nPanel*nThd.
   integer, intent(out) :: Panels(ldPanels,nThd)

   ! .. Local Scalars ..
   integer :: iThd, j, nP, jWS, i

   ! .. Parameters ..
   real*8, parameter :: negOne=-1.0_8, zero=0.0_8, one=1.0_8

   ! .. External Functions ..
   !integer, external :: omp_get_thread_num, omp_get_max_threads

   ! ===========================
   ! .. Executable Statements ..

   ! Form triangular factor of reflectors.
   ! Since everyone needs to own a copy of Y and T anyway, we might as well use LAPACK.
   ! Is there a benefit to forming T from qrf of leading block?
   !call dlarft('F','C',m-jCmpt,nBlock,A(jCmpt+1,jCmpt+1),ldA,tau(jCmpt+1),T,ldT)

!$--NOT OMP PARALLEL PRIVATE(iThd,j,nP,jWS,i)-- SHARED(style,m,n,jCmpt,nBlock,nPanel,nThd,A,ldA,T,ldT,W,ldW,WS,ldWS,Panels,ldPanels)
   !iThd=omp_get_thread_num()
   iThd = 0
   jWS = iThd*sdWS+1
   ! T nBlock*nBlock. Y (m-jCmpt)*nBlock. WS nBlock*nP. A (m-jCmpt)*nP. USAGE = nP*(nBlock+m-jCmpt) + nBlock*(nBlock+m-jCmpt)
   call getPanels(style,nBlock+m-jCmpt,n,jCmpt+nBlock,nPanel,nBlock*(nBlock+m-jCmpt),iThd,nThd,Panels,ldPanels)

   call firstPanel(Panels,ldPanels,iThd,j,nP)
   do while (nP>0)
      ! PROCESS A(jCmpt+1:m,JP) where JP=j:j+nP-1 covers jCmpt+nBlock+1:n
      ! I. Y'*A -> WS
      ! I.1 Copy WS(1:nBlock,jWS:jWS+nP-1) := A(jCmpt+1:jCmpt+nBlock,j:j+nP-1) prepares Y1' triangle multiply.
      call dlacpy('A',nBlock,nP,A(jCmpt+1,j),ldA,WS(1,jWS),ldWS)
      ! I.2 Triangle multiply WS(1:nBlock,jWS:jWS+nP-1) := Y1'*WS(1:nBlock,jWS:jWS+nP-1)
      !  Y1 = A(jCmpt+1:jCmpt+nBlock,jCmpt+1:jCmpt+nBlock)
      call dtrmm('L','L','T','U',nBlock,nP,one,A(jCmpt+1,jCmpt+1),ldA,WS(1,jWS),ldWS)
      ! I.2 WS(1:nBlock,jWS:jWS+nP-1) += Y2'*A(jCmpt+nBlock+1:m,j:j+nP-1)
      !  Y2 = A(jCmpt+nBlock+1:m,jCmpt+1:jCmpt+nBlock)
      call dgemm('T','N',nBlock,nP,m-jCmpt-nBlock,one,A(jCmpt+nBlock+1,jCmpt+1),ldA,A(jCmpt+nBlock+1,j),ldA,one,WS(1,jWS),ldWS)
      ! II. WS(1:nBlock,jWS:jWS+nP-1) := T'*WS(1:nBlock,jWS:jWS+nP-1)
      call dtrmm('L','U','T','N',nBlock,nP,one,T,ldT,WS(1,jWS),ldWS)
      ! III. Save WS to W=(TtYtA)'
      do i=1,nBlock
         call dcopy(nP,WS(i,jWS),ldWS,W(j,i),1)
      enddo
      !call dlacpy('A',nBlock,nP,WS(1,jWS),ldWS,W(1,j),ldW)
      ! IV. A -= Y*WS
      ! IV.1 A2 -= Y2*WS
      !  A(jCmpt+nBlock+1:m,j:j+nP-1) -= A(jCmpt+nBlock+1:m,jCmpt+1:jCmpt+nBlock)*WS(1:nBlock,jWS:jWS+nP-1)
      call dgemm('N','N',m-jCmpt-nBlock,nP,nBlock,negOne,A(jCmpt+nBlock+1,jCmpt+1),ldA,WS(1,jWS),ldWS,one,A(jCmpt+nBlock+1,j),ldA)
      ! IV.2 WS(1:nBlock,jWS:jWS+nP-1) := Y1*WS(1:nBlock,jWS:jWS+nP-1)
      !  Y1=A(jCmpt+1:jCmpt+nBlock,jCmpt+1:jCmpt+nBlock)
      call dtrmm('L','L','N','U',nBlock,nP,one,A(jCmpt+1,jCmpt+1),ldA,WS(1,jWS),ldWS)
      ! IV.3 A(jCmpt+1:jCmpt+nBlock,j:j+nP-1) -= WS(1:nBlock,jWS:jWS+nP-1)
      do i=0,nP-1
         call daxpy(nBlock,negOne,WS(1,jWS+i),1,A(jCmpt+1,j+i),1)
      enddo

      call nextPanel(Panels,ldPanels,iThd,j,nP)
   enddo
!$--NOT OMP END -- PARALLEL

   ! ===========================
end subroutine updateA
