! --- CARRQR --- ==============================================================
! Communication Avoiding Rank-Revealing QR.
subroutine carrqr(m,n,blkMax,A,ldA,jpvt,tau,info)
   implicit none

   ! .. Scalar Arguments ..
   integer, intent(in) :: m, n, blkMax, ldA
   ! m: rows in A
   ! n: cols in A
   ! blkMax: max columns to permute per sample block
   integer, intent(out) :: info
   ! info: error info

   ! .. Array Arguments ..
   real*8, intent(inout) :: A(ldA,n)
   ! A[in]: matrix to factor
   ! A[out]=[Y \ R]
   integer, intent(out) :: jpvt(n)
   ! jpvt: column pivots;
   real*8, intent(out) :: tau(min(m,n))
   ! tau: reflector coefficients

   ! .. Local Scalars ..
   integer :: k, blk, j, jc
   integer :: nThd, maxThd, nPanel
   real*8 :: dtmp

   ! .. Local work arrays, leading dimensions ..
   real*8, allocatable, dimension(:,:) :: Y, R, Wt, W, WS, A0
   real*8, allocatable, dimension(:) :: nA
   integer, allocatable, dimension(:,:) :: Panels
   ! Wt: T'*Y'*A
   ! WS: workspace with column-based processor affinity.
   integer :: ldY, ldR, ldWt, ldW, ldWS, sdWS, ldPanels, nY, nR
   
   ! .. Parameters ..
   integer, parameter :: numaMax=12, style=1
   real*8, parameter :: negOne=-1.0, zero=0.0, one=1.0
   logical, parameter :: checkError=.false., timeSubroutines=.false.

   ! .. Functions ..
   !integer, external :: omp_get_max_threads
   real*8, external :: qrpError, QR_Error

   ! ===========================
   ! .. Executable Statements ..
   k = min(m,n)

   ! Check args: carrqr(m,n,blkMax,A,ldA,jpvt,tau,info)
   info=0
   if (m<0) then
      info=-1
   elseif (n<0) then
      info=-2
   elseif (blkMax<=0.or.blkMax>k) then
      info=-3
   elseif (ldA<min(m,1)) then
      info=-5
   endif
   if (info/=0) then
      call xerbla('carrqr',-info)
      return
   endif

   ! Quick return if possible.
   if (m==0.or.n==0.or.k==0) return

   ! cores, threads
   !nThd=omp_get_max_threads()
   nThd = 1
   maxThd=nThd
   !write(*, '(a,i)') 'maxThd=', maxThd
   ! array dimensions
   ldY=m
   ldR=blkMax
   ldWt=blkMax ! Wt=T'Y'A
   ldW=n
   ldWS=2*blkMax
   nPanel=blkMax
   sdWS=max(blkMax,nPanel)
   ldPanels=16
   nY=blkMax*maxThd
   nR=2*blkMax*maxThd
   
   ! Allocate local arrays.
   allocate(Y(ldY,nY),R(ldR,nR),Wt(ldWt,nR),W(ldW,blkMax),WS(ldWS,maxThd*sdWS),&
            nA(nR),Panels(ldPanels,maxThd))
   if (checkError) then
      !write(*, '(a)') 'Copying A->A0'
      allocate(A0(ldA,n))
      call dlacpy('A',m,n,A,ldA,A0,ldA)
   endif
   
   ! Construct original permutation order, zero R.
   do j=1,n
      jpvt(j)=j
   enddo

   jc=0 ! completed columns

   do
      blk=min(blkMax,k-jc) ! number of columns to move in this block.
      
      !if (checkError) write(*, '(a)') 'reduceQrcp'
      call reduceQrcp(m,n,jc,blkMax,nThd,A,ldA,jpvt,Y,ldY,R,ldR,Wt,ldWt,nA)

      ! QR factorize new leading bcmp columns of A
      !if (checkError) write(*,'(a)')'qrfA'
      do j=1,blk
         tau(jc+j)=Y(jc+j,j)
      enddo
      call dlacpy('U',blk,blk,R(1,1),ldR,A(jc+1,jc+1),ldA)
      call dlacpy('L',m-jc-1,blk,Y(jc+2,1),ldY,A(jc+2,jc+1),ldA)
      !call dlarft('F','C',m-jc,blk,Y(jc+1,1),ldY,tau(jc+1),Y(jc+1,1),ldY)

      !if (checkError) write(*,'(a)')'updateA'
      call updateA(style,m,n,jc,blk,nPanel,nThd,A,ldA,Y(jc+1,1),ldY,W,ldW,WS,ldWS,sdWS,Panels,ldPanels)

      if (checkError) then
         dtmp=qrpError(m,n,jc+blk,A0,ldA,A,ldA,tau,jpvt)
         !write (*,'(a,es9.2)') 'qr(A) err:', dtmp
      endif

      jc=jc+blk
      if (jc==k) exit

   enddo

   deallocate(Y,R,Wt,W,WS,nA,Panels)
   if (checkError) then
      deallocate(A0)
   endif

end subroutine carrqr
! =============================================================================




! --- REDUCEQRCP --- ==========================================================
subroutine reduceQrcp(m,n,jc,blkMax,nThd,A,ldA,jpvt,Y,ldY,R,ldR,Wt,ldWt,nA)
   implicit none

   ! .. Scalar Arguments ..
   integer, intent(in) :: m, n, jc, blkMax, nThd, ldA, ldY, ldR, ldWt
   ! m, n: rows, cols in A
   ! jc: cols in A complete
   ! blkMax: maximum reflectors to assemble in binary reduction.
   ! Note: This determines the column block size, so we don't want to choose
   !    a smaller blocking factor. The internal subroutine is smart enough to
   !    only return the number of reflectors that are possible.
   ! nThd: number of threads
   
   ! .. Array Arguments ..
   real*8, intent(inout) :: A(ldA,n)
   real*8, intent(out) :: Y(ldY,blkMax*nThd), R(ldR,2*blkMax*nThd),&
                          Wt(ldWt,2*blkMax*nThd), nA(2*blkMax*nThd)

   integer, intent(inout) :: jpvt(n)

   ! .. Local Scalars ..
   integer :: iThd, reduceDist, j1, j2, n1, n2, nSplit, i
   
   !integer, external :: omp_get_thread_num
   ! ===========================
   ! .. Executable Statements ..

   reduceDist=blkMax
   do
      nSplit=1+(n-jc-1)/(2*reduceDist)
      !write(*,'(a,i)') 'nSplit=',nSplit
!$--NOT OMP PARALLEL--- DO SHARED(blkMax,m,n,jc,nSplit,reduceDist,ldA,ldY,ldR,ldWt) PRIVATE(j1,j2,n1,n2,iThd) SCHEDULE(STATIC)
      do i=0,nSplit-1
         j1 = jc + 1 + 2*reduceDist*i
         n1 = min(blkMax,n-j1+1)
         j2 = j1 + reduceDist
         n2 = max(0,min(blkMax,n-j2+1))
         !iThd=omp_get_thread_num()
         iThd = 0
         !write(*, '(a1,i2,a4,i4,a4,i4,a4,i4,a4,i4)')'t',iThd,' j1=',j1,' n1=',n1,' j2=',j2,' n2=',n2
         call splitQrcp(m,n1,n2,jc,A(1,j1),A(1,j2),ldA,&
            jpvt(j1),jpvt(j2),&
            Y(1,1+iThd*blkMax),ldY,&
            R(1,1+(2*iThd)*blkMax),R(1,1+(2*iThd+1)*blkMax),ldR,&
            Wt(1,1+(2*iThd)*blkMax),Wt(1,1+(2*iThd+1)*blkMax),ldWt,&
            nA(1+(2*iThd)*blkMax),nA(1+(2*iThd+1)*blkMax))
      enddo
!$--NOT OMP END PARALLEL --- DO
      reduceDist=reduceDist*2
      if (reduceDist>=n-jc) exit

   enddo

end subroutine reduceQrcp
! =============================================================================






! --- SPLITQRCP --- ===========================================================
! This subroutine is intended to be called with one thread.
! Several calls should be made concurrently and input arguments should be adjusted
! to avoid any simultaneous access to input arrays.
! A1, A2, R1, R2, W1t, W2t, should be separate column blocks. A1 will be filled with pivots.
! Y should have thread affinity for good reuse. Y must be thick enough to store
! size(A1,2) reflectors. The upper triangle of Y will store T for fast future computation.
! To avoid data races, it is important to have leading dimensions rounded up to cachelines.
! We store Wt in column major for better data locality among processes: Wt = [Wt1;Wt2;...]
subroutine splitQrcp(m,n1,n2,jc,A1,A2,ldA,jpvt1,jpvt2,&
                     Y,ldY,R1,R2,ldR,W1t,W2t,ldWt,nA1,nA2)
   implicit none

   ! .. Scalar Arguments ..
   integer, intent(in) :: m, n1, n2, jc, ldA, ldY, ldR, ldWt
   ! m: rows in A
   ! n1, n2: cols in A1, A2
   ! jc: cols in A complete
   
   ! .. Array Arguments ..
   real*8, intent(inout) :: A1(ldA,n1), A2(ldA,n2)
   real*8, intent(out) :: Y(ldY,n1), R1(ldR,n1), R2(ldR,n2),&
                          W1t(ldWt,n1), W2t(ldWt,n2),&
                          nA1(n1), nA2(n2)
   ! nA1, nA2: column norms in A1, A2

   integer, intent(inout) :: jpvt1(n1), jpvt2(n2)
   ! jpvt1, jpvt2: column pivot indices corresponding to A1, A2

   ! .. Local Scalars ..
   integer :: rc, ip, jp, maxInd, maxBlk, k
   ! rc: reflectors complete in Y / rows complete in R
   ! ip: row in process
   ! jp: col in process in A1
   
   ! paneling variables.
   integer :: i
   ! i: temp
   real*8 :: maxNrm, dtmp

   ! .. Parameters ..
   real*8, parameter :: negOne=-1.0_8, zero=0.0_8, one=1.0_8
   logical, parameter :: isCheckSplit=.true.
   ! .. External Functions ..
   real*8 :: dnrm2
   ! ===========================
   ! .. Executable Statements ..

   !write(*, '(a1,i2,a4,i4,a4,i4)')'t',iThd,' n1=',n1,' n2=',n2
   if (n1==0) return

   ! Reset max tracking to invalid state.
   maxInd=0
   maxBlk=0
   maxNrm=negOne
   do i=1,n1
      nA1(i)=dnrm2(m-jc,A1(jc+1,i),1)
      if (nA1(i)>maxNrm) then
         maxInd=i
         maxBlk=1
         maxNrm=nA1(i)
      endif
   enddo
   do i=1,n2
      nA2(i)=dnrm2(m-jc,A2(jc+1,i),1)
      if (nA2(i)>maxNrm) then
         maxInd=i
         maxBlk=2
         maxNrm=nA2(i)
      endif
   enddo

   rc = 0 ! Reflectors completed.
   k = min(n1, m-jc)
   do while (rc<k)
      jp=rc+1 ! Leading column in A1, Y, and row in R1, R2, W1t, W2t in process.
      ip=jc+rc+1 !Leading row in A1, A2, and Y in process.
      !write(*,'(a1,i2,a4,i2,a3,i2,a4,i4,a4,i2)') 't',iThd,' rc=',rc,' k=', k,' jp=',jp, ' ip=',ip
      
      if (maxBlk==1) then
         if (maxInd/=jp) then
            call dswap(m,A1(1,jp),1,A1(1,maxInd),1)
            call dswap(rc,R1(1,jp),1,R1(1,maxInd),1)
            call dswap(rc,W1t(1,jp),1,W1t(1,maxInd),1)
            nA1(maxInd)=nA1(jp)
            i=jpvt1(jp)
            jpvt1(jp)=jpvt1(maxInd)
            jpvt1(maxInd)=i
         endif
      elseif (maxBlk==2) then
         call dswap(m,A1(1,jp),1,A2(1,maxInd),1)
         call dswap(rc,R1(1,jp),1,R2(1,maxInd),1)
         call dswap(rc,W1t(1,jp),1,W2t(1,maxInd),1)
         nA2(maxInd)=nA1(jp)
         i=jpvt1(jp)
         jpvt1(jp)=jpvt2(maxInd)
         jpvt2(maxInd)=i
      else
         write(*,'(a)')'Invalid pivot.'
         stop
      endif

      ! Copy column to Y to construct next reflector.
      call dcopy(m-ip+1,A1(ip,jp),1,Y(ip,jp),1)

      ! Apply prior reflections.
      if (rc>0) then
         call dgemv('N',m-ip+1,rc,negOne,Y(ip,1),ldY,W1t(1,jp),1,one,Y(ip,jp),1)
      endif

      ! Construct reflector. dtmp = tau
      call dlarfg(m-ip+1,Y(ip,jp),Y(ip+1,jp),1,dtmp)
      R1(jp,jp)=Y(ip,jp)
      Y(ip,jp)=one
      ! Store -tau*Y(ip:m,1:rc)'*Y(ip:m,jp) inside Y(1:rc,jp)
      if (rc>0) then
         call dgemv('T',m-ip+1,rc,-dtmp,Y(ip,1),ldY,Y(ip,jp),1,zero,Y(jc+1,jp),1)
      endif

      ! Reflector inner products.
      call dgemv('T',m-ip+1,n1-jp,dtmp,A1(ip,jp+1),ldA,Y(ip,jp),1,zero,W1t(jp,jp+1),ldWt)
      call dgemv('T',m-ip+1,n2,   dtmp,A2(ip,1),   ldA,Y(ip,jp),1,zero,W2t(jp,1),   ldWt)
      if (rc>0) then
         call dgemv('T',rc,n1-jp,one,W1t(1,jp+1),ldWt,Y(jc+1,jp),1,one,W1t(jp,jp+1),ldWt)
         call dgemv('T',rc,n2,   one,W2t(1,1),   ldWt,Y(jc+1,jp),1,one,W2t(jp,1),   ldWt)
      endif

      ! Update next row in R(rPrc,:)=A(iPrc,:)-Y(iPrc,1:rPrc)*Wt(1:rPrc,:).
      call dcopy(n1-jp,A1(ip,jp+1),ldA,R1(jp,jp+1),ldR)
      call dcopy(n2,   A2(ip,1),   ldA,R2(jp,1),   ldR)
      call dgemv('T',jp,n1-jp,negOne,W1t(1,jp+1),ldWt,Y(ip,1),ldY,one,R1(jp,jp+1),ldR)
      call dgemv('T',jp,n2,   negOne,W2t(1,1),   ldWt,Y(ip,1),ldY,one,R2(jp,1),   ldR)

      Y(ip,jp)=dtmp ! Save tau to diagonal for now.
      if (rc>0) then
         call dtrmv('U','N','N',rc,Y(jc+1,1),ldY,Y(jc+1,jp),1) ! Finish processing T.
      endif

      maxInd=0
      maxBlk=0
      maxNrm=negOne
      do i=jp+1,n1
         dtmp=abs(R1(jp,i))/nA1(i)
         nA1(i)=nA1(i)*sqrt(max(zero, (one+dtmp)*(one-dtmp)))
         if (nA1(i)>maxNrm) then
            maxInd=i
            maxBlk=1
            maxNrm=nA1(i)
         endif
      enddo
      do i=1,n2
         dtmp=abs(R2(jp,i))/nA2(i)
         nA2(i)=nA2(i)*sqrt(max(zero, (one+dtmp)*(one-dtmp)))
         if (nA2(i)>maxNrm) then
            maxInd=i
            maxBlk=2
            maxNrm=nA2(i)
         endif
      enddo
      
      rc = rc+1
   enddo 

   !call checkSplit(m,n1,n2,jc,A1,A2,ldA,Y,ldY,R1,R2,ldR,W1t,W2t,ldWt)

end subroutine splitQrcp
! =============================================================================




subroutine checkSplit(m,n1,n2,jc,A1,A2,ldA,Y,ldY,R1,R2,ldR)
   implicit none

   ! .. Scalar Arguments ..
   integer, intent(in) :: m, n1, n2, jc, ldA, ldY, ldR
   ! m: rows in A
   ! n1, n2: cols in A1, A2
   ! jc: cols in A complete
   
   ! .. Array Arguments ..
   real*8, intent(in) :: A1(ldA,n1), A2(ldA,n2)
   real*8, intent(in) :: Y(ldY,n1), R1(ldR,n1), R2(ldR,n2)

   integer :: nE, ldErr, mE, lwork, k, info, i, j
   real*8 :: dtmp, relErr
   real*8, allocatable, dimension(:) :: work, tau
   real*8, allocatable, dimension(:,:) :: Err
   real*8, parameter :: negOne=-1.0_8
   real*8 :: dlange

   ! ===========================
   ! .. Executable Statements ..

   ! Get problem dimensions
   nE = max(n1,n2)
   mE = m-jc
   k = min(n1,mE)
   ldErr = mE   
   call dormqr('L','T',mE,nE,k,Y(jc+1,1),ldY,tau,Err,ldErr,dtmp,-1,info)
   lwork=dtmp

   allocate(Err(ldErr,nE),work(lwork),tau(k))

   ! Copy tau
   do i=1,k
      tau(i)=Y(jc+i,i)
   enddo

   ! Check A1:
   ! Copy A1 to Err
   call dlacpy('A',mE,n1,A1(jc+1,1),ldA,Err,ldErr)
   ! Apply Q'*A1 = R1 =: Err
   call dormqr('L','T',mE,n1,k,Y(jc+1,1),ldY,tau,Err,ldErr,work,lwork,info)
   ! Remove R1
   do j=1,n1
      do i=1,min(j,k)
         Err(i,j) = Err(i,j) - R1(i,j)
      enddo
   enddo

   relErr=dlange('F',mE,n1,Err,ldErr,dtmp)/dlange('F',mE,n1,A1(jc+1,1),ldA,dtmp)
   !write(*, '(a,es)', advance='no') 'R1_re=', relErr

   ! Check A2:
   ! Copy A2 to Err
   call dlacpy('A',mE,n2,A2(jc+1,1),ldA,Err,ldErr)
   ! Apply Q'*A2 = R2 =: Err
   call dormqr('L','T',mE,n2,k,Y(jc+1,1),ldY,tau,Err,ldErr,work,lwork,info)
   ! Remove R2
   do j=1,n2
      do i=1,k
         Err(i,j) = Err(i,j) - R2(i,j)
      enddo
   enddo
   
   relErr=dlange('F',k,n2,Err,ldErr,dtmp)/dlange('F',mE,n2,A2(jc+1,1),ldA,dtmp)
   !write(*, '(a,es)', advance='no') 'R2_re=', relErr

   deallocate(Err,work,tau)

end subroutine checkSplit







