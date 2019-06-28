***********************************************************************
*     This program fits a2F to the optical memory function
*     Energy is in units of eV. 
*     The chi-squared function is weighted with sma**2 factors for 
*     each data point.
*     From the file param_in.dat it reads:
*     temperature in kelvin,  lower and upper limit of Pi(w) in eV
*     max number of blocks, first number of blocks, ifitab, ilump
*     the integer parameter ifitab:
*     with ifitab=0 the program fits Pi(w)
*     with ifitab=1 the program fits Pi(w) + correction to epsinf
*     with ifitab=2 Pi(w) is fitted + corrections to epsinf and wp
*     ilump is the number of x-values over which the input is averaged
***********************************************************************
      program fitdat
*          
      parameter(nw=2000,ne=200,mp=99)
      double precision x(nw)
      complex*16 y(nw),s(nw)
      integer imx
      common/data/x,y,s,imx
*
      double precision fe(nw,ne),few(nw,ne)
      complex*16 ll(nw,ne,mp),llw(nw,ne,mp)      
      integer ifunc,nfunc
      common/funcarr/fe,few,ll,llw,ifunc,nfunc
*                   
      double precision pm(mp)
      real ch0
*      
      complex*16 cmem,cdmdp(mp),yy,sg,cm
      character*3 datin
*
      integer m,mmx,i,ilump,iter
      double precision pmold(mp),wmx,t,kb,m1,m2,h2,a,b,xx,xs,y1,y2,s1,s2
*      
      double precision w0,dw
      common/inblock/w0,dw
*      
      integer ifitab
      double precision kbt,tpi,a0,b0
      common/memoir/kbt,tpi,a0,b0,ifitab           
*----------------------------------------------------------------------                                                   
      tpi=8*atan(1.)
***** 1 kelvin = 86.18 micro eV ******      
      kb=86.18e-6
*
***** load initial parameters
      open(10,file='param_in.txt')
      read(10,*) t,w0,wmx
      read(10,*) mmx,m,ifitab,ilump
      read(10,*) a0,b0
      dw=(wmx-w0)/m
      do 20 j=1,m
       read(10,*) h2
       pm(j)=dabs(h2)
20    continue
*     ifitab=0 fit a2F
*     ifitab=1 fit a2F and correct epsinf
*     ifitab=2 fit a2F; correct epsinf and wp
      close(10)
      kbt=kb*t
*   
      a=a0 
      b=b0  
*
      write(*,*) 'initial parameters have been read'
*
***** load experimental data
      write(*,*) 'temperature in kelvin (3 digits)'
      read(*,*) datin
      write(*,*) 'temperature '
      read(*,*) it
      kbt=kb*it
      open(15,file='mem'//datin//'.dat')
      do 5 i=1,100000                        
       read(15,*,end=6) xx,y1,y2,s1,s2
       if (xx.gt.0.7) goto 6  
5     continue
6     imx=i-1
      close(15)            
      open(25,file='mem'//datin//'.dat')
      imx=imx/ilump      
      do 15 i=1,imx
       xs=0.
       yy=(0.,0.)
       sg=(0.,0.)
       do 13 j=1,ilump                       
        read(25,*,end=14) xx,y1,y2,s1,s2
        xs=xs+xx
        yy=yy+cmplx(y1,y2)
        sg=sg+cmplx(s1,s2)
c        xs=xs+xx
c        yy=yy+cmplx(y1,y2)
c        sg=sg+cmplx(s1,s2)        
13     continue     
14     x(i)=xs/ilump
       y(i)=yy/ilump
       s(i)=0.0001*sg/real(ilump)**1.5
c       write(*,*) i,x(i),y(i),s(i)
15    continue
      close(25)     
*
***** start the Levenberg-Marquardt fitting
80    write(*,*) 'fit with ',m,' blocks.'       
      mtot=m+ifitab 
      pm(m+2)=a
      pm(m+1)=b 
      iter=m*8
      ifunc=0
      nfunc=imx
      call dofit(imx,pm,ch0,mtot,iter)
      a=pm(m+2)
      b=pm(m+1)
c
      write(*,*) 'writing output '
*
      open(15,file='fit'//datin//'.dat')
      do 45 i=1,imx  
       call myfunc(i,pm,cmem,cdmdp,mtot)
       xx=x(i)
       yy = cmem
*       yy=1/(a/(xx+y(i))+b*xx)-xx
*       cm=1/(a/(xx+cmem)+b*xx)-xx
       write(15,*) xx,dble(yy),dimag(yy)
45    continue       
      close(15)
*
      open(16,file='glue'//datin//'.dat')
      write(16,*) w0,0
      do 55, j=1,m
       w1=w0+dw*(j-1)
       w2=w1+dw
       write(16,*) w1,dabs(pm(j))
       write(16,*) w2,dabs(pm(j))       
55    continue  
      write(16,*) w2,0
      close(16)
*
      open(17,file='param'//datin//'out.txt')
      write(17,*) t,w0,wmx
      write(17,*) mmx,m,ifitab,ilump
      write(17,*) a,b,ch0     
      do 65, j=1,m
       write(17,*) dabs(pm(j))
65    continue
      close(17)
*   
      do 90 j=1,m
       pmold(j)=pm(j)
90    continue
      a=pm(m+2)
      b=pm(m+1)      
      do 92 j=1,m
       pm(2*j-1)=pmold(j)
       pm(2*j)=pmold(j)
92    continue
      pm(2*m+2)=a
      pm(2*m+1)=b
      if ((m*2).gt.mmx) goto 100
      m=m*2
      dw=(wmx-w0)/m
      goto 80        
*   
100   end
***********************************************************************

***********************************************************************
***** Fit the data
***********************************************************************
      subroutine dofit(imx,pm,ch0,m,iter)
*           
      parameter(mp=99)
      integer imx,iter
      double precision pm(mp),chs
      real ch0,chlm
      double precision cvr(mp,mp),alf(mp,mp),bet(mp)            
*      
      integer i,flag
      double precision almd,alamol,ochs
*---------------------------------------------------------------------- 
      almd=-1.
      flag=0
      alamol=almd
      chlm=0.7
*
      do 45 i=1,iter  
       call mrqmin(pm,cvr,alf,bet,m,chs,ochs,almd)  
       ch0=(real(chs)/imx)**0.5       
       write(*,*) i,'  almd=',almd
     * ,' <chi(i)^2/sigma(i)^2)^2>^1/2=',ch0
      if (almd.gt.alamol) then
       flag=flag+1
      else
       flag=0
      endif
      alamol=almd
      if (ch0.lt.chlm) then
       write(*,*) 'Iteration interrupted. Reason: '
       write(*,*) ch0,' = (chs/imx)^0.5 dropped below limit = ',chlm
       goto 46 
      endif
      if (flag.gt.6) then
       write(*,*) 'Iteration interrupted. Reason: '
       write(*,*) almd,' = almd increased 6 times '
       goto 46
      endif
45    continue
*
      write(*,*) 'Iteration interrupted. Reason: '
      write(*,*) 'number of iterations equals ', i
46    almd=0
      call mrqmin(pm,cvr,alf,bet,m,chs,ochs,almd)
60    continue
      return
      end
*********************************************************************

*********************************************************************
      SUBROUTINE MRQMIN(pm,cvr,alf,bet,m,chs,ochs,almd)     
      PARAMETER(mp=99)
      double precision pm(mp),cvr(mp,mp),alf(mp,mp),bet(mp)
      double precision chs,ochs,almd      
      integer m      
      integer j,k
      double precision pmtry(mp),dpm(mp)
*
      integer ifitab
      double precision kbt,tpi,a0,b0
      common/memoir/kbt,tpi,a0,b0,ifitab           
*
      IF(almd.LT.0.)THEN
        almd=0.001
        CALL MRQCOF(pm,alf,bet,m,chs)
        ochs=chs
      ENDIF
      DO 11 J=1,m
        DO 10 K=1,m
          cvr(J,K)=alf(J,K)
10      CONTINUE
        cvr(J,J)=alf(J,J)*(1.+almd)
        dpm(J)=bet(J)
11    CONTINUE
      CALL GAUSSJ(cvr,dpm,m)
* here we make sure that the next try stays within the boundary 
* value that blockheight > 0      
      DO 15 J=1,m-ifitab
       if ((pm(j)+dpm(j)).lt.-0.01) then
        pmtry(j)=0.01
        dpm(j)=pmtry(j)-pm(j)
       else
        pmtry(j)=pm(j)+dpm(j)     
       endif
15    continue 
*
       DO 16 J=1+m-ifitab,m
        pmtry(j)=pm(j)+dpm(j) 
16    CONTINUE
*
      CALL MRQCOF(pmtry,cvr,dpm,m,chs) 
      IF(chs.LT.ochs)THEN
        almd=0.1*almd
        ochs=chs
        DO 18 J=1,m
          DO 17 K=1,m
            alf(J,K)=cvr(J,K)
17        CONTINUE
          bet(J)=dpm(J)
          pm(j)=pmtry(j)
18      CONTINUE
      ELSE
        almd=10.*almd
        chs=ochs
      ENDIF            
      RETURN
      END
*********************************************************************

*********************************************************************
      SUBROUTINE MRQCOF(pm,alf,bet,m,chs)
*      
      parameter(nw=2000)
      double precision x(nw)
      complex*16 y(nw),s(nw)
      integer n
      common/data/x,y,s,n 
*           
      parameter(mp=99)
      double precision pm(mp),alf(mp,mp),bet(mp),chs
      integer m
*     
      integer i,j,k
      double precision sinv1,sinv2,dy1,dy2,wt1,wt2
      complex*16 dydpm(mp),ymod 
*
      DO 12 J=1,m
        DO 11 K=1,J
          alf(J,K)=0.
11      CONTINUE
        bet(J)=0.
12    CONTINUE
      chs=0.
*      
      DO 15 i=1,n
       CALL myfunc(i,pm,YMOD,dydpm,m)   
        sinv1=1./(dble(s(i))**2)
        sinv2=1./(dimag(s(i))**2)                      
        DY1=dble(Y(I)-YMOD)
        DY2=dimag(Y(I)-YMOD)        
        DO 14 J=1,M
          WT1=dble(dydpm(j))*sinv1    
          WT2=dimag(dydpm(j))*sinv2  
          DO 13 K=1,J
           alf(J,K)=alf(J,K)+WT1*dble(dydpm(k))+WT2*dimag(dydpm(k))
13        CONTINUE
          bet(J)=bet(J)+DY1*WT1+DY2*WT2
14      CONTINUE
        chs=chs+DY1*DY1*sinv1+DY2*DY2*sinv2
15    CONTINUE
*
      DO 17 J=2,m
        DO 16 K=1,J-1
          alf(K,J)=alf(J,K)
16      CONTINUE
17    CONTINUE
      RETURN
      END
*********************************************************************


***********************************************************************
* Added 2 february 2012 by D. van der Marel
***********************************************************************
* Calculation of mem and d(mem)/d(parameter) 
* Called by SUBROUTINE MRQCOF
*---------------------------------------------------------------------
* Label       || Meaning                               || input/output
*---------------------------------------------------------------------
* x           || Frequency                             || in
* pm(m)       || parameters for mem                    || in
* mem         || mem                                   || out
* dmdp(nfix)  || derivate of mem with resp to pm(i)    || out
* m           || no of blocks in histogram             || in
*---------------------------------------------------------------------
      subroutine myfunc(iw,pm,cmem,cdmdp,mtot)
*      
      parameter(nw=2000,ne=200,mp=99)
      double precision x(nw)
      complex*16 cy(nw),cs(nw)
      integer n
      common/data/x,cy,cs,n 
*
      double precision fe(nw,ne),few(nw,ne)
      complex*16 ll(nw,ne,mp),llw(nw,ne,mp)      
      integer ifunc,nfunc
      common/funcarr/fe,few,ll,llw,ifunc,nfunc     
* 
      double precision pm(mp)
      complex*16 cmem,cdmdp(mp)
      integer iw,mtot,m
*  
      integer ie,j,intmx   
      double precision w,nfe,eintmn,eintmx,e,de,y,nfew,w1,w2,ew,a,b
      complex*16 ee,eew,lj,lwj,dum,dsum(mp),cz,cdzdp(mp),cdmda,cdmdb
*      
      double precision w0,dw
      common/inblock/w0,dw
      integer ifitab
      double precision kbt,tpi,a0,b0
      common/memoir/kbt,tpi,a0,b0,ifitab
* 
      w=x(iw)
*
      ifunc=ifunc+1
*      write(*,*) 'myfunc with iw, ifunc, nfunc = ',iw,ifunc,nfunc
      m=mtot-ifitab
      if (ifitab.eq.0) then
       b=b0
       a=a0
      endif 
      if (ifitab.eq.1) then
       a=a0
       b=pm(m+1)
      endif
      if (ifitab.eq.2) then
       a=pm(m+2)
       b=pm(m+1)
      endif 
      eintmn=-w-21*kbt
      eintmx=21*kbt
      intmx=ne
      de=(eintmx-eintmn)/intmx
***** here the loop begins for the integration of the memory function
      cz=(0.,0.)
      do 5 j=1,m
       cdzdp(j)=(0.,0.)
c       if(pm(j).lt.0) write(*,*) 'in myfunc pm(',j,')=',pm(j),'<0'
5     continue
      do 60 ie=1,intmx
       e=eintmn+ie*de
       ew=e+w
*
****   Fermi function at energy e
       if (ifunc.le.nfunc) then
        y=e/kbt 
        if (abs(y).lt.60) then
         nfe=1./(exp(y)+1.)
         else
         if (y.gt.0) then
          nfe=0.
          else
          nfe=1.
         endif 
        endif
****    Fermi function at energy e+w 
        y=(e+w)/kbt
        if (abs(y).lt.60) then
         nfew=1./(exp(y)+1.)
        else
         if (y.gt.0) then
          nfew=0.
         else
          nfew=1.
         endif 
        endif        
*
        fe(iw,ie)=nfe        
        few(iw,ie)=nfew
*        
       else
        nfe=fe(iw,ie)
        nfew=few(iw,ie)
       endif  
*
****   selfenergy at energy e and ew
       ee=(0.,0.)
       eew=(0.,0.)       
       do 33 j=1,m
        w1=w0+(j-1)*dw
        w2=w1+dw 
        if (ifunc.le.nfunc) then
         call kernint(e,w1,w2,lj)
         ll(iw,ie,j)=lj
        else
         lj=ll(iw,ie,j)
        endif                     
        ee=ee+lj*pm(j)
        if (ifunc.le.nfunc) then
         call kernint(ew,w1,w2,lwj)                   
         llw(iw,ie,j)=lwj
        else
         lwj=llw(iw,ie,j)
        endif  
        eew=eew+lwj*pm(j)  
33     continue 
       cz=cz+de*(nfe-nfew)/(w-eew+conjg(ee))  
*            
       do 35 j=1,m
        w1=w0+(j-1)*dw
        w2=w1+dw 
        lj=ll(iw,ie,j)        
        lwj=llw(iw,ie,j)       
        dum=(lwj-conjg(lj))*de        
        cdzdp(j)=cdzdp(j)+dum*(nfe-nfew)/(w-eew+conjg(ee))**2        
35     continue  
       
60    continue   
*
***   memory function: mem
      cmem=-w+a*w/(cz-b*w**2)       
*
***   memory function derivatives: dmdp
      do 65 j=1,m
      cdmdp(j)=-cdzdp(j)*(cmem+w)**2/(a*w)
65    continue        
      cdmda=(cmem+w)/a
      cdmdb=w*(cmem+w)**2/a      
      cdmdp(m+2)=cdmda
      cdmdp(m+1)=cdmdb
c      write(*,*) 'in myfunc m(x)=',cmem
      return
      end
***********************************************************************

***********************************************************************
      subroutine kernint(w,w1,w2,l12)
      double precision w,w1,w2
      complex*16 l12
*      
      double precision a1,a2
      complex*16 ca1,ca2,cb1,cb2,cca1,cca2,ccb1,ccb2,gammln
*
      integer ifitab
      double precision kbt,tpi,a0,b0
      common/memoir/kbt,tpi,a0,b0,ifitab
*       
      a1=0.5*w1/kbt
      a2=0.5*w2/kbt
      if (a1.lt.89) a1=log(sinh(a1))+log(2.)
      if (a2.lt.89) a2=log(sinh(a2))+log(2.)
      ca1=0.5+(0.,1.)*(w1-w)/(tpi*kbt)
      ca2=0.5+(0.,1.)*(w2-w)/(tpi*kbt)
      cb1=0.5-(0.,1.)*(w1+w)/(tpi*kbt)
      cb2=0.5-(0.,1.)*(w2+w)/(tpi*kbt)
      cca1=gammln(ca1)
      cca2=gammln(ca2)
      ccb1=gammln(cb1)
      ccb2=gammln(cb2)
      l12=a2-a1+cca2-cca1+ccb2-ccb1
      l12=l12*(0.,-1.)*tpi*kbt
*      
      return
      end
***********************************************************************


***********************************************************************
      complex*16 FUNCTION GAMMLN(XX)
      complex*16 xx
*      
      complex*16 x,tmp,ser
      double precision COF(6),STP,HALF,ONE,FPF
*
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
*
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
***********************************************************************


*********************************************************************
      SUBROUTINE GAUSSJ(A,B,m)
      PARAMETER(mp=99)
      double precision A(mp,mp),B(mp)
      integer m
*      
      integer IPIV(mp),INDXR(mp),INDXC(mp)
      integer j,i,k,irow,icol,l,ll
      double precision big,dum,pivinv
*
      DO 11 J=1,m
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,m
        BIG=0.
        DO 13 J=1,m
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,m
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                write(*,*) 'Singular matrix, ipiv'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,m
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DUM=B(IROW)
          B(IROW)=B(ICOL)
          B(ICOL)=DUM
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) write(*,*) 'Singular matrix, A'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,m
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        B(ICOL)=B(ICOL)*PIVINV
        DO 21 LL=1,m
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,m
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            B(LL)=B(LL)-B(ICOL)*DUM
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=m,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,m
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
*********************************************************************





