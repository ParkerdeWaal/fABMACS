cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine is hard coded for alanine dipeptide, given in as transparent
! a representation as possible. This is written for readability, and to 
! demonstrate simple differences between WTmetaD and mABP. This is not 
! written for speed or elegance. -BMD 2015
!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hellof(istep,xx,ff)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)      
!declare some stuff 
      real*8 dpop(2,300,300),pop(300,300)
      real*8 dela(3,300,300),ddela(300,300)
      real*8 rfree(300,300),decon(300,300),ocog(3)
      real*8 alp,alp2,bee,cee,omega,DelT,pi,pi2
      real*8 scal,bolt,alpe,bah,bah2,are
      real*8 rcom1,rcom2,rcom3
      integer nbin,imabp,ialp
      save !and save the stuff
      real ff(8*3),xx(8*3)
      real*8 jaco(2,8*3),cog(3)
      real*8 deid(3,4),qa(12)
      real*8 ave1(2),freest,s,sf
      real*8 orml,sumd,at1,at2,bang,arg

      if(istep.eq.0)then
         bah=1d0
         bah2 = 1d0
         bah3=0d0
         bolt = 300d0/11604.467d0 *96.48d0 !hard code TEMPERATURE!!!!!!
         nbin = 300 !number of bins in phi-psi
         open(11,file='./reffreeE') !reference free energy
         do k=1,nbin
            do l=1,nbin
               read(11,*) unk, unk, rfree(k,l) !read it in
            enddo
            read(11,*)
         enddo
         close(11) !close it

         open(11,file='params.in') !user inputs
         read(11,*) bee !b parameter in mABP
         read(11,*) cee !c parameter in mABP (is taken as c*dt)
         read(11,*) alp !alpha
         read(11,*) imabp !0=mABP, 1=WTmetaD, 2=muTmetaD, 3=SHUS
         read(11,*) scal !SHUS power== gamma/(1-gamma) or m for muTmetaD
         read(11,*) are !the r parameter for muT
         close(11) !close it

         
         pi=acos(-1d0) !have some pi
         pi2=pi*2d0    !twice
         alpe=alpe*pi/180d0
c         alp=8d0
c         alp = alp*pi/180d0 !convert to rads
c         alp2 = alp*alp !save the alpha^2
         ialp=int(alp)

         dx = (pi2-0d0)/dble(nbin) !grid width
         tol = alp*dx
         tol2=4d0*dx

         do j=1,nbin !pre-compute delta_alpha "hill functions"
            xatj = (j-1)*dx
            sum=0d0
            sum2=0d0
            do i=1,nbin
               xati = (i-1)*dx
               arg2=xati-xatj
               if(arg2.gt.pi)then
                  arg2=arg2-pi*2d0
               endif
               if(arg2.lt.-pi)then
                  arg2 = arg2 +pi*2d0
               endif         
!-----MOLLIIIIIII
               ex = arg2/tol
               aex=ex*ex
               aex=sqrt(aex)
               aex2=aex*aex
               func=0d0
               dermoli=0d0
               if(aex.lt.1d0)then
                  func = exp(-scal/(1d0-aex2))/exp(-scal) !so height is 1 at x=0
                  aa=scal*func/((1d0-aex2)*(1d0-aex2))!should be negative?
                  dermoli=aa*2d0*ex/tol !right?
               endif
               ex = arg2/(tol/2d0)
               aex=ex*ex
               aex=sqrt(aex)
               aex2=aex*aex
               funca=0d0
               dermolia=0d0
               if(aex.lt.1d0)then
                  funca = exp(-scal/(1d0-aex2))/exp(-scal) !so height is 1 at x=0
                  aa=scal*funca/((1d0-aex2)*(1d0-aex2))!should be negative?
                  dermolia=aa*2d0*ex
               endif
c               dela(1,i,j)=exp(-arg2*arg2/alp2)
               ex = arg2/(tol/3d0)
               aex=ex*ex
               aex=sqrt(aex)
               aex2=aex*aex
               funcb=0d0
               dermolib=0d0
               if(aex.lt.1d0)then
                  funcb = exp(-scal/(1d0-aex2))/exp(-scal) !so height is 1 at x=0
                  aa=scal*funcb/((1d0-aex2)*(1d0-aex2))!should be negative?
                  dermolib=aa*2d0*ex
               endif
c               dela(1,i,j)=exp(-arg2*arg2/alp2)
               dela(1,i,j)=func!(func+funca+funcb)/3d0
               ddela(i,j) = dermoli!(dermoli +dermolia+dermolib)/3d0 !this is for gradient

               sum=sum+dx*dela(1,i,j)
               sum2=sum2+dx*exp(-arg2*arg2/(10d0*pi/180d0)**2)
               if(j.eq.150)then
                  write(90,*) i, dela(1,i,j), 
     .                 exp(-arg2*arg2/(10d0*pi/180d0)**2), sum, sum2
               endif
!--------------------------------------------------------------------
!              this is used for the "vanishingly small" alpha,
!              and uses a width of 2 degrees
               ex2 = arg2/tol2
               aex2=ex2*ex2
               aex2=sqrt(aex2)
               aex22=aex2*aex2
               func2=0d0
               if(aex2.lt.1d0)then
                  func2 = exp(-scal/(1d0-aex22))/exp(-scal) !so height is 1 at x=0
               endif
c               dela(2,i,j)=exp(-arg2*arg2/2d0**2)  !for mabp stuff and muTmetaD
               dela(2,i,j)=func2!exp(-arg2*arg2/0.0523d0**2)  !for mabp stuff and muTmetaD
!--------------------------------------------------------------------
c               ddela(i,j) = dela(1,i,j)*2d0*arg2/alp2 !this is for gradient
            enddo

         enddo
         close(90)
         omega = bolt*bee*cee !WTmetaD conversion, and muTmetaD
         DelT = bolt*bee/(1d0-bee) !ditto

         do j=1,nbin
            do i=1,nbin                  
               pop(i,j) = 0d0   !initialize arrays
               if(imabp.eq.0)then
               decon(i,j)=1d0
               else
               decon(i,j)=0d0
               endif
               dpop(1,i,j) = 0d0
               dpop(2,i,j) = 0d0
            enddo
         enddo
      com1=0d0
      com2=0d0
      com3=0d0
      do i=1,8!22
         j=1
         com1=com1+xx(3*i-(3-j))
         j=2
         com2=com2+xx(3*i-(3-j))
         j=3
         com3=com3+xx(3*i-(3-j))
      enddo
      rcom1=com1/8d0!22d0
      rcom2=com2/8d0!22d0
      rcom3=com3/8d0!22d0 !set refs

      endif                     !end the initial stuff
!do basic PBC here with a cubic box, careful to hardcode the box
!hold COM at center too since gromacs is a...
      com1=0d0
      com2=0d0
      com3=0d0
      do i=1,8!22
         j=1
         com1=com1+xx(3*i-(3-j))
         j=2
         com2=com2+xx(3*i-(3-j))
         j=3
         com3=com3+xx(3*i-(3-j))
      enddo
      com1=com1/8d0!22d0
      com2=com2/8d0!22d0
      com3=com3/8d0!22d0

c      width=2.69741d0 !this shouldnt be active...
c      do i=1,22
c         ip=i
c         do j=1,3
c            pos=xx(3*ip-(3-j))
c            if(pos.lt.0d0)then
c               xx(3*ip-(3-j))=pos+width
c            elseif(pos.gt.width)then
c               xx(3*ip-(3-j))=pos-width
c            endif
c         enddo
c      enddo
      np=8
      width=2.69297
      width2=width*width
      width24=width2/4d0
      do l=1,1
      do i=1,np-1
         do j=i+1,i+1!np
c            if(ipair(i,j).eq.1)then
               do k=1,3
                  xxi=xx(i*3-(3-k))
                  xxj=xx(j*3-(3-k))
c                  d0i=xxi-x0(i*3-(3-k))
c                  d0j=xxj-x0(j*3-(3-k))
                  dr=xxi-xxj
                  dr2=dr*dr                  
                  if(dr2.gt.width24)then !somebody need unwrappin, just roll j
                     if((dr+width)*(dr+width).lt.width24)then
                        xx(j*3-(3-k))=xx(j*3-(3-k))-width
                     elseif((dr-width)*(dr-width).lt.width24)then
                        do ii=1,i
                           xx(ii*3-(3-k))=xx(ii*3-(3-k))-width !should loop backward on i?
                        enddo
                     else
                        write(*,*) "oh shit... i cant map this guy! "
                        write(*,*) i, j, xxi, xxj, dr, sqrt(dr2)
                        stop
                     endif
                  endif
               enddo
c            endif
         enddo
      enddo
      enddo!l
      do i=1,np
         do j=1,3
            cog(j)=cog(j)+xx(i*3-(3-j))
         enddo
      enddo
      cog=cog/dble(np)
      if(istep.eq.0)ocog=cog

      do j=1,3
         dc=cog(j)-ocog(j)
         dc=dc*dc
         if(dc.gt.width24)then  !jumped
            do i=1,np
               xx(i*3-(3-j))=xx(i*3-(3-j))-
     .              width*(cog(j)-ocog(j))/sqrt(dc)
            enddo
            cog(j)=cog(j)-width*(cog(j)-ocog(j))/sqrt(dc)
         endif
      enddo
      ocog = cog

     
      jaco = 0d0 !this will hold the phi-psi derivatives
      natom1=1!5 
      natom2=2!7
      natom3=3!9
      natom4=4!15
      qa(1) = xx(3*natom1-(3-1))
      qa(2) = xx(3*natom1-(3-2))
      qa(3) = xx(3*natom1-(3-3))
      qa(4) = xx(3*natom2-(3-1))
      qa(5) = xx(3*natom2-(3-2))
      qa(6) = xx(3*natom2-(3-3))
      qa(7) = xx(3*natom3-(3-1))
      qa(8) = xx(3*natom3-(3-2))
      qa(9) = xx(3*natom3-(3-3))
      qa(10) = xx(3*natom4-(3-1))
      qa(11) = xx(3*natom4-(3-2))
      qa(12) = xx(3*natom4-(3-3))
      do i=1,3
         do j=1,4
            deid(i,j) = 0d0
         enddo
      enddo
      angle = 0d0
      call bmdDI(qa,deid,angle)
      angl1 = angle*pi/180d0
      do idm=1,3
         ip = natom1*3-(3-idm)
         jaco(1,ip)=jaco(1,ip)+deid(idm,1)
      enddo
      do idm=1,3
         ip = natom2*3-(3-idm)
         jaco(1,ip)=jaco(1,ip)+deid(idm,2)
      enddo
      do idm=1,3
         ip = natom3*3-(3-idm)
         jaco(1,ip)=jaco(1,ip)+deid(idm,3)
      enddo
      do idm=1,3
         ip = natom4*3-(3-idm)
         jaco(1,ip)=jaco(1,ip)+deid(idm,4)
      enddo
      
      !second angle
      natom1=5!7
      natom2=6!9
      natom3=7!15
      natom4=8!17
      qa(1) = xx(3*natom1-(3-1))
      qa(2) = xx(3*natom1-(3-2))
      qa(3) = xx(3*natom1-(3-3))
      qa(4) = xx(3*natom2-(3-1))
      qa(5) = xx(3*natom2-(3-2))
      qa(6) = xx(3*natom2-(3-3))
      qa(7) = xx(3*natom3-(3-1))
      qa(8) = xx(3*natom3-(3-2))
      qa(9) = xx(3*natom3-(3-3))
      qa(10) = xx(3*natom4-(3-1))
      qa(11) = xx(3*natom4-(3-2))
      qa(12) = xx(3*natom4-(3-3))
      do i=1,3
         do j=1,4
            deid(i,j) = 0d0
         enddo
      enddo
      angle = 0d0
      call bmdDI(qa,deid,angle)
      angl2 = angle*pi/180d0

      do idm=1,3
         ip = natom1*3-(3-idm)
         jaco(2,ip)=jaco(2,ip)+deid(idm,1)
      enddo
      do idm=1,3
         ip = natom2*3-(3-idm)
         jaco(2,ip)=jaco(2,ip)+deid(idm,2)
      enddo
      do idm=1,3
         ip = natom3*3-(3-idm)
         jaco(2,ip)=jaco(2,ip)+deid(idm,3)
      enddo
      do idm=1,3
         ip = natom4*3-(3-idm)
         jaco(2,ip)=jaco(2,ip)+deid(idm,4)
      enddo

      !find the current position in bins
      dx = (pi2-0d0)/dble(nbin)
      ibin = int(angl1/dx)+1      
      if(ibin.gt.nbin)ibin=nbin
      if(ibin.lt.1)ibin=1     
      ibin1=ibin
      ibin = int(angl2/dx)+1      
      if(ibin.gt.nbin)ibin=nbin
      if(ibin.lt.1)ibin=1     
      ibin2=ibin

!-------------------------------------------------------------------------
!      UPDATE THE BIAS POTENTIAL HERE ---
!-------------------------------------------------------------------------

      if(imabp.eq.0)then !if mABP, then do this to get forces
c      denom = 1d0+cee*pop(ibin1,ibin2)            
c      ave1(1) = -bolt*cee*dpop(1,ibin1,ibin2)/denom 
c      denom = 1d0+cee*pop(ibin1,ibin2)      
c      ave1(2) = -bolt*cee*dpop(2,ibin1,ibin2)/denom 
      denom = 1d0+cee*(1d0-bee)*pop(ibin1,ibin2)            
      ave1(1) = -cee*bee*bolt*dpop(1,ibin1,ibin2)/denom 
      denom = 1d0+cee*(1d0-bee)*pop(ibin1,ibin2)      
      ave1(2) = -cee*bee*bolt*dpop(2,ibin1,ibin2)/denom 

      elseif(imabp.eq.1)then !do WTmetaD
      ave1(1) = -dpop(1,ibin1,ibin2)
      ave1(2) = -dpop(2,ibin1,ibin2)            
      s=exp(-pop(ibin1,ibin2)/DelT)*omega
      elseif(imabp.eq.2)then !do muTmetaD
      ave1(1) = -dpop(1,ibin1,ibin2)
      ave1(2) = -dpop(2,ibin1,ibin2)            
c                  alpno = 2d0*pi/180d0
c                  relate = alp/alpno !should be squared
c                  relate = relate *relate
c                  relate = relate * omega/DelT
c         fac=(relate*bah2+1d0)**scal
c         s=omega*fac/(relate*decon(ibin1,ibin2)+1d0)         

         fac=(are*bah2+1d0)**scal
cc         s=omega*fac/decon(ibin1,ibin2)
         s=omega*fac/(are*decon(ibin1,ibin2)+1d0)         
      else !shus-ish
      ave1(1) = -dpop(1,ibin1,ibin2)
      ave1(2) = -dpop(2,ibin1,ibin2)      
      s=omega/(log(pop(ibin1,ibin2)+1d0)**scal+1d0)
      endif
      !update grids
      bah = 1d80
      bah2=1d0
      bah3=0d0
      snorm = 0d0
      jst=ibin2-(ialp+4) !11 because used 10 in def of tol
c      jst=max(jst,1)
      jnd=ibin2+(ialp+4)
c      jnd=min(nbin,jnd)
      do j =jst,jnd!1,nbin!jst,jnd!1,nbin!jst,jnd
         jb=j
         if(j.lt.1)jb=j+nbin
         if(j.gt.nbin)jb=j-nbin         

c      do j = 1,nbin!jst,jnd
         atj = (j-1)*dx
         jbin2=jb !why do this?
         arg2=angl2-atj
         if(arg2.gt.pi)then
            arg2=arg2-pi2
         endif
         if(arg2.lt.-pi)then
            arg2 = arg2 +pi2
         endif         
         boomc= dela(2,ibin2,jb) !this is hill with "vanishingly small" width
         booma=dela(1,ibin2,jb)  !this is hill for biasing
         ist=ibin1-(ialp+4)           !11 because used 10 in def of tol
c         ist=max(ist,1)
         ind=ibin1+(ialp+4)
c         ind=min(nbin,ind)
         do i=ist,ind
c         do i = 1,nbin
            ib=i
            if(i.lt.1)ib=i+nbin
            if(i.gt.nbin)ib=i-nbin         
            ati = (i-1)*dx
            jbin = ib
            arg1=angl1-ati
            if(arg1.gt.pi)then
               arg1=arg1-pi2
            endif
            if(arg1.lt.-pi)then
               arg1 = arg1 +pi2
            endif         
            boomcc= dela(2,ibin1,ib)!second part of "vanishingly small" width
            boom = dela(1,ibin1,ib)*booma ! product of two hills for bias
            if(imabp.eq.0)then !mabp
            pop(jbin,jbin2) = pop(jbin,jbin2) + boom
!deconvoluted free energy estimate ------------------------------
            decon(jbin,jbin2)=decon(jbin,jbin2)+boomcc*boomc
!----------------------------------------------------------------
            dpop(1,jbin,jbin2) = dpop(1,jbin,jbin2)+ddela(ibin1,ib)
     .      *dela(1,ibin2,jb)
            dpop(2,jbin,jbin2) = dpop(2,jbin,jbin2)+ddela(ibin2,jb)
     .          *dela(1,ibin1,ib)
            else !metaD, or SHUS, depending only on how "s" is defined
!---------- Notice the array updates are different only by a factor of "s"
            decon(jbin,jbin2)=decon(jbin,jbin2)+boomcc*boomc!for the PMF
            pop(jbin,jbin2) = pop(jbin,jbin2) + boom*s
            if(decon(jbin,jbin2).lt.bah)bah=decon(jbin,jbin2)
            if(decon(jbin,jbin2).gt.bah2)bah2=decon(jbin,jbin2)
            if(pop(jbin,jbin2).gt.bah3)bah3=pop(jbin,jbin2)
            dpop(1,jbin,jbin2) = dpop(1,jbin,jbin2)+ddela(ibin1,ib)
     .      *dela(1,ibin2,jb)*s
            dpop(2,jbin,jbin2) = dpop(2,jbin,jbin2)+ddela(ibin2,jb)
     .          *dela(1,ibin1,ib)*s
c
c            dpop(1,jbin,jbin2) = dpop(1,jbin,jbin2)+ddela(ibin1,i)*s
c     .           *dela(1,ibin2,j)
c            dpop(2,jbin,jbin2) = dpop(2,jbin,jbin2)+ddela(ibin2,j)*s
c     .           *dela(1,ibin1,i)
            endif
         enddo
      enddo
!----------------------------------------------------------------------
!Add the bias forces back to system force
      do i=1,8*3!22*3 
         ff(i) = ff(i)+(ave1(1))*jaco(1,i)+(ave1(2))*jaco(2,i)
      enddo

c      do i=1,8!22
c         j=1
c         ip=3*i-(3-j)
c         ff(ip) = ff(ip) - 100d0*bolt*(com1-rcom1)/22d0
c         j=2
c         ip=3*i-(3-j)
c         ff(ip) = ff(ip) - 100d0*bolt*(com2-rcom2)/22d0
c         j=3
c         ip=3*i-(3-j)
c         ff(ip) = ff(ip) - 100d0*bolt*(com3-rcom3)/22d0
c      enddo
!----------------------------------------------------------------------
!     Write a restart file and get convergence curve
      if(mod(istep,50000).eq.0)then 
         !Get the zero-of energy
         bang=0d0    
         do j=1,nbin
            do i=1,nbin
               if(imabp.eq.0)then !mABP
               !mABP sans-MOLLY
c                  freest=-bolt*log(decon(i,j)*pop(i,j)**(bee/(1d0-bee)))
                  freest=-bolt/(1d0-bee)
     .                 *log(cee*(1d0-bee)*pop(i,j)+1d0)

               elseif(imabp.eq.1)then !WTmetaD
                  sf=exp(-pop(i,j)/DelT)*omega
                  freest = bolt*log(sf)-pop(i,j)
               elseif(imabp.eq.2)then
c                  alpno = 2d0*pi/180d0
c                  relate = alp/alpno !should be squared
c                  relate = relate *relate
c                  relate = relate * omega/DelT
c         fac=(relate*bah2+1d0)**scal
c         sf=omega*fac/(relate*decon(i,j)+1d0)         
         fac=(are*bah2+1d0)**scal
         sf=omega*fac/(are*decon(i,j)+1d0)         
                  freest = bolt*log(sf)-pop(i,j)
               else
                  freest=bolt*log(1d0/(log(pop(i,j)+1d0)**scal+1d0))
     .                 -pop(i,j)                  
               endif
               if(freest.lt.bang)bang=freest !Sets the zero
            enddo
         enddo
         er = 0d0
         open(89,file='freeE')
         sumd=0d0
         orml=0d0
         do j=1,nbin
            at2 = dx * (dble(j) - 0.5d0)
            do i=1,nbin           
               at1 = dx * (dble(i) - 0.5d0)
               if(imabp.eq.0)then !mABP
               !mABP sans-MOLLY
c                  freest=-bolt*log(decon(i,j)*pop(i,j)**(bee/(1d0-bee)))  
                  freest=-bolt/(1d0-bee)
     .                 *log(cee*(1d0-bee)*pop(i,j)+1d0)

               elseif(imabp.eq.1)then !WTmetaD
                  sf=exp(-pop(i,j)/DelT)*omega
                  freest = bolt*log(sf)-pop(i,j)
               elseif(imabp.eq.2)then
c                  alpno = 2d0*pi/180d0
c                  relate = alp/alpno !should be squared
c                  relate = relate *relate
c                  relate = relate * omega/DelT
c         fac=(relate*bah2+1d0)**scal
c         sf=omega*fac/(relate*decon(i,j)+1d0)         

         fac=(are*bah2+1d0)**scal
c                  fac=bah2**scal
c                  sf=omega*fac/decon(i,j)
         sf=omega*fac/(are*decon(i,j)+1d0)         
                  freest = bolt*log(sf)-pop(i,j)
               else
                  freest=bolt*log(1d0/(log(pop(i,j)+1d0)**scal+1d0))
     .                 -pop(i,j)                  
               endif

               freest=freest-bang !put min to zero
               if(rfree(i,j).lt.30d0)then
c                  write(*,*) rfree(i,j), freest
               sumd=sumd+dabs(freest-rfree(i,j))
               orml = orml+1d0
               endif
               write(89,*) at1, at2, freest, decon(i,j)/bah2 !write free energy to fort.89
            enddo
            write(89,*)
         enddo
         close(89)
!                 TimeStep, Angle1, Angle2, Convergence metric, hill height (not mabp)
         write(88,*) istep, angl1, angl2, sumd/orml,  s
         call flush(88) !write to fort.88
c         stop
      endif


      return
      end !end the ABP routine(s)-------------------------
!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------
      subroutine bmdDI(qa,deid,angle)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      real*8 deid(3,4),qa(12),angle,radian
      
      radian = 180d0/acos(-1d0)
      do i=1,4
         do  j=1,3
            deid(j,i)=0d0
         enddo
      enddo
      
      xia = qa(1)
      yia = qa(2)
      zia = qa(3)
      xib = qa(4)
      yib = qa(5)
      zib = qa(6)
      xic = qa(7)
      yic = qa(8)
      zic = qa(9)
      xid = qa(10)
      yid = qa(11)
      zid = qa(12)
      xba = xib - xia
      yba = yib - yia
      zba = zib - zia
      xcb = xic - xib
      ycb = yic - yib
      zcb = zic - zib      
      xdc = xid - xic
      ydc = yid - yic
      zdc = zid - zic
      xt = yba*zcb - ycb*zba
      yt = zba*xcb - zcb*xba
      zt = xba*ycb - xcb*yba
      xu = ycb*zdc - ydc*zcb
      yu = zcb*xdc - zdc*xcb
      zu = xcb*ydc - xdc*ycb
      xtu = yt*zu - yu*zt
      ytu = zt*xu - zu*xt
      ztu = xt*yu - xu*yt
      rt2 = xt*xt + yt*yt + zt*zt
      ru2 = xu*xu + yu*yu + zu*zu
      if(rt2.lt. 0.000000001)then
       rt2=0.000000001
      endif
      if(ru2.lt. 0.000000001)then
       ru2 = 0.000000001
      endif
      rtru = sqrt(rt2 * ru2)
      if(rtru.ne.0d0)then
      rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
      cosine = (xt*xu + yt*yu + zt*zu) / rtru
      sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
      cosine = min(1.0d0,max(-1.0d0,cosine))
      angle = radian * acos(cosine)
      if (sine .lt. 0.0d0)  angle = -angle+360d0
      if(angle.ge.360d0)angle=angle-360d0
      dedphi = 1d0
      xca = xic - xia
      yca = yic - yia
      zca = zic - zia      
      xdb = xid - xib
      ydb = yid - yib
      zdb = zid - zib
      dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
      dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
      dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
      dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
      dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
      dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
      dedxia = zcb*dedyt - ycb*dedzt
      dedyia = xcb*dedzt - zcb*dedxt
      dedzia = ycb*dedxt - xcb*dedyt
      dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
      dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
      dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
      dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
      dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
      dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
      dedxid = zcb*dedyu - ycb*dedzu
      dedyid = xcb*dedzu - zcb*dedxu
      dedzid = ycb*dedxu - xcb*dedyu
      deid(1,1) = deid(1,1) + dedxia
      deid(2,1) = deid(2,1) + dedyia 
      deid(3,1) = deid(3,1) + dedzia
      deid(1,2) = deid(1,2) + dedxib
      deid(2,2) = deid(2,2) + dedyib
      deid(3,2) = deid(3,2) + dedzib
      deid(1,3) = deid(1,3) + dedxic
      deid(2,3) = deid(2,3) + dedyic
      deid(3,3) = deid(3,3) + dedzic
      deid(1,4) = deid(1,4) + dedxid
      deid(2,4) = deid(2,4) + dedyid
      deid(3,4) = deid(3,4) + dedzid
      else
       write(*,*) 'failed dihed'
       stop
      endif
c      deid = deid !* radian
      RETURN
      END
