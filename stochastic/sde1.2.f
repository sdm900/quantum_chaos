*************************************************************************
C      PROGRAM: SDE
C       ----------------
C       SIMULATES SDE
*************************************************************************
C
C       Reference: Unpublished/ Peter D. Drummond, 16 Nov 1994
C       Algorithm:  Stochastic phase space
C
***************************************************************************
C
C       ps      program switches:
C                ps(1) = algorithm switch, to alter the sde
C                ps(2) = sde type; =1,2 (for dimension-doubling)
C                ps(3) = integration option ; not currently used
C                ps(4) = modes/lattice point
C                ps(5) = variables/lattice point (ie. =ps(2)*ps(4))
C                ps(6) = unused
C                ps(7) = real noises per point   - depends on the sde
C                ps(8) = iterations of the implicit sde solver
C
C  PARAMETER VALUES
C
        integer MSEED,MPATH,MPT,MSTEP
        integer MM,MV,MN,MPA,MPRINTS
        parameter (MSEED = 100000000)
        parameter (MPATH = 1000000)
        parameter (MPT = 1024)
        parameter (MSTEP = 1024)
        parameter (MM = 10)
        parameter (MV = 10)
        parameter (MN = 10)
        parameter (MPA = 10)
        parameter (MPRINTS = 10)
C
C  DATA DECLARATIONS
C
C       include 'param.f'
        complex*16 z_i(MV),z_c(MV),z_f(MV)
        complex*16 ps_c, ps_f,ci,ci2
        complex*16  A(MV),pa(MPA),lpa(MPA)
        real*8 delw_1(MN), delw_2(MN), delw_t(MN)
        real*8 time,deltc,deltf,sigma,pi,pi2
        real*8 gauss,x,i(2*MV),s(MN)
        real*8 nc,n2c,nf,n2f,error,time1,nsigma
        real*8 errormax,errpoint,timec,timef
        real*8 sigmamax,sigpoint
        real*8 out_c(2,0:MM,0:MPT),out_f(2,0:MM,0:MPT)
        integer seed,paths,points,steps,moms,n
        integer path,point,step,mom,fn
        integer var,noise,ps(10),j(10)
        integer n_print
C
C  CONSTANTS
C
        ci=(0.,1.d0)
        pi=4.*atan(1.d0)
        pi2=2.*pi
        fn=9
C
C  INPUT SIMULATION DATA TO START
C
        open(fn, file='SDE1.dat',status='unknown',form='formatted')
        write(*,*) ' _________________________________________________ '
        write(*,*) '|                                                 |'
        write(*,*) '| SDE VERSION 1.0                                 |'
        write(*,*) '|                                                 |'
        write(*,*) '| SIMULATES SDE - COPYRIGHT PETER D DRUMMOND 1996 |'
                                  |'
        write(*,*) '|                                                 |'
        write(*,*) '| Estimates of sampling error, step error made.   |'
        write(*,*) '|                                                 |'
        write(*,*) '| Output is re-usable as an input OR gnuplot file.|'
                                  |'
        write(*,*) '|                                                 |'
        write(*,*) '|_________________________________________________|'
        write(fn,*) '# PROGRAM: SDE VERSION 1.0'
        call input(1,'seed                ',1.d0,MSEED*1.d0,x,seed,fn)
        call input(1,'paths               ',1.d0,MPATH*1.d0,x,paths,fn)
        call input(1,'points              ',1.d0,MPT*1.d0,x,points,fn)
        call input(1,'steps/pnt           ',1.d0,MSTEP*1.d0,x,steps,fn)
        call input(1,'iterations          ',1.d0,MSTEP*1.d0,x,ps(8),fn)
        call input(1,'max moment          ',1.d0,MM*1.d0,x,moms,fn)
        call input(1,'max time            ',0.0d0,1.D100,time,n,fn)
        deltc = time/(steps*points)
        deltf=deltc/2.
        sigma=1./sqrt(deltf)                            !sd = sqrt(1/ timestep)
        n_print=max0(1,paths/MPRINTS)
C        n_fractf=0
C        n_fractc=0
C
D       write(*,*) 'deltc,deltf,sigma = ', deltc,deltf,sigma
C  INITIALISE RANDOM NUMBERS
C
        call in_rand (seed)                             !initialize random
numbers
C
C  INITIALISE EQUATION PARAMETERS
C
        call equ_par (ps,pa,fn)
        call equ_in (ps,lpa,z_i,fn)
        write(*,*) ' _________________________________________________'
        write(*,*) '|                                                 |'
        write(*,*) '| INPUT std deviation, sigma=sigma_x=sigma_y      |'
        write(*,*) '| IE assumed equal real and imaginary fluctuations|'
        write(*,*) '|_________________________________________________|'
        call input(ps(4),'input std deviation ',0.0d0,1.D100,i,j,fn)
        call input(ps(7),'stochastic sd: [b_n]',0.0d0,1.D100,s,j,fn)
        do var=1,ps(7)
                s(var)=sigma*s(var)
        end do
D       write(*,*) 'deltc,deltf,sigma = ', deltc,deltf,sigma
C
C  LOOP AROUND ENSEMBLE
C
         do  path=1,paths
                if ((path/n_print)*n_print.eq.path) then
                        write(*,*) 'Starting path ',path,' / ',paths
                end if
                call gauss(ps(4),i,delw_1,delw_2,delw_t)
                ps_c=(0.,0.)
                ps_f=(0.,0.)
                do var=1,ps(4)
                        z_c(var) = z_i(var)+delw_1(var)+ci*delw_2(var)
                        z_f(var) = z_c(var)
                end do
                if (ps(2).eq.2) then
                do var=1,ps(4)
                        z_c(ps(4)+var) = z_c(var)
                        z_f(ps(4)+var) = z_c(var)
                end do
                end if
                call calc(ps,ps_c,z_c,out_c,moms,0)
                call calc(ps,ps_f,z_f,out_f,moms,0)
                do point=1,points
                  do step=1,steps
                     call gauss(ps(7),s,delw_1,delw_2,delw_t)
C
C  COARSE AND FINE INTEGRATION STEPS
C
D               write(*,*) 'ps_c,z_c,delw_t,deltc'
D               write(*,*) ps_c,z_c,delw_t,deltc
D               write(*,*) 'pa,ps(5)'
D               write(*,*) pa,ps(5)
D               write(*,*) 'COARSE:STEP=' ,step
                call intg(ps,ps_c,z_c,delw_t,deltc,pa,lpa)
D               write(*,*) 'FINE:STEP=' ,step
                call intg(ps,ps_f,z_f,delw_1,deltf,pa,lpa)
                     call intg(ps,ps_f,z_f,delw_2,deltf,pa,lpa)
                  end do
                  call calc(ps,ps_c,z_c,out_c,moms,point)
                  call calc(ps,ps_f,z_f,out_f,moms,point)
                end do
C
C     END PATH
C
        end do
C
C     END ENSEMBLE
C
        sigmamax=0.
        errormax=0.
        errpoint=0
        sigpoint=0
        do  mom=1,moms
             write(*,*) 'MOMENT = ' , mom
             write(*,*) 'Time     n     error(tot)     error(step)'
             write(*,*) '*******************************************'
             do  point=0,points
                nc =out_c(1,mom,point)/paths
                n2c=out_c(2,mom,point)/paths
                nf =out_f(1,mom,point)/paths
                n2f=out_f(2,mom,point)/paths
                error=abs(nc-nf)
                nsigma=0.
                if (paths.gt.1)  then
                        nsigma = sqrt(abs(n2f-nf**2)/(paths -1))
                end if
                if (nsigma.gt.sigmamax) then
                        sigmamax=nsigma
                        sigpoint=point
                endif
                nsigma = error+nsigma
                if (error.gt.errormax) then
                        errormax=error
                        errpoint=point
                endif
                time1=point*deltc*steps
                write(*,700) time1, nf, nsigma, error
                write(fn,700) time1, nf, nsigma, error
             end do
        write(fn,*)
        end do
700             format (4G14.6)
             write(*,*) '*******************************************'
C             timec =  n_fractc*time/(paths*steps*points)
C             timef =  n_fractf*time/(2*paths*steps*points)
             write(*,*) '*SUCCESSFUL FINISH'
             write(*,*) '*COMPLETED              ', paths, ' PATHS'
             write(*,*) '*TOTAL TIME          =  ', time
C             write(*,*) '*TIME ON SURFACE(C)  =  ', timec
C            write(*,*) '*TIME ON SURFACE(F)  =  ', timef
             write(*,*) '*WORST STEP ERROR    =  ', errormax
             write(*,*) '*OCCURRED AT TIME    =  ',errpoint*deltc*steps
             write(*,*) '*WORST STD ERROR    =   ', sigmamax
             write(*,*) '*OCCURRED AT TIME    =  ',sigpoint*deltc*steps
             write(*,*) '*******************************************'
C             write(fn,*) '#TIME ON SURFACE(C)  =  ', timec
C             write(fn,*) '#TIME ON SURFACE(F)  =  ', timef
             write(fn,*) '#WORST STEP ERROR    =  ', errormax
             write(fn,*) '#OCCURRED AT TIME    =  ',errpoint*deltc*steps
             write(fn,*) '#WORST STD ERROR    =   ', sigmamax
             write(fn,*) '#OCCURRED AT TIME    =  ',sigpoint*deltc*steps
             stop
             end
C
C     END PROGRAM
C
*************************************************************************
C       SUBROUTINE: EQU_IN
C       ----------------
C       RETURNS EQUATION PARAMETERS
*************************************************************************
        subroutine  equ_in (ps,lpa,z_i,fn)
C
C       Reference: None
C       Algorithm:  Simple data input routine
C
        complex*16 x,z_i(*),lpa(*)
        integer var,fn,n,ps(10)
        do var=1,ps(5)
        write(*,*) ' _________________________________________________ '
        write(*,*) '|                                                 |'
        write(*,*) '| Variable,',var,'                           |'
        write(*,*) '|                                                 |'
        write(*,*) '|_________________________________________________|'
        call input(2,'complex gain        ',-1.d9,1.d9,lpa(var),n,fn)
        call input(2,'initial amplitude   ',0.d0,1.d9,x,n,fn)
        z_i(var)=x
        end do
        return
        end
*************************************************************************
C       SUBROUTINE: INTG
C       ----------------
C       INTEGRATES ONE TIME STEP
*************************************************************************
        subroutine intg(ps,psi,z,delw,delt,pa,lpa)
C
C       Reference:   Drummond & Mortimer
C       Algorithm:   Iterated central partial difference
C
C
C  PARAMETER VALUES
C
               integer MSEED,MPATH,MPT,MSTEP
               integer MM,MV,MN,MPA
               parameter (MSEED = 100000000)
               parameter (MPATH = 1000000)
               parameter (MPT = 1024)
               parameter (MSTEP = 1024)
               parameter (MM = 10)
               parameter (MV = 10)
               parameter (MN = 10)
               parameter (MPA = 10)
C
C  DATA DECLARATIONS
C
        complex*16 z(MV),z_m(MV),del_z(MV),pa(*),lpa(*),psi
        real*8 delw(MN)
        real*8 delt
        integer ps(10),iter,var
        do var=1,ps(5)
                z_m(var)=z(var)
        end do
        do iter=1,ps(8)
                   call deriv(ps,z_m,delw,del_z,pa)!get deriv
                        do var=1,ps(5)
                                del_z(var)=del_z(var)+lpa(var)*z_m(var)
                                z_m(var)=z(var)+del_z(var)*delt/2.d0
                        end do
        end do
        do var=1,ps(5)
                z(var)=z(var)+del_z(var)*delt
        end do
        return
        end
*************************************************************************
C       SUBROUTINE: CALC
C       ----------------
C       CALCULATES MOMENTS
*************************************************************************
        subroutine  calc(ps,psi,z,out,moms,point)
C
C       Reference: none
C       Algorithm: calculates  moments
C
C
C  PARAMETER VALUES
C
               integer MSEED,MPATH,MPT,MSTEP
               integer MM,MV,MN,MPA
               parameter (MSEED = 100000000)
               parameter (MPATH = 1000000)
               parameter (MPT = 1024)
               parameter (MSTEP = 1024)
               parameter (MM = 10)
               parameter (MV = 10)
               parameter (MN = 10)
               parameter (MPA = 10)
        complex*16 z(MV),psi,zc
        real*8 n1
        real*8 out(2,0:MM,0:MPT)
        integer mom,moms,point,ps(*),var,exp
        do  mom = 1,moms
                exp=1+(mom-1)/ps(4)
                var=mom-(exp-1)*ps(4)
              zc =z(var)
              if (ps(2).eq.2) zc=z(var+ps(4))
              zc=dconjg(zc)
              n1=       ((z(var))*zc )**exp             !compute variable
              out(1,mom,point) = out(1,mom,point)+n1    !store variable
              out(2,mom,point) = out(2,mom,point)+n1**2 !store variance
        end do
        return
        end
*************************************************************************
C       SUBROUTINE: INPUT
C       ----------------
C       RETURNS INPUT PARAMETERS
C       COPYRIGHT:  ALL RIGHTS RESERVED - Peter D. Drummond, 27/04/96
*************************************************************************
        subroutine  input(nmax,name,min,max,r_v,i_v,f)
C
C       Reference: None
C       Algorithm:  Simple data input routine
C
C       WARNING -  KEEP    ML > = MAXD+1  TO AVOID OUT OF BOUNDS
C
        integer len_trim        !FUNCTION TO RETURN TRIMMED STRING LENGTH
        integer n,nmax,ML
        parameter (ML = 10)
        real*8 min,max,r_v(nmax)
        real*4 r4(ML)
        character*20 name
        character*180 data
        integer i_v(nmax),f,iflag,inerr
        iflag = -1
        inerr = 0
        if (nmax.le.0) return
        if (nmax.gt.ML) then
                write(*,*) nmax,' values wanted, only ' ,ML, ' allowed.'
                stop
        end if
        do while (iflag.lt.0.and.inerr.ge.0)
        write(*,*)
     1  ' ____________________________________________________________'
        write(*,*) '| Please input',nmax,' number(s) for ', name,'|'
        write(*,*) '| In range: ',min,'-',max           , '|'
        write(*,*)
     1  '|____________________________________________________________|'
            read(*,'(a180)',iostat=inerr) data
            iflag=1
            n=1
            if (data(1:1).eq.' ') n=2
            if (data(n:n).eq.'#') then
                if (data(n+1:n+1).eq.'t') then
                        write(f,'(a60)',iostat=inerr) data
                        iflag=-1
                end if
                data(n:n) = ' '
            end if
            if (inerr.ge.0.and.iflag.eq.1) then
                read(data,*,iostat=iflag) (r_v(n),n=1,nmax)
                if (iflag.ge.0) then
                  do n=1,nmax
                        if (r_v(n).lt.min) then
                        write(*,*) 'Input too small, you must try again'
                            iflag = -1
                        end if
                        if (r_v(n).gt.max) then
                        write(*,*) 'Input too big, you must try again'
                            iflag = -1
                        end if
                  end do
                end if
            end if
        end do
        do n=1,nmax
                i_v(n)=int(r_v(n))
                r4(n)=r_v(n)
        end do
        write(*,*) 'OK: ',(r4(n),n=1,nmax)
        write(data,*) (r4(n),n=1,nmax)
        write(f,'(1x, 4a)') '#',data(:len_trim(data)),'       ',name
        return
        end
*************************************************************************
C       FUNCTION: LEN_TRIM
C       ------------------
C       RETURN THE LENGTH OF A STRING LESS TRAILING SPACES
*************************************************************************
        integer function len_trim(str)
C
        character*(*) str
        integer l
        if (str .eq. ' ') then
           l = 1
        else
           l = len(str)
           do while (str(l:l) .eq. ' ')
              l = l-1
           end do
        end if
        len_trim = l
        return
        end
*************************************************************************
C       FUNCTION: GAUSS
C       ----------------
C       GAUSSIAN RANDOM NUMBER GENERATOR
C       COPYRIGHT:  ALL RIGHTS RESERVED - Peter D. Drummond, 27/04/96
*************************************************************************
        subroutine gauss(nl,s,b1,b2,bc)
C
C       Reference: Numerical Recipes
C       Algorithm: Box Mueller
C       Function:
C               - returns nl real gaussian random numbers in b1,b2,bc
C               - variance = <b^2> = s(n)^2, mean = 0, in  b1(n),b2(n)
C               - b1,b2 arrays are for the fine lattice calculation!
C               - bc is the average for the coarse lattice calculation!
        integer nl,n
        real*8 s(nl)
        real*8 randmzt,pi2,r1,r2
        real*8 bc(nl)
C               real*8 test1,test2,test12       !DEBUG
        real*8 b1(nl),b2(nl)
        complex*16 c
        pi2=8.*atan(1.d0)
        do  n=1,nl                      !lattice
                r1=randMZT()
                r2=randMZT()
                c=dsqrt(-2.d0*dlog(r1))*cdexp((0.d0,1.d0)*pi2*r2)*s(n)
                b1(n)=dreal(c)
                b2(n)=dimag(c)
                bc(n)=0.5D0*(dreal(c)+dimag(c))
        end do
C       write(*,*) 'Gauss',b1(1,1),b2(1,1),bc(1,1)
C       test1=0.d0
C       test2=0.d0
C       test12=0.d0
C       do n=1,nl
C       test1=test1+b1(n)**2
C       test2=test2+b2(n)**2
C       test12=test12+b1(n)*b2(n)
C       end do
C       write(*,*) ' variance =  ', test1/(nl),test2/(nl),test12/(nl)
C       write(*,*) ' expected =  ', s(1)**2,s(1)**2,'  0.0'
        return
        end
*************************************************************************
C       SUBROUTINE: INPUT
C       ----------------
C       RETURNS INPUT PARAMETERS
*************************************************************************
        subroutine  oinput (name,min,max,r_var,i_var,fn)
C
C       Reference: None
C       Algorithm:  Simple data input routine
C
        real*8 min,max,r_var
        character*12 name, data
        integer i_var,fn,iflag,inerr
        iflag = -1
        inerr = 0
        do while (iflag.lt.0.and.inerr.ge.0)
            write(*,*) 'Input: ',name,' Range:',min,'-',max
c                   write(*,*) '1'
            read(*,'(a12)',iostat=inerr) data
c                   write(*,*) '2',data
            if (data(1:1).eq.'#') data(1:1) = ' '
            if (data(2:2).eq.'#') data(2:2) = ' '
            if (inerr.ge.0) then
                read(data,*,iostat=iflag) r_var
c                       write(*,*) '3',data,r_var
c                       write(*,*) 'input',r_var,min,max
                if (iflag.ge.0) then
                        if (r_var.lt.min) then
                            write(*,*) 'Input too small'
                            iflag = -1
                        end if
                        if (r_var.gt.max) then
                            write(*,*) 'Input too large'
                            iflag = -1
                        end if
                end if
            end if
        end do
        i_var=int(r_var)
        write(*,*) name ,' = ', r_var
        write(fn,*) '# ',  r_var,'       ',name
        return
        end
*************************************************************************
C       FUNCTION: GAUSS
C       ----------------
C       GAUSSIAN RANDOM NUMBER GENERATOR
*************************************************************************
        subroutine ogauss(delw1,delw2,sigma,pi2,noises)
C
C       Reference: Numerical Recipes
C       Algorithm: Box Mueller
C
        real*8 sigma,randmzt,pi2,r1,r2,delw1(*),delw2(*)
        complex*16 c
        integer noise,noises
        do 100 noise=1,noises
        r1=randMZT()
        r2=randMZT()
        c = dsqrt(-2.d0*dlog(r1))*cdexp((0.,1.)*pi2*r2)*sigma
        delw1(noise)=dreal(c)
100     delw2(noise)=dimag(c)
        return
        end
*************************************************************************
C       FUNCTION: RANDOM
C       ----------------
C       PORTABLE PSEUDORANDOM NUMBER GENERATOR
*************************************************************************
        real*8 function randMZT ()
C
C       Reference: F.James, Comput. Phys. Commun. 60 (1990) 329.
C       Algorithm: G.Marsaglia, A.Zaman and W.-W. Tsang,
C               Stat. Prob. Lett. 9 (1990) 35.
