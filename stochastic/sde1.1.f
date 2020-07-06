C
C       COMMON BLOCK
C
            real*8 u(97)
            real*8 c, cd, cm
            integer i97, j97
            common /ranset/ u, c, cd, cm, i97, j97
C
C       VARIABLES
C
            real*8 uni
C
        uni = u(i97) - u(j97)
        if (uni .lt. 0.0) uni = uni + 1.0
        u(i97) = uni
        i97 = i97 - 1
        if (i97 .eq. 0) i97 = 97
        j97 = j97 - 1
        if (j97 .eq. 0) j97 = 97
        c = c-cd
        if (c .lt. 0.0) c = c+cm
        uni = uni-c
        if (uni .lt. 0.0) uni = uni + 1.0
C
C       CODE TO AVOID ZEROES
C
        if (uni .eq. 0.0) then
            uni = u(j97) * 2.0**-24
             if (uni .eq. 0.0) uni = 2.0**-48
        end if
C
        randmzt = uni
        return
        end ! random
************************************************************
C       SUBROUTINE: IN_RAND
C       ------------------
C       INITIALIZING ROUTINE FOR RANDOM
*************************************************************
        subroutine in_rand (ijkl)
C
C       ARGUMENT
C
C         In
            integer ijkl
C
C       Input: 0 <= ijkl <= 900 000 000
C
C       To get standard values in Marsaglia-Zaman paper,
C               (i=12, j=34, k=56, l=78) put ijkl = 54217137
C
C       COMMON BLOCK
C
            real*8 u(97)
            real*8 c, cd, cm
            integer i97, j97
            common /ranset/ u, c, cd, cm, i97, j97
C
C       VARIABLES
C
            integer ij, kl, i, j, k, l, m, ii, jj
            real*8 s, t
C
        ij = ijkl/30082
        kl = ijkl - 30082*ij
        i = mod(ij/177, 177) + 2
        j = mod(ij, 177) + 2
        k = mod(kl/169, 178) + 1
        l = mod(kl, 169)
        do ii = 1, 97
            s = 0.0
            t = 0.5
            do jj = 1, 24
                m = mod(mod(i*j, 179) * k, 179)
                i = j
                j = k
                k = m
                l = mod(53*l + 1, 169)
                if (mod(l*m, 64) .ge. 32) s = s+t
                t = 0.5*t
            end do
            u(ii) = s
        end do
        c = 362436.0/16777216.0
        cd = 7654321.0/16777216.0
        cm = 16777213.0/16777216.0
        i97 = 97
        j97 = 33
        return
        end
*********************************************************************
C
********************************************************************************
C       SUBROUTINE: EQU_PAR
C       ----------------
C       RETURNS EQUATION PARAMETERS
C       COPYRIGHT:  ALL RIGHTS RESERVED - Peter D. Drummond, 27/04/96
********************************************************************************
        subroutine  equ_par (ps,pa,f)
C
C       Reference: None
C       Algorithm:  Simple data input routine
C       ps(1) allows switching the algorithm if necessary
C
        real*8 x(10)
        integer j(10),m
        complex*16 pa(*)
        integer ps(*),f
        write(*,*) ' _________________________________________________ '
        write(*,*) '|                                                 |'
        write(*,*) '| STOCHASTIC EQUATIONS CURRENTLY AVAILABLE:       |'
        write(*,*) '|                                                 |'
        write(*,*) '|                                                 |'
        write(*,*) '|  0 -      Kubo oscillator                       |'
        write(*,*) '|  1 -      Landau-Ginzburg                       |'
        write(*,*) '|  2 -      Degenerate Paramp                     |'
        write(*,*) '|  3 -      Nondegenerate Paramp                  |'
        write(*,*) '|  4 -      Degenerate Gap Paramp                 |'
        write(*,*) '|_________________________________________________|'
        call input(1,'equation            ',0.d0,4.d0,x,ps(1),f)
        write(*,*) ' _________________________________________________ '
        write(*,*) '|                                                 |'
        write(*,*) '| PHASE-SPACE TYPES CURRENTLY AVAILABLE:          |'
        write(*,*) '|                                                 |'
        write(*,*) '|  1 -      Standard (Wigner) space               |'
        write(*,*) '|                                                 |'
        write(*,*) '|  2 -      Doubled (+P)  space                   |'
        write(*,*) '|                                                 |'
        write(*,*) '|_________________________________________________|'
        call input(1,'type                ',1.0d0,2.d0,x,ps(2),f)
        write(*,*) ' _________________________________________________ '
        write(*,*) '|                                                 |'
        write(*,*) '| INTEGRATION OPTIONS CURRENTLY AVAILABLE:        |'
        write(*,*) '|                                                 |'
        write(*,*) '|  0 -      Standard (Stratonovic) sde            |'
        write(*,*) '|                                                 |'
C        write(*,*) '|  1 -     Hybrid (Graham) path integral          |'
        write(*,*) '|                                                 |'
        write(*,*) '|_________________________________________________|'
        call input(1,'    option          ',0.0d0,0.d0,x,ps(3),f)
        call input(1,'    modes           ',1.0d0,4.0d0,x,ps(4),f)
        write(*,*) ' _________________________________________________ '
        write(*,*) '|                                                 |'
        write(*,*) '| STRATONOVIC STOCHASTIC DIFFERENTIAL EQUATION:   |'
        write(*,*) '|                                                 |'
        write(*,*) '|  da_n/dt  = D_n +  L_n[a] + da_n/dt|(NL)        |'
        write(*,*) '|                                                 |'
        write(*,*) '|  L_n[a] represents a generalized linear term.   |'
        write(*,*) '|                                                 |'
        write(*,*) '|                                                 |'
        write(*,*) '| NONLINEAR EQUATION CHOSEN:                      |'
        write(*,*) '|                                                 |'
     |'
        write(*,*) '|   a_n      =  complex variable                  |'
        write(*,*) '|   B_n      =  noise term                        |'
        write(*,*) '|   C_n      =  nonlinear term                    |'
        write(*,*) '|   D_n      =  driving term                      |'
        write(*,*) '|                                                 |'
        write(*,*) '| <B_n1(t1)B_n2(t2)> = b_n1.d_(n1,n2).d(t1-t2) are|'
        write(*,*) '| real delta-correlated stochastic fields.        |'
        write(*,*) '|                                                 |'
        write(*,*) '|  *    =  complex conjg.                         |'
        write(*,*) '|                                                 |'
       if (ps(1).eq.0) then
        ps(4)=1
        ps(7)=2
        write(*,*) '|  0 -      KuBo osCillator:                      |'
        write(*,*) '|                                                 |'
        write(*,*) '|  da1/dt  = C1.a1^2 + a1.(B1 + iB2)              |'
        write(*,*) '|                                                 |'
        if (ps(2).eq.2) then
                write(*,*) 'WARNING: setting  ps(2)=1'
                ps(2)=1
        end if
        else if (ps(1).eq.1.and.ps(2).eq.1) then
         if (ps(4).eq.1) then
        ps(7)=2
        write(*,*) '|  1.1W -      Landau-GinzBurg: (Wigner)          |'
        write(*,*) '|                                                 |'
        write(*,*) '|  da1/dt  = C1.(a1^2.a1*) + B1 +iB2              |'
        write(*,*) '|                                                 |'
         else
                ps(4)=2
                ps(7)=4
        write(*,*) '|  1.2W -      Landau-GinzBurg: (Wigner)          |'
        write(*,*) '|                                                 |'
        write(*,*) '|  da1/dt = 1.1W +i[C2r.a1.|a2|^2 + C2i.a2^2.a1*] |'
        write(*,*) '|                                                 |'
        write(*,*) '|  da2/dt = 1.1W +i[C2r.a2.|a1|^2 + C2i.a1^2.a2*] |'
        write(*,*) '|                                                 |'
         end if
        else if (ps(1).eq.1.and.ps(2).eq.2) then
        ps(4)=1
        ps(7)=2
        write(*,*) '|  1.1P -      Landau-GinzBurg: (+P)              |'
        write(*,*) '|                                                 |'
        write(*,*) '|  da1/dt  = D1 + C1.a2*.a1^2 + a1.Rt(C1).B1      |'
        write(*,*) '|  da2/dt  = D1 + C1.a1*.a2^2 + a2.Rt(C1).B2      |'
        write(*,*) '|                                                 |'
        else if (ps(1).eq.2.and.ps(2).eq.1) then
        write(*,*) '|  2 -      Degenerate Paramp: (Wigner)           |'
        ps(4)=2
        ps(7)=4
        write(*,*) '|                                                 |'
        write(*,*) '|  da1/dt = D1 +C1.a2.a1*  +  B1+iB2              |'
        write(*,*) '|  da2/dt = D2 -C2*.a1^2/2+  B3+iB4               |'
        write(*,*) '|                                                 |'
        else if (ps(1).eq.2.and.ps(2).eq.2) then
        ps(4)=2
        ps(7)=2
        write(*,*) '|  2 -      Degenerate Paramp: (+P)               |'
        write(*,*) '|                                                 |'
        write(*,*) '|  da1/dt = D1 +C1.a2.a3* +  Rt(C1.a2).B1         |'
        write(*,*) '|  da2/dt = D2 -C2*.a1^2/2                        |'
        write(*,*) '|  da3/dt = D1 +C1.a4.a1* + Rt(C1.a4).B2          |'
        write(*,*) '|  da4/dt = D2 -C2*.a3^2/2                        |'
        write(*,*) '|                                                 |'
        else if (ps(1).eq.3.and.ps(2).eq.1) then
        ps(4)=3
        ps(7)=6
        write(*,*) '|  3 -      Nondegenerate Paramp: (Wigner)        |'
        write(*,*) '|                                                 |'
        write(*,*) '|  da1/dt = D1 + C1.a3.a2*  +  B1+iB2             |'
        write(*,*) '|  da2/dt = D2 + C1.a3.a1*  +  B3+iB4             |'
        write(*,*) '|  da3/dt = D3 - C2*.a1.a2 +  B5+iB6              |'
        write(*,*) '|                                                 |'
        else if (ps(1).eq.3.and.ps(2).eq.2) then
        ps(4)=3
        ps(7)=4
        write(*,*) '|  3 -      Nondegenerate Paramp: (+P)            |'
        write(*,*) '|                                                 |'
        write(*,*) '|  da1/dt = D1 +C1.a3.a5* +  Rt(C1.a3).(B1+iB2)   |'
        write(*,*) '|  da2/dt = D2 +C1.a3.a4* +  Rt(C1.a3).(B1-iB2)   |'
        write(*,*) '|  da3/dt = D3 -C2*.a1.a2                         |'
        write(*,*) '|  da4/dt = D1 +C1.a6.a2* +  Rt(C1.a6).(B3+iB4)   |'
        write(*,*) '|  da5/dt = D2 +C1.a6.a1* +  Rt(C1.a6).(B3-iB4)   |'
        write(*,*) '|  da6/dt = D3 -C2*.a4.a5                         |'
        write(*,*) '|                                                 |'
        else if (ps(1).eq.4.and.ps(2).eq.1) then
        ps(4)=4
        ps(7)=8
        write(*,*) '|  4 -      Degenerate Gap Paramp: (Wigner)       |'
        write(*,*) '|                                                 |'
        write(*,*) '|  da1/dt = D1 + C1.a2.a1*   + C3.a3  + B1+iB2    |'
        write(*,*) '|  da2/dt = D2 -C2*.a1^2/2 + C4.a4  + B3+iB4      |'
        write(*,*) '|  da3/dt = D3 + C1.a4.a3*   - C3*.a1  + B5+iB6   |'
        write(*,*) '|  da4/dt = D4 -C2*.a3^2/2 - C4*.a2  + B7+iB8     |'
        write(*,*) '|                                                 |'
        else if (ps(1).eq.4.and.ps(2).eq.2) then
        ps(4)=4
        ps(7)=4
        write(*,*) '|  4 -      Degenerate Gap Paramp: (+P)           |'
        write(*,*) '|                                                 |'
        write(*,*) '|  da1/dt = D1 + C1.a2.a5*  + C3.a3 + Rt(C1.a2).B1|'
        write(*,*) '|  da2/dt = D2 -C2*.a1^2/2 + C4.a4                |'
        write(*,*) '|  da3/dt = D3 + C1.a4.a7*  - C3*.a1 +Rt(C1.a4).B2|'
        write(*,*) '|  da4/dt = D4 -C2*.a3^2/2 - C4*.a2               |'
        write(*,*) '|  da5/dt = D1 + C1.a6.a1*  + C3.a7 + Rt(C1.a6).B3|'
        write(*,*) '|  da6/dt = D2 -C2*.a5^2/2 + C4.a8                |'
        write(*,*) '|  da7/dt = D3 + C1.a4.a3*  - C3*.a5 +Rt(C1.a8).B4|'
        write(*,*) '|  da8/dt = D4 -C2*.a7^2/2 - C4*.a6               |'
        else
        write(*,*) '|  OPTIONS CHOSEN NOT CURRENTLY AVAILABLE         |'
        stop
        end if
        write(*,*) '|_________________________________________________|'
        ps(5)=ps(2)*ps(4)
        do m=1,ps(4)
        call input(2,'driving  [D_n]      ',-1.d9,1.d9,pa(m),j,f)
        end do
        do m=1,ps(4)
        call input(2,'coupling [C_n]      ',-1.d9,1.d9,pa(m+ps(4)),j,f)
        end do
        return
        end
*************************************************************************
C       SUBROUTINE: DERIV
C       ----------------
C       RETURNS STOCHASTIC DIFFERENTIAL EQUATION COEFFICIENTS
C
C       Reference: Drummond, Chaturvedi, Walls
C       Algorithm:   Two-photon absorption, Stratonovich
C       COPYRIGHT:  ALL RIGHTS RESERVED - Peter D. Drummond, 27/04/96
*******************************************************************************
C
C  DOCUMENTATION
C
********************************************************************************
C
C       IMPORT:
C
C       ps      set of program switches, defining algorithm
C       a       complex field vector, at center of interval in time
C       bc      real noise field vector for SDE
C       pa      set of program parameters for equation
C
C       EXPORT:
C
C       da      complex derivative for the a-field (excluding linear terms)
C
C***************************************************************************
*****
        subroutine deriv(ps,a,bc,da,pa)
        complex*16 pa(*),c,dc,ci,z1,z2,c2i,c2r
        integer ps(*),v
        complex*16  a(*),da(*)
        real*8     bc(*)
C                       write(*,*) 'PAUSE:DERIV',(ps(v), v=1,9)
C       read (*,*)
        ci=(0.d0,1.d0)
        c=pa(1+ps(4))
        dc=dconjg(pa(2+ps(4)))
        IF (ps(1).eq.0) then                    !KUBO SINGLE MODE
           da(1)=a(1)*(bc(1)+(0.d0,1.d0)*bc(2))! stochastic input
           da(1)=pa(1)+c*a(1)**2+da(1)
         ELSE IF (ps(1).eq.1) then              !LANDAU-GINZBURG SINGLE MODE
         if (ps(2).lt.2) then                           !Wigner
           if (ps(4).eq.1) then
           da(1) = pa(1)+c*a(1)**2*dconjg(a(1))
           da(1)= da(1)+ bc(1)+ci*bc(2)! stochastic input
           else                                 !LANDAU-GINZBURG TWIN MODE
           c2r=ci*dreal(pa(2+ps(4)))
           c2i=ci*dimag(pa(2+ps(4)))
                                                !Wigner
           da(1) = pa(1)+c*a(1)**2*dconjg(a(1))
           da(1) = da(1)+c2r*a(1)*a(2)*dconjg(a(2))
           da(1) = da(1)+c2i*a(2)**2*dconjg(a(1))
           da(1) = da(1)+ bc(1)+ci*bc(2)! stochastic input
           da(2) = pa(2)+c*a(2)**2*dconjg(a(2))
           da(2) = da(2)+c2r*a(2)*a(1)*dconjg(a(1))
           da(2) = da(2)+c2i*a(1)**2*dconjg(a(2))
           da(2)= da(2)+ bc(3)+ci*bc(4)! stochastic input
           end if
          else                                  !LANDAU-GINZBURG SINGLE MODE
                                        !+P
           da(1) = pa(1)+c*a(1)**2*dconjg(a(2))
           da(2) = pa(1)+c*a(2)**2*dconjg(a(1))
           da(1) = da(1)+cdsqrt(pa(1))*a(1)*bc(1)       !
           da(2) = da(2)+cdsqrt(pa(1))*a(2)*bc(2)
          end if
        ELSE IF (ps(1).eq.2) then               !DEGENERATE PARAMETRIC
OSCILLATOR
          if (ps(2).lt.2) then                          !Wigner
                da(1) = pa(1)+c*a(2)*dconjg(a(1))
                da(2) = pa(2)-0.5d0*dc*a(1)**2
                da(1)= da(1)+bc(1)+ci*bc(2)! stochastic term
                da(2)= da(2)+bc(3)+ci*bc(4)! stochastic term
          else                                          !+P
                da(1) = pa(1)+c*a(2)*dconjg(a(3))
                da(1) = da(1)+sqrt(c*a(2))*bc(1)
                da(2) = pa(2)-dc*a(1)**2/2.d0
                da(3) = pa(1)+c*a(4)*dconjg(a(1))
                da(3) = da(3)+sqrt(c*a(4))*bc(2)
                da(4) = pa(2)-dc*a(3)**2/2.d0
         end if
         ELSE IF        (ps(1).eq.3) then       !NON-DEGENERATE PARAMETRIC
OSCILLATOR
          if (ps(2).lt.2) then                          !Wigner
                da(1) = pa(1)+c*a(3)*dconjg(a(2))
                da(1) = da(1)+bc(1)+ci*bc(2)! stochastic term
                da(2) = pa(2)+c*a(3)*dconjg(a(1))
                da(2) = da(2)+bc(3)+ci*bc(4)! stochastic term
                da(3) = pa(3)-dc*a(1)*a(2)
                da(3) = da(3)+bc(5)+ci*bc(6)! stochastic term
          else                                          !+P
                z1=bc(1)+ci*bc(2)
                z2=bc(3)+ci*bc(4)
                da(1) = pa(1)+c*a(3)*dconjg(a(5))
                da(1) = da(1)+sqrt(c*a(3))*z1
                da(2) = pa(2)+c*a(3)*dconjg(a(4))
                da(2) = da(2)+sqrt(c*a(3))*dconjg(z1)
                da(3) = pa(3)-dc*a(1)*a(2)
                da(4) = pa(1)+c*a(6)*dconjg(a(2))
                da(4) = da(4)+sqrt(c*a(6))*z2
                da(5) = pa(2)+c*a(6)*dconjg(a(1))
                da(5) = da(5)+sqrt(c*a(6))*dconjg(z2)
                da(6) = pa(3)-dc*a(4)*a(5)
         end if
         ELSE IF        (ps(1).eq.4) then               !Gap
          if (ps(2).lt.2) then                          !Wigner
            da(1) = pa(1)+ pa(7)*a(3)+c*a(2)*dconjg(a(1))
            da(1) = da(1)+bc(1)+ci*bc(2)! stochastic term
            da(2) = pa(2)+pa(8)*a(4)-0.5d0*dc*a(1)**2
            da(2) = da(2)+bc(3)+ci*bc(4)! stochastic term
            da(3) = pa(3)-dconjg(pa(7))*a(1)+a(4)*c*dconjg(a(3))
            da(3) = da(3)+bc(5)+ci*bc(6)! stochastic term
            da(4) = pa(4)-dconjg(pa(8))*a(2)-0.5d0*dc*a(3)**2
            da(4) = da(4)+bc(7)+ci*bc(8)! stochastic term
          else                                          !+P
            da(1)= pa(1)+pa(7)*a(3)+c*a(2)*dconjg(a(5))
            da(1)= da(1)+sqrt(pa(1)*a(2))*bc(1)
            da(2)= pa(2)+pa(8)*a(4)-dc*a(1)**2/2.d0
            da(3)= pa(3)-dconjg(pa(7))*a(1)+c*a(4)*dconjg(a(7))
            da(3)= da(3)+sqrt(pa(2)*a(4))*bc(2)
            da(4)= pa(4)-dconjg(pa(8))*a(2)-0.5d0*dc*a(3)**2
            da(5)= pa(1)+pa(7)*a(7)+c*a(6)*dconjg(a(1))
            da(5)= da(5)+sqrt(pa(1)*a(6))*bc(3)
            da(6)= pa(2)+pa(8)*a(8)-dc*a(5)**2/2.d0
            da(7)= pa(3)-dconjg(pa(7))*a(5)+c*a(8)*dconjg(a(3))
            da(7)= da(7)+sqrt(pa(2)*a(8))*bc(4)
            da(8)= pa(4)-dconjg(pa(8))*a(6)-0.5d0*dc*a(7)**2
          end if
        ELSE IF (ps(1).eq.11) then                      ! Degen. Paramp+subtract
          if (ps(2).lt.2) then                          !Wigner subtraction
                da(1) = pa(1)+c*a(2)*dconjg(a(1))
                da(2) = pa(2)-0.5d0*dc*a(1)**2
                da(1)= da(1)+bc(1)+ci*bc(2)! stochastic term
                da(2)= da(2)+bc(3)+ci*bc(4)! stochastic term
                da(3) = pa(1)+a(4)*c*dconjg(a(3))!simplified
                da(3)= da(3)+bc(1)+ci*bc(2)! stochastic term
                da(4)= pa(2)+bc(3)+ci*bc(4)! stochastic term
C        write(*,*) 'parameter= ', pa(1),pa(2),pa(3)
C        write(*,*) 'z= ', a(1),a(2),a(3),a(4)
C        write(*,*) 'ad= ', da(1),da(2),da(3)
          else                                          !+P subtraction
                da(1) = pa(1)+c*a(2)*dconjg(a(5))
                da(1) = da(1)+sqrt(c*a(2))*bc(1)
                da(2) = pa(2)-dc*a(1)**2/2.d0
                da(3) = pa(1)+c*a(4)*dconjg(a(7))
                da(3) = da(3)+sqrt(c*a(4))*bc(1)
                da(4) = pa(2)
                da(5) = pa(1)+c*a(6)*dconjg(a(1))
                da(5) = da(5)+sqrt(c*a(6))*bc(2)
                da(6) = pa(2)-dc*a(5)**2/2.d0
                da(7) = pa(1)+c*a(8)*dconjg(a(3))
                da(7) = da(7)+sqrt(c*a(8))*bc(2)
                da(8) = pa(2)
C        write(*,*) 'parameter= ', pa(1),pa(2),pa(3)
C        write(*,*) 'z= ', a(1),a(2),a(3)
C        write(*,*) 'ad= ', da(1),da(2),da(3)
         end if
        END IF
D        write(*,*) 'DERIV: PROGRAM SWITCHES =',ps(1),ps(2),ps(3),ps(4),ps(5)
D        write(*,*) 'parameter= ', pa(1),pa(2),pa(3)                    !DEBUG
D        write(*,*) 'a1-4= ', a(1),a(2),a(3),a(4)               !DEBUG
D        write(*,*) 'bc1-2= ', bc(1),bc(2)                              !DEBUG
D        write(*,*) 'da1-4= ', da(1),da(2),da(3),da(4)  !DEBUG
        return
        end
