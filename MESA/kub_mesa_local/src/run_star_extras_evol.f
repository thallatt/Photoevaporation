! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      
      implicit none
      
      integer :: time0, time1, clock_rate

      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% max_timestep = 10**(s% x_ctrl(7))*secyer
         s% dt_next = 10**(s% x_ctrl(7))*secyer   !set as thermal timescale  s% cgrav(1)*s% mstar/s% r(1)*s% L(1)]

         !turn on other mass loss routine and other energy
         s% other_adjust_mdot=> <<escape_model>>_mdot !EL_mdot !ELR_mdot !HD_mdot ! HD_highmass_mdot
         s% other_energy => varTEQ_energy !energy_rogers !
      end subroutine extras_controls

      !!!!!!!!mass loss!!!!!!!do not forget to change the pointer!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine EL_mdot(id, ierr)
         !use star_def
         !use const_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: eps_EUV, a, REUV, FEUV
         double precision :: Rhill
         double precision :: Ktide
         double precision :: sigma
         double precision :: Tphoto
         double precision :: mH,muH
         double precision :: Pphoto,g,H,Peuv,Pratio, zHeight, Rp, Rs,rhoS,rhoB0,rhoBP,eLim,rrLim,hv,alphaREC,Vs

         type (star_info), pointer :: s

         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         
         eps_EUV = s% x_ctrl(2)  !mass loss efficiency factor
         a = s% x_ctrl(1)        !orbital separation (AU)
         !REUV = s% x_ctrl(3) * s% photosphere_r * Rsun   !REUV (in cm)

         ! calcuate the radius at which EUV is absorbed
         sigma=2e-18 !sigma0*(20/13.6)^-3, sigma0 = 6e-18
         Tphoto= 10**(s% log_surface_temperature)
         mH= 1.67e-24  
	 muH=1.00794 !g/mol molar mass!!! why?
         Pphoto= 10**(s% log_surface_pressure)  !dyn/cm2
        
        
         hv= 3.2e-11 ! ergs, is 20 eV
         alphaREC=2.7e-13 !cm^3/s

         g=(s% cgrav(1)*s% mstar)/((Rsun*(10**(s% log_surface_radius)))**2)
         Peuv=(mH*g)/(sigma)  !convert from dyn/cm2 to bars !? it is not a convertion
         H=(kerg*Tphoto)/(2*mH*g)    !in cm
         Pratio=LOG(Peuv/Pphoto)
         zHeight=-H*Pratio
         Rp=Rsun*(10**(s% log_surface_radius))
         Reuv=Rp+zHeight

         !Energy Limited 
         !Feuv = <<escape_Cpow>> * (s% star_age * 1.0e-9)**(<<escape_betapow>>)*(s% x_ctrl(1))**(-2)  !29.12*((s% star_age)/1d9)**(-1.24)*(s% x_ctrl(1))**(-2) !+1d6 the shift
         !if(FEUV > <<escape_sat>>) FEUV = <<escape_sat>>
         !if(Feuv>8.85d4) Feuv = 8.85d4
         if((s% star_age) < <<t_sat>>) then
            Feuv=3585.743422231086*(s% x_ctrl(1))**(-2) + 1048.7394601572157*(s% x_ctrl(1))**(-2)
         else
            Feuv= (s% x_ctrl(1))**(-2)*(<<a_x>>*((s% star_age)/<<t_sat>>)**(<<b_x>>) + <<a_euv>>*((s% star_age)/<<t_sat>>)**(<<b_euv>>))
            !Feuv= (s% x_ctrl(2))**(-2)*(3522.3247240527944*((s% star_age)/4d6)**(-0.9) + 1027.20279*((s% star_age)/4d6)**(-1.1)) !3685.1546645010812*(s% x_ctrl(2))**(-2)*((s% star_age)/226d6)**(-2.15) + 1059.4107025353496*(s% x_ctrl(2))**(-2)*((s% star_age)/226d6)**(-2.5)
         endif

         Rhill=(s% x_ctrl(1)*1.495978921d13*((s% mstar)/(3*Msun))**(1/3))
         Ktide=1-(3*Reuv)/(2*Rhill)+ 1/(2*(Rhill/Reuv)**3)
         eLim = -(eps_EUV*pi*Feuv*(Reuv)**3)/(s% cgrav(1)*s% mstar*Ktide)

         
         !Set stellar mass loss rate with E-limited mass loss (in g/s)
         s% mstar_dot = eLim 
         

         

      end subroutine EL_mdot

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine HD_mdot(id, ierr)
         !use star_def
         !use const_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: a, Rpl_RE, Mpl_ME, Teq_K, FEUV, lb, bt0, lhy !a==d0 lb==Lambda critical
         real(dp) :: zeta, eta, beta, alp1, alp2, alp3, K, C !coefficients of approximation
         type (star_info), pointer :: s

         ierr = 0

         ! get the star_info pointer using id
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !xctrl inputs 
         a = s% x_ctrl(1)        !orbital separation (AU)
         
         Rpl_RE = (10.**(s% log_surface_radius)) * 6.96d10 / 6.378d8
         Mpl_ME = s% star_mass * 1.9885d33 / 5.9722d27
         Teq_K = 10.**(s% log_surface_temperature)

         !FEUV flux from star (incident on the planet)
         ! Tu et al. 2015 
         !Feuv = <<escape_Cpow>> * (s% star_age * 1.0e-9)**(<<escape_betapow>>)*(s% x_ctrl(1))**(-2) ! erg cm^-2 s^-1
         if((s% star_age) < <<t_sat>>) then
            Feuv=3585.743422231086*(s% x_ctrl(1))**(-2) + 1048.7394601572157*(s% x_ctrl(1))**(-2)
         else
            Feuv= (s% x_ctrl(1))**(-2)*(<<a_x>>*((s% star_age)/<<t_sat>>)**(<<b_x>>) + <<a_euv>>*((s% star_age)/<<t_sat>>)**(<<b_euv>>))
            !Feuv= (s% x_ctrl(2))**(-2)*(3522.3247240527944*((s% star_age)/4d6)**(-0.9) + 1027.20279*((s% star_age)/4d6)**(-1.1)) !3685.1546645010812*(s% x_ctrl(2))**(-2)*((s% star_age)/226d6)**(-2.15) + 1059.4107025353496*(s% x_ctrl(2))**(-2)*((s% star_age)/226d6)**(-2.5)
         endif
         !29.12 * (s% star_age * 1.0e-9)**(-1.24) / a**2 /2. 
         !if(FEUV > <<escape_sat>>) FEUV = <<escape_sat>>

         !here we implement approximation
         
         zeta = - 1.297796148718774 + 6.861843637445744
         eta  = 0.884595403184073 + 0.009459807206476

         beta = 32.019929264625155 - 16.408393523348366
         alp1 = 0.422232254541188 - 1.0
         alp2 = - 1.748858849270155 + 3.286179370395197
         alp3 = 3.767941293231585 - 2.75

         K = (zeta + eta * LOG(a))
         C = beta + alp1 * LOG(FEUV) + alp2 * LOG(a) + alp3 * LOG(Rpl_RE)

         lb = EXP(C / K);

         !BETA0
         bt0 = 6.6726d-8*Mpl_ME * 5.9722d27/ (Rpl_RE * 6.378d8* (1.3807d-16* Teq_K/ 1.6726d-24))

         !!! TestCF         
         !!! TestCFJ      


         IF (bt0 .LE. lb) THEN
            zeta = -6.861843637445744
            eta  = -0.009459807206476

            beta = 32.019929264625155
            alp1 = 0.422232254541188
            alp2 = -1.748858849270155
            alp3 = 3.767941293231585

            K = (zeta + eta * LOG(a))
            C = beta + alp1 * LOG(FEUV) + alp2 * LOG(a) + alp3 * LOG(Rpl_RE)
            lhy = exp(C + K * LOG(bt0));!testCFJ
         ELSE
            zeta = -1.297796148718774
            eta  = 0.884595403184073

            beta = 16.408393523348366
            alp1 = 1.0
            alp2 = -3.286179370395197
            alp3 = 2.75

            K = (zeta + eta * LOG(a)) 
            C = beta + LOG(FEUV) + alp2 * LOG(a) + alp3 * LOG(Rpl_RE)

            lhy = exp(C + K * LOG(bt0)) !testCF
         END IF
         

         !Set stellar mass loss rate with HBA mass loss (in g/s)
         s% mstar_dot = - lhy

        !write(*,*) 'REUV (cm), FEUV (erg/cm/s), mstar_dot (g/s)', REUV, FEUV, s% mstar_dot

      end subroutine HD_mdot


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! for >40 Mearth planets; low-Lambda planets become energy-limited
      subroutine HD_highmass_mdot(id, ierr)
         !use star_def
         !use const_def
         ! from energy limited
         double precision :: Rhill
         double precision :: Ktide
         double precision :: sigma
         double precision :: Tphoto
         double precision :: mH,muH
         double precision :: Pphoto,g,H,Peuv,Pratio, zHeight, Rp, Rs,rhoS,rhoB0,rhoBP,eLim,rrLim,hv,alphaREC,Vs

         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: a, Rpl_RE, Mpl_ME, Teq_K, FEUV, lb, bt0, lhy, eps_EUV, REUV !a==d0 lb==Lambda critical
         real(dp) :: zeta, eta, beta, alp1, alp2, alp3, K, C, thresh !coefficients of approximation ! thresh=mass threshold
         type (star_info), pointer :: s

         ierr = 0

         ! get the star_info pointer using id
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !xctrl inputs 
         a = s% x_ctrl(1)        !orbital separation (AU)
         eps_EUV = s% x_ctrl(2)  !mass loss efficiency factor
         
         Rpl_RE = (10.**(s% log_surface_radius)) * 6.96d10 / 6.378d8
         Mpl_ME = s% star_mass * 1.9885d33 / 5.9722d27
         Teq_K = 10.**(s% log_surface_temperature)

         ! calcuate the radius at which EUV is absorbed
         sigma=2e-18 !sigma0*(20/13.6)^-3, sigma0 = 6e-18
         Tphoto= 10**(s% log_surface_temperature)
         mH= 1.67e-24  
	 muH=1.00794 !g/mol molar mass!!! why?
         Pphoto= 10**(s% log_surface_pressure)  !dyn/cm2
        
        
         hv= 3.2e-11 ! ergs, is 20 eV
         alphaREC=2.7e-13 !cm^3/s

         g=(s% cgrav(1)*s% mstar)/((Rsun*(10**(s% log_surface_radius)))**2)
         Peuv=(mH*g)/(sigma)  !convert from dyn/cm2 to bars !? it is not a convertion
         H=(kerg*Tphoto)/(2*mH*g)    !in cm
         Pratio=LOG(Peuv/Pphoto)
         zHeight=-H*Pratio
         Rp=Rsun*(10**(s% log_surface_radius))
         Reuv=Rp+zHeight

         Rhill=(s% x_ctrl(1)*1.495978921d13*((s% mstar)/(3*Msun))**(1/3))
         Ktide=1-(3*Reuv)/(2*Rhill)+ 1/(2*(Rhill/Reuv)**3)

         !FEUV flux from star (incident on the planet)
         ! Tu et al. 2015 
         !Feuv = <<escape_Cpow>> * (s% star_age * 1.0e-9)**(<<escape_betapow>>)*(s% x_ctrl(1))**(-2) ! erg cm^-2 s^-1
         if((s% star_age) < <<t_sat>>) then
            Feuv=3585.743422231086*(s% x_ctrl(1))**(-2) + 1048.7394601572157*(s% x_ctrl(1))**(-2)
         else
            Feuv= (s% x_ctrl(1))**(-2)*(<<a_x>>*((s% star_age)/<<t_sat>>)**(<<b_x>>) + <<a_euv>>*((s% star_age)/<<t_sat>>)**(<<b_euv>>))
            !Feuv= (s% x_ctrl(2))**(-2)*(3522.3247240527944*((s% star_age)/4d6)**(-0.9) + 1027.20279*((s% star_age)/4d6)**(-1.1)) !3685.1546645010812*(s% x_ctrl(2))**(-2)*((s% star_age)/226d6)**(-2.15) + 1059.4107025353496*(s% x_ctrl(2))**(-2)*((s% star_age)/226d6)**(-2.5)
         endif
         !29.12 * (s% star_age * 1.0e-9)**(-1.24) / a**2 /2. 
         !if(FEUV > <<escape_sat>>) FEUV = <<escape_sat>>

         !here we implement approximation
         
         zeta = - 1.297796148718774 + 6.861843637445744
         eta  = 0.884595403184073 + 0.009459807206476

         beta = 32.019929264625155 - 16.408393523348366
         alp1 = 0.422232254541188 - 1.0
         alp2 = - 1.748858849270155 + 3.286179370395197
         alp3 = 3.767941293231585 - 2.75

         K = (zeta + eta * LOG(a))
         C = beta + alp1 * LOG(FEUV) + alp2 * LOG(a) + alp3 * LOG(Rpl_RE)

         lb = EXP(C / K);

         !BETA0
         bt0 = 6.6726d-8*Mpl_ME * 5.9722d27/ (Rpl_RE * 6.378d8* (1.3807d-16* Teq_K/ 1.6726d-24))
         thresh=4d1
         !!! TestCF         
         !!! TestCFJ      

         IF (Mpl_ME .GT. thresh) THEN
            print *, 'GREATER'
         ELSE
            print *, 'LESSEQ'
         END IF
         IF (Mpl_ME .GT. thresh) THEN
            IF (bt0 .LE. lb) THEN
               print *, 'Lambda is less than critical lambda, lambda is', bt0
               zeta = -6.861843637445744
               eta  = -0.009459807206476
            
               beta = 32.019929264625155
               alp1 = 0.422232254541188
               alp2 = -1.748858849270155
               alp3 = 3.767941293231585
            
               K = (zeta + eta * LOG(a))
               C = beta + alp1 * LOG(FEUV) + alp2 * LOG(a) + alp3 * LOG(Rpl_RE)
               lhy = exp(C + K * LOG(bt0));!testCFJ
            ELSE
               print *, 'Lambda is greater than critical lambda, lambda is', bt0
               print *, 'critical lambda is ', lb
               eLim = (eps_EUV*pi*Feuv*(Reuv)**3)/(s% cgrav(1)*s% mstar*Ktide) ! positive sign so negative below
               lhy = eLim
            END IF
         ELSE
            IF (bt0 .LE. lb) THEN
               print *, 'Lambda is less than critical lambda, lambda is', bt0
               zeta = -6.861843637445744
               eta  = -0.009459807206476
            
               beta = 32.019929264625155
               alp1 = 0.422232254541188
               alp2 = -1.748858849270155
               alp3 = 3.767941293231585
            
               K = (zeta + eta * LOG(a))
               C = beta + alp1 * LOG(FEUV) + alp2 * LOG(a) + alp3 * LOG(Rpl_RE)
               lhy = exp(C + K * LOG(bt0));!testCFJ
            ELSE
               print *, 'Lambda is greater than critical lambda, lambda is', bt0
               zeta = -1.297796148718774
               eta  = 0.884595403184073
            
               beta = 16.408393523348366
               alp1 = 1.0
               alp2 = -3.286179370395197
               alp3 = 2.75
            
               K = (zeta + eta * LOG(a)) 
               C = beta + LOG(FEUV) + alp2 * LOG(a) + alp3 * LOG(Rpl_RE)
            
               lhy = exp(C + K * LOG(bt0)) !testCF
            END IF
         END IF
         !Set stellar mass loss rate with HBA mass loss (in g/s)
         s% mstar_dot = - lhy

        !write(*,*) 'REUV (cm), FEUV (erg/cm/s), mstar_dot (g/s)', REUV, FEUV, s% mstar_dot

      end subroutine HD_highmass_mdot

      !!!!!!!!!!!!!!!!!other heating!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine varTEQ_energy(id, ierr)
         use const_def, only: Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer :: k
         real(dp) :: Rpl_RE, Mpl_ME
         real(dp) :: ei_time, ei_lbol, ei_teff, ei_rstar, sage, Lstar

         !use star_lib, only: star_ptr
	 double precision :: extraheat,junk,diag_heat
         integer :: i,n,counter,z,p,numOfAcc,zeta,jafter,jbefore,indexI
	 double precision :: core_epsrad,Rpl,pressureDep,pressureHere,random_dp
	 real(dp) :: tauHere,Vesc, KE,massTot,Cd,phi,Mpl,Rhopl,H,g,mH,Tacc
	 !real(dp), DIMENSION(30000000) :: readR,readM,readT
         !REAL(dp), DIMENSION(:), ALLOCATABLE :: readR,readM,readT
         real(dp), DIMENSION(700) :: arrayKE
         real(dp), DIMENSION(5000) :: arrayI

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
               
         s% extra_heat(:) = s% extra_power_source

         open(unit = 12, file = <<heating_input>>, status = 'old', action = "read")
         !log10(age [years]) log10(Lbol [Lsun]) log10(Teff [K]) log10(Rstar [Rsun])

         ei_time = 0.0
         sage = s% star_age !in years

         DO WHILE (sage .gt. 10.0 ** ei_time)
            read(12, *) ei_time, ei_lbol, ei_teff, ei_rstar
            
         END DO

110      format(4(3X, F10.6))
         

         Lstar = (10.0**ei_lbol) * 3.9d33/4.0/s% x_ctrl(1)/s% x_ctrl(1)/2.238d26/3.1416
         
         
         close(unit = 12)

         	

         do k = 1, 100 !s% nz            
            s% irradiation_heat(k) = 1.* 3.1416* Lstar * ((exp(s% lnr(k)))**2) / sum(s% dm(1:100)) !1.36d8 !*2.0 !Lstar * 90. !
                        
         end do

         ! EXTRA HEATING CONSTANT IN SPACE AND TIME
         ! Heat due to Thermal inertia of the core
         ! x_ctrl(4) = cv in erg/K/g from inlist

         k = s% nz
         extraheat = -s% x_ctrl(4) * s% M_center * s% T(s% nz) * s% dlnT_dt(s% nz) / s% dm(s% nz) ! erg/g/sec
           !assuming dlnT_dt is in s^-1

         ! EXTRA HEATING CONSTANT IN SPACE AND TIME
         ! Heat produced by radioactive decay due to 40K, 235U, 238U, 232Th, respectively

         k = s% nz
         core_epsrad = 36.7d-8 * exp(-5.543d-10 * s% star_age) ! erg/g/sec   40K
         core_epsrad = core_epsrad + 4.6d-8 * exp(-9.8485d-10 * s% star_age) ! erg/g/sec  235U
         core_epsrad = core_epsrad + 2.3d-8 * exp( -1.5513d-10 * s% star_age)! erg/g/sec  238U
         core_epsrad = core_epsrad + 1.3d-8 * exp( -0.4948d-10 * s% star_age)! erg/g/sec  232Th

         s% extra_heat(k) = (extraheat+ s% x_ctrl(5) * s% M_center * core_epsrad / s% dm(k)) ! erg/g/sec, core heat flux density


         return
      end subroutine varTEQ_energy


      subroutine energy_rogers(id, ierr)
         type (star_info), pointer :: s
         !use const_def, only: Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
	 
	 !use star_lib, only: star_ptr
	 double precision :: extraheat,junk,diag_heat
         integer :: k,i,n,counter,z,p,numOfAcc,zeta,jafter,jbefore,indexI
	 double precision :: core_epsrad,Rpl,pressureDep,pressureHere,random_dp
	 real(dp) :: tauHere,Vesc, KE,massTot,Cd,phi,Mpl,Rhopl,H,g,mH,Tacc
	 !real(dp), DIMENSION(30000000) :: readR,readM,readT
         !REAL(dp), DIMENSION(:), ALLOCATABLE :: readR,readM,readT
         real(dp), DIMENSION(700) :: arrayKE
         real(dp), DIMENSION(5000) :: arrayI


         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         ierr = 0


	 ! INITALIZE
         s% extra_heat = 0.d0
         s% d_extra_heat_dlndm1(:) = 0d0
         s% d_extra_heat_dlntm1(:) = 0d0

         ! EXTRA HEATING CONSTANT IN SPACE AND TIME
         ! Heat due to Thermal inertia of the core
         ! x_ctrl(3) = cv in erg/K/g from inlist

         k = s% nz
         extraheat = -s% x_ctrl(4) * s% M_center * s% T(s% nz) * s% dlnT_dt(s% nz) / s% dm(s% nz) ! erg/g/sec
           !assuming dlnT_dt is in s^-1

         ! EXTRA HEATING CONSTANT IN SPACE AND TIME
         ! Heat produced by radioactive decay due to 40K, 235U, 238U, 232Th, respectively

         k = s% nz
         core_epsrad = 36.7d-8 * exp(-5.543d-10 * s% star_age) ! erg/g/sec   40K
         core_epsrad = core_epsrad + 4.6d-8 * exp(-9.8485d-10 * s% star_age) ! erg/g/sec  235U
         core_epsrad = core_epsrad + 2.3d-8 * exp( -1.5513d-10 * s% star_age)! erg/g/sec  238U
         core_epsrad = core_epsrad + 1.3d-8 * exp( -0.4948d-10 * s% star_age)! erg/g/sec  232Th

         s% extra_heat(k) = (extraheat+ s% x_ctrl(5) * s% M_center * core_epsrad / s% dm(k)) ! erg/g/sec, core heat flux density




      end subroutine energy_rogers



      subroutine ELR_mdot(id, ierr)
         use star_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         real(dp) :: rsol
         real(dp) :: msol

         double precision :: eps0
         double precision :: Rhill
         double precision :: Ktide
         double precision :: mH,muH
         double precision :: sigma
         double precision :: Tphoto
         double precision :: Pphoto,g,H,Peuv,Pratio,Reuv,zHeight, Rp, Rs,rhoS,rhoB0,rhoBP,eLim,rrLim,hv, alphaREC,Feuv,Vs

         integer :: k,i,j,numTass
         real(dp) :: tauHere
         type(star_info), pointer :: s

         call star_ptr(id, s, ierr)

         if ( s% x_ctrl(2) > 0.d0 ) then
           eps0=s% x_ctrl(2) 
	 endif
         do k = 1, s% nz, +1
           tauHere= s% tau(k)
           if (tauHere .ge. 2/3) exit
  	 end do
         ierr = 0
         ! calcuate the radius at which EUV is absorbed
         sigma=2e-18 !sigma0*(20/13.6)^-3, sigma0 = 6e-18
         Tphoto= 10**(s% log_surface_temperature)
         mH= 1.67e-24  
	 muH=1.00794 !g/mol molar mass!!! why?
         Pphoto= 10**(s% log_surface_pressure)  !dyn/cm2
        
        
         hv= 3.2e-11 ! ergs, is 20 eV
         alphaREC=2.7e-13 !cm^3/s

         g=(s% cgrav(1)*s% mstar)/((Rsun*(10**(s% log_surface_radius)))**2)
         Peuv=(mH*g)/(sigma)  !convert from dyn/cm2 to bars !?? it is not a convertion
         H=(kerg*Tphoto)/(2*mH*g)    !in cm
         Pratio=LOG(Peuv/Pphoto)
         zHeight=-H*Pratio
         Rp=Rsun*(10**(s% log_surface_radius))
         Reuv=Rp+zHeight

         !Energy Limited 
         !Feuv = <<escape_Cpow>> * (s% star_age * 1.0e-9)**(<<escape_betapow>>)*(s% x_ctrl(1))**(-2) 
         !29.12*((s% star_age)/1d9)**(-1.24)*(s% x_ctrl(1))**(-2) !+1d6 the shift
         !if(FEUV > <<escape_sat>>) FEUV = <<escape_sat>>
         !Feuv=29.12*((s% star_age)/1d9)**(-1.24)*(s% x_ctrl(1))**(-2) !+1d6 the shift
         !if(Feuv>8.85d4) Feuv = 8.85d4

         Rhill=(s% x_ctrl(1)*1.495978921d13*((s% mstar)/(3*Msun))**(1/3))
         Ktide=1-(3*Reuv)/(2*Rhill)+ 1/(2*(Rhill/Reuv)**3)
         eLim = -(eps0*pi*Feuv*(Reuv)**3)/(s% cgrav(1)*s% mstar*Ktide)

         !Recombination dominated mass loss (full ionization)  
         Vs= ((kerg*10000)/(mH/2))**0.5   !isothermal speed sound Ct of ionized hydrogen gas at 1d4 K  ! mH/2
         !!Rs=(((s% cgrav(1)*s% mstar)-3*s% cgrav(1)*Msun*Reuv**3)/(s% x_ctrl(1)*1.496e13))/(2*Vs**2)  !for highly inflated planets 
         !Rs=((s% cgrav(1))*(s% mstar))/(2*Vs**2)
         !rhoB0=(muH*g/(2*sigma*kerg*10000))
         !rhoBP=(((Feuv/hv)*sigma*rhoB0)/alphaREC)**0.5   !the denisty of the base is the balance between neutral and ionized H
         !rhoS=rhoBP*muH*exp((s% cgrav(1)*s% mstar)/(Rp*Vs**2)*(Rp/Reuv-1))
         !rrLim =  -(4*pi*rhoS*Vs*Rs**2) 
         
         rrLim= -pi*(((s% cgrav(1)*s% mstar)/(Vs**2))**2)*Vs*mH*((Feuv*s% cgrav(1)*s% mstar)/(hv*alphaREC*Reuv**2*Vs**2))**(0.5)*&
         exp(2- (s% cgrav(1)*s% mstar)/(Vs**2*Reuv))

         
         if (ABS(eLim)<ABS(rrLim)) then  !condition for radiative or energy lim
            s% mstar_dot = eLim
         else 
            s% mstar_dot = rrLim
         endif


	 !s% accreted_material_j= 
	 !s% star_mdot = +PLmassperyear

         !s% accreted_material_j = &
              !s% x_ctrl(5)*sqrt(s% cgrav(1) * s% mstar * s% photosphere_r*Rsun)

         !write(*,*) "debug", s% mstar_dot/Msun*secyer, 10**(s% x_ctrl(6))
         !s% mstar_dot = s% mstar_dot + 10**(s% x_ctrl(6))*Msun/secyer

            
      end subroutine ELR_mdot

      
      
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         call system_clock(time0,clock_rate)
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call system_clock(time1,clock_rate)
         dt = dble(time1 - time0) / clock_rate / 60
         write(*,'(/,a50,f12.2,99i10/)') 'runtime (minutes), retries, backups, steps', &
            dt, s% num_retries, s% num_backups, s% model_number
         ierr = 0
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 4
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         !use const_def
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: mp_m


         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         ! column 1 planet mass in Earth masses
         names(1) = "mp_mearth"
         vals(1) = s% m(1)/m_earth

         ! column 2 planet radius in Earth radii
         names(2) = "rp_rearth"
         vals(2) = s% r(1)/r_earth

         ! column 3: envelope mass fraction
         names(3) = "f_env"
         vals(3) = (s% m(1) - s% M_center) / s% m(1)

         ! column 4: mean planet density (in g/cm^3)
         names(4) = "rho_ave"
         vals(4) =  s% m(1) / (pi4/3d0 * (s% r(1))**3 )

         ierr = 0
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)
      end function extras_finish_step
      
      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg    
         num_ints = i
         
         i = 0
         ! call move_dbl       
         
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info

      end module run_star_extras
      
