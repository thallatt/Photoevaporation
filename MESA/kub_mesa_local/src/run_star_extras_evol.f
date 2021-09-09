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
      real(dp), allocatable, dimension(:), save     ::  Rp_evap,Mp_evap
      real(dp), allocatable, dimension(:,:,:), save ::  Mdot_evap
      real(dp), allocatable, dimension(:,:), save   ::  Mdot_table
      integer      :: NR,NM
      real(dp), dimension(6)     :: tableflux=(/ 3.18831e3, 3.5368e3, &
            8.8419e3, 3.5368e4, 7.9577e4, 3.5368e5 /)
      real(dp)                   :: Tflux=3.53677651e6 ! flux in evaporation table

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
         s% other_adjust_mdot=> <<escape_model>>_mdot !EL_mdot !ELR_mdot !HD_mdot ! HD_highmass_mdot ! owen12_mdot
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
         double precision :: age
         double precision :: mH,muH
         double precision :: Pphoto,g,H,Peuv,Pratio, zHeight, Rp, Rs,rhoS,rhoB0,rhoBP,eLim,rrLim,hv,alphaREC,Vs
         character(len=40) :: age_str
         character(len=40) :: period_
         character(len=100) :: preamble_str="python3 /Applications/kub_mesa_local/retrieve_lxuv.py "
         character(len=2) :: space_str=" "
         double precision :: from_file

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

         print *,'R_XUV/R_P=',Reuv/Rp

         !Energy Limited 
         !if((s% star_age) < <<t_sat>>) then
         !   Feuv=3585.743422231086*(s% x_ctrl(1))**(-2) + 1048.7394601572157*(s% x_ctrl(1))**(-2)
         !else
         !   Feuv= (s% x_ctrl(1))**(-2)*(<<a_x>>*((s% star_age)/<<t_sat>>)**(<<b_x>>) + <<a_euv>>*((s% star_age)/<<t_sat>>)**(<<b_euv>>))
         !endif

         write(age_str,*) (s% star_age)/1d6
         write(period_,*) <<pstar>>

         CALL SYSTEM(preamble_str//age_str//period_)
         OPEN(unit=2,file="temp_lxuv.txt",status='old', action='read')
         READ(2,*) from_file
         Feuv=from_file/4./3.14159/(s% x_ctrl(1))/(s% x_ctrl(1))/2.238d26

         Rhill=(s% x_ctrl(1)*1.495978921d13*((s% mstar)/(3*Msun))**(1/3))
         Ktide=1-(3*Reuv)/(2*Rhill)+ 1/(2*(Rhill/Reuv)**3)
         eLim = -(eps_EUV*pi*Feuv*(Reuv)**3)/(s% cgrav(1)*s% mstar*Ktide)

         !Set stellar mass loss rate with E-limited mass loss (in g/s)
         s% mstar_dot = eLim 
        
         if(<<a_x>>==0.) then
             s% mstar_dot = 0.
         endif

      end subroutine EL_mdot

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EL+Parker wind analytic model
      subroutine EL_Parker_mdot(id, ierr) ! planet parameters: /Applications/mesa-r11554/star/public/star_data.inc
         !use star_def
         !use const_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: eps_EUV, a, REUV, FEUV
         double precision :: Rhill
         double precision :: Ktide
         double precision :: sigma
         double precision :: Tphoto
         double precision :: age
         double precision :: mH,muH
         double precision :: Pphoto,g,H,Peuv,Pratio, zHeight, Rp, Rs,rhoS,rhoB0,rhoBP,eLim,rrLim,hv,alphaREC,Vs
         real(dp) :: Rpl_RE, Mpl_ME, Teq_K, CSound, Pres, RBondi
         character(len=40) :: age_str
         character(len=40) :: period_
         character(len=40) :: radius_str
         character(len=40) :: pres_str
         character(len=40) :: cs_str
         character(len=40) :: f_rr_rb_str
         character(len=100) :: preamble_str="python3 /Applications/kub_mesa_local/retrieve_lxuv.py "
         character(len=100) :: preamble_str2="python3 /Applications/kub_mesa_local/retrieve_mdot_parker.py "
         character(len=2) :: space_str=" "
         double precision :: from_file
         double precision :: from_file2

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
         !if((s% star_age) < <<t_sat>>) then
         !   Feuv=3585.743422231086*(s% x_ctrl(1))**(-2) + 1048.7394601572157*(s% x_ctrl(1))**(-2)
         !else
         !   Feuv= (s% x_ctrl(1))**(-2)*(<<a_x>>*((s% star_age)/<<t_sat>>)**(<<b_x>>) + <<a_euv>>*((s% star_age)/<<t_sat>>)**(<<b_euv>>))
         !endif

         write(age_str,*) (s% star_age)/1d6
         write(period_,*) <<pstar>>

         CALL SYSTEM(preamble_str//age_str//period_)
         OPEN(unit=2,file="temp_lxuv.txt",status='old', action='read')
         READ(2,*) from_file
         Feuv=from_file/4./3.14159/(s% x_ctrl(1))/(s% x_ctrl(1))/2.238d26

         Rhill=(s% x_ctrl(1)*1.495978921d13*((s% mstar)/(3*Msun))**(1/3))
         Ktide=1-(3*Reuv)/(2*Rhill)+ 1/(2*(Rhill/Reuv)**3)
         eLim = (eps_EUV*pi*Feuv*(Reuv)**3)/(s% cgrav(1)*s% mstar*Ktide) ! +: - applied at end

         Rpl_RE = (10.**(s% log_surface_radius)) * 6.96d10 / 6.378d8
         Mpl_ME = s% star_mass * 1.9885d33 / 5.9722d27
         Teq_K =10.**(s% log_surface_temperature)
         
         CSound = (s% photosphere_csound)
         Pres = 10.**(s% log_surface_pressure)
         RBondi = 0.5 * s% cgrav(1) * Mpl_ME * 5.9722d27 / CSound**2.

         write(radius_str,*) Rpl_RE*6.378d8
         write(pres_str,*) Pres
         write(cs_str,*) CSound
         write(f_rr_rb_str,*) (Rpl_RE*6.378d8/RBondi)**-4.*exp(3.-4./(Rpl_RE*6.378d8/RBondi))

         CALL SYSTEM(preamble_str2//radius_str//pres_str//cs_str//f_rr_rb_str)
         OPEN(unit=2,file="temp_mdot_parker.txt",status='old', action='read')
         READ(2,*) from_file2

         IF (from_file2 .ge. eLim) THEN
            s% mstar_dot = -from_file2
            print *,'parker'
         ELSE
            s % mstar_dot = -eLim
         ENDIF

         if(<<a_x>>==0.) then
            s% mstar_dot = 0.
         else
            if((s% star_age)/1d6 > 1000.) then
                s% mstar_dot = 0.
             else
                s% mstar_dot = s % mstar_dot
             endif
            !s% mstar_dot = - lhy
         endif

      end subroutine EL_Parker_mdot

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine HD_mdot(id, ierr)
         !use star_def
         !use const_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: a, Rpl_RE, Mpl_ME, Teq_K, FEUV, lb, bt0, lhy !a==d0 lb==Lambda critical
         real(dp) :: zeta, eta, beta, alp1, alp2, alp3, K, C !coefficients of approximation
         character(len=40) :: age_str
         character(len=40) :: period_
         character(len=100) :: preamble_str="python3 /Applications/kub_mesa_local/retrieve_lxuv.py "
         character(len=2) :: space_str=" "
         double precision :: from_file
         type (star_info), pointer :: s

         ierr = 0

         ! get the star_info pointer using id
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !xctrl inputs 
         a = s% x_ctrl(1)        !orbital separation (AU)
         
         Rpl_RE = (10.**(s% log_surface_radius)) * 6.96d10 / 6.378d8
         Mpl_ME = s% star_mass * 1.9885d33 / 5.9722d27
         Teq_K =10.**(s% log_surface_temperature)
         
         !FEUV flux from star (incident on the planet)
         ! Tu et al. 2015 
         !if((s% star_age) < <<t_sat>>) then
         !   Feuv=3585.743422231086*(s% x_ctrl(1))**(-2) + 1048.7394601572157*(s% x_ctrl(1))**(-2)
         !else
         !   Feuv= (s% x_ctrl(1))**(-2)*(<<a_x>>*((s% star_age)/<<t_sat>>)**(<<b_x>>) + <<a_euv>>*((s% star_age)/<<t_sat>>)**(<<b_euv>>))
         !endif
         !here we implement approximation
         
         write(age_str,*) (s% star_age)/1d6
         write(period_,*) <<pstar>>

         CALL SYSTEM(preamble_str//age_str//period_)
         OPEN(unit=2,file="temp_lxuv.txt",status='old', action='read')
         READ(2,*) from_file
         Feuv=from_file/4./3.14159/(s% x_ctrl(1))/(s% x_ctrl(1))/2.238d26
  
         print *, 'feuv', Feuv
  
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
         if(<<a_x>>==0.) then
            s% mstar_dot = 0.
         else
            s% mstar_dot = - lhy
         endif
         print *, 'mdot =.1e12', lhy/1d12
        !write(*,*) 'REUV (cm), FEUV (erg/cm/s), mstar_dot (g/s)', REUV, FEUV, s% mstar_dot
         
      end subroutine HD_mdot


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! for >40 Mearth planets use different interpolation
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
         double precision :: from_file ! mass loss from interpolation

         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: a, Rpl_RE, Mpl_ME, Teq_K, FEUV, lb, bt0, lhy, eps_EUV, REUV !a==d0 lb==Lambda critical
         real(dp) :: zeta, eta, beta, alp1, alp2, alp3, K, C, thresh !coefficients of approximation ! thresh=mass threshold
         type (star_info), pointer :: s
         character(len=40) :: Rpl_RE_str
         character(len=40) :: Mpl_ME_str
         character(len=40) :: Teq_str
         character(len=40) :: Feuv_str
         character(len=100) :: preamble_str="python3 /Applications/kub_mesa_local/interpol_new.py 1. "
         character(len=2) :: space_str=" "
         character(len=40) :: age_str
         character(len=40) :: period_
         character(len=100) :: preamble_str_lxuv="python3 /Applications/kub_mesa_local/retrieve_lxuv.py "
         double precision :: from_file_lxuv
         
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
         sigma=2e-18
         Tphoto= 10**(s% log_surface_temperature)
         mH= 1.67e-24  
         muH=1.00794 !g/mol molar mass
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
         !if((s% star_age) < <<t_sat>>) then
         !   Feuv=3585.743422231086*(s% x_ctrl(1))**(-2) + 1048.7394601572157*(s% x_ctrl(1))**(-2)
         !else
         !   Feuv= (s% x_ctrl(1))**(-2)*(<<a_x>>*((s% star_age)/<<t_sat>>)**(<<b_x>>) + <<a_euv>>*((s% star_age)/<<t_sat>>)**(<<b_euv>>))
         !endif
         
         write(age_str,*) (s% star_age)/1d6
         write(period_,*) <<pstar>>

         CALL SYSTEM(preamble_str_lxuv//age_str//period_)
         OPEN(unit=2,file="temp_lxuv.txt",status='old', action='read')
         READ(2,*) from_file_lxuv
         Feuv=from_file_lxuv/4./3.14159/(s% x_ctrl(1))/(s% x_ctrl(1))/2.238d26
         
         thresh=0.

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

         write(Rpl_RE_str,*) Rpl_RE
         write(Mpl_ME_str,*) Mpl_ME
         write(Teq_str,*) Teq_K
         write(Feuv_str,*) Feuv

         bt0 = 6.6726d-8*Mpl_ME * 5.9722d27/ (Rpl_RE * 6.378d8* (1.3807d-16* Teq_K/ 1.6726d-24))

         if(<<a_x>>==0.) then
             s% mstar_dot = 0.
         else
             if(Mpl_ME .GT. thresh)then
                 CALL SYSTEM(preamble_str//Feuv_str//space_str//Teq_str//space_str//Rpl_RE_str//space_str//Mpl_ME_str)
                 ! retrieve mdot from temporary file
                 OPEN(unit=2,file="temp_mdot.txt",status='old', action='read')
                 READ(2,*) from_file
                 !Set stellar mass loss rate with interpolated mass loss (in g/s)
                 s% mstar_dot = - from_file
                 print *,'rpl',Rpl_RE_str
                 print *,'mpl',Mpl_ME_str
                 print *,'teq',Teq_str
                 print *,'xuv',Feuv_str
                 print *,'lambda',bt0
                 close(2)
             else
                 if(bt0 .LE. lb)then
                     zeta = -6.861843637445744
                     eta  = -0.009459807206476
         
                     beta = 32.019929264625155
                     alp1 = 0.422232254541188
                     alp2 = -1.748858849270155
                     alp3 = 3.767941293231585
                    
                     K = (zeta + eta * LOG(a))
                     C = beta + alp1 * LOG(FEUV) + alp2 * LOG(a) + alp3 * LOG(Rpl_RE)
                     lhy = exp(C + K * LOG(bt0));!testCFJ
                     s% mstar_dot = - lhy
                 else
                     zeta = -1.297796148718774
                     eta  = 0.884595403184073
                     
                     beta = 16.408393523348366
                     alp1 = 1.0
                     alp2 = -3.286179370395197
                     alp3 = 2.75
                     
                     K = (zeta + eta * LOG(a)) 
                     C = beta + LOG(FEUV) + alp2 * LOG(a) + alp3 * LOG(Rpl_RE)
                     
                     lhy = exp(C + K * LOG(bt0)) !testCF
                     s% mstar_dot = - lhy
                 end if
            end if
         end if
        
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
         character(len=40) :: age_str
         character(len=40) :: period_
         character(len=100) :: preamble_str="python3 /Applications/kub_mesa_local/retrieve_lbol.py "
         character(len=2) :: space_str=" "
         double precision :: from_file

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
               
         s% extra_heat(:) = s% extra_power_source
    
         open(unit = 12, file = <<heating_input>>, status = 'old', action = "read")
         !log10(age [years]) log10(Lbol [Lsun]) log10(Teff [K]) log10(Rstar [Rsun])

         ei_time = 0.0
         sage = s% star_age !in years

         write(age_str,*) (s% star_age)/1d6
         write(period_,*) <<pstar>>

         DO WHILE (sage .gt. 10.0 ** ei_time)
            read(12, *) ei_time, ei_lbol, ei_teff, ei_rstar
         END DO

110      format(4(3X, F10.6))

         CALL SYSTEM(preamble_str//age_str//space_str//period_)
         OPEN(unit=2,file="temp_lbol.txt",status='old', action='read')
         READ(2,*) from_file
         
         Lstar = from_file /4.0/s% x_ctrl(1)/s% x_ctrl(1)/2.238d26/3.1416 !from_file /4.0/s% x_ctrl(1)/s% x_ctrl(1)/2.238d26/3.1416 !(10.0**ei_lbol) * 3.9d33/4.0/s% x_ctrl(1)/s% x_ctrl(1)/2.238d26/3.1416 ! flux = Lstar_bol / (4*pi*(a_cm**2))
         
         close(unit = 12)


         do k = 1, 100 !s% nz            
            s% irradiation_heat(k) = 1.* 3.1416* Lstar * ((exp(s% lnr(k)))**2) / sum(s% dm(1:100)) !1.36d8 !*2.0 !Lstar * 90. !
            !print *,'sigma',sum(s% dm(1:100))/(4.*3.14159*(exp(s% lnr(k)))**2)
         end do

         ! EXTRA HEATING CONSTANT IN SPACE AND TIME
         ! Heat due to Thermal inertia of the core
         ! x_ctrl(4) = cv in erg/K/g from inlist

         k = s% nz
         extraheat = -s% x_ctrl(4) * s% M_center * s% T(s% nz) * s% dlnT_dt(s% nz) / s% dm(s% nz) ! erg/g/sec
           !assuming dlnT_dt is in s^-1

         ! EXTRA HEATING CONSTANT IN SPACE AND TIME
         ! Heat produced by radioactive decay due to 40K, 235U, 238U, 232Th, respectively
         ! x_ctrl(5) = fraction of core mass in "chondrite-like" rocky material

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
         double precision :: Thresh,elx,ele,Fx,Fe

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
        muH=1.00794 !g/mol molar mass
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
         if((s% star_age) < <<t_sat>>) then
            Feuv=3585.743422231086*(s% x_ctrl(1))**(-2) + 1048.7394601572157*(s% x_ctrl(1))**(-2)
            Fx=1048.7394601572157*(s% x_ctrl(1))**(-2)
            Fe=3585.743422231086*(s% x_ctrl(1))**(-2)
         else
            Feuv= (s% x_ctrl(1))**(-2)*(<<a_x>>*((s% star_age)/<<t_sat>>)**(<<b_x>>) + <<a_euv>>*((s% star_age)/<<t_sat>>)**(<<b_euv>>))
            Fx= (s% x_ctrl(1))**(-2) * <<a_x>>*((s% star_age)/<<t_sat>>)**(<<b_x>>)
            Fe=(s% x_ctrl(1))**(-2) * <<a_euv>>*((s% star_age)/<<t_sat>>)**(<<b_euv>>)
            !Feuv= (s% x_ctrl(2))**(-2)*(3522.3247240527944*((s% star_age)/4d6)**(-0.9) + 1027.20279*((s% star_age)/4d6)**(-1.1)) !3685.1546645010812*(s% x_ctrl(2))**(-2)*((s% star_age)/226d6)**(-2.15) + 1059.4107025353496*(s% x_ctrl(2))**(-2)*((s% star_age)/226d6)**(-2.5)
         endif

         Rhill=(s% x_ctrl(1)*1.495978921d13*((s% mstar)/(3*Msun))**(1/3))
         Ktide=1-(3*Reuv)/(2*Rhill)+ 1/(2*(Rhill/Reuv)**3)
         eLim = -(eps0*pi*Feuv*(Reuv)**3)/(s% cgrav(1)*s% mstar*Ktide)

         elx = -(eps0*pi*Fx*(Rp)**3)/(s% cgrav(1)*s% mstar*Ktide)

         ele = -(eps0*pi*Fe*(Reuv)**3)/(s% cgrav(1)*s% mstar*Ktide)

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

         Thresh=1d40*((s% x_ctrl(1))/0.1)**2*(eLx/1d12)**2*(Rp/(6.371d8*10))**(-3)

         if ((Fe * 4.*3.14159*((s% x_ctrl(1)) * 1.5d13)**2 * 3.2d11) > (Thresh)) then
            if ((Fe) > (1d4)) then
                s% mstar_dot = rrLim
            else
                s% mstar_dot = ele
            endif
         else
              s% mstar_dot = elx
         endif

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


      subroutine owen12_mdot(id, ierr)
         use star_def
         implicit none
         integer, intent(in)  :: id
         real(dp)             :: L10R,L10M,L10Mdot_u, L10Mdot_b,L10Mdot
         real(dp)             :: Tstat, PLfoff, Lx_Lbol
         real(dp)             :: mdot,mdot_EUV
         real(dp)             :: HEflux
         real(dp)             :: EUVflux
         real(dp)             :: a_au
         real(dp)             :: Xrayflux
         real(dp) :: Rpl_RE, Mpl_ME
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         integer, intent(out) :: ierr
         character(len=40) :: age_str
         character(len=40) :: period_
         character(len=100) :: preamble_str="python3 /Applications/kub_mesa_local/retrieve_lxuv.py "
         character(len=2) :: space_str=" "
         double precision :: from_file
         type (star_info), pointer     :: s
         
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         ! USED to calculate photoevaporation rates of planets

          call read_evap

          !use seperation and Teq to calculate Lbol
          !Lbol=16.0*pi*((s% x_ctrl(12))*(au))**2.0*boltz_sigma*(s% x_ctrl(19))**(4.0)

         Rpl_RE = (10.**(s% log_surface_radius)) * 6.96d10 / 6.378d8
         Mpl_ME = s% star_mass * 1.9885d33 / 5.9722d27

         a_au = (s% x_ctrl(1))

          !Tstat=7d7
          !Lx_Lbol=0.00025118864315095795
          !PLfoff=1.19

          !compute X-ray flux
          !HEflux= Lx_Lbol*1414710.6052612918*(s% x_ctrl(1))**(-2)

        ! x ray only 
        !   if((s% star_age) < <<t_sat>>) then
        !    HEflux=1048.7394601572157*(s% x_ctrl(1))**(-2)
        ! else
        !    HEflux= (s% x_ctrl(1))**(-2)*(<<a_x>>*((s% star_age)/<<t_sat>>)**(<<b_x>>))
        ! endif
         
        ! XUV
        !if((s% star_age) < <<t_sat>>) then
        !    HEflux=3585.743422231086*(s% x_ctrl(1))**(-2) + 1048.7394601572157*(s% x_ctrl(1))**(-2)
        !    EUVflux=3585.743422231086*(s% x_ctrl(1))**(-2)
        ! else
        !    HEflux= a_au**(-2)*(<<a_x>>*((s% star_age)/<<t_sat>>)**(<<b_x>>) + <<a_euv>>*((s% star_age)/<<t_sat>>)**(<<b_euv>>))
        !     EUVflux=<<a_euv>>*((s% star_age)/<<t_sat>>)**(<<b_euv>>)*(s% x_ctrl(1))**(-2)
        ! endif

         write(age_str,*) (s% star_age)/1d6
         write(period_,*) <<pstar>>

         CALL SYSTEM(preamble_str//age_str//period_)
         OPEN(unit=2,file="temp_lxuv.txt",status='old', action='read')
         READ(2,*) from_file
         HEflux=from_file/4./3.14159/(s% x_ctrl(1))/(s% x_ctrl(1))/2.238d26         

        !if (s% star_age .gt.  Tstat) HEflux=HEflux*((s% star_age)/Tstat)**(-PLfoff)
          
          L10R=log10(10.**(s% log_surface_radius)*Rsun)
          L10M=log10(s% star_mass * 1.9885d33)
          call interpolate_evap(L10R,L10M,L10Mdot)
          
          mdot=10.0**L10Mdot*(HEflux/Tflux)
          
          print*, 'mdot',mdot/1d12
          
          call get_EUV_rate(HEflux,mdot_EUV,(s% star_mass)*1.9885d33,10.**(s% log_surface_radius)*Rsun)
          if (mdot_EUV .gt. mdot) mdot=mdot_EUV

         s% mstar_dot = -mdot

         if(<<a_x>>==0.) then
             s% mstar_dot = 0.
         endif

       end subroutine owen12_mdot

       subroutine get_EUV_rate(HEflux,mdot,Mp,Rp)

        implicit none

        real(dp), intent(in) :: HEflux,Mp,Rp !cgs units
        real(dp), intent(out):: mdot !cgs units
        real(dp), parameter  :: eff=0.15 ! efficiency

        ! calculates a lower bound to the mdot due to pure EUV mass-loss

        mdot=eff*pi*HEflux*Rp**3.0/(4.0*standard_cgrav*Mp) ! low EUV rate Equation 1 of Owen & Jackson (2012)

        return
       end subroutine get_EUV_rate

       subroutine read_evap

       !reads in photoevaporation tables

       character*32          ::      evap_table
       integer               ::      table_length
       integer               ::      i,j,k,n
       real(dp), allocatable, dimension(:) :: evap_M, evap_R, evap_Mdot

       ! Hard code this for now
       table_length=5600
       NR=140
       NM=40

       allocate(Rp_evap(NR))
       allocate(Mp_evap(NM))
       allocate(Mdot_evap(NM,NR,6))
       allocate(Mdot_table(NM,NR))

       allocate(evap_M(table_length))
       allocate(evap_R(table_length))
       allocate(evap_Mdot(table_length))

       evap_M=0.0
       evap_R=0.0
       evap_Mdot=0.0

       do n=1,6
        write(evap_table,'(a,I1.1,a)') "Owen_Rates/Mdot",n,".out"      
        !read in the file
        open(unit=73,file=evap_table)

        do i=1,table_length
           read(73,*) evap_R(i), evap_M(i), evap_Mdot(i)
        end do

        close(73)
        !put in store tables
        k=1
        do i=1,NR
          do j=1,NM
             Rp_evap(i)=evap_R(k)
             Mp_evap(j)=evap_M(k)
             Mdot_evap(j,i,n)=evap_Mdot(k) !reversal to be consistent with MATLAB style indexing
             k=k+1
          end do
        end do
  
       end do

       ! read in new table
       evap_table='Owen_Rates/Mdot_out_Mar15.out'
       open(unit=74, file=evap_table)
       do i=1,table_length
       	read(74,*) evap_R(i), evap_M(i), evap_Mdot(i)
       end do
       close(74)

       !put in store tables
       k=1
       do i=1,NR
          do j=1,NM
             Rp_evap(i)=evap_R(k)
             Mp_evap(j)=evap_M(k)
             Mdot_table(j,i)=evap_Mdot(k) !reversal to be consistent with MATLAB style indexing
             k=k+1
          end do
       end do

     
       end subroutine read_evap

      subroutine interpolate_evap(LR,LM,LMD)
      
      !performs bi-linear interpolation of the evaporation tables
      
      real(dp), intent(in)     :: LR,LM !logarithmic values of Radius and Mass [cgs]
      real(dp), intent(out)    :: LMD   !lograithmic values of the mass-loss rate [cgs]

      integer                  :: NRad,NMass !lower index of evaporation table
      integer                  :: i,j,n
      
      !variables for linear interpolation
      real(dp)                 :: int1, int2, int3

      !variables for nearest neighbour interpolation
      real(dp)                 :: store_dist
      real(dp), dimension(4,2) :: distance !distance to point, mdot at that point

      logical                  :: use_nearest

      use_nearest=.false.

      Nrad=0
      NMass=0
      do i=1,NR
         if (LR .lt. Rp_evap(i)) then
            NRad=i-1
            exit
         endif
      end do

      do j=1,NM
         if (LM .lt. Mp_evap(j)) then
            NMass=j-1
            exit
         endif
      end do

!       if (Nrad==0) then
!          print*, "Planet Radius not on evaporation table!"
!          print*, 10**LR
!          return
!       endif
!       if (Nmass==0) then
!          print*, "Planet Mass not on evaporation table!"
!          print*, 10**LM
!          return
!       endif
      
      !now check not near boundary

      if (Mdot_table(NMass  ,Nrad  ) .lt. 1) use_nearest=.true.
      if (Mdot_table(NMass+1,Nrad  ) .lt. 1) use_nearest=.true.
      if (Mdot_table(NMass  ,Nrad+1) .lt. 1) use_nearest=.true.
      if (Mdot_table(NMass+1,Nrad+1) .lt. 1) use_nearest=.true.

      !Now perform interpolation
      if (use_nearest) then
         !use closetest non-zero value
         !Also throw warning to terminal
         distance(1,1)=((LR-Rp_evap(Nrad  ))**2+(LM-Mp_evap(NMass  ))**2)
         distance(2,1)=((LR-Rp_evap(Nrad+1))**2+(LM-Mp_evap(NMass  ))**2)
         distance(3,1)=((LR-Rp_evap(Nrad  ))**2+(LM-Mp_evap(NMass+1))**2)
         distance(4,1)=((LR-Rp_evap(Nrad+1))**2+(LM-Mp_evap(NMass+1))**2)

         distance(1,2)=Mdot_table(NMass  ,NRad  )
         distance(2,2)=Mdot_table(NMass+1,NRad  )
         distance(3,2)=Mdot_table(NMass  ,NRad+1)
         distance(4,2)=Mdot_table(NMass+1,NRad+1)
         
         store_dist=1d100
         LMD=0.
         do i=1,4
            if (distance(i,1).lt.store_dist) then
               if (distance(i,2) .gt. 1) then
                  LMD=distance(i,2)
                  store_dist=distance(i,1)
               endif
            endif
         end do   
!          if (LMD .gt. 0) then
!             print*, "Warning Using Nearest Neighbour Evaporation Interpolation"
!          endif
      else
         !linear interpolataion

         !Mass first then Radius
         int1=Mdot_table(Nmass,NRad)+(Mdot_table(Nmass+1,NRad)-Mdot_table(Nmass,Nrad)) &
         &  *((LM-Mp_evap(Nmass))/(Mp_evap(Nmass+1)-Mp_evap(Nmass)))
      
         int2=Mdot_table(Nmass,NRad+1)+(Mdot_table(Nmass+1,NRad+1)-Mdot_table(Nmass,Nrad+1)) &
         &  *((LM-Mp_evap(Nmass))/(Mp_evap(Nmass+1)-Mp_evap(Nmass)))
      
         int3=int1+(int2-int1)/(Rp_evap(NRad+1)-Rp_evap(Nrad))*(LR-Rp_evap(Nrad))
         LMD=int3
      endif

      ! need to deallocate saved arrays before next time step
      deallocate(Rp_evap)
      deallocate(Mp_evap)
      deallocate(Mdot_evap)
      deallocate(Mdot_table)

      end subroutine interpolate_evap

      end module run_star_extras
      
