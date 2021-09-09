     !!!!!!!!!!!!!!!!!other heating!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine varTEQ_energy(id, ierr)
         use const_def, only: Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer :: k
         real(dp) :: ei_time, ei_lbol, ei_teff, ei_rstar, sage, Lstar

         !use star_lib, only: star_ptr
	     double precision :: extraheat,junk,diag_heat
         integer :: i,n,counter,z,p,numOfAcc,zeta,jafter,jbefore,indexI
	     double precision :: core_epsrad,pressureDep,pressureHere,random_dp,L_tid,Q_tid,a_orbit,Rpl
	     real(dp) :: tauHere,Vesc, KE,massTot,Cd,phi,Mpl,Rhopl,H,g,mH,Tacc,oblty
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
               
         ! IS THIS PURELY HEAT DISSIPATION FROM CORE? 
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
         
         !Lstar = from_file /4.0/s% x_ctrl(1)/s% x_ctrl(1)/2.238d26/3.1416 !from_file /4.0/s% x_ctrl(1)/s% x_ctrl(1)/2.238d26/3.1416 !(10.0**ei_lbol) * 3.9d33/4.0/s% x_ctrl(1)/s% x_ctrl(1)/2.238d26/3.1416 ! flux = Lstar_bol / (4*pi*(a_cm**2))
         Lstar = 4d33 /4.0/s% x_ctrl(1)/s% x_ctrl(1)/2.238d26/3.1416
         close(unit = 12)

         ! zone 1 is at photosphere
         do k = 1, 100 !s% nz            
            s% irradiation_heat(k) = 1.* 3.1416* Lstar * ((exp(s% lnr(k)))**2) / sum(s% dm(1:100)) !1.36d8 !*2.0 !Lstar * 90. !
            !print *,'sigma',sum(s% dm(1:100))/(4.*3.14159*(exp(s% lnr(k)))**2)
         end do

         ! EXTRA HEATING CONSTANT IN SPACE AND TIME
         ! Heat due to Thermal inertia of the core
         ! x_ctrl(4) = cv in erg/K/g from inlist

         ! nz is innermost zone, i.e. zone 0 is at the photosphere
         k = s% nz
         extraheat = -s% x_ctrl(4) * s% M_center * s% T(s% nz) * s% dlnT_dt(s% nz) / s% dm(s% nz) ! erg/g/sec
           !assuming dlnT_dt is in s^-1

         ! EXTRA HEATING CONSTANT IN SPACE AND TIME
         ! Heat produced by radioactive decay due to 40K, 235U, 238U, 232Th, respectively
         ! x_ctrl(5) = fraction of core mass in "chondrite-like" rocky material i.e. ==1 Earth-like core compositions

         k = s% nz
         core_epsrad = 36.7d-8 * exp(-5.543d-10 * s% star_age) ! erg/g/sec   40K
         core_epsrad = core_epsrad + 4.6d-8 * exp(-9.8485d-10 * s% star_age) ! erg/g/sec  235U
         core_epsrad = core_epsrad + 2.3d-8 * exp( -1.5513d-10 * s% star_age)! erg/g/sec  238U
         core_epsrad = core_epsrad + 1.3d-8 * exp( -0.4948d-10 * s% star_age)! erg/g/sec  232Th

         ! obliquity tides (eccentricity=0)
         ! NOTE: Msun might change!
         oblty=pi/4.
         Q_tid=1d6
         Rpl=Rsun*(10**(s% log_surface_radius))
         a_orbit = s% x_ctrl(1) * 1.5d13
         L_tid = 4.5 * (s% cgrav(1) * Msun / a_orbit**3.)**0.5 * sin(oblty)**2. * (Q_tid * (1.+ cos(oblty)**2.))**-1. & 
          * (s% cgrav(1) * Msun**2. * Rpl**-1.) * (Rpl / a_orbit)**6.

         s% extra_heat(k) = (extraheat+ s% x_ctrl(5) * s% M_center * core_epsrad / s% dm(k)) ! erg/g/sec, core heat flux density ! deposited at zone above core (k=nz)
         !s% extra_heat(k) = (L_tid / s% dm(k) + extraheat + s% x_ctrl(5) * s% M_center * core_epsrad / s% dm(k))
         


         return
      end subroutine varTEQ_energy
