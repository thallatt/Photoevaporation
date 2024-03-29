! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************
 
      module other_adjust_mdot

      ! consult star/other/README for general usage instructions
      ! control name: use_other_adjust_mdot = .true.
      ! procedure pointer: s% other_adjust_mdot => my_routine


      implicit none
      
      contains

      ! your routine will be called after winds and before mass adjustment
   
      subroutine zero_mdot(id, ierr)
         use star_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = 0
         !this one is the NULL escape for tests
      end subroutine null_other_adjust_mdot



      subroutine energy_limited_mdot(id, ierr)
         !use star_def
         !use const_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: eps_EUV, a, REUV, FEUV
         type (star_info), pointer :: s

         ierr = 0

         ! get the star_info pointer using id
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !Give xctrl inputs more interpretable names
         eps_EUV = s% x_ctrl(2)  !mass loss efficiency factor
         a = s% x_ctrl(1)        !orbital separation (AU)
         REUV = s% x_ctrl(3) * s% photosphere_r * Rsun   !REUV (in cm)

         !FEUV flux from star (incident on the planet)
         FEUV = 29.7 * (s% star_age * 1.0e-9)**(-1.23) / a**2  ! erg cm^-2 s^-1

         !Set stellar mass loss rate with E-limited mass loss (in g/s)
         s% mstar_dot = - eps_EUV * pi * FEUV * REUV**3 / (standard_cgrav * s% mstar)

        !write(*,*) 'REUV (cm), FEUV (erg/cm/s), mstar_dot (g/s)', REUV, FEUV, s% mstar_dot

      end subroutine energy_limited_mdot





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

         !FEUV flux from star (incident on the planet) TEMPORARY
         FEUV = 29.7 * (s% star_age * 1.0e-9)**(-1.23) / a**2 /2.0 ! erg cm^-2 s^-1

         !here we implement approximation
         !!! Lborder BRING CONSTANTS TO X_CTRL???
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
         

         !Set stellar mass loss rate with E-limited mass loss (in g/s)
         s% mstar_dot = - lhy

        !write(*,*) 'REUV (cm), FEUV (erg/cm/s), mstar_dot (g/s)', REUV, FEUV, s% mstar_dot

      end subroutine HD_mdot


      end module other_adjust_mdot
      
      
      
      
