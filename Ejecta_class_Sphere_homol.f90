!!!############################### Module containing quantities and routines associated with the specific environment/ejecta ##################################################################################################
!!!######### The code should be written in a way that allows staraightforward replacement of this module with a different one (characterizing a different environment), without affecting the rest of the code ################

!! Copied from '/media/indrek/Samsung_T5/work/Fortran_sandbox/OOP/Test2' as 'Module_OOP_multiple_ejecta_classdef.f90'

!!!###############################################################################################################################################################################################
!!!################################################################### Specific ejecta objects ###################################################################################################
!!!###############################################################################################################################################################################################


!!!##########################################################################################################################################################################
!!!############################################################ Homologous sphere ###########################################################################################
!!!##########################################################################################################################################################################

Module class_Sphere_homol
  use RealKind
  use constants
  implicit none
  private
!   real(kind=rk) :: pi = 3.1415926535897931d0
  
  !!!######### Class for homologously expanding sphere ############
  type, public :: Sphere_homol   !!!! Assumes v = r/t and R_sph -> 0 at t -> 0
  
    !!!########### Primary parameters ############
      real(kind=rk) :: v_ej  !! outer edge velocity; dimensionless (i.e. actually velocity/c)
      real(kind=rk) :: M_ej
      real(kind=rk) :: xi_min  ! r_min/R_sph
      real(kind=rk) :: alpha_r  !!! rho ~ r^{-alpha_r}
      
      real(kind=rk) :: Zav
      real(kind=rk) :: Zz
      
      real(kind=rk) :: mu_mol
      real(kind=rk) :: c_Tambov = 1.0_rk  !!! To increase the heat capacity of plasma to stabilize temperature. Must be careful not to increase cooling times to comparable levels to t_dyn, of for the heat capacity not to rival radiation. 
      
      
      real(kind=rk) :: R_sph_0
    !!!###########################################
    
    !!!############# Grid ########################
      integer :: n_r
      real(kind=rk), dimension(:), allocatable :: xi_cell_vec, xi_cell_bnd_vec    !! vector of normalized radii xi = 1 <=> r = R_sph
      real(kind=rk), dimension(:), allocatable :: kT_vec                          !! kT_vec - dimensionless temperature (i.e. theta). Treated as in the local COMOVING frame (EXCEPT, currently at 08.05.23, in the initialization phase)
      real(kind=rk), dimension(:), allocatable :: M_cell_vec
    !!!###########################################
    
    !!!############ Choose initial radiation #####
      integer :: i_ch_init_rad   !!! 1 - Planck
    !!!###########################################
    
    
! ! !     real(kind=rk) :: w_fiducial    ! Kludge (or maybe not); NB: Can change as the code runs, e.g. by MC_photon_split(ct) if photons are split to improve statistics 
        
    contains
      procedure :: rho => get_density         !!! Treated as in the EXTERNAL frame; NBNB: rho set in Get_ejecta_params is in the local COMOVING frame and used as such in Evolve_photon_diffusion
      procedure :: v_rad => get_velocity
      procedure :: R_sph => get_Rsph
      procedure :: chk_esc
      procedure :: initialize
      procedure :: initialize_radiation
!       procedure :: initialize_radiation => Initialize_Radiation_Passive_fireball_Planck   !!! Could be pointing to some other prescription for photon initialization
      
      procedure :: Emit_therm_brems => Emit_thermal_brems_Sphere_homol
      procedure :: tcool_therm_brems => Cooling_time_thermal_brems
      
      procedure, NOPASS :: Print_spectrum => Printout_timedep_spectrum
      procedure :: Print_Sphere => Printout_Sphere_qties
      procedure :: Printout_Energy_conserv
      
      procedure :: chk_within   !! 10.05.23
      
  end type Sphere_homol
  !!!###############################################################
  
  Contains
  
    !!!########################## Initialize ######################################################
  
    Subroutine initialize(self)
      use auxiliary, ONLY: Setup_grid
      use Init_printout, ONLY: FF
      implicit none
        class(Sphere_homol) :: self
        real(kind=rk) :: d_xi
        real(kind=rk) :: kT_in
        integer :: i
        
        real(kind=rk), dimension(:), allocatable :: M_aux
        real(kind=rk) :: M_ej_norm, M_ej, v_ej, xi_min, alpha_r, Zav, Zz, mu_mol, c_Tambov, R_sph_0, dM !! , kT_vec
        integer :: n_r
        integer :: i_ch_init_rad    
        
        10 format(e14.7)
        11 format(i4)
            
        !!!################
          open(100,ACTION='READ',FILE = 'data/Passive_fireball/Passive_fireball_params.dat')
            do i = 1,5
            read(100,*)
            end do
            read(100,*) M_ej_norm
            read(100,*)
            read(100,*) v_ej
            read(100,*)
            read(100,*) xi_min
            read(100,*)
            read(100,*) alpha_r
            read(100,*)
            read(100,*) Zav
            read(100,*)
            read(100,*) Zz
            read(100,*)
            read(100,*) mu_mol
            read(100,*)
            read(100,*) c_Tambov
            read(100,*)
            read(100,*) n_r
            read(100,*)
            read(100,*) kT_in
            read(100,*)
            read(100,*) R_sph_0
            read(100,*)
            read(100,*) i_ch_init_rad
          close(100)
        !!!##############################################
        
        if (c_Tambov.lt.1.0_rk) c_Tambov = 1.0_rk    !! If c_Tambov<1, set to 1
        
        !!!######## Write back to file with filename extenstion FF ##################
          open(110,ACTION='WRITE',FILE = 'data/Passive_fireball/write_back_params/Passive_fireball_params_' // trim(adjustl(FF)) // '.dat')
            write(110,*) 'Parameters repeated back from Passive_fireball_params.dat for reproducibility'
            write(110,*) 'Parameters for a homologously expanding enjecta, initially for computing radiative transfer in a passively expanding fireball'
            write(110,*) 'The following quantities are passed to the ejecta object (Ejecta_obj)'
            write(110,*) ''
            
            write(110,*) 'M_ej_norm: ejecta mass in solar units'
            write(110,10) M_ej_norm
            write(110,*) 'v_ej: ejecta velocity at the outer edge (c units)'
            write(110,10) v_ej
            write(110,*) 'xi_min: inner (relative) radius of ejecta'
            write(110,10) xi_min
            write(110,*) 'alpha_r: ejecta density profile (rho ~ r^(-alpha_r))'
            write(110,10) alpha_r
            write(110,*) 'Zav'
            write(110,10) Zav
            write(110,*) 'Zz'
            write(110,10) Zz
            write(110,*) 'mu_mol'
            write(110,10) mu_mol
            write(110,*) 'c_Tambov'
            write(110,10) c_Tambov            
            write(110,*) 'n_r: number of radial zones'
            write(110,11) n_r
            write(110,*) 'kT_in: initial fixed temperature (can/should be overridden by a more physical prescription)'
            write(110,10) kT_in
            write(110,*) 'R_sph_0: initial radius of the sphere'
            write(110,10) R_sph_0
            write(110,*) 'i_ch_init_rad: Choose initial radiation field: 1 - Planck, 2 - Wien, 3 - Powerlaw'
            write(110,11) i_ch_init_rad
          close(110)
        !!!##########################################################################

          M_ej = M_ej_norm*Msun
        
          self%M_ej = M_ej
          self%v_ej = v_ej
          self%xi_min = xi_min
          self%alpha_r = alpha_r
          self%Zav = Zav
          self%Zz = Zz
          self%mu_mol = mu_mol
          self%c_Tambov = c_Tambov
          self%n_r = n_r
          
          self%i_ch_init_rad = i_ch_init_rad
          
          allocate(self%kT_vec(n_r))
  !           self%kT_vec(:) = kT_vec(:)  !! can/should be set to something more physical (rather than constant) (from module)
          self%kT_vec(:) = kT_in  !! can/should be set to something more physical (rather than constant)  (from file)
          
          self%R_sph_0 = R_sph_0
        !!!#############################################
          
        allocate(self%xi_cell_vec(n_r), self%xi_cell_bnd_vec(n_r+1))
        allocate(self%M_cell_vec(n_r))

        
        if (1.eq.0) then  !! uniform grid in xi
          call Setup_grid(self%xi_cell_bnd_vec,d_xi,self%xi_min,1.0_rk)  ! xi_max = 1 (by definition)
          self%xi_cell_vec(1:n_r) = (self%xi_cell_bnd_vec(1:n_r) + self%xi_cell_bnd_vec(2:n_r+1))/2.0_rk
          self%M_cell_vec(1:n_r) = M_ej*(self%xi_cell_bnd_vec(2:n_r+1)**(3.0_rk - alpha_r) - self%xi_cell_bnd_vec(1:n_r)**(3.0_rk - alpha_r))/(1.0_rk - xi_min**(3.0_rk - alpha_r))   !!! Need to be extended to handle alpha_r = 3
        else   !!! uniform grid in mass
          allocate(M_aux(n_r+1))
          call Setup_grid(M_aux,dM,0.0_rk,M_ej)
          self%xi_cell_bnd_vec = (M_aux*(1.0_rk - xi_min**(3.0_rk - alpha_r))/M_ej + xi_min**(3.0_rk - alpha_r))**(1.0/(3.0_rk - alpha_r))
          self%xi_cell_vec(1:n_r) = (self%xi_cell_bnd_vec(1:n_r) + self%xi_cell_bnd_vec(2:n_r+1))/2.0_rk
          self%M_cell_vec(1:n_r) = M_aux(2:n_r+1) - M_aux(1:n_r)
          deallocate(M_aux)
        end if
      
    End Subroutine initialize
    
    
    Subroutine initialize_radiation(self)
      use Init_printout, ONLY: FF
      implicit none
      
        class(Sphere_homol) :: self
        real(kind=rk) :: l_rad, alpha_x, x_min, x_max
        integer :: i
        
        10 format(e14.7)
        11 format(i4)
        
        if (self%i_ch_init_rad.eq.1) then
          call Initialize_Radiation_Passive_fireball_Planck(self)
        else if (self%i_ch_init_rad.eq.2) then  
          !!!################ Read data ###################
            open(100,ACTION='READ',FILE = 'data/Passive_fireball/Init_rad_Wien.dat')
              do i = 1,5
              read(100,*)
              end do
              read(100,*) l_rad
            close(100)
          !!!######## Write back to file with filename extenstion FF ##################
            open(110,ACTION='WRITE',FILE = 'data/Passive_fireball/write_back_params/Init_rad_Wien_' // trim(adjustl(FF)) // '.dat')
            write(110,*) 'Parameters repeated back from Passive_fireball_params.dat for reproducibility'
            write(110,*) 'Parameters for initial radiation field (NB: temperature already set in Passive_fireball_params.dat)'
            write(110,*) 'Read by class_Sphere_homol/initialize_radiation'
            write(110,*) ''
            write(110,*) 'l_rad: Compactness, defined as sigma_T*\int u_ph*dr'
            write(110,10) l_rad
            close(110)
          !!!##############################################
          call Initialize_Radiation_Passive_fireball_Wien(self,l_rad)
        else if (self%i_ch_init_rad.eq.3) then  
          !!!################ Read data ###################
            open(100,ACTION='READ',FILE = 'data/Passive_fireball/Init_rad_Powerlaw.dat')
              do i = 1,5
              read(100,*)
              end do
              read(100,*) l_rad
              read(100,*)
              read(100,*) alpha_x
              read(100,*)
              read(100,*) x_min
              read(100,*)
              read(100,*) x_max
            close(100)
          !!!######## Write back to file with filename extenstion FF ##################
            open(120,ACTION='WRITE',FILE = 'data/Passive_fireball/write_back_params/Init_rad_Powerlaw_' // trim(adjustl(FF)) // '.dat')
            write(120,*) 'Parameters repeated back from Passive_fireball_params.dat for reproducibility'
            write(120,*) 'Parameters for initial radiation field (NB: temperature already set in Passive_fireball_params.dat)'
            write(120,*) 'Read by class_Sphere_homol/initialize_radiation'
            write(120,*) ''
            write(120,*) 'l_rad: Compactness, defined as sigma_T*\int u_ph*dr'
            write(120,10) l_rad
            write(120,*) 'alpha_x: Power-law slope, dN/dlnx ~ x^{-alpha_x}'
            write(120,10) alpha_x
            write(120,*) 'x_min: Minimum photon energy (mc^2 units)'
            write(120,10) x_min
            write(120,*) 'x_max: Maximum photon energy (mc^2 units)'
            write(120,10) x_max
            close(120)
          !!!##############################################
          call Initialize_Radiation_Passive_fireball_Powerlaw(self,l_rad,alpha_x,x_min,x_max)
        else
          stop 'Initialize_Radiation: cannot handle i_ch_init_rad'
        end if
    
    End Subroutine initialize_radiation
    
        
    Subroutine Initialize_Radiation_Passive_fireball_Planck(self)
      !!! Initialization specific to passife fireball, using type(Sphere_homol) :: Sphere_object
    
      use Photon_distribution, ONLY: Photons_vec
      use Photon_distribution_escape, ONLY: Photons_escape, i_exist_max_escape
      
! ! !       use Passive_fireball_params, ONLY: R_sph_0
!       use Ejecta_generic, ONLY: Ejecta_obj      !!! Only used here to obtain Ejecta_obj%i_choose_ejecta for a consistency check
!       use Ejecta_objects, ONLY: Sphere_object   !! Use specific object (passive fireball), which requires i_choose_ejecta = 1
      
!       use Ejecta_choice_param, ONLY: i_choose_ejecta  !! Only used for double checking Ejecta_obj%i_choose_ejecta
      
      use MC_photon_params, ONLY: n_MC_ph, w_fiducial
      
      use Ph_thermal, ONLY: f_Planck, Planck_moments
      use auxiliary, ONLY: Setup_grid, Findvalue_1dim_v3
      use Geometry, ONLY: get_v_vec_from_angles_in_r_system

      implicit none
      
        class(Sphere_homol) :: self
        real(kind=rk) :: v_ej, xi_min, ct, d_Vol
        real(kind=rk), dimension(:), allocatable :: kT_vec, xi_cell_bnd_vec
        real(kind=rk) :: theta, n_ph, u_ph, N_ph_tot
        real(kind=rk) :: x_min, x_max, lnx_min, lnx_max, d_lnx, x, lnx_val, mu, phi, mu_loc, phi_loc, r
        real(kind=rk) :: r_cell_min, r_cell_max
        real(kind=rk), dimension(:,:), allocatable :: P_cumul_r_x, lnx
        real(kind=rk) :: w, rnd1, rnd2, rnd3, rnd4, rnd5, rnd6
        real(kind=rk) :: n_MC_dV_real
        real(kind=rk) :: R_sph_0  !!! NB: locally-defined (no longer rea from Module Passive_fireball_params)
        real(kind=rk) :: alpha, N_particles, n_tot
        integer :: n_r, i_r, n_x, i_x, i_MC, i_MC_dV, n_MC_dV, i_dummy
      
!         if ((Ejecta_obj%i_choose_ejecta.ne.1).or.(i_choose_ejecta.ne.1)) then
!           print *,'Initialize_Radiation_Passive_fireball: Ejecta_obj%i_choose_ejecta.ne.1 or i_choose_ejecta.ne.1', Ejecta_obj%i_choose_ejecta, i_choose_ejecta
!           print *,'This routine is specific to passive fireball, which should have i_choose_ejecta = 1' 
!           stop
!         end if
            
        if (allocated(Photons_vec)) deallocate(Photons_vec)
        
        if (allocated(Photons_escape)) deallocate(Photons_escape)
        i_exist_max_escape = 0
        
        n_r = self%n_r
        v_ej = self%v_ej      !!! seems unused in this routine
        xi_min = self%xi_min
        R_sph_0 = self%R_sph_0
        
        allocate(kT_vec(n_r),xi_cell_bnd_vec(n_r+1))
        kT_vec = self%kT_vec
        xi_cell_bnd_vec = self%xi_cell_bnd_vec
        
        ct = R_sph_0/v_ej
        
        n_x = 1000
        allocate(lnx(n_r,n_x))   !! Different grid for each location (allowing for different temperatures)
        allocate(P_cumul_r_x(n_r,n_x))
        P_cumul_r_x(:,:) = 0.0_rk
        
!         Vol_tot = 4.0_rk*pi/3.0_rk*(1.0_rk - xi_min**3.0_rk)*R_sph_0**3.0_rk   !!! Total volume of the ejecta 
        
        N_ph_tot = 0.0_rk  !!! Total number of physical photons
        do i_r = 1,n_r
!           print *,'i_r', i_r
!           d_r = (xi_cell_bnd_vec(i_r+1) - xi_cell_bnd_vec(i_r))*R_sph_0
          d_Vol = 4.0_rk*pi/3.0_rk*(xi_cell_bnd_vec(i_r+1)**3.0_rk - xi_cell_bnd_vec(i_r)**3.0_rk)*R_sph_0**3.0_rk
          
          theta = kT_vec(i_r)
          call Planck_moments(theta,n_ph,u_ph)
          
          
          
          !!!############# Check that the adjusted heat content is not too high ############################################## 
            alpha = 3.0_rk/2.0_rk !! wrong if kT relativistic (then should be 3)
            N_particles = self%M_cell_vec(i_r)/(self%mu_mol*m_p)
            n_tot = N_particles/d_Vol
            
            print *,'n_ph/n_tot, c_Tambov', n_ph/n_tot, self%c_Tambov
            if (self%c_Tambov.ge.max(1.0_rk,(0.1_rk*n_ph/n_tot))) then  !! Require the adjusted matterheat content to be less than 0.1 of the radiation density
              print *,'self%c_Tambov.ge.(0.1_rk*n_ph/n_tot); self%c_Tambov, n_ph/n_tot:', self%c_Tambov, n_ph/n_tot
              stop
            end if
          !!!#################################################################################################################

          
          !!!############## Grids for cumulative distr. in ph. energy #################
            x_min = theta*1.0e-7_rk
            x_max = theta*1.0e2_rk
            lnx_min = log(x_min)
            lnx_max = log(x_max)
            call Setup_grid(lnx(i_r,:),d_lnx,lnx_min,lnx_max)
          !!!#########################################################################
          
          
          N_ph_tot = N_ph_tot + n_ph*d_Vol
          
          P_cumul_r_x(i_r,:) = 0.0_rk   ! redundant
          do i_x = 2,n_x
            x = exp((lnx(i_r,i_x) + lnx(i_r,i_x-1))/2.0_rk)
            P_cumul_r_x(i_r,i_x) = P_cumul_r_x(i_r,i_x-1) + f_Planck(x,theta)*d_lnx
            
!             print *,'x, d_lnx, theta, f_Planck(x,theta)', x, d_lnx, theta, f_Planck(x,theta)
!             print *, 'P_cumul_r_x(i_r,i_x)', P_cumul_r_x(i_r,i_x)
!             stop
            
          end do  !!!### do i_x = 1,n_x
          
          P_cumul_r_x(i_r,:) = P_cumul_r_x(i_r,:)/P_cumul_r_x(i_r,n_x) 
                              
        end do  !!!### do i_r = 1,n_r
        
        print *,'N_ph_tot', N_ph_tot
        
!         print *,'P_cumul_r_x(1,:)', P_cumul_r_x(1,:)
!         stop
        
        !!!############# Draw photons ##################
        !!########### Strategy 3 (in notebook) ####
        
        allocate(Photons_vec(n_MC_ph))
        Photons_vec(:)%exist = .FALSE.
        
        w = N_ph_tot/n_MC_ph
        
! ! !         self%w_fiducial = w  !! Moved w_fiducial out of Ejecta class
        w_fiducial = w
        
        i_MC = 0
        do i_r = 1,n_r
!           print *,'i_r', i_r
          d_Vol = 4.0_rk*pi/3.0_rk*(xi_cell_bnd_vec(i_r+1)**3.0_rk - xi_cell_bnd_vec(i_r)**3.0_rk)*R_sph_0**3.0_rk
          theta = kT_vec(i_r)
          call Planck_moments(theta,n_ph,u_ph)
          
          w = N_ph_tot/n_MC_ph
          
          n_MC_dV_real = n_ph*d_Vol/w
          n_MC_dV = n_MC_dV_real
          w = w*n_MC_dV_real/n_MC_dV
          
!           print *,'n_MC_dV', n_MC_dV
 
          do i_MC_dV = 1,n_MC_dV
            i_MC = i_MC + 1
          
            call RANDOM_NUMBER(rnd1)
            call RANDOM_NUMBER(rnd2)
            call RANDOM_NUMBER(rnd3)
            
            call Findvalue_1dim_v3(rnd1,P_cumul_r_x(i_r,:),lnx(i_r,:),lnx_val,i_dummy)
            
!             print *,'P_cumul_r_x(i_r,:)', P_cumul_r_x(i_r,:)
!             print *,'lnx(i_r,:)', lnx(i_r,:)
!             print *,'lnx_val', lnx_val
!             print *,'rnd1', rnd1
!             stop
            
            mu = 2.0_rk*rnd2 - 1.0_rk
            phi = 2.0_rk*pi*rnd3
            
            call RANDOM_NUMBER(rnd6)
            r_cell_min = xi_cell_bnd_vec(i_r)*R_sph_0
            r_cell_max = xi_cell_bnd_vec(i_r+1)*R_sph_0
            r = (rnd6*r_cell_max**3.0_rk + (1.0_rk - rnd6)*r_cell_min**3.0_rk)**(1.0_rk/3.0_rk)      !! uniform distribution in shperical geometry
! !             r = (rnd6*r_cell_max + (1.0_rk - rnd6)*r_cell_min)  !! WRONG, for testing
            
            Photons_vec(i_MC)%ct = ct
            Photons_vec(i_MC)%r = r
            Photons_vec(i_MC)%x = exp(lnx_val)
            Photons_vec(i_MC)%theta = acos(mu)
            Photons_vec(i_MC)%phi = phi
            Photons_vec(i_MC)%weight = w
            
            Photons_vec(i_MC)%ct_aux = 0.0_rk !!! unused
            
            Photons_vec(i_MC)%id = i_MC   !!! Not unique
            Photons_vec(i_MC)%exist = .TRUE.
            
            !!!#### Draw location ####
              call RANDOM_NUMBER(rnd4)
              call RANDOM_NUMBER(rnd5)
              mu_loc = 2.0_rk*rnd4 - 1.0_rk
              phi_loc = 2.0_rk*pi*rnd5
                      
              Photons_vec(i_MC)%r_loc_vec(1) = r*(1.0_rk - mu_loc**2.0_rk)**0.5_rk*cos(phi_loc)
              Photons_vec(i_MC)%r_loc_vec(2) = r*(1.0_rk - mu_loc**2.0_rk)**0.5_rk*sin(phi_loc)
              Photons_vec(i_MC)%r_loc_vec(3) = r*mu_loc
            !!!#######################
            
            call get_v_vec_from_angles_in_r_system(Photons_vec(i_MC)%theta,Photons_vec(i_MC)%phi,Photons_vec(i_MC)%r_loc_vec(:),Photons_vec(i_MC)%v_vec(:))   !! Calculates Photons_vec(i_MC)%v_vec
                        
            Photons_vec(i_MC)%ct_sourced = Photons_vec(i_MC)%ct
            Photons_vec(i_MC)%r_sourced = Photons_vec(i_MC)%r
            
            Photons_vec(i_MC)%count_scatter = 0
          
          end do
        
        end do
        
        if (i_MC.ne.n_MC_ph) then
          print *,'Initialize_Radiation_Passive_fireball: i_MC.ne.n_MC_ph', i_MC, n_MC_ph
!           stop
        end if
        
        
        print *,'Initial number of photons (from MC)', sum(Photons_vec%weight, mask = Photons_vec%exist)
        print *,'Initial energy in MC photons', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist)
        print *,'Initial energy per photon (from MC)', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist)/sum(Photons_vec%weight, mask = Photons_vec%exist)
        
        deallocate(kT_vec,xi_cell_bnd_vec)
        deallocate(P_cumul_r_x)
        deallocate(lnx)

    End Subroutine Initialize_Radiation_Passive_fireball_Planck
    
    
    
    Subroutine Initialize_Radiation_Passive_fireball_Wien(self,l_rad)
      !!! Initialization specific to passife fireball, using type(Sphere_homol) :: Sphere_object
    
      use Photon_distribution, ONLY: Photons_vec
      use Photon_distribution_escape, ONLY: Photons_escape, i_exist_max_escape
      
! ! !       use Passive_fireball_params, ONLY: R_sph_0
!       use Ejecta_generic, ONLY: Ejecta_obj      !!! Only used here to obtain Ejecta_obj%i_choose_ejecta for a consistency check
!       use Ejecta_objects, ONLY: Sphere_object   !! Use specific object (passive fireball), which requires i_choose_ejecta = 1
      
!       use Ejecta_choice_param, ONLY: i_choose_ejecta  !! Only used for double checking Ejecta_obj%i_choose_ejecta
      
      use MC_photon_params, ONLY: n_MC_ph, w_fiducial
      
      use Ph_thermal, ONLY: f_Wien
      use auxiliary, ONLY: Setup_grid, Findvalue_1dim_v3
      use Geometry, ONLY: get_v_vec_from_angles_in_r_system

      implicit none
      
        class(Sphere_homol) :: self
        real(kind=rk), intent(in) :: l_rad
        real(kind=rk) :: v_ej, M_ej, xi_min, alpha_r, Zav, ct, d_Vol
        real(kind=rk), dimension(:), allocatable :: kT_vec, xi_cell_bnd_vec
        real(kind=rk) :: theta, n_ph, u_ph, N_ph_tot, n_tot
        real(kind=rk) :: x_min, x_max, lnx_min, lnx_max, d_lnx, x, lnx_val, mu, phi, mu_loc, phi_loc, r
        real(kind=rk) :: r_cell_min, r_cell_max
        real(kind=rk), dimension(:,:), allocatable :: P_cumul_r_x, lnx
        real(kind=rk) :: w, rnd1, rnd2, rnd3, rnd4, rnd5, rnd6
        real(kind=rk) :: n_MC_dV_real
        real(kind=rk) :: R_sph_0  !!! NB: locally-defined (no longer rea from Module Passive_fireball_params)
        real(kind=rk) :: tau_T, rho
        
        real(kind=rk) :: E_matter_therm, E_matter_therm_eff, N_particles, N_particles_eff, alpha    !! runtime printout only (E_matter_therm ignores relativistic effects)
        integer :: n_r, i_r, n_x, i_x, i_MC, i_MC_dV, n_MC_dV, i_dummy
      
!         if ((Ejecta_obj%i_choose_ejecta.ne.1).or.(i_choose_ejecta.ne.1)) then
!           print *,'Initialize_Radiation_Passive_fireball: Ejecta_obj%i_choose_ejecta.ne.1 or i_choose_ejecta.ne.1', Ejecta_obj%i_choose_ejecta, i_choose_ejecta
!           print *,'This routine is specific to passive fireball, which should have i_choose_ejecta = 1' 
!           stop
!         end if
            
        if (allocated(Photons_vec)) deallocate(Photons_vec)
        
        if (allocated(Photons_escape)) deallocate(Photons_escape)
        i_exist_max_escape = 0
        
        n_r = self%n_r
        v_ej = self%v_ej
        M_ej = self%M_ej
        xi_min = self%xi_min
        R_sph_0 = self%R_sph_0
        
        alpha_r = self%alpha_r
        Zav = self%Zav
        
        allocate(kT_vec(n_r),xi_cell_bnd_vec(n_r+1))
        kT_vec = self%kT_vec
        xi_cell_bnd_vec = self%xi_cell_bnd_vec
        
        ct = R_sph_0/v_ej
        
        n_x = 1000
        allocate(lnx(n_r,n_x))   !! Different grid for each location (allowing for different temperatures)
        allocate(P_cumul_r_x(n_r,n_x))
        P_cumul_r_x(:,:) = 0.0_rk
        
!         Vol_tot = 4.0_rk*pi/3.0_rk*(1.0_rk - xi_min**3.0_rk)*R_sph_0**3.0_rk   !!! Total volume of the ejecta

        !!! Initial opacity through the ejecta
        if (alpha_r.eq.1.0_rk) then
          tau_T = 2.0_rk*M_ej/(4.0_rk*pi*R_sph_0**2.0_rk)*Zav*sigma_T/m_p*(-log(xi_min))/(1.0_rk - xi_min**2.0_rk)
        else
          tau_T = (3.0_rk - alpha_r)/(1.0_rk - alpha_r)*M_ej/(4.0_rk*pi*R_sph_0**2.0_rk)*Zav*sigma_T/m_p*(1.0_rk - xi_min**(1.0_rk - alpha_r))/(1.0_rk - xi_min**(3.0_rk - alpha_r))
        end if
        
        print *,'alpha_r, M_ej, R_sph_0, Zav, sigma_T, m_p, xi_min', alpha_r, M_ej, R_sph_0, Zav, sigma_T, m_p, xi_min
        print *,'tau_T', tau_T
        
        E_matter_therm = 0.0_rk
        E_matter_therm_eff = 0.0_rk
        N_ph_tot = 0.0_rk  !!! Total number of physical photons
        do i_r = 1,n_r
!           print *,'i_r', i_r
!           d_r = (xi_cell_bnd_vec(i_r+1) - xi_cell_bnd_vec(i_r))*R_sph_0
          d_Vol = 4.0_rk*pi/3.0_rk*(xi_cell_bnd_vec(i_r+1)**3.0_rk - xi_cell_bnd_vec(i_r)**3.0_rk)*R_sph_0**3.0_rk
          
          
          theta = kT_vec(i_r)
          
          r = R_sph_0*self%xi_cell_vec(i_r)
          rho = self%rho(ct,r)
          u_ph = l_rad*Zav*rho/(m_p*tau_T)
          n_ph = u_ph/(3.0_rk*theta)
                    
          print *,'n_ph', n_ph
          
!           call Planck_moments(theta,n_ph,u_ph)

          
            alpha = 3.0_rk/2.0_rk !! wrong if kT relativistic (then should be 3)
            N_particles = self%M_cell_vec(i_r)/(self%mu_mol*m_p)
            N_particles_eff = N_particles*self%c_Tambov    !!! Artificially increase the number of particles/heat capacity
            
            n_tot = N_particles/d_Vol
            
            !!!############# Check that the adjusted heat content is not too high ############################################## 
              print *,'n_ph/n_tot, c_Tambov', n_ph/n_tot, self%c_Tambov
              if (self%c_Tambov.ge.max(1.0_rk,(0.1_rk*n_ph/n_tot))) then  !! Require the adjusted matterheat content to be less than 0.1 of the radiation density
                print *,'Initialize_Radiation_Passive_fireball_Wien: self%c_Tambov.ge.(0.1_rk*n_ph/n_tot); self%c_Tambov, n_ph/n_tot:', self%c_Tambov, n_ph/n_tot
                print *,'i_r', i_r
                stop
              end if
            !!!#################################################################################################################
            
          !!!### Just for printout: thermal energy in particles (NB: E_matter_therm ignores relativistic effects) ####  
            E_matter_therm = E_matter_therm + N_particles*theta*alpha
            E_matter_therm_eff = E_matter_therm_eff + N_particles_eff*theta*alpha
          !!!#########################################################################################################
          
          
          
          !!!############## Grids for cumulative distr. in ph. energy #################
            x_min = theta*1.0e-7_rk
            x_max = theta*1.0e2_rk
            lnx_min = log(x_min)
            lnx_max = log(x_max)
            call Setup_grid(lnx(i_r,:),d_lnx,lnx_min,lnx_max)
          !!!#########################################################################
          
          
          N_ph_tot = N_ph_tot + n_ph*d_Vol
          
          P_cumul_r_x(i_r,:) = 0.0_rk   ! redundant
          do i_x = 2,n_x
            x = exp((lnx(i_r,i_x) + lnx(i_r,i_x-1))/2.0_rk)
!             P_cumul_r_x(i_r,i_x) = P_cumul_r_x(i_r,i_x-1) + f_Planck(x,theta)*d_lnx
            P_cumul_r_x(i_r,i_x) = P_cumul_r_x(i_r,i_x-1) + f_Wien(x,theta,u_ph)*d_lnx
            
!             print *,'x, d_lnx, theta, f_Planck(x,theta)', x, d_lnx, theta, f_Planck(x,theta)
!             print *, 'P_cumul_r_x(i_r,i_x)', P_cumul_r_x(i_r,i_x)
!             stop
            
          end do  !!!### do i_x = 1,n_x
          
          P_cumul_r_x(i_r,:) = P_cumul_r_x(i_r,:)/P_cumul_r_x(i_r,n_x) 
                              
        end do  !!!### do i_r = 1,n_r
        
        print *,'N_ph_tot', N_ph_tot
        
!         print *,'P_cumul_r_x(1,:)', P_cumul_r_x(1,:)
!         stop
        
        !!!############# Draw photons ##################
        !!########### Strategy 3 (in notebook) ####
        
        allocate(Photons_vec(n_MC_ph))
        Photons_vec(:)%exist = .FALSE.
        
        w = N_ph_tot/n_MC_ph
        
! !         self%w_fiducial = w  !! Moved w_fiducial out of Ejecta class
        w_fiducial = w
        
        i_MC = 0
        do i_r = 1,n_r
!           print *,'i_r', i_r
          d_Vol = 4.0_rk*pi/3.0_rk*(xi_cell_bnd_vec(i_r+1)**3.0_rk - xi_cell_bnd_vec(i_r)**3.0_rk)*R_sph_0**3.0_rk
          theta = kT_vec(i_r)
!           call Planck_moments(theta,n_ph,u_ph)

          !!!### Construction ###
          r = R_sph_0*self%xi_cell_vec(i_r)
          rho = self%rho(ct,r)
          u_ph = l_rad*Zav*rho/(m_p*tau_T)   !!! 
          n_ph = u_ph/(3.0_rk*theta)
          !!!####################
          
          w = N_ph_tot/n_MC_ph
          
          n_MC_dV_real = n_ph*d_Vol/w
          n_MC_dV = n_MC_dV_real
          w = w*n_MC_dV_real/n_MC_dV
          
!           print *,'n_MC_dV', n_MC_dV
 
          do i_MC_dV = 1,n_MC_dV
            i_MC = i_MC + 1
          
            call RANDOM_NUMBER(rnd1)
            call RANDOM_NUMBER(rnd2)
            call RANDOM_NUMBER(rnd3)
            
            call Findvalue_1dim_v3(rnd1,P_cumul_r_x(i_r,:),lnx(i_r,:),lnx_val,i_dummy)
            
!             print *,'P_cumul_r_x(i_r,:)', P_cumul_r_x(i_r,:)
!             print *,'lnx(i_r,:)', lnx(i_r,:)
!             print *,'lnx_val', lnx_val
!             print *,'rnd1', rnd1
!             stop
            
            mu = 2.0_rk*rnd2 - 1.0_rk
            phi = 2.0_rk*pi*rnd3
            
            call RANDOM_NUMBER(rnd6)
            r_cell_min = xi_cell_bnd_vec(i_r)*R_sph_0
            r_cell_max = xi_cell_bnd_vec(i_r+1)*R_sph_0
            r = (rnd6*r_cell_max**3.0_rk + (1.0_rk - rnd6)*r_cell_min**3.0_rk)**(1.0_rk/3.0_rk)      !! uniform distribution in shperical geometry
! !             r = (rnd6*r_cell_max + (1.0_rk - rnd6)*r_cell_min)  !! WRONG, for testing
            
            Photons_vec(i_MC)%ct = ct
            Photons_vec(i_MC)%r = r
            Photons_vec(i_MC)%x = exp(lnx_val)
            Photons_vec(i_MC)%theta = acos(mu)
            Photons_vec(i_MC)%phi = phi
            Photons_vec(i_MC)%weight = w
            
            Photons_vec(i_MC)%ct_aux = 0.0_rk !!! unused
            
            Photons_vec(i_MC)%id = i_MC   !!! Not unique
            Photons_vec(i_MC)%exist = .TRUE.
            
            !!!#### Draw location ####
              call RANDOM_NUMBER(rnd4)
              call RANDOM_NUMBER(rnd5)
              mu_loc = 2.0_rk*rnd4 - 1.0_rk
              phi_loc = 2.0_rk*pi*rnd5
                      
              Photons_vec(i_MC)%r_loc_vec(1) = r*(1.0_rk - mu_loc**2.0_rk)**0.5_rk*cos(phi_loc)
              Photons_vec(i_MC)%r_loc_vec(2) = r*(1.0_rk - mu_loc**2.0_rk)**0.5_rk*sin(phi_loc)
              Photons_vec(i_MC)%r_loc_vec(3) = r*mu_loc
            !!!#######################
            
            call get_v_vec_from_angles_in_r_system(Photons_vec(i_MC)%theta,Photons_vec(i_MC)%phi,Photons_vec(i_MC)%r_loc_vec(:),Photons_vec(i_MC)%v_vec(:))   !! Calculates Photons_vec(i_MC)%v_vec
                        
            Photons_vec(i_MC)%ct_sourced = Photons_vec(i_MC)%ct
            Photons_vec(i_MC)%r_sourced = Photons_vec(i_MC)%r
            
            Photons_vec(i_MC)%count_scatter = 0
          
          end do
        
        end do
        
        if (i_MC.ne.n_MC_ph) then
          print *,'Initialize_Radiation_Passive_fireball: i_MC.ne.n_MC_ph', i_MC, n_MC_ph
!           stop
        end if
        
        print *,'Initial number of photons (from MC)', sum(Photons_vec%weight, mask = Photons_vec%exist)
        print *,'Initial energy in MC photons', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist)
        print *,'Initial energy per photon (from MC)', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist)/sum(Photons_vec%weight, mask = Photons_vec%exist)
        
        print *,'Initial thermal energy in matter; E_matter_therm, E_matter_therm_eff:', E_matter_therm, E_matter_therm_eff
        
        deallocate(kT_vec,xi_cell_bnd_vec)
        deallocate(P_cumul_r_x)
        deallocate(lnx)

    End Subroutine Initialize_Radiation_Passive_fireball_Wien    
    
    
    
    
    
    Subroutine Initialize_Radiation_Passive_fireball_Powerlaw(self,l_rad,alpha_x,x_min,x_max)
    
      use MC_photon_params, ONLY: n_MC_ph, w_fiducial
    
      use Photon_distribution, ONLY: Photons_vec
      use Photon_distribution_escape, ONLY: Photons_escape, i_exist_max_escape
      
      use auxiliary, ONLY: Setup_grid, Findvalue_1dim_v3
      use Geometry, ONLY: get_v_vec_from_angles_in_r_system
    
      implicit none
    
      class(Sphere_homol) :: self
      real(kind=rk), intent(in) :: l_rad, alpha_x, x_min, x_max
      real(kind=rk) :: u_ph_ov_n_ph
      real(kind=rk) :: v_ej, M_ej, xi_min, alpha_r, Zav, ct, d_Vol
      real(kind=rk), dimension(:), allocatable :: xi_cell_bnd_vec
      real(kind=rk) :: n_ph, u_ph, N_ph_tot
      real(kind=rk) :: lnx_min, lnx_max, d_lnx, x, lnx_val, mu, phi, mu_loc, phi_loc, r
      real(kind=rk) :: r_cell_min, r_cell_max
      real(kind=rk), dimension(:), allocatable :: P_cumul_r_x, lnx   !!! Was a matrix in Plack and Wien versions with varying temperature within the ejecta
      real(kind=rk) :: w, rnd1, rnd2, rnd3, rnd4, rnd5, rnd6
      real(kind=rk) :: n_MC_dV_real
      real(kind=rk) :: R_sph_0  !!! NB: locally-defined (no longer rea from Module Passive_fireball_params)
      real(kind=rk) :: tau_T, rho
      real(kind=rk) :: alpha, N_particles, n_tot
      integer :: n_r, i_r, n_x, i_x, i_MC, i_MC_dV, n_MC_dV, i_dummy
      
        if (allocated(Photons_vec)) deallocate(Photons_vec)
        
        if (allocated(Photons_escape)) deallocate(Photons_escape)
        i_exist_max_escape = 0
        
        if (alpha_x.eq.0.0_rk) then
          u_ph_ov_n_ph = (x_max - x_min)/log(x_max/x_min)
        else if (alpha_x.eq.1.0_rk) then
          u_ph_ov_n_ph = -log(x_max/x_min)/(x_max**(-1.0_rk) - x_min**(-1.0_rk))
        else
          u_ph_ov_n_ph = -alpha_x/(1.0_rk - alpha_x)*(x_max**(1.0_rk - alpha_x) - x_min**(1.0_rk - alpha_x))/(x_max**(-alpha_x) - x_min**(-alpha_x))    !!! Average energy per photon
        end if
        
        print *,'u_ph_ov_n_ph', u_ph_ov_n_ph
!         stop
        
        
        n_r = self%n_r
        v_ej = self%v_ej
        M_ej = self%M_ej
        xi_min = self%xi_min
        R_sph_0 = self%R_sph_0
        
        alpha_r = self%alpha_r
        Zav = self%Zav
        
        allocate(xi_cell_bnd_vec(n_r+1))
        xi_cell_bnd_vec = self%xi_cell_bnd_vec
        
        ct = R_sph_0/v_ej
        
        n_x = 1000
        allocate(lnx(n_x))   !! Different grid for each location (allowing for different temperatures)
        allocate(P_cumul_r_x(n_x))
        P_cumul_r_x(:) = 0.0_rk
        
        !!!############## Grids for cumulative distr. in ph. energy #################
          lnx_min = log(x_min)
          lnx_max = log(x_max)
          call Setup_grid(lnx(:),d_lnx,lnx_min,lnx_max)
        !!!##########################################################################
        
        !!! Initial opacity through the ejecta
        if (alpha_r.eq.1.0_rk) then
          tau_T = 2.0_rk*M_ej/(4.0_rk*pi*R_sph_0**2.0_rk)*Zav*sigma_T/m_p*(-log(xi_min))/(1.0_rk - xi_min**2.0_rk)
        else
          tau_T = (3.0_rk - alpha_r)/(1.0_rk - alpha_r)*M_ej/(4.0_rk*pi*R_sph_0**2.0_rk)*Zav*sigma_T/m_p*(1.0_rk - xi_min**(1.0_rk - alpha_r))/(1.0_rk - xi_min**(3.0_rk - alpha_r))
        end if
        
        print *,'alpha_r, M_ej, R_sph_0, Zav, sigma_T, m_p, xi_min', alpha_r, M_ej, R_sph_0, Zav, sigma_T, m_p, xi_min
        print *,'tau_T', tau_T
!         stop
        
        N_ph_tot = 0.0_rk  !!! Total number of physical photons
        do i_r = 1,n_r
          d_Vol = 4.0_rk*pi/3.0_rk*(xi_cell_bnd_vec(i_r+1)**3.0_rk - xi_cell_bnd_vec(i_r)**3.0_rk)*R_sph_0**3.0_rk
          
          r = R_sph_0*self%xi_cell_vec(i_r)
          rho = self%rho(ct,r)
          u_ph = l_rad*Zav*rho/(m_p*tau_T)
          n_ph = u_ph/u_ph_ov_n_ph
          
          !!!############# Check that the adjusted heat content is not too high ############################################## 
            alpha = 3.0_rk/2.0_rk !! wrong if kT relativistic (then should be 3)
            N_particles = self%M_cell_vec(i_r)/(self%mu_mol*m_p)
            n_tot = N_particles/d_Vol
            
            print *,'n_ph/n_tot, c_Tambov', n_ph/n_tot, self%c_Tambov
            if (self%c_Tambov.ge.max(1.0_rk,(0.1_rk*n_ph/n_tot))) then  !! Require the adjusted matterheat content to be less than 0.1 of the radiation density
              print *,'self%c_Tambov.ge.(0.1_rk*n_ph/n_tot); self%c_Tambov, n_ph/n_tot:', self%c_Tambov, n_ph/n_tot
              stop
            end if
          !!!#################################################################################################################
          
          
          N_ph_tot = N_ph_tot + n_ph*d_Vol
        end do !!! do i_r = 1,n_r
        
        print *,'N_ph_tot', N_ph_tot
        
        P_cumul_r_x(:) = 0.0_rk   ! redundant
        do i_x = 2,n_x
          x = exp((lnx(i_x) + lnx(i_x-1))/2.0_rk)
          P_cumul_r_x(i_x) = P_cumul_r_x(i_x-1) + x**(-alpha_x)*d_lnx
        end do  !!!### do i_x = 1,n_x
        P_cumul_r_x(:) = P_cumul_r_x(:)/P_cumul_r_x(n_x) 
        
        
      
        !!!############# Draw photons ##################
        
        allocate(Photons_vec(n_MC_ph))
        Photons_vec(:)%exist = .FALSE.
        
        w = N_ph_tot/n_MC_ph
        
!         self%w_fiducial = w  !!! Moved w_fiducial out of Ejecta class
        w_fiducial = w
        
        i_MC = 0
        do i_r = 1,n_r
!           print *,'i_r', i_r
          d_Vol = 4.0_rk*pi/3.0_rk*(xi_cell_bnd_vec(i_r+1)**3.0_rk - xi_cell_bnd_vec(i_r)**3.0_rk)*R_sph_0**3.0_rk

          !!!### Construction ###
          r = R_sph_0*self%xi_cell_vec(i_r)
          rho = self%rho(ct,r)
          u_ph = l_rad*Zav*rho/(m_p*tau_T)   !!!
          n_ph = u_ph/u_ph_ov_n_ph
          
          w = N_ph_tot/n_MC_ph
          
          n_MC_dV_real = n_ph*d_Vol/w
          n_MC_dV = n_MC_dV_real
          w = w*n_MC_dV_real/n_MC_dV
          
          
          do i_MC_dV = 1,n_MC_dV
            i_MC = i_MC + 1
          
            call RANDOM_NUMBER(rnd1)
            call RANDOM_NUMBER(rnd2)
            call RANDOM_NUMBER(rnd3)
            
            call Findvalue_1dim_v3(rnd1,P_cumul_r_x(:),lnx(:),lnx_val,i_dummy)
            
!             print *,'P_cumul_r_x(i_r,:)', P_cumul_r_x(i_r,:)
!             print *,'lnx(i_r,:)', lnx(i_r,:)
!             print *,'lnx_val', lnx_val
!             print *,'rnd1', rnd1
!             stop
            
            mu = 2.0_rk*rnd2 - 1.0_rk
            phi = 2.0_rk*pi*rnd3
            
            call RANDOM_NUMBER(rnd6)
            r_cell_min = xi_cell_bnd_vec(i_r)*R_sph_0
            r_cell_max = xi_cell_bnd_vec(i_r+1)*R_sph_0
            r = (rnd6*r_cell_max**3.0_rk + (1.0_rk - rnd6)*r_cell_min**3.0_rk)**(1.0_rk/3.0_rk)      !! uniform distribution in shperical geometry
! !             r = (rnd6*r_cell_max + (1.0_rk - rnd6)*r_cell_min)  !! WRONG, for testing
            
            Photons_vec(i_MC)%ct = ct
            Photons_vec(i_MC)%r = r
            Photons_vec(i_MC)%x = exp(lnx_val)
            Photons_vec(i_MC)%theta = acos(mu)
            Photons_vec(i_MC)%phi = phi
            Photons_vec(i_MC)%weight = w
            
            Photons_vec(i_MC)%ct_aux = 0.0_rk !!! unused
            
            Photons_vec(i_MC)%id = i_MC   !!! Not unique
            Photons_vec(i_MC)%exist = .TRUE.
            
            !!!#### Draw location ####
              call RANDOM_NUMBER(rnd4)
              call RANDOM_NUMBER(rnd5)
              mu_loc = 2.0_rk*rnd4 - 1.0_rk
              phi_loc = 2.0_rk*pi*rnd5
                      
              Photons_vec(i_MC)%r_loc_vec(1) = r*(1.0_rk - mu_loc**2.0_rk)**0.5_rk*cos(phi_loc)
              Photons_vec(i_MC)%r_loc_vec(2) = r*(1.0_rk - mu_loc**2.0_rk)**0.5_rk*sin(phi_loc)
              Photons_vec(i_MC)%r_loc_vec(3) = r*mu_loc
            !!!#######################
            
            call get_v_vec_from_angles_in_r_system(Photons_vec(i_MC)%theta,Photons_vec(i_MC)%phi,Photons_vec(i_MC)%r_loc_vec(:),Photons_vec(i_MC)%v_vec(:))   !! Calculates Photons_vec(i_MC)%v_vec
                        
            Photons_vec(i_MC)%ct_sourced = Photons_vec(i_MC)%ct
            Photons_vec(i_MC)%r_sourced = Photons_vec(i_MC)%r
            
            Photons_vec(i_MC)%count_scatter = 0
          
          end do  !!! do i_MC_dV = 1,n_MC_dV
          
        end do  !!!### do i_r = 1,n_r
        
        
        if (i_MC.ne.n_MC_ph) then
          print *,'Initialize_Radiation_Passive_fireball: i_MC.ne.n_MC_ph', i_MC, n_MC_ph
!           stop
        end if
        
        print *,'Initial number of photons (from MC)', sum(Photons_vec%weight, mask = Photons_vec%exist)
        print *,'Initial energy in MC photons', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist)
        print *,'Initial energy per photon (from MC)', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist)/sum(Photons_vec%weight, mask = Photons_vec%exist)
        
        
        deallocate(xi_cell_bnd_vec)
        deallocate(P_cumul_r_x)
        deallocate(lnx)
    
    End Subroutine Initialize_Radiation_Passive_fireball_Powerlaw
!     
!     
!     
!      
    
    
    
    
        
    !!!######################## End: Initialize ###################################################
    
    
    
    Function get_density(self,t,r) result(rho)
      !!! t is actually ct
      implicit none
      class(Sphere_homol), intent(in) :: self
      real(kind=rk), intent(in) :: t, r
      real(kind=rk) :: v_ej, M_ej, r_min, alpha_r, R_sph, xi
      real(kind=rk) :: rho
      
        v_ej = self%v_ej
        M_ej = self%M_ej
        alpha_r = self%alpha_r
        
        R_sph = v_ej*t
        xi = r/R_sph
        r_min = R_sph*self%xi_min
        
        if ((r.le.R_sph).and.(r.ge.r_min)) then
          rho = (3.0d0 - alpha_r)*M_ej/(4.0d0*pi*R_sph**3.0d0*(1.0d0 - (self%xi_min)**(3.0d0 - alpha_r)))*xi**(-alpha_r)
        else
          rho = 0.0d0
        end if
    End Function get_density
    
    
    Function get_velocity(self,t,r) result(v_rad)
      !!! t is actually ct
      implicit none
      class(Sphere_homol), intent(in) :: self
      real(kind=rk), intent(in) :: t, r
      real(kind=rk) :: v_ej, r_min, R_sph
      real(kind=rk) :: v_rad
      
        v_ej = self%v_ej      
        R_sph = v_ej*t
        r_min = R_sph*self%xi_min
      
        if ((r.le.R_sph).and.(r.ge.r_min)) then
          v_rad = r/t
        else
          v_rad = 0.0_rk
        end if

    End Function get_velocity
    
    
    Function get_Rsph(self,t) result(R_sph)
      !!! t is actually ct
      implicit none
      class(Sphere_homol) :: self
      real(kind=rk), intent(in) :: t
      real(kind=rk) :: R_sph
        R_sph = self%v_ej*t
    End Function get_Rsph
    
    
    
    Logical Function chk_esc(self,t,r)
      !!! t is actually ct
      implicit none
      class(Sphere_homol) :: self
      real(kind=rk), intent(in) :: t, r
      real(kind=rk) :: R_sph
        R_sph = self%v_ej*t
        
        if (r.gt.R_sph) then
          chk_esc = .TRUE.
        else
          chk_esc = .FALSE.
        end if
    End Function chk_esc
    
    
    Logical Function chk_within(self,t,r)
      !!! t is actually ct
      implicit none
      class(Sphere_homol) :: self
      real(kind=rk), intent(in) :: t, r
      real(kind=rk) :: R_sph, r_min
        R_sph = self%v_ej*t
        r_min = R_sph*self%xi_min
                
        if ((r.le.R_sph).and.(r.ge.r_min)) then
          chk_within = .TRUE.
        else
          chk_within = .FALSE.
        end if
    End Function chk_within
    
    
    
    !!!##############################################################################################################################################################
    !!!##################################################### Emission from homologous sphere ########################################################################
    !!!##############################################################################################################################################################
    
      !!!#### STILL APPROXIMATE (Gaunt factor =1)
      Subroutine Emit_thermal_brems_Sphere_homol(self,ct,d_ct)
      
        !!! ct before updating by d_ct
      
        use data_type_def, ONLY: Photon_param
        use Cumulative_distributions, ONLY: P_cumul_Therm_brems, y_min_brems, y_max_brems, n_y_brems, gff, kT_floor_brems, &
                                            x_pr_min
                                            !!#### NBNBNB: x_pr_min is from Compton grids, used here to avoid emitting photons below what can be handled by Compton
            
        use Photon_distribution, ONLY: Photons_vec
        use MC_photon_params, ONLY: w_fiducial
            
        use auxiliary, ONLY: Setup_grid, Findvalue_1dim_v3
        use Geometry, ONLY: get_v_vec_from_angles_in_r_system
      
        implicit none
        
        class(Sphere_homol) :: self
        real(kind=rk), intent(in) :: d_ct, ct
        real(kind=rk) :: v_ej, M_ej, R_sph, CF0, Cf, d_Vol, rho_av_LAB, rho_av_RF, n_e, Z2_nI, kT
        real(kind=rk) :: N_dot_range, N_ph_emit, n_MC_dV_real, w, lny, y, x_pr, mu_pr, x, mu, phi, r_cell_min, r_cell_max, r, mu_loc, phi_loc
        real(kind=rk) :: y_min_loc
        real(kind=rk) :: rnd1, rnd2, rnd3, rnd4, rnd5, rnd6, rnd7, rnd8
        integer :: n_r, i_r, i_MC, n_MC_dV, i_MC_dV, n_MC_curr, i_dummy
                
        real(kind=rk), dimension(:), allocatable :: P_cumul_loc, lny_loc
        real(kind=rk) :: Z_bulk, Beta_bulk, Gamma_bulk, ct_emit, d_En_matter, d_En_matter_MC
        real(kind=rk) :: d_kT, kT_new, N_particles, N_particles_eff, alpha
        real(kind=rk) :: E_dot
        real(kind=rk) :: dummy, E_matter_therm
        integer :: n_y_loc, i_y_min
        
        type(Photon_param), dimension(:), allocatable :: Photons_emit, Photons_tmp
        
!         real(kind=rk), save :: TEST_d_En_matter = 0.0_rk, TEST_d_En_matter_MC = 0.0_rk  !!! Local TEST only, can be removed
        
          n_r = self%n_r
          v_ej = self%v_ej
          M_ej = self%M_ej
          R_sph = self%R_sph(ct)
          
        
! ! !           gff = 1.0_rk  !!! NBNB: TEMPORARY (gaunt factor) -> put in Variables.f90/Cumulative_distributions
          Cf0 = (8.0_rk/(3.0_rk*pi))**0.5_rk*c0*sigma_T*afs*gff
          
! ! !           w_fiducial = self%w_fiducial  ! Moved w_fiducial out of Ejecta_class

! ! !           n_y_brems = size(P_cumul_Therm_brems%lny(:))  !! Size of the entire cumulative distribution: already set in Variables
                
          E_matter_therm = 0.0_rk       
                
          i_MC = 0
          do i_r = 1,n_r
            d_Vol = 4.0_rk*pi/3.0_rk*(self%xi_cell_bnd_vec(i_r+1)**3.0_rk - self%xi_cell_bnd_vec(i_r)**3.0_rk)*R_sph**3.0_rk
            
            Beta_bulk = v_ej*self%xi_cell_vec(i_r)                              !!! Introduced 08.05.23
            Gamma_bulk = 1.0_rk/(1.0_rk - Beta_bulk**2.0_rk)**0.5_rk            !!! Introduced 08.05.23
            
            rho_av_LAB = self%M_cell_vec(i_r)/d_Vol   !!! Lab frame
            rho_av_RF = rho_av_LAB/Gamma_bulk
            n_e = rho_av_RF*self%Zav/m_p             !!! Update 08.05.23: introduced rho_av --> rho_av_RF, i.e. n_e now treated as in the local rest frame
            Z2_nI = rho_av_RF*self%Zz/m_p      !!    !!! Reminder: implement Zz properly!! Update 08.05.23: introduced rho_av --> rho_av_RF, i.e. Z2_nI now treated as in the local rest frame
            
            kT = self%kT_vec(i_r)        !! kT treated as in the local comoving frame
            
            if (kT.gt.kT_floor_brems) then
            
              y_min_loc = x_pr_min/kT
            
              Cf = Cf0*kT**(-0.5_rk)*n_e*Z2_nI

            
              !!!############ Construction #############
                call Findvalue_1dim_v3(log(y_min_loc),P_cumul_Therm_brems%lny(:),P_cumul_Therm_brems%lny(:),dummy,i_y_min)  !! P_cumul_Therm_brems%P_y_unnorm(i) = (E_1_ymin - E_1) = I(y) - I(y_min) in the notes
                n_y_loc = n_y_brems - i_y_min
                i_y_min = i_y_min + 1  !! 1st grid point above y_min_loc 
                allocate(P_cumul_loc(n_y_loc),lny_loc(n_y_loc))
                P_cumul_loc(1:n_y_loc) = P_cumul_Therm_brems%P_y_unnorm(i_y_min:n_y_brems) - P_cumul_Therm_brems%P_y_unnorm(i_y_min)  !! = I(i_y) - I(i_y_min)
                lny_loc(1:n_y_loc) = P_cumul_Therm_brems%lny(i_y_min:n_y_brems)
                
                N_dot_range = Cf*P_cumul_loc(n_y_loc)    !!! In the local comoving frame
                N_ph_emit = N_dot_range*d_ct*d_Vol/c0    !!! In the lab frame (using d_ct*d_Vol = invariant)
                
                P_cumul_loc(1:n_y_loc) = P_cumul_loc(1:n_y_loc)/P_cumul_loc(n_y_loc)  ! normalize
              !!!#######################################
              
!               E_dot = -Cf*kT                    !!! OUT 10.05.23
              E_dot = -Cf*kT *(exp(-exp(lny_loc(1))) - exp(-exp(lny_loc(n_y_loc))))   !!! IN 10.05.23: exponents account for power em. in a finite range (above x_pr_min), for energy balance with the emitted photons.
                                                                                     !!! Physically, the photons below x_pr_min are assumed to be reabsorbed immediately
                                                                                     !!! NB: the double exponent is correct
              
              if (E_dot.gt.0.0_rk) then
                print *,'Emit_thermal_brems_Sphere_homol ERROR: E_dot > 0', E_dot
                stop
              end if
              
  !             N_dot = Cf*P_cumul_Therm_brems%Ndot_norm
  !             N_ph_emit = N_dot*d_ct*d_Vol/c0
              
              n_MC_dV_real = N_ph_emit/w_fiducial
              
              if (n_MC_dV_real.lt.1.0_rk) then   !!! 10.05.23: Added option for n_MC_dV_real.lt.1.0
                call random_number(rnd8)
                if (rnd8.le.n_MC_dV_real) then
                  n_MC_dV = 1
                else
                  n_MC_dV = 0
                end if
                w = w_fiducial
              else
                n_MC_dV = n_MC_dV_real
                w = w_fiducial*n_MC_dV_real/dble(n_MC_dV)
              end if
              
  !             print *,'######### Number of emitted brems MC photons: n_MC_dV', n_MC_dV
  !             print *,'kT', kT
              
              if (.not.(allocated(Photons_emit)).and.(n_MC_dV.gt.0)) then   !!! 10.05.23: added '.and.(n_MC_dV.gt.0)'
                allocate(Photons_emit(n_MC_dV*n_r))
                Photons_emit%exist = .FALSE.            !!! added 10.05.23
              end if 
              
! !               print *, allocated(Photons_emit), i_MC, n_MC_dV
              
              !!#### Extend array if necesary #############
!               if ((i_MC + n_MC_dV).gt.(size(Photons_emit))) then
              if (((i_MC + n_MC_dV).gt.(size(Photons_emit))).and.(allocated(Photons_emit))) then
                allocate(Photons_tmp(i_MC))
                Photons_tmp(1:i_MC) = Photons_emit(1:i_MC)
                deallocate(Photons_emit)
                allocate(Photons_emit(i_MC + n_MC_dV))
                Photons_emit%exist = .FALSE.
                Photons_emit(1:i_MC) = Photons_tmp(1:i_MC)
                deallocate(Photons_tmp)
              end if
              
              d_En_matter_MC = 0.0_rk
              
              do i_MC_dV = 1,n_MC_dV
                i_MC = i_MC + 1
                        
                call RANDOM_NUMBER(rnd1)
                call Findvalue_1dim_v3(rnd1,P_cumul_loc(:),lny_loc(:),lny,i_dummy)
  !               call Findvalue_1dim_v3(rnd1,P_cumul_Therm_brems%P_y(:),P_cumul_Therm_brems%lny(:),lny,i_dummy)
                
                y = exp(lny)
                x_pr = y*kT    !!! Fluid frame 
                
                call RANDOM_NUMBER(rnd2)
                call RANDOM_NUMBER(rnd3)
                mu_pr = 2.0_rk*rnd2 - 1.0_rk
                phi = 2.0_rk*pi*rnd3
                
                !!#### TEST #####
                  if (x_pr.lt.1.0e-10_rk) then
                    print *,'rnd1', rnd1
                    print *,'kT', kT
                    print *,'P_cumul_loc', P_cumul_loc
                    print *,'lny_loc', lny_loc
                    stop
                  end if
                !!!################

                !!!################ Draw time and location #######################
                  call RANDOM_NUMBER(rnd7)
! ! !                   rnd7 = 0.0_rk  !!! NBNBNBNB: WRONG:TEST
                  ct_emit = ct + d_ct*rnd7
                  R_sph = self%R_sph(ct_emit)
                
                  call RANDOM_NUMBER(rnd6)
                  r_cell_min = self%xi_cell_bnd_vec(i_r)*R_sph
                  r_cell_max = self%xi_cell_bnd_vec(i_r+1)*R_sph
                  r = (rnd6*r_cell_max**3.0_rk + (1.0_rk - rnd6)*r_cell_min**3.0_rk)**(1.0_rk/3.0_rk)      !! uniform distribution in spherical geometry
                !!!###############################################################  
                
                !!!########## Local ejecta velocity #########################################
                  Beta_bulk = self%v_rad(ct_emit,r)    !!! NBNBNB: PUT IN
                  Gamma_bulk = 1.0_rk/(1.0_rk - Beta_bulk**2.0_rk)**0.5_rk
                  Z_bulk = Beta_bulk*Gamma_bulk 
                !!!##########################################################################
                
                !!!####### Lorentz transformation #############
                  x = x_pr*(Gamma_bulk + Z_bulk*mu_pr)
                  mu = (mu_pr + Beta_bulk)/(1.0_rk + Beta_bulk*mu_pr)
                !!#############################################
                
                Photons_emit(i_MC)%ct = ct_emit !! Should distribute in time, but carefil with changing R_sph
                Photons_emit(i_MC)%r = r
                Photons_emit(i_MC)%x = x
                Photons_emit(i_MC)%theta = acos(mu)
                Photons_emit(i_MC)%phi = phi
                Photons_emit(i_MC)%weight = w
                
                Photons_emit(i_MC)%ct_aux = 0.0_rk !!! unused
                
                Photons_emit(i_MC)%id = i_MC   !!! Not unique
                
                if (x.gt.x_pr_min) then   !! Discard if below Compton grid threshold
                  Photons_emit(i_MC)%exist = .TRUE.
                else
                  Photons_emit(i_MC)%exist = .FALSE.  ! Kludge
                end if
                
                !!!#### Draw location ####
                  call RANDOM_NUMBER(rnd4)
                  call RANDOM_NUMBER(rnd5)
                  mu_loc = 2.0_rk*rnd4 - 1.0_rk
                  phi_loc = 2.0_rk*pi*rnd5
                          
                  Photons_emit(i_MC)%r_loc_vec(1) = r*(1.0_rk - mu_loc**2.0_rk)**0.5_rk*cos(phi_loc)
                  Photons_emit(i_MC)%r_loc_vec(2) = r*(1.0_rk - mu_loc**2.0_rk)**0.5_rk*sin(phi_loc)
                  Photons_emit(i_MC)%r_loc_vec(3) = r*mu_loc
                !!!#######################
                
                call get_v_vec_from_angles_in_r_system(Photons_emit(i_MC)%theta,Photons_emit(i_MC)%phi,Photons_emit(i_MC)%r_loc_vec(:),Photons_emit(i_MC)%v_vec(:))   !! Calculates Photons_vec(i_MC)%v_vec
                          
                Photons_emit(i_MC)%ct_sourced = Photons_emit(i_MC)%ct
                Photons_emit(i_MC)%r_sourced = Photons_emit(i_MC)%r
              
                Photons_emit(i_MC)%count_scatter = 0
                
                d_En_matter_MC = d_En_matter_MC - w*x_pr  !!! Minus is right, means that matter enegy 'gain' is negative
                
  !               print *,'x_pr', x_pr
                
              end do  !!! do i_MC_dV = 1,n_MC_dV
              deallocate(P_cumul_loc,lny_loc)
              
              !!!####### Update temperature in the shell NO LONGER based on energy loss due to actually emitted photons ###########################################
              !!!####### NBNB: based on emitted energy in the fluid frame, energy conserv NOT perfect if photons not sampled perfectly isotropically ####
              
                alpha = 3.0_rk/2.0_rk !! wrong if kT relativistic (then should be 3)
                N_particles = self%M_cell_vec(i_r)/(self%mu_mol*m_p)
                N_particles_eff = N_particles*self%c_Tambov    !!! Artificially increase the number of particles/heat capacity
            
                d_En_matter = E_dot*d_ct*d_Vol/c0
            
                d_kT = d_En_matter/(alpha*N_particles_eff)
                kT_new = kT + d_kT
                self%kT_vec(i_r) = kT_new
                
!                 print *,'d_En_matter_MC/d_En_matter', d_En_matter_MC/d_En_matter
                
                E_matter_therm = E_matter_therm + alpha*N_particles*kT_new  !!! NBNBNB: Just for runtime printouts; NB: E_matter_therm ignores relativistic effects

                if (kT_new.lt.0.0_rk) then
                  print *,'Emit_thermal_brems_Sphere_homol ERROR: kT_new < 0'
                  print *,'i_r', i_r
                  print *,'kT_new, kT, d_kT ', kT_new, kT, d_kT 
                  print *,'d_En_matter, N_particles, N_particles_eff', d_En_matter, N_particles, N_particles_eff
                  print *,'self%M_cell_vec(i_r)', self%M_cell_vec(i_r)
                  print *,'n_MC_dV, N_ph_emit', n_MC_dV, N_ph_emit
                  stop
                end if
              !!!########################################################################################
              
!               TEST_d_En_matter_MC = TEST_d_En_matter_MC + d_En_matter_MC
!               TEST_d_En_matter = TEST_d_En_matter + d_En_matter
!               print *,'d_En_matter, d_En_matter_MC', d_En_matter, d_En_matter_MC
!               print *,'TEST_d_En_matter, TEST_d_En_matter_MC', TEST_d_En_matter, TEST_d_En_matter_MC

            end if !! if (kT.gt.kT_floor_brems) then
        
          end do  !! do i_r = 1,n_r
                    
!           print *,'Thermal energy in matter', E_matter_therm        
                    
          !!!# Construction, KLUDGE ###
          
          if (i_MC.gt.0) then   !! Were any photons emitted? (added 10.05.23)
            n_MC_curr = size(Photons_vec)
            allocate(Photons_tmp(n_MC_curr))
            Photons_tmp(1:n_MC_curr) = Photons_vec(1:n_MC_curr)
            deallocate(Photons_vec)
            allocate(Photons_vec(n_MC_curr + i_MC))
            Photons_vec%exist = .FALSE.
            Photons_vec(1:n_MC_curr) = Photons_tmp(1:n_MC_curr)
            deallocate(Photons_tmp)
            Photons_vec(n_MC_curr+1:n_MC_curr+i_MC) = Photons_emit(1:i_MC)
            deallocate(Photons_emit)
          end if
          
      End Subroutine Emit_thermal_brems_Sphere_homol  
      
      
      
      Function Cooling_time_thermal_brems(self,ct) result(ct_cool_min)
      
        !!! ct_cool_min in the lab frame
      
        use Cumulative_distributions, ONLY: gff
        
        implicit none
        
          class(Sphere_homol) :: self
          real(kind=rk), intent(in) :: ct
          real(kind=rk) :: ct_cool_min
          real(kind=rk), dimension(:), allocatable :: t_cool
          real(kind=rk) :: R_sph, Cf0, d_Vol, rho_av_LAB, rho_av_RF, n_e, Z2_nI, kT, E_dot, alpha, n_tot
          real(kind=rk) :: Beta_bulk, Gamma_bulk, v_ej
          integer :: i_r, n_r
        
          n_r = self%n_r
          R_sph = self%R_sph(ct)
          v_ej = self%v_ej
          
!           gff = 1.0_rk  !!! NBNB: TEMPORARY (gaunt factor) -> put in Variables.f90/Cumulative_distributions
          Cf0 = (8.0_rk/(3.0_rk*pi))**0.5*c0*sigma_T*afs*gff
          
          allocate(t_cool(n_r))
          
          do i_r = 1,n_r
            d_Vol = 4.0_rk*pi/3.0_rk*(self%xi_cell_bnd_vec(i_r+1)**3.0_rk - self%xi_cell_bnd_vec(i_r)**3.0_rk)*R_sph**3.0_rk
            
            Beta_bulk = v_ej*self%xi_cell_vec(i_r)                              !!! Introduced 08.05.23
            Gamma_bulk = 1.0_rk/(1.0_rk - Beta_bulk**2.0_rk)**0.5_rk            !!! Introduced 08.05.23
                    
            rho_av_LAB = self%M_cell_vec(i_r)/d_Vol
            rho_av_RF = rho_av_LAB/Gamma_bulk
            n_e = rho_av_RF*self%Zav/m_p
            Z2_nI = rho_av_RF*self%Zz/m_p   !!    !!! REminder: implement Zz properly!!
            
            kT = self%kT_vec(i_r)
            E_dot = Cf0*kT**(0.5)*n_e*Z2_nI
            
            alpha = 3.0_rk/2.0_rk
            n_tot = rho_av_RF/(self%mu_mol*m_p)
            
            t_cool(i_r) = alpha*n_tot*kT/E_dot*Gamma_bulk  !!! Introduced Gamma_bulk on 08.05.23: t_cool in the lab frame
          end do
          
          ct_cool_min = c0*minval(t_cool)
          
          deallocate(t_cool)
          
          if (1.eq.0) then
            print *,'Cooling_time_thermal_brems: kT(1)', self%kT_vec(1)
            print *,'Cooling_time_thermal_brems: min(kT), max(kT)', minval(self%kT_vec(:)), maxval(self%kT_vec(:))
            print *,'Cooling_time_thermal_brems: minloc(kT), maxloc(kT)', minloc(self%kT_vec(:)), maxloc(self%kT_vec(:))
          end if
      
      End Function Cooling_time_thermal_brems
    
    !!!##############################################################################################################################################################
    !!!#################################### End: emission from a homologous sphere ##################################################################################
    !!!##############################################################################################################################################################
    
    
    !!!##############################################################################################################################################################
    !!!##################################################### Printouts from homologous sphere #######################################################################
    !!!##############################################################################################################################################################
    
    
      Subroutine Printout_timedep_spectrum(ct,sw_reset_Photons_escape)
            
        use Photon_distribution_escape, ONLY: Photons_escape, i_exist_max_escape
        use Photon_distribution, ONLY: Photons_vec
        use auxiliary, ONLY: Setup_grid, bin
        use Init_printout, ONLY: FF
        
        implicit none
        
!         class(Sphere_homol), intent(in) :: self
        real(kind=rk), intent(in) :: ct
        real(kind=rk), dimension(:), allocatable :: x_vec, lnx_vec, Nx, Nx_MC
        real(kind=rk), dimension(:), allocatable, save :: Nx_esc, Nx_esc_MC
        real(kind=rk) :: lnx_min, lnx_max, d_lnx, x
        integer :: n_x_bin, i_ph, n_ph, n_ph_esc, i
        integer, save :: i_save = 0
        integer, save :: i_exist_max_escape_chk = 0
        logical, intent(in) :: sw_reset_Photons_escape
        
        3 format(e14.7,' ',e14.7,' ',e14.7)
        4 format(e14.7,' ',e14.7,' ',e14.7,' ',e14.7)
        
          n_x_bin = 40 ! 20
          lnx_min = log(1.0e-6_rk)
          !lnx_max = log(1.0e-2_rk) ! Should be increased
          lnx_max = log(1.0e2_rk)
          allocate(x_vec(n_x_bin),lnx_vec(n_x_bin))        
          call Setup_grid(lnx_vec,d_lnx,lnx_min,lnx_max)
          x_vec(:) = exp(lnx_vec)

          if (.NOT.(allocated(Nx))) then  !!! 01.10.20
            allocate(Nx(n_x_bin), Nx_MC(n_x_bin))
            Nx(:) = 0.0_rk
            Nx_MC(:) = 0.0_rk
          end if
          
          if (allocated(Photons_vec)) then          
            n_ph = size(Photons_vec)
            do i_ph = 1,n_ph
              x = Photons_vec(i_ph)%x
              
              if (Photons_vec(i_ph)%exist) then
                call Bin(x,x_vec,Photons_vec(i_ph)%weight,Nx)   
                call Bin(x,x_vec,1.0_rk,Nx_MC) 
              end if
              
            end do
          end if
          
          if (i_save.eq.0) then
            open(100,FILE='results/Spectra/Spec_photons_' // trim(adjustl(FF)) // '.dat')
  !           open(100,FILE='results/Spectra/Spec_photons.dat')
          else
            open(100,POSITION='APPEND',FILE='results/Spectra/Spec_photons_' // trim(adjustl(FF)) // '.dat')
  !           open(100,POSITION='APPEND',FILE='results/Spectra/Spec_photons.dat')
          end if

          do i = 1,n_x_bin
            x = exp(log(x_vec(i)) + d_lnx/2.0_rk)
            write(100,4) ct, x, Nx(i)/d_lnx, Nx_MC(i)
          end do
          close(100)
          deallocate(Nx,Nx_MC)
          
          

          !!!############ Escaping photons #################
            
            if (.NOT.(allocated(Nx_esc))) then  !!! 01.10.20
              allocate(Nx_esc(n_x_bin),Nx_esc_MC(n_x_bin))
              Nx_esc(:) = 0.0_rk
              Nx_esc_MC(:) = 0.0_rk
            end if
            
            if (sw_reset_Photons_escape) then  !! Is Photons_escape reset at each timestep? (01.10.20)
              Continue  !! DO NOT RESET Nx_esc, Nx_esc_MC, i_exist_max_escape_chk
            else
              i_exist_max_escape_chk = 0
              Nx_esc(:) = 0.0_rk
              Nx_esc_MC(:) = 0.0_rk
            end if
                      
            if (allocated(Photons_escape)) then          
              n_ph_esc = size(Photons_escape)
              do i = 1,n_ph_esc
                x = Photons_escape(i)%x
                if (Photons_escape(i)%exist) then
                
                  if (i.le.i_exist_max_escape) then   !!! SHOULD (!) be redundant
                    call Bin(x,x_vec,Photons_escape(i)%weight,Nx_esc)  
                    call Bin(x,x_vec,1.0_rk,Nx_esc_MC)
                  else
                    open(600,FILE='results/Misc/Photons_escape_ERROR_' // trim(adjustl(FF)) // '.dat')
                    write(600,*) 'Printout_timedep_spectrum'
                    write(600,*) i, i_exist_max_escape
                    close(600)
                  end if
                  i_exist_max_escape_chk = i_exist_max_escape_chk + 1
    
                end if
              end do
            end if
            
            print *,'i_exist_max_escape, i_exist_max_escape_chk', i_exist_max_escape, i_exist_max_escape_chk
            
            if (i_save.eq.0) then
              open(110,FILE='results/Spectra/Spec_esc_cumul_' // trim(adjustl(FF)) // '.dat')
  !             open(110,FILE='results/Spectra/Spec_esc_cumul.dat')
            else
              open(110,POSITION='APPEND',FILE='results/Spectra/Spec_esc_cumul_' // trim(adjustl(FF)) // '.dat')
  !             open(110,POSITION='APPEND',FILE='results/Spectra/Spec_esc_cumul.dat')
            end if

            do i = 1,n_x_bin
              x = exp(log(x_vec(i)) + d_lnx/2.0_rk)
              write(110,4) ct, x, Nx_esc(i)/d_lnx, Nx_esc_MC(i)
            end do
            close(110)
            
            
            if (.NOT.(sw_reset_Photons_escape)) then
              deallocate(Nx_esc,Nx_esc_MC)
            end if
          
          !!!################################################
          
          i_save = i_save + 1
          
          print *,'Printout_timedep_spectrum DONE'
          
      End Subroutine Printout_timedep_spectrum
      
      
      
      
      Subroutine Printout_Sphere_qties(self,ct)
      
!         use Ejecta_generic, ONLY: Ejecta_obj
        use Init_printout, ONLY: FF
        
        use Photon_distribution, ONLY: Photons_vec
        use auxiliary, ONLY: bin
        
        implicit none
        
        class(Sphere_homol), intent(in) :: self
        real(kind=rk), intent(in) :: ct
        real(kind=rk), dimension(:), allocatable :: E_rad_r, N_rad_r, N_MC_r
        real(kind=rk) :: R_sph, x, r, xi, d_Vol, u_rad, n_rad
        real(kind=rk) :: alpha, N_particles, n_tot, n_e, Z2_nI
        integer :: n_r, i_r, n_ph, i_ph
        integer, save :: i_save = 0
        
        4 format(e14.7,' ',e14.7,' ',e14.7,' ',e14.7)
        14 format(e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7)
        
        
!           if (Ejecta_obj%i_choose_ejecta.eq.1) then   !!! Printouts specific to a passive homologous ejecta
        
            R_sph = self%R_sph(ct)
            n_r = size(self%kT_vec)
          
            allocate(E_rad_r(n_r),N_rad_r(n_r),N_MC_r(n_r))
            E_rad_r(:) = 0.0_rk
            N_rad_r(:) = 0.0_rk
            N_MC_r(:) = 0.0_rk
            
            !!!############################ Get radiation energy density #######################################
              if (allocated(Photons_vec)) then          
                n_ph = size(Photons_vec)
                do i_ph = 1,n_ph
                  x = Photons_vec(i_ph)%x
                  r = Photons_vec(i_ph)%r
                  xi = r/R_sph
!                   call Findvalue_1dim_v4(xi,self%xi_cell_bnd_vec(:),self%xi_cell_bnd_vec(:),dummy,i_cell)
                  
                  if (Photons_vec(i_ph)%exist) then
                    call Bin(xi,self%xi_cell_bnd_vec(:),Photons_vec(i_ph)%weight,N_rad_r(:))
                    call Bin(xi,self%xi_cell_bnd_vec(:),Photons_vec(i_ph)%weight*Photons_vec(i_ph)%x,E_rad_r(:))
                    call Bin(xi,self%xi_cell_bnd_vec(:),1.0_rk,N_MC_r(:))
                  end if
                end do
              end if
            !!!#################################################################################################
        


            if (i_save.eq.0) then
              open(200,FILE='results/Ejecta/Ejecta_qties_' // trim(adjustl(FF)) // '.dat')
    !           open(200,FILE='results/Ejecta/Ejecta_qties.dat')
            else
              open(200,POSITION='APPEND',FILE='results/Ejecta/Ejecta_qties_' // trim(adjustl(FF)) // '.dat')
    !           open(200,POSITION='APPEND',FILE='results/Ejecta/Ejecta_qties.dat')
            end if
            
            do i_r = 1,n_r
              d_Vol = 4.0_rk*pi/3.0_rk*(self%xi_cell_bnd_vec(i_r+1)**3.0_rk - self%xi_cell_bnd_vec(i_r)**3.0_rk)*R_sph**3.0_rk
              
              u_rad = E_rad_r(i_r)/d_Vol
              n_rad = N_rad_r(i_r)/d_Vol
              
              alpha = 3.0_rk/2.0_rk !! wrong if kT relativistic (then should be 3)
              N_particles = self%M_cell_vec(i_r)/(self%mu_mol*m_p)
              n_tot = N_particles/d_Vol
              
              n_e = self%M_cell_vec(i_r)/d_Vol*self%Zav/m_p   !!! LAB FRAME
              
              Z2_nI = self%M_cell_vec(i_r)/d_Vol*self%Zz/m_p  !!! rho*Zz/m_p

!               write(200,14) ct, self%xi_cell_vec(i_r), self%kT_vec(i_r), self%M_cell_vec(i_r), R_sph, u_rad, n_rad, n_tot, n_e, N_MC_r(i_r), d_Vol         !!! added 09.05.23, after test_12
!              write(200,14) ct, self%xi_cell_vec(i_r), self%kT_vec(i_r), self%M_cell_vec(i_r), R_sph, u_rad, n_rad, n_tot, n_e, N_MC_r(i_r), d_Vol, Z2_nI  !!! added 10.05.23, after test_13
              write(200,14) ct, self%xi_cell_vec(i_r), self%kT_vec(i_r), self%M_cell_vec(i_r), R_sph, u_rad, n_rad, n_tot, n_e, N_MC_r(i_r), d_Vol, Z2_nI, self%xi_cell_bnd_vec(i_r)  !!! added 11.05.23, after test_13a
                    !!!     0         1                      2                   3               4      5      6      7     8       9          10     11             12
              
                            !!! Eventually, should create a more logical order of params.
            
            end do
            close(200)
            
  !           if (alpha_r.eq.1.0_rk) then
  !             tau_T = 2.0_rk*M_ej/(4.0_rk*pi*R_sph_0**2.0_rk)*Zav*sigma_T/m_p*(-log(xi_min))/(1.0_rk - xi_min**2.0_rk)
  !           else
  !             tau_T = (3.0_rk - alpha_r)/(1.0_rk - alpha_r)*M_ej/(4.0_rk*pi*R_sph_0**2.0_rk)*Zav*sigma_T/m_p*(1.0_rk - xi_min**(1.0_rk - alpha_r))/(1.0_rk - xi_min**(3.0_rk - alpha_r))
  !           end if
  
            deallocate(E_rad_r,N_rad_r,N_MC_r)
            
!           else
!             !!! do nothing
!           end if
        
        
          i_save = i_save + 1 
        
      End Subroutine Printout_Sphere_qties
      
      
      
      Subroutine Printout_Energy_conserv(self,ct,sw_reset_Photons_escape)
      
        use Energy_conserv, ONLY: En_Compt_ph_loss, En_Compt_bulk_gain, En_Compt_bulk_heat_gain, En_brems_ph_loss, En_brems_bulk_gain, En_brems_bulk_heat_gain
      
        use Photon_distribution, ONLY: Photons_vec
        use Photon_distribution_escape, ONLY: Photons_escape
        use Init_printout, ONLY: FF
      
        implicit none
        
        class(Sphere_homol), intent(in) :: self
        real(kind=rk), intent(in) :: ct
        real(kind=rk) :: En_rad_tot, En_rad_esc_step, R_sph
        real(kind=rk), save :: En_rad_esc = 0.0_rk
        integer, save :: i_save = 0
        logical, intent(in) :: sw_reset_Photons_escape
        
          4 format(e14.7,' ',e14.7,' ',e14.7,' ',e14.7)
          14 format(e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7,' ',e14.7)
          
          R_sph = self%R_sph(ct)
          
          if (allocated(Photons_vec)) then
            En_rad_tot = sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist)
          else
            En_rad_tot = 0.0_rk
          end if
            
          if (allocated(Photons_escape)) then  
            En_rad_esc_step = sum(Photons_escape%x*Photons_escape%weight, mask = Photons_escape%exist)
          else
            En_rad_esc_step = 0.0_rk
          end if
          
          if (sw_reset_Photons_escape) then
            En_rad_esc = En_rad_esc + En_rad_esc_step
          else
            En_rad_esc = En_rad_esc_step
          end if
          
            if (i_save.eq.0) then
              open(200,FILE='results/Tests/Energy_conserv_' // trim(adjustl(FF)) // '.dat')
            else
              open(200,POSITION='APPEND',FILE='results/Tests/Energy_conserv_' // trim(adjustl(FF)) // '.dat')
            end if
            
            write(200,14) ct, R_sph, En_rad_tot, En_rad_esc, En_Compt_ph_loss, En_Compt_bulk_gain, En_Compt_bulk_heat_gain, En_brems_ph_loss, En_brems_bulk_gain, En_brems_bulk_heat_gain
            
            close(200)
            
            i_save = i_save + 1 
            
            print *,'Printout_Energy_conserv DONE'
        
      End Subroutine Printout_Energy_conserv
      
    
    
    !!!##############################################################################################################################################################
    !!!####################################################End: Printouts from homologous sphere ####################################################################
    !!!##############################################################################################################################################################    
    
    
!!!##########################################################################################################################################################################
!!!################################################# End: Homologous sphere #################################################################################################
!!!##########################################################################################################################################################################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
End Module class_Sphere_homol



!!!#######################################################################################################################################
!!!################################ Store specific ejecta objects ########################################################################

Module class_Some_yet_undefined_ejecta_type
  implicit none
  type, public :: Some_yet_undefined_ejecta_type
  end type Some_yet_undefined_ejecta_type
End Module class_Some_yet_undefined_ejecta_type


Module Ejecta_objects
  use class_Sphere_homol
  use class_Some_yet_undefined_ejecta_type
  implicit none
  type(Sphere_homol) :: Sphere_object
  type(Some_yet_undefined_ejecta_type) :: Some_yet_undefined_ejecta_quantity
End Module Ejecta_objects

!!!############################### End: Store specific ejecta objects ####################################################################
!!!#######################################################################################################################################


!!!###############################################################################################################################################################################################
!!!########################################################### End: Specific ejecta objects ######################################################################################################
!!!###############################################################################################################################################################################################
