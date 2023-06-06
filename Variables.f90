
!!!############################################## Type definitions ##################################################################################

Module data_type_def  

  use Realkind
 
    type Photon_param
  ! 	sequence		!! gives error if in
      real(kind=rk) :: ct, r, x, theta, phi, weight !, mu
      real(kind=rk) :: ct_aux
      integer :: id                 !!! NBNBNB: NOT UNIQUE
      logical :: exist
      
      real(kind=rk), dimension(3) :: r_loc_vec, v_vec      !!! As of 16.08.18, used only when calculating diffusion in a spherical shell with holes (i_medium = 4)
                                                          !!! NBNBNB: v_vec is normalized to unity (it is essentially beta)	

      real(kind=rk) :: ct_sourced  !!! time the original particle giving rise to it was created (could be generations ago)  (introduced 03.07.19)
      real(kind=rk) :: r_sourced                                                    !! Introduced 22.04.20 (diagnostics only); set in Emit_el_cool_delta_t_and_append (exclusively at present)
                                   
      !!!######### Commented out for passive diffusion model; was in for SLSNe simul ############                             
!       !!! diagnostics: not needed for code to run ###
      integer :: count_scatter
!       integer :: count_shock_crossings_ph = 0                                       !! Introduced 20.04.20 (diagnostics only)
!       integer :: count_shock_crossings_el = 0                                       !! Introduced 20.04.20 (diagnostics only)
!       integer :: count_shock_crossings_ph_current = 0                               !! Introduced 21.04.20 (diagnostics only)
!       integer :: count_generation = 0                                               !! Introduced 21.04.20 (diagnostics only)
!       integer :: i_zone_born = 0                                                    !! Introduced 22.04.20 (diagnostics only); set in Emit_el_cool_delta_t_and_append (exclusively at present)
!       integer :: i_zone = 0                 !!! 05.11.20: introduced explicit zone label, for tackling nonthermal pair prod. in two zones; initialized to 0 ("no zone")  
      !!!########################################################################################                   
                         
    end type Photon_param
    
    
    !!!############### Compton scattering of photons with a thermal electron distribution ###############################
      type grid_MC_Compton_type
        real(kind=rk), dimension(:), allocatable :: x_pr_vec, z_vec, mu_rel_pr_vec, x_erf_vec, mu_sc_erf_vec, kT_vec ! , theta_vec
      end type grid_MC_Compton_type
      
      type P_cumul_scatt
        real(kind=rk), dimension(:,:), allocatable :: kappa_pr_norm, P_mu_erf_ph
        real(kind=rk), dimension(:,:,:), allocatable :: P_mu_pr_e
        real(kind=rk), dimension(:,:,:,:), allocatable :: P_z 
      end type P_cumul_scatt
    !!!##################################################################################################################
    
    !!!############ Bremsstrahlung emission from thermal distribution of electrons ###########################
      type P_cumul_Therm_brems_type
        real(kind=rk), dimension(:), allocatable :: P_y, P_y_unnorm, lny   !!! Cuurently (02.05.23), only P_y_unnorm = I(y) - I(y_min) is used (prop. to number of photons emitted between y_min and y, where y = x/theta)
        real(kind=rk) :: Ndot_norm
      end type P_cumul_Therm_brems_type
    !!!#######################################################################################################
    
End Module data_type_def

!!!##################################################################################################################################################


!!!##################################### Quantities associated with MC particles ############################################################

Module Photon_distribution
  use data_type_def, ONLY: Photon_param
  implicit none
    type(Photon_param), dimension(:), allocatable ::  Photons_vec
!     type(Photon_param), dimension(:), allocatable ::  Photons_secondary  !! removed
End Module Photon_distribution


Module Photon_distribution_escape
  use data_type_def, ONLY: Photon_param
  implicit none
    type(Photon_param), dimension(:), allocatable :: Photons_escape
    integer :: i_exist_max_escape                !!! Index of last existing photon; reset at each timestep IF sw_reset_Photons_escape = .TRUE.
End Module Photon_distribution_escape

!!!###########################################################################################################################################


!!!################################# Physical processes: grids, rates, cumul distr. ###################################

Module Cumulative_distributions
  !! NB: Should probably separate modules according to processes; also more than cumulative distributions, so should probably be renamed

  use Realkind
  use data_type_def

  implicit none
  
  !!!############ Compton for photons interacting with thermal electrons #################################################################
  type(P_cumul_scatt) :: P_cumul_MC_Compton  !!! P_cumul -> P_cumul_Compton -> P_cumul_MC_Compton;
  type(grid_MC_Compton_type) :: MC_Compton_grids		!!! MC_Compton_grids% [ x_pr_vec, z_vec, mu_rel_pr_vec, x_erf_vec, mu_sc_erf_vec, kT_vec ]   !! moved here on 26.10.21
  
  real(kind=rk) :: x_pr_min = 1.0e-8_rk       !! 1.0e-6_rk in SLSN code
  real(kind=rk) :: x_pr_max = 1.0e0_rk        !! 2.0e5_rk in SLSN code
  real(kind=rk) :: mu_rel_min = -1.0_rk
  real(kind=rk) :: mu_rel_max = 1.0_rk
  
  real(kind=rk) :: z_min = 1.0e-5_rk          !! 
  real(kind=rk) :: z_max = 1.0e0_rk           !! 1.0e-1 in SLSN code
  real(kind=rk) :: kT_min = 1.0e-7_rk         !!
  real(kind=rk) :: kT_max = 1.0e0_rk         !! 1.0e-4 in SLSN code; kT_max = 1.0e-1_rk up to test_15 (included)
  real(kind=rk) :: x_erf_min = 1.0e-8_rk      !! 1.0e-6 in SLSN code
  real(kind=rk) :: x_erf_max = 1.0e0_rk       !! 1.0e6 in SLSN code; CAREFUL, no MeV GeV photons if x_erf_max too low
  real(kind=rk) :: mu_sc_erf_min = -1.0_rk
  real(kind=rk) :: mu_sc_erf_max = 1.0_rk
  
  integer :: n_z = 100
  integer :: n_x = 200
  integer :: n_mu = 400
  integer :: n_x_erf = 400      !! 200 ! hurry up (for testing)
  integer :: n_mu_sc_erf = 400  !! 200 ! hurry up (for testing)
  integer :: n_kT = 60          !!! 30 up to test_15 (included)
  !!!#####################################################################################################################################
  
  
  !!!############ Bremsstrahlung emission from thermal distribution of electrons #########################################################
    type(P_cumul_Therm_brems_type) :: P_cumul_Therm_brems
    real(kind=rk) :: y_min_brems = 1.0e-8_rk !! (10.05.23) !!! 1.0e-5_rk   !!! minimum y = x/theta   
    real(kind=rk) :: y_max_brems = 20.0_rk !! (10.05.23) !! 10.0_rk     !!! maximum y = x/theta
    integer :: n_y_brems = 1000
    
    real(kind=rk), parameter :: gff = 1.0_rk !! preliminary
    real(kind=rk), parameter :: kT_floor_brems = 1.0e-6_rk
  !!!#####################################################################################################################################
  
End Module Cumulative_distributions  



Module Energy_conserv
  use realkind
  implicit none
    
    !!!####### Compton ######
      real(kind=rk) :: En_Compt_ph_loss = 0.0_rk                  !!! Total energy gained by the flow (in the lab frame);               ... + weight*(x - x_sc)
      real(kind=rk) :: En_Compt_bulk_gain = 0.0_rk                !!! Energy given to the bulk flow;                                    ... + weight*(x_pr*cos_theta_pr - x_sc_pr*cos_theta_sc_pr)*Z_bulk
      real(kind=rk) :: En_Compt_bulk_heat_gain = 0.0_rk           !!! Thermal energy gain by flow, transformed into the lab frame;      ... + weight*(x_pr - x_sc_pr)*Gamma_bulk
    !!!######################
    
    !!!####### Thermal brems ##########
      real(kind=rk) :: En_brems_ph_loss = 0.0_rk
      real(kind=rk) :: En_brems_bulk_gain = 0.0_rk
      real(kind=rk) :: En_brems_bulk_heat_gain = 0.0_rk
    !!!################################
                                                        
End Module Energy_conserv

!!!####################################################################################################################



!!!#################################### Environment parameters ##############################################

!!######################## NBNBNBNBNBNBNB: Superseded by a data file data/Passive_fireball_params.dat #################################
! Module Passive_fireball_params
!   use realkind
!   implicit none
!   
!     !!!####################### These quantities will be passed to the ejecta object (Ejecta_obj) ##########################
!     real(kind=rk), parameter :: M_ej = 1.0_rk
!     real(kind=rk), parameter :: v_ej = 0.1_rk
!     real(kind=rk), parameter :: xi_min = 0.0_rk ! 0.5
!     real(kind=rk), parameter :: alpha_r = 1.0_rk
!         
!     real(kind=rk), parameter :: Zav = 1.0_rk
!     
!     integer, parameter :: n_r = 20 ! 10 ! 30
!     real(kind=rk), dimension(n_r), parameter :: kT_vec = 1.0e-3_rk   !! Can be used to set a fixed an/or initial value (i.e. at a time when R_sph = R_sph_0); MC radiation initialization is dependent on this
!     
!     real(kind=rk), parameter :: R_sph_0 = 1.0_rk  !! Initial radius of the sphere
!     !!!####################################################################################################################
!     
! End Module Passive_fireball_params


Module Ejecta_choice_param
  use realkind
  implicit none
    integer, parameter :: i_choose_ejecta = 1  !!! 1 - Passive fireball; uses Sphere_homol class for ejecta (type(Ejecta) :: Ejecta_obj --> type(Sphere_homol) :: Sphere_object)
End Module Ejecta_choice_param

!!!##########################################################################################################



!!!################################## Parameters related to MC particles #########################################################

Module MC_photon_params
  use realkind
  implicit none
    integer, parameter :: n_MC_ph_0 = 1000 ! 10000 ! 30000   !!! Initial target number of MC photons in the system. CAREFUL not to mix up with a locally defined number of MC photons in subroutines;
                                            
    integer :: n_MC_ph = n_MC_ph_0            !!! Used at initialization, subsequently as the target MC photons density during the code run (i.e. can change as the code runs)
                                            
    integer :: n_ph_cull = 100*n_MC_ph_0 ! 5*n_MC_ph
    
    real(kind=rk) :: w_fiducial  !!! set at the initialization stage, using numbre of physical photons and n_MC_ph above; NB: Can change further as the code runs, e.g. by MC_photon_split(ct) if photons are split to improve statistics 
    
    
    
End Module MC_photon_params

!!!###############################################################################################################################



!!!################################### General parameters/switces related code running #####################################

Module Switches
  implicit none
    logical, parameter :: sw_Compton = .TRUE.
    logical, parameter :: sw_Therm_brems = .TRUE. ! .FALSE. ! .TRUE.  !!! NBNBNBNBNB: NOT relativistically correct (03.05.23)!!!
    
    logical, parameter :: sw_reset_Photons_escape = .TRUE.  !!! .FALSE. !!!### if TRUE, Photons_escape is RESET at the end of each timestep
    
End Module Switches

!!!#########################################################################################################################




