Program Main
  
  use RealKind
  use Constants
  use Initialize_physics, ONLY: Setup_grids_Compton_MC, P_cumul_distributions_Compton_MC, P_cumul_Thermal_brems_MC
  use Ejecta_init, ONLY: Initialize_Ejecta
  use Propagator, ONLY: Evolve_photon_diffusion
  use Printouts   !!! , ONLY: Printout_timedep_spectrum, Printout_Sphere_qties
! !   use MC_init, ONLY: Initialize_Radiation
  use Init_printout, ONLY: Read_Filename
  use Init_printout, ONLY: FF
  
  use Photon_distribution_escape, ONLY: Photons_escape

  use auxiliary, ONLY: MC_photon_split, MC_target_number
  
  use Ejecta_generic, ONLY: Ejecta_obj
    
  use Switches, ONLY: sw_Compton, sw_Therm_Brems, sw_reset_Photons_escape
  
    implicit none
    
!     real(kind=rk), dimension(:), allocatable ::
    real(kind=rk) :: ct_max, ct_in, d_ln_ct_nominal, ct, ct_next
    integer :: n_t_nominal
    integer :: i_t
        
!     sw_Compton = .TRUE.
!     sw_Therm_Brems = .TRUE.
    
      if (sw_Compton) then
        call Setup_grids_Compton_MC       !!! USES: Cumulative_distributions (Compton section)
        call P_cumul_distributions_Compton_MC   !!! FILLS: Cumulative_distributions, ONLY: P_cumul_MC_Compton
                                                !!! USES: Cumulative_distributions
                                                !!! CALLS: f_MW_norm, KN_cross_tot
      end if
      
      if (sw_Therm_Brems) then
        call P_cumul_Thermal_brems_MC
      end if
    
    call Read_Filename()   !!! Sets the filename extension for printouts in Init_printout/FF
    print *,'Running simulation ', FF
    
    call Initialize_Ejecta()    !!!! NBNBNB: ejecta object could be generalized further if the type-bound routines accepted r_loc_vec instead of r; ALSO includes initialization of MC photons
!     call Initialize_Radiation()   ! no longer a separate routine
    
    

    ct_in = Ejecta_obj%ct_0
    print *,'t_in (days)', ct_in/c0/3600.0/24.0
    
!     call Printout_timedep_spectrum(ct_in,.FALSE.)  !! Initial spectrum to file
    call Save_output(ct_in,sw_reset_Photons_escape)   !!! 2. parameter sw_reset_Photons_escape!!! Should be made .TRUE. to save memory and time !!! ct --> ct_in (08.05.23)
    
    !!!################## Preliminary #####################
    ct_max = ct_in*100.0_rk ! 300.0_rk ! 10.0_rk ! 3.0_rk ! 10.0_rk
    n_t_nominal = 100 ! 300 ! 30 ! 100
    d_ln_ct_nominal = (log(ct_max) - log(ct_in))/(dble(n_t_nominal) - 1.0_rk)
    !!!####################################################

    ct = ct_in
    i_t = 0
    do while (ct.lt.ct_max)
      i_t = i_t + 1
      print *,'i_t, t (days)', i_t, ct/c0/3600.0/24.0
      print *,'Ejecta_obj%domain_size(ct)', Ejecta_obj%domain_size(ct)

      ct_next = exp(log(ct) + d_ln_ct_nominal)
      call Evolve_photon_diffusion(ct,ct_next)
      ct = ct_next
      
      call Save_output(ct,sw_reset_Photons_escape)   !!! 2. parameter sw_reset_Photons_escape !!  Should be made .TRUE. to save memory and time
      
      !!!######## Adjust the target number of MC particles, also adjust n_ph_cull accordingly #################
        call MC_target_number(Ejecta_obj%domain_size(ct),Ejecta_obj%domain_size(ct_in))
      !!!######################################################################################################
      
      call MC_photon_split(ct,.FALSE.,2)  !! Don't force split, use splitting factor 2
      
      if (sw_reset_Photons_escape.AND.allocated(Photons_escape)) deallocate(Photons_escape)                     !!!### if TRUE, Photons_escape is RESET at the end of each timestep
            
    end do
    
    print *, FF
    print *,'Kahv on valmist!'

     
End Program Main
