Module Propagator

  use realkind
  use constants

  implicit none

  !!!########### Runtime quantities, accessible to all subroutines in the module ##############
  private
  
!   real(kind=rk) ::
  
  !!!##########################################################################################
  
  public Evolve_photon_diffusion
  
  Contains


!     Subroutine Get_ejecta_params(r,ct,rho,Beta_bulk,kT,R_sph,Zav)        !!! rho_loc, Beta_bulk_ph, kT (out): at the photon's location; R_sphere - radius of the sphere at ct
    Subroutine Get_ejecta_params(r,ct,rho,Beta_bulk,kT,Zav,Zz)              !! removed R_sph to maintain generality

      !! In the SLSN code kT is determined in Aux_routines.f90/Get_thermal_quantites and is essentially the blackbody radiation temperature
      
      !!! 08.05.23: rho now treated as in the local rest frame
      
      use Ejecta_generic, ONLY: Ejecta_obj
    
      implicit none
      
      real(kind=rk), intent(in) :: r, ct
      real(kind=rk), intent(out) :: rho, Beta_bulk, kT, Zav, Zz
      real(kind=rk) :: Gamma_bulk
            
        
        Beta_bulk = Ejecta_obj%v_rad(ct,r)
        Gamma_bulk = 1.0_rk/(1.0_rk - Beta_bulk**2.0_rk)**0.5_rk      !!! Introduced 08.05.23
        
        rho = Ejecta_obj%rho(ct,r)/Gamma_bulk                         !!! Introduced Gamma_bulk on 08.05.23: rho in the local rest frame
        
        kT = Ejecta_obj%kT(ct,r)
        Zav = Ejecta_obj%Zav()
        Zz = Ejecta_obj%Zz()

    End Subroutine Get_ejecta_params
    
    
    
    
    Subroutine Get_photon_max_timestep(eps,ct,ct_max_step)
      !!! Should one also have criteria based on opacities (Compt, and especially brems since it destroys the photon)???
      
      use Ejecta_generic, ONLY: Ejecta_obj
      implicit none
      real(kind=rk), intent(in) :: eps, ct
      real(kind=rk), intent(out) :: ct_max_step
      real(kind=rk) :: ct_lightcross, ct_cool_min
        !ct_max_step = Ejecta_obj%R_sph(ct)*eps  !!! Rudimentary   
        
        ct_lightcross = Ejecta_obj%domain_size(ct)  !!! Rudimentary  
        ct_cool_min = Ejecta_obj%Min_cool_time(ct)  !!! Minimal cooling time in the domain
        
        
        
        ct_max_step = eps*min(ct_lightcross,ct_cool_min)
        
        if (1.eq.0) print *,'Get_photon_max_timestep: ct_lightcross, ct_cool_min, ct_max_step', ct_lightcross, ct_cool_min, ct_max_step
        
    End Subroutine Get_photon_max_timestep
    
    
 
    Subroutine Check_escape(Photon,chk_esc)
      use data_type_def, ONLY: Photon_param
      use Ejecta_generic, ONLY: Ejecta_obj
      implicit none
      type(Photon_param) :: Photon
      logical, intent(out) ::  chk_esc
        chk_esc = Ejecta_obj%chk_esc(Photon%ct,Photon%r)
    End Subroutine Check_escape
    
    
    
    Subroutine Update_Temp(r,ct,d_En_matter,chk_err)
      !!! d_En_matter - dim.-less energy gained by matter due to Compton scattering, in the local fluid frame
      use Ejecta_generic, ONLY: Ejecta_obj
      implicit none
      
        real(kind=rk), intent(in) :: r, ct, d_En_matter
        real(kind=rk) :: N_particles, N_particles_eff, kT, d_kT, kT_new, alpha
        logical, intent(out) :: chk_err
        
        chk_err = .FALSE.
        
        alpha = 3.0_rk/2.0_rk !! wrong if kT relativistic (then should be 3)
        N_particles = Ejecta_obj%M_cell(ct,r)/(Ejecta_obj%mu_mol()*m_p)
        N_particles_eff = N_particles*Ejecta_obj%c_Tambov()
        
        kT = Ejecta_obj%kT(ct,r)
        d_kT = d_En_matter/(alpha*N_particles_eff)
        kT_new = kT + d_kT
        call Ejecta_obj%set_kT(ct,r,kT_new)
        
        if (kT_new.lt.0.0_rk) then
          print *,'Update_Temp ERROR: kT_new < 0'
          print *,'kT_new, kT, d_kT ', kT_new, kT, d_kT 
          print *,'d_En_matter, N_particles, N_particles_eff', d_En_matter, N_particles, N_particles_eff
          print *,'Ejecta_obj%M_cell(ct,r)', Ejecta_obj%M_cell(ct,r)
          print *,'r/Ejecta_obj%domain_size(ct)', r/Ejecta_obj%domain_size(ct)
          chk_err = .TRUE.
!           stop
        end if
              
    End Subroutine Update_Temp
 
 
  
    Subroutine Evolve_photon_diffusion(ct_in,ct_fin)
    
      !!! Based on bits from /media/indrek/Samsung_T5/work/Proged_T5/Nonlinear_MC/Continuous_enviroment/devel/thermal pair prod/Testbed_removing_chi_etc_271221/Non_lin_MC.f90: Evolve_spherical_diffusion_twozone_continuous_ejecta

      !! ct_max -> ct_fin
      !! Removed R_sphere_loc; R_sphere has a new meaning of CURRENT radius (previously was the radius at the beginning of the time step)
      !! eps2 -> eps
      !! n_e_loc, rho_loc -> n_e, rho
      
      use Photon_distribution, ONLY: Photons_vec
      use Ejecta_generic, ONLY: Ejecta_obj
      use Geometry, ONLY: get_v_vec_from_angles_in_r_system, get_photon_angles_in_rz_system
      use auxiliary, ONLY: MC_append_escape, MC_photon_cleanup, MC_photon_cull, Findvalue_1dim_v3_nearestindex
      use Phys_procs, ONLY: Scatter_photon

      use Cumulative_distributions, ONLY: P_cumul_MC_Compton,     &		!! P_cumul_MC_Compton% [kappa_pr_norm, P_mu_erf_ph,P_mu_pr_e,P_z]
                                            MC_Compton_grids,     &
                                            gff
                                            
                                            
      use Photon_distribution_escape, ONLY: Photons_escape, i_exist_max_escape          !!! NBNBNB: Just for runtime printouts to terminal, doesn't affect the code run    
      
      use Energy_conserv, ONLY: En_Compt_ph_loss, En_brems_ph_loss, En_brems_bulk_gain, En_brems_bulk_heat_gain
      
      !!!###### Temporary, only for runtime tests #####
!          use Energy_conserv, ONLY: En_Compt_ph_loss
      !!!##############################################   
      
      use Switches, ONLY: sw_Therm_brems  ! sw_Compton,   (CCurrently Compton computed anyway)

      implicit none
      
      real(kind=rk), intent(in) :: ct_in, ct_fin
      real(kind=rk) :: ct, d_ct
      real(kind=rk) :: eps, ct_max_step
      real(kind=rk) :: rho, n_e, Z2_nI, kT, Zav, Zz
      real(kind=rk) :: Beta_bulk_ph, Gamma_bulk_ph, Z_bulk_ph
      real(kind=rk) :: mu, Doppler, x_pr, mu_pr
      real(kind=rk) :: theta, phi
      real(kind=rk), dimension(3) :: r_loc_vec, v_vec
      real(kind=rk) :: x_sc, x_sc_pr, theta_sc, phi_sc, z_drawn
      real(kind=rk) :: alpha_Compton, alpha_photoion, alpha_brems, alpha_tot, lambda
      real(kind=rk) :: rnd1, rnd3
      real(kind=rk) :: d_s_max, d_s_inter, d_s_step
      real(kind=rk) :: ratio1, ratio2, dummy
      real(kind=rk) :: d_En_matter
      real(kind=rk) :: Cf0_brems_abs
      integer :: i_steps
      integer :: i, n_ph_bin
      integer :: i_x, i_kT
      logical :: chk_inter, chk_esc, chk_within
      logical :: chk_Compton_event, chk_photoion_event, chk_brems_event
      logical :: chk_err

      
      ct = ct_in
      
      eps = 0.1_rk ! 0.001_rk ! 0.5_rk  ! 0.1_rk   !! preliminary  
      
      
      Cf0_brems_abs = (2.0_rk*pi/3.0_rk)**0.5_rk/(4.0_rk*pi**2.0_rk)*sigma_T*afs*lambda_C**3.0_rk*gff
      print *,'Cf0_brems_abs', Cf0_brems_abs
      
      
      i_steps = 0
      do while (ct.lt.ct_fin)
      
        i_steps = i_steps + 1
        
        if (i_steps/100*100.eq.i_steps) print *,'i_steps, ct, ct_fin', i_steps, ct, ct_fin
        
        !!!###### How many photons ######
          n_ph_bin = 0
          if (allocated(Photons_vec)) then
            n_ph_bin = size(Photons_vec)    !! For some reason size(Photons_vec) if UNALLOCATED
          else
            n_ph_bin = 0
          end if    
        !!!##############################

        !!!##### Maximal timestep ###############  
          call Get_photon_max_timestep(eps,ct,ct_max_step) !! limit should be based on the timescale over which the target quantities (or geometry) changes. In SLSN simul, this included the photon fileld itself.
        !!!######################################
        
        !!!################################################ Main loop through photons; emission #####################################################################################################
          d_ct = min(ct_max_step,ct_fin-ct)

          if (1.eq.0) print *,'BEFORE: minval(Photons_vec%x), maxval(Photons_vec%x)', minval(Photons_vec%x, mask = Photons_vec%exist), maxval(Photons_vec%x, mask = Photons_vec%exist)
          
          if (sw_Therm_brems) then
            call Ejecta_obj%Emit_thermal(ct,d_ct)
            n_ph_bin = size(Photons_vec)  !! Should recheck if allocated
          end if
          
          if (1.eq.0) then
            print *,'AFTER: minval(Photons_vec%x), maxval(Photons_vec%x)', minval(Photons_vec%x, mask = Photons_vec%exist), maxval(Photons_vec%x, mask = Photons_vec%exist)
            print *,'i_steps, n_ph_bin', i_steps, n_ph_bin
            print *,'ct, ct_fin', ct, ct_fin
          end if
          
          ct = ct + d_ct
          
          do i = 1,n_ph_bin
          
!             print *,'i', i
!             if (i/100*100.eq.i) print *,'i', i
            
            do while ((Photons_vec(i)%ct.lt.ct).and.(Photons_vec(i)%exist))

              d_s_max = ct - Photons_vec(i)%ct

              !!!################################## Get ejecta qties ####################################
                call Get_ejecta_params(Photons_vec(i)%r,Photons_vec(i)%ct,rho,Beta_bulk_ph,kT,Zav,Zz)        !!! rho,Beta_bulk_ph, kT (out): at the photon's location
                                                                                                             !!! NB: rho, kT in the fluid frame (as ooposed to e.g. Ejecta_obj%rho)
                                                
                Gamma_bulk_ph = 1.0_rk/(1.0_rk - Beta_bulk_ph**2.0_rk)**0.5_rk                                     !!! NBNB: kT was not here in the SLSN code
                Z_bulk_ph = Beta_bulk_ph*Gamma_bulk_ph  
              !!!########################################################################################
              
              !!!############### Define quantities for computing free paths, no longer separately for each process (05.03.21) ################
                mu = cos(Photons_vec(i)%theta)
                Doppler = 1.0_rk/(Gamma_bulk_ph - Z_bulk_ph*mu)
                x_pr = Photons_vec(i)%x/Doppler
              !!!#############################################################################################################################
              
              !!!#### Currently (08.50.23) densities treated as in the LAB frame ##################
                n_e = rho*Zav/m_p                   !!! NBNBNBNB: Zav should be redefined if there are a significant number of pairs)!!! 08.05.23: rho and n_e in the local REST FRAME
                Z2_nI = rho*Zz/m_p                  !!! 08.05.23: rho and n_e in the local REST FRAME
              !!!##################################################################################
              
              !!!################################################# Opacities ########################################################################
                !!!########## Photoionization ###########
                  alpha_photoion = 0.0_rk
                !!!######################################
                
                !!!########## Thermal bremsstrahlung ######
                  if (sw_Therm_brems) then
                    alpha_brems = Cf0_brems_abs*n_e*Z2_nI*kT**(-0.5_rk)*x_pr**(-3.0_rk)*(1.0_rk - exp(-x_pr/kT))*(1.0_rk - Beta_bulk_ph*mu)*Gamma_bulk_ph  !!! 08.05.23: introduced Gamma_bulk_ph, so n_e, Z2_nI now in the RF 
                  else
                    alpha_brems = 0.0_rk
                  end if
                !!!########################################
                
                !!!########## Compton scattering #########
                  call Findvalue_1dim_v3_nearestindex(x_pr,MC_Compton_grids%x_pr_vec(:),MC_Compton_grids%x_pr_vec(:),dummy,i_x)         !!!! NBNBNBNB: i_x used in multiple locations independently
                  call Findvalue_1dim_v3_nearestindex(kT,MC_Compton_grids%kT_vec(:),MC_Compton_grids%kT_vec(:),dummy,i_kT)      !!!! NBNBNB: kT was from Sphere_params in SLSN code, NOT so here
                  alpha_Compton = sigma_t*n_e*(1.0_rk - Beta_bulk_ph*mu)*Gamma_bulk_ph*P_cumul_MC_Compton%kappa_pr_norm(i_x,i_kT)             !!! sigma --> sigma_t (18.04.23); 08.05.23: introduced Gamma_bulk_ph, so n_e now in the RF 
                !!!#######################################
                
                alpha_tot = alpha_Compton + alpha_photoion + alpha_brems
              !!!################################################# End: Opacities ###################################################################
              
              !!!##### Draw photon propagation distance ####################################################### 
                if (alpha_tot.gt.0.0_rk) then
                  lambda = 1.0/alpha_tot
                  call RANDOM_NUMBER(rnd1)
                  d_s_inter = -log(rnd1)*lambda		!! distance photon would propagate before scattering in uniform medium
                else
                  lambda = 1.0e90_rk
                  d_s_inter = 1.0e90_rk
                end if
              !!!##############################################################################################
              
              
              !!!########## Determine whether an interaction occurs ##############
                d_s_step = min(d_s_inter,d_s_max)
                if (d_s_inter.lt.d_s_max) then
                  chk_inter = .TRUE.
                else
                  chk_inter = .FALSE.
                end if
              !!!#################################################################
              
              !!!######## Update photon location ########################
                Photons_vec(i)%ct = Photons_vec(i)%ct + d_s_step
                Photons_vec(i)%r_loc_vec(:) = Photons_vec(i)%r_loc_vec(:) + d_s_step*Photons_vec(i)%v_vec(:)
                Photons_vec(i)%r = norm2(Photons_vec(i)%r_loc_vec(1:3))
                call get_photon_angles_in_rz_system(Photons_vec(i)%r_loc_vec,Photons_vec(i)%v_vec,Photons_vec(i)%theta,Photons_vec(i)%phi)  !! ## Updatingpdate angles (which are rel. to radial direction) after propagating the photon
              !!!########################################################
              
              !!!########## Update ejecta qties after propagating the photon ##############
                call Get_ejecta_params(Photons_vec(i)%r,Photons_vec(i)%ct,rho,Beta_bulk_ph,kT,Zav,Zz)        !!! rho_loc, Beta_bulk_ph (out): at the photon's location
                                                                                                             !!! NB: rho, kT in the fluid frame (as ooposed to e.g. Ejecta_obj%rho)
                Gamma_bulk_ph = 1.0_rk/(1.0_rk - Beta_bulk_ph**2.0_rk)**0.5_rk                                     !!! NBNB: kT was not here in the SLSN code
                Z_bulk_ph = Beta_bulk_ph*Gamma_bulk_ph 
              !!!##########################################################################
              
              !!!############### Added updated mu, Doppler, x_pr here (21.06.21) #############################################################
                mu = cos(Photons_vec(i)%theta)
                Doppler = 1.0_rk/(Gamma_bulk_ph - Z_bulk_ph*mu)
                x_pr = Photons_vec(i)%x/Doppler   !! In fact, x_pr is not used below
                mu_pr = (mu - Beta_bulk_ph)/(1.0_rk - Beta_bulk_ph*mu)
              !!!#############################################################################################################################
              
              !!!###### Check whether photon escapes ############
                chk_esc = Ejecta_obj%chk_esc(Photons_vec(i)%ct,Photons_vec(i)%r)
                chk_within = Ejecta_obj%chk_within(Photons_vec(i)%ct,Photons_vec(i)%r)
                
!                 if (Photons_vec(i)%r.gt.R_sphere) then   !!! if (Photons_vec(i)%r.gt.R_sphere_loc) then  !! Photons escape if r > R_sphere

!                 !!! TEST, WRONG!!! #########
!                   if (chk_esc) then
!                     if (Photons_vec(i)%theta.gt.(pi/2.0_rk)) Photons_vec(i)%theta = pi - Photons_vec(i)%theta
!                     chk_inter = .FALSE.
!                     chk_esc = .FALSE.
!                   end if
!                 !!!############################


                if (chk_esc) then
                  call MC_append_escape(Photons_vec(i))   !! append photon to escape array (Photons_escape), adjust i_exist_max_escape
                  Photons_vec(i)%exist = .FALSE.
                  chk_inter = .FALSE.
                end if
                
                if (.NOT.(chk_within)) then  !! Outside ejecta, but not necessarily escaping (i.e. within a cavity)
                  chk_inter = .FALSE.
                end if
              !!!################################################
              
              !!!################# Perform interaction #############################
                if (chk_inter) then
                  ratio1 = alpha_Compton/(alpha_Compton + alpha_brems + alpha_photoion)
                  ratio2 = (alpha_Compton + alpha_brems)/(alpha_Compton + alpha_brems + alpha_photoion)
                  call RANDOM_NUMBER(rnd3)
                  if (rnd3.le.ratio1) then
                    chk_Compton_event = .TRUE.
                    chk_brems_event = .FALSE.
                    chk_photoion_event = .FALSE.
                  else if (rnd3.le.ratio2) then
                    chk_Compton_event = .FALSE.
                    chk_brems_event = .TRUE.
                    chk_photoion_event = .FALSE.
                  else
                    chk_Compton_event = .FALSE.
                    chk_brems_event = .FALSE.
                    chk_photoion_event = .TRUE.
                  end if
!                 end if
                
                  if (chk_Compton_event.and.(rho.gt.0.0_rk)) then   !!! rho.gt.0.0_rk ensures that the photon hasn't travelled into a cavity after the integaction was determined to happen (chk_inter)
                                                                    !!! Should improve the propagation into and out of the cavity so that the photon travels exactly to the boundary
                                                
! ! !                     Z_bulk_ph = 0.0_rk !! WRONG, TEST
!                     call Scatter_photon(x_sc,x_sc_pr,theta_sc,phi_sc,Photons_vec(i)%x,Photons_vec(i)%theta,Photons_vec(i)%phi,Photons_vec(i)%weight,Z_bulk_ph,kT)    !!! Check how kT is determined!
                    call Scatter_photon(x_sc,x_sc_pr,theta_sc,phi_sc,z_drawn,Photons_vec(i)%x,Photons_vec(i)%theta,Photons_vec(i)%phi,Photons_vec(i)%weight,Z_bulk_ph,kT)    !!! Check how kT is determined!
                    
                    !!!######## Update local temperature ############
                      d_En_matter = (x_pr - x_sc_pr)*Photons_vec(i)%weight               !!! Matter gain/photon loss in the fluid frame
                      call Update_Temp(Photons_vec(i)%r,Photons_vec(i)%ct,d_En_matter,chk_err)
                      if (chk_err) then
                        print *,'Evolve_photon_diffusion/Update_Temp ERROR: kT < 0'
                        print *,' x_pr, x_sc_pr, Photons_vec(i)%weight', x_pr, x_sc_pr, Photons_vec(i)%weight
                        print *,'Z_bulk_ph, kT, z_drawn', Z_bulk_ph, kT, z_drawn
                        stop
                      end if
                    !!!##############################################
                    
                    !!!################## Update photon ###################
                      Photons_vec(i)%x = x_sc
                      Photons_vec(i)%theta = theta_sc
                      Photons_vec(i)%phi = phi_sc                        
                      Photons_vec(i)%count_scatter = Photons_vec(i)%count_scatter + 1
                      
                      !!!#### Need to update v_vec after altering the angles #################  
                        theta = Photons_vec(i)%theta
                        phi = Photons_vec(i)%phi
                        r_loc_vec(:) = Photons_vec(i)%r_loc_vec(:)
                        call get_v_vec_from_angles_in_r_system(theta,phi,r_loc_vec(:),v_vec(:))
                        Photons_vec(i)%v_vec(:) = v_vec(:) 
                      !!!#####################################################################
                    !!!####################################################

                  else if (chk_brems_event.and.(rho.gt.0.0_rk)) then
                  
                    d_En_matter = x_pr*Photons_vec(i)%weight
                    
                    !!!################### Energy deposition in the flow ##########################
                      En_brems_ph_loss = En_brems_ph_loss + Photons_vec(i)%x*Photons_vec(i)%weight
                      En_brems_bulk_gain = En_brems_bulk_gain + Z_bulk_ph*x_pr*mu_pr*Photons_vec(i)%weight
                      En_brems_bulk_heat_gain = En_brems_bulk_heat_gain + Gamma_bulk_ph*x_pr*Photons_vec(i)%weight
                    !!!############################################################################
                    
                    call Update_Temp(Photons_vec(i)%r,Photons_vec(i)%ct,d_En_matter,chk_err)
                    if (chk_err) then
                      print *,'Evolve_photon_diffusion/Update_Temp ERROR (from bremsstrahlung absorption)'
                      stop
                    end if
                    
                    Photons_vec(i)%exist = .FALSE.
                    
                    
                  else if (chk_photoion_event.and.(rho.gt.0.0_rk)) then !!! rho.gt.0.0_rk ensures that the photon hasn't travelled into a cavity after the integaction was determined to happen (chk_inter)
                                                                        !!! Should improve the propagation into and out of the cavity so that the photon travels exactly to the boundary
                    print *,'Cannot handle photoionization yet'
                    
                    print *,'rnd3', rnd3
                    print *,'ratio1', ratio1
                    print *,'ratio2', ratio2
                    print *,'alpha_Compton', alpha_Compton
                    print *,'alpha_photoion', alpha_photoion
                    print *,'alpha_brems', alpha_brems
                    print *,'d_s_inter, d_s_max, lambda', d_s_inter, d_s_max, lambda
                    
                    stop
                  end if
                  
                end if  !!! if (chk_inter) then
                
              !!!###################################################################
                            
            end do !!! do while ((Photons_vec(i)%ct.lt.ct).and.(Photons_vec(i)%exist))
          end do    !! do i = 1,n_ph_bin
        !!!################################################ End: Main loop through photons ##########################################################################################################  
        
        
        !!!######### Photon emission from thermal plasma ###########
          ! ct (updated)
          ! d_ct (time interval)
          ! call Emit_thermal
        !!!#########################################################
        
        
        call MC_photon_cleanup      !!! 08.05.23
        call MC_photon_cull         !!! 08.05.23
        
        if (1.eq.0) then
        !!!#### Temporary ####
          print *,'Energy in photons', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist)
          if (allocated(Photons_escape)) then
            print *,'Energy in escaped photons, i_exist_max_escape', sum(Photons_escape%x*Photons_escape%weight, mask = Photons_escape%exist), i_exist_max_escape
            print *,'Total energy in system', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist) + sum(Photons_escape%x*Photons_escape%weight, mask = Photons_escape%exist)
            print *,'Number of escaped photons', sum(Photons_escape%weight, mask = Photons_escape%exist)
          end if
          print *,'Number of photons', sum(Photons_vec%weight, mask = Photons_vec%exist)
          print *,'Number of MC photons', sum(Photons_vec%weight/Photons_vec%weight, mask = Photons_vec%exist)
        !!!###################
        end if
        
!         stop
      end do  !!! do while (ct.lt.ct_fin)
      
      print *,'##############################################################'
      
      print *,'Average energy per photon', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist)/sum(Photons_vec%weight, mask = Photons_vec%exist)
      print *,'Average energy per escaping photon', sum(Photons_escape%x*Photons_escape%weight, mask = Photons_escape%exist)/sum(Photons_escape%weight, mask = Photons_escape%exist)
      
      print *,'Energy in photons', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist)
      print *,'Energy in escaped photons, i_exist_max_escape', sum(Photons_escape%x*Photons_escape%weight, mask = Photons_escape%exist), i_exist_max_escape
      print *,'Total energy in system', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist) + sum(Photons_escape%x*Photons_escape%weight, mask = Photons_escape%exist)
      
      print *,'Photon energy loss by Compton', En_Compt_ph_loss
      print *, 'Total energy budget (shold be constant)', sum(Photons_vec%x*Photons_vec%weight, mask = Photons_vec%exist) + sum(Photons_escape%x*Photons_escape%weight, mask = Photons_escape%exist) + En_Compt_ph_loss
      
      print *,'Number of photons', sum(Photons_vec%weight, mask = Photons_vec%exist)
      print *,'Number of escaped photons', sum(Photons_escape%weight, mask = Photons_escape%exist)
      
      print *,'##############################################################'


    End Subroutine Evolve_photon_diffusion
  
  
End Module Propagator  
