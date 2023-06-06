Module Phys_procs

  use realkind
  use constants
  
  implicit none
  
    Contains
  

!!!##############################################################################################################################################################################     
!!!################################################## Scattering routines #######################################################################################################   
              
!       Subroutine Scatter_photon(x_sc,x_sc_pr,theta_sc,phi_sc,x,theta,phi,weight,Z_bulk,kT)
      Subroutine Scatter_photon(x_sc,x_sc_pr,theta_sc,phi_sc,z_val,x,theta,phi,weight,Z_bulk,kT)
      
      
  !     Subroutine Scatter_photon(x_sc,x_sc_pr,x_pr,theta_sc,phi_sc,x,theta,phi,Z_bulk,kT)    !!! added x_pr, Removed again
      
  !       use grids_MC_Compton, ONLY: MC_Compton_grids				!! MCMC_Compton_grids% [x_pr_vec, z_vec, mu_rel_pr_vec, x_erf_vec, mu_sc_erf_vec, kT_vec]
        use Cumulative_distributions, ONLY: P_cumul_MC_Compton,     &		!! P_cumul_MC_Compton% [kappa_pr_norm, P_mu_erf_ph,P_mu_pr_e,P_z]
                                            MC_Compton_grids              !! moved MC_Compton_grids to Cumulative_distributions on 26.10.21
        use auxiliary, ONLY: Findvalue_1dim_v3, Findvalue_1dim_v3_nearestindex
  !       use Test_qties, ONLY: Test_qties_MC			!! only for test

        use Energy_conserv, ONLY: En_Compt_ph_loss, En_Compt_bulk_gain, En_Compt_bulk_heat_gain
      
      
        implicit none
        
        real(kind=rk), intent(in) :: x, theta, phi, Z_bulk, kT, weight
        real(kind=rk), intent(out) :: x_sc, theta_sc, phi_sc, x_sc_pr, z_val ! , x_pr
        real(kind=rk) :: rnd2, rnd3, rnd4, rnd5, rnd6
        real(kind=rk) :: cos_theta_pr, sin_theta_pr, mu_rel_val, gamma_e, beta_e, Doppler_pr_erf, x_erf, x_sc_erf ! , x_sc_pr
        real(kind=rk) :: cos_sc_erf, sin_sc_erf, phi_sc_erf, cos_pk_erf, sin_pk_erf, cos_pk_sc_erf, mu_sc
        real(kind=rk) :: cos_kk_pr, sin_kk_pr, phi_sc_pr, cos_theta_sc_pr, sin_theta_sc_pr, cos_phi_kk_sc, phi_kk_sc, dummy
        real(kind=rk) :: mu, Gamma_bulk, Beta_bulk, Doppler, x_pr
        integer :: i_dummy, i_x, i_kT, i_mu_val, i_x_erf_val
        integer :: i_mu_sc_erf				!! for testing only
  !       integer, parameter :: tol = 1.0e-12_rk        !!! NBNBNBNBNB: ERROR: corrected 17.12.18
        real(kind=rk), parameter :: tol = 1.0e-7_rk ! tol = 1.0e-8_rk ! tol = 1.0e-12_rk   !!! relaxed tolerance a bit on 18.01.19  !!! relaxed to 1e-7 on 08.05.19
        logical :: sw_crash
        
        mu = cos(theta)
        Gamma_bulk = (1.0_rk + Z_bulk**2.0_rk)**0.5_rk
        Beta_bulk = Z_bulk/Gamma_bulk
        Doppler = 1.0_rk/(Gamma_bulk - Z_bulk*mu)
        x_pr = x/Doppler
        
        !!!#### Angle of the incident photon in the flow frame ########
        cos_theta_pr = (mu - Beta_bulk)/(1.0_rk - Beta_bulk*mu)
        sin_theta_pr = (1.0_rk - cos_theta_pr**2.0_rk)**0.5_rk
              
        !!!##### Get index of the flow frame energy and temperature in the flow ##########################
        call Findvalue_1dim_v3_nearestindex(x_pr,MC_Compton_grids%x_pr_vec(:),MC_Compton_grids%x_pr_vec(:),dummy,i_x)     !!! If x_pr > MC_Compton_grids%x_pr_vec(n_x) use the last gridpoint! Dangerous??? Is in deep KN though
        call Findvalue_1dim_v3_nearestindex(kT,MC_Compton_grids%kT_vec(:),MC_Compton_grids%kT_vec(:),dummy,i_kT)
        !!!######################################################################################
        
        !!!##### Draw electron angle relative to the photon direction in the flow frame #########
        call RANDOM_NUMBER(rnd2)
        call Findvalue_1dim_v3_nearestindex(rnd2,P_cumul_MC_Compton%P_mu_pr_e(i_x,i_kT,:),MC_Compton_grids%mu_rel_pr_vec(:),mu_rel_val,i_mu_val)
        !!!######################################################################################
        
        !!!################# Draw electron energy in the flow frame #############################
        call RANDOM_NUMBER(rnd3)
  !       call Findvalue_1dim_v3(rnd3,P_cumul_MC_Compton%P_z(i_x,i_kT,i_mu_val,:),MC_Compton_grids%z_vec(:),z_val,i_dummy)
        call Findvalue_1dim_v3(rnd3,P_cumul_MC_Compton%P_z(:,i_mu_val,i_kT,i_x),MC_Compton_grids%z_vec(:),z_val,i_dummy)   !! inverted (20.10.20)
        !!!######################################################################################
        
        gamma_e = (z_val**2.0_rk + 1.0_rk)**0.5_rk
        beta_e = z_val/gamma_e
        Doppler_pr_erf = 1.0_rk/(gamma_e*(1.0_rk - beta_e*mu_rel_val))
        x_erf = x_pr/Doppler_pr_erf
        
        !!!!######### Get scattering angle in the ERF: cos_sc_erf, phi_sc_erf, x_sc_erf ##########
        call Findvalue_1dim_v3_nearestindex(x_erf,MC_Compton_grids%x_erf_vec(:),MC_Compton_grids%x_erf_vec(:),dummy,i_x_erf_val)
        call RANDOM_NUMBER(rnd4)
        

        call Findvalue_1dim_v3(rnd4,P_cumul_MC_Compton%P_mu_erf_ph(i_x_erf_val,:),MC_Compton_grids%mu_sc_erf_vec(:),cos_sc_erf,i_dummy)
        
        if (1.eq.0) then   !!! TEST with isotropic scattering; WRONG!!!!!!!
          cos_sc_erf = 2.0*rnd4 - 1.0_rk
        end if
      
        x_sc_erf = x_erf/(1.0_rk + x_erf*(1.0_rk - cos_sc_erf))
        
  !   !     x_sc_erf = x_erf			!!! NBNBNBNB: WRONG, coherent in the ERF; COMMENT OUT
        
        call RANDOM_NUMBER(rnd5)		!!! drawing phi_sc (azimuth of the scattering angle rel. to el.directon around k)
        phi_sc_erf = 2.0_rk*pi*rnd5
        !!!######################################################################################
        
    !     z_val, Z_bulk, mu_rel_val, cos_sc_erf, phi_sc_erf, x_sc_erf
        
        sin_sc_erf = (1.0_rk - cos_sc_erf**2.0_rk)**0.5_rk
        cos_pk_erf = (mu_rel_val - beta_e)/(1.0_rk - beta_e*mu_rel_val)
        sin_pk_erf = (1.0_rk - cos_pk_erf**2.0_rk)**0.5_rk
        
        !!!## angle between the scattered photon and the electron in the ERF ##
        cos_pk_sc_erf = cos_pk_erf*cos_sc_erf + sin_pk_erf*sin_sc_erf*cos(phi_sc_erf)
        
        !!!## angle between the scattered and incident photons in the flow frame ##
        cos_kk_pr = 1.0_rk - (1.0_rk - cos_sc_erf)*(1.0_rk - beta_e*mu_rel_val)/(1.0_rk + beta_e*cos_pk_sc_erf)
        sin_kk_pr = (1.0_rk - cos_kk_pr**2.0_rk)**0.5_rk
        
        !!!######## Assuming isotropic electrons in the flow frame, draw from a uniform distribution the azimuth ############
        !!!######## of the scattered photon relative to the flow firection around the direction of the incident photon ######
        call RANDOM_NUMBER(rnd6)
        phi_sc_pr = 2.0_rk*pi*rnd6
        !!!##################################################################################################################
        
        !!!## angle of the scattered photon relative to the flow direction in the flow frame (cos_theta_sc_pr,phi_sc) ############
        cos_theta_sc_pr = cos_theta_pr*cos_kk_pr + sin_theta_pr*sin_kk_pr*cos(phi_sc_pr)
        sin_theta_sc_pr = (1.0_rk - cos_theta_sc_pr**2.0_rk)**0.5_rk
        if (abs(sin_theta_pr).eq.0.0_rk) then
          cos_phi_kk_sc = -cos(phi_sc_pr)	!!!### Relative azimuth of [k_pr,k_sc_pr] around Gamma_bulk is random (TAKEN as minus the azimuth of k_sc_pr around k_pr) IF momentum of k parallel to Gamma_bulk ###
        else if (abs(sin_theta_sc_pr).eq.0.0_rk) then
          cos_phi_kk_sc = 1.0_rk			!!!#### Value irrelevant if momentum of k_sc_pr parallel to Gamma_bulk ####
        else
          cos_phi_kk_sc = (cos_kk_pr - cos_theta_pr*cos_theta_sc_pr)/(sin_theta_pr*sin_theta_sc_pr)	! relative azimuth of the photons around Gamma_bulk
    !       print *,'abs(sin_theta_pr)', abs(sin_theta_pr)
        end if
        
        sw_crash = .FALSE.
        if (cos_phi_kk_sc.gt.1.0_rk) then
          if (cos_phi_kk_sc.lt.(1.0_rk + tol)) then
            cos_phi_kk_sc = 1.0_rk
            phi_kk_sc = 0.0_rk
          else
            print *,'Scatter_photon ERROR: cos_phi_kk_sc > 1.0', cos_phi_kk_sc
            sw_crash = .TRUE.
          end if
        else if (cos_phi_kk_sc.lt.-1.0_rk) then
          if (cos_phi_kk_sc.gt.(-1.0_rk - tol)) then
            cos_phi_kk_sc = -1.0_rk
            phi_kk_sc = pi
          else
            print *,'Scatter_photon ERROR: cos_phi_kk_sc < -1.0', cos_phi_kk_sc
            print *,'tol', tol, -1.0_rk - tol, cos_phi_kk_sc.gt.(-1.0_rk - tol)
            sw_crash = .TRUE.
          end if
        else
          phi_kk_sc = acos(cos_phi_kk_sc)		!! gives a value between [0,pi]
        end if
        
  !       phi_kk_sc = acos(cos_phi_kk_sc)		!! gives a value between [0,pi]
        if (phi_sc_pr.lt.pi) then
          phi_sc = phi - phi_kk_sc
        else
          phi_sc = phi + phi_kk_sc
        end if
    !     print *,'phi, phi_kk_sc, cos_phi_kk_sc', phi, phi_kk_sc, cos_phi_kk_sc
    !     print *,'a, b', (cos_kk_pr - cos_theta_pr*cos_theta_sc_pr), (sin_theta_pr*sin_theta_sc_pr)
    !     pause
        !!!#######################################################################################################################
        
        !!!##### scattered photon energy and direction in the lab frame #######
          x_sc_pr = x_sc_erf*(gamma_e + z_val*cos_pk_sc_erf) 
          x_sc = x_sc_pr*(Gamma_bulk + Z_bulk*cos_theta_sc_pr)
          mu_sc = (cos_theta_sc_pr + Beta_bulk)/(1.0_rk + Beta_bulk*cos_theta_sc_pr)
          
          if (mu_sc.gt.1.0_rk) then
            if (mu_sc.lt.(1.0_rk + tol)) then
              mu_sc = 1.0_rk
              theta_sc = 0.0_rk
            else
              print *,'Scatter_photon ERROR: mu_sc > 1.0', mu_sc
              sw_crash = .TRUE.
            end if
          else if (mu_sc.lt.-1.0_rk) then
            if (mu_sc.gt.(-1.0_rk - tol)) then
              mu_sc = -1.0_rk
              theta_sc = pi
            else
              print *,'Scatter_photon ERROR: mu_sc < -1.0', mu_sc
              sw_crash = .TRUE.
            end if
          else
            theta_sc = acos(mu_sc)		!! gives a value between [0,pi]
          end if

  !       theta_sc = acos(mu_sc)
        !!!####################################################################
        
        
        !!!####################### Construction (06.01.22) ###########################################################
          !!!################### Energy deposition into the flow ################################
            En_Compt_ph_loss = En_Compt_ph_loss + weight*(x - x_sc)
            En_Compt_bulk_gain = En_Compt_bulk_gain + weight*(x_pr*cos_theta_pr - x_sc_pr*cos_theta_sc_pr)*Z_bulk
            En_Compt_bulk_heat_gain = En_Compt_bulk_heat_gain + weight*(x_pr - x_sc_pr)*Gamma_bulk         !!! Thermal energy gain by flow, transformed into the lab frame
          !!!####################################################################################
        !!!###########################################################################################################
        
        !!!! OUT: x_sc, x_sc_pr, theta_sc, phi_sc
        
  !       !!!#### Test quantities: not used by the code itself #######
  !       !!!#### NBNBNBNBNB: probably substantially wastes time #####
  !       call Findvalue_1dim_v3_nearestindex(cos_theta_sc_pr,Test_qties_MC%mu_grid(:),Test_qties_MC%mu_grid(:),dummy,i_mu_sc_erf)
  !       Test_qties_MC%N_mu_sc_erf(i_mu_sc_erf) = Test_qties_MC%N_mu_sc_erf(i_mu_sc_erf) + 1.0_rk
  !       !!!#########################################################

  !       print *,'x, x_sc, gamma_e-1, kT', x, x_sc, gamma_e-1, kT


        if ((phi_sc.ge.0.0_rk).or.(phi_sc.lt.0.0_rk)) then
          Continue
        else
          print *,'############ Scatter_photon ERROR: phi_sc ###############', phi_sc
          sw_crash = .TRUE.
        end if
        
        if (sw_crash) then
          print *,'Gamma_bulk, Z_bulk', Gamma_bulk, Z_bulk
          print *,'mu, Doppler', mu, Doppler
          
          print *,'z_val, gamma_e, beta_e', z_val, gamma_e, beta_e
          print *,'Doppler_pr_erf', Doppler_pr_erf
          print *,'x_erf', x_erf
          print *,'x_sc_erf', x_sc_erf
          
          print *,'x, x_pr', x, x_pr
          print *,'x_sc, x_sc_pr', x_sc, x_sc_pr
          print *,'kT', kT
          
          print *,'mu_rel_val', mu_rel_val
          print *,'cos_pk_erf', cos_pk_erf
          
          print *,'cos_sc_erf, sin_sc_erf', cos_sc_erf, sin_sc_erf
          print *,'cos_pk_sc_erf', cos_pk_sc_erf
          print *,'phi_sc_erf', phi_sc_erf
          
          print *,'cos_theta_sc_pr, sin_theta_sc_pr', cos_theta_sc_pr, sin_theta_sc_pr
          print *,'cos_theta_pr, sin_theta_pr', cos_theta_pr, sin_theta_pr
          print *,'cos_kk_pr, sin_kk_pr', cos_kk_pr, sin_kk_pr
          print *,'phi_sc_pr', phi_sc_pr
          
          print *,'cos_phi_kk_sc, phi_kk_sc', cos_phi_kk_sc, phi_kk_sc
          
          print *,'phi', phi
          print *,'##########################################################'
          
          stop
          
        end if

      End Subroutine Scatter_photon
      
!!!################################################## End: scattering routines ##################################################################################################   
!!!##############################################################################################################################################################################       

         
    
!!!##############################################################################################################################################################################     
!!!################################################## Photon emission from thermal plasma #######################################################################################       
      
! ! !       POOOOOOLELI: should this be part of the Ejecta_obj object?
      
!       Subroutine Emit_thermal(d_ct)
!         implicit none
!         
!         real(kind=rk), intent(in) :: d_ct
!       
!       End Subroutine Emit_thermal
      

        
  
  
!!!############################################### End: photon emission from thermal plasma #####################################################################################   
!!!##############################################################################################################################################################################     
 
End Module Phys_procs
