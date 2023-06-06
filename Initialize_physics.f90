Module Initialize_physics

  use realkind
  use constants
  
  implicit none
  
    Contains


!!!################################################## Precalculation: Photon interactions ##############################################################################


!!!#####################################################################################################################################################      
!!!################################### Compton scattering ##############################################################################################
!!!#####################################################################################################################################################  
        
    Subroutine Setup_grids_Compton_MC
      
      !! Routine sets up x_pr_vec, z_vec, mu_rel_pr_vec, kT_vec (in grids)
      
!       use grids_MC_Compton, ONLY: x_pr_min, x_pr_max, z_min, z_max, mu_rel_min, mu_rel_max, kT_min, kT_max, x_erf_min, x_erf_max, mu_sc_erf_min, mu_sc_erf_max,   &
!                                   n_z, n_x, n_mu, n_x_erf, n_mu_sc_erf, n_kT
      use Cumulative_distributions, ONLY: MC_Compton_grids,     &  !! Fills MC_Compton_grids% [x_pr_vec, z_vec, mu_rel_pr_vec, x_erf_vec, mu_sc_erf_vec, kT_vec]
                                          x_pr_min, x_pr_max, z_min, z_max, mu_rel_min, mu_rel_max, kT_min, kT_max, x_erf_min, x_erf_max, mu_sc_erf_min, mu_sc_erf_max,   &
                                          n_z, n_x, n_mu, n_x_erf, n_mu_sc_erf, n_kT
      use auxiliary
      
      real(kind=rk) :: dlnx_pr, lnx_pr_min, lnx_pr_max, dlnz_pr, lnz_min, lnz_max, d_mu_rel, ln_kT_min, ln_kT_max, d_ln_kT, lnx_erf_min, lnx_erf_max, d_ln_x_erf, d_mu_sc_erf
      real(kind=rk), dimension(:), allocatable :: lnx_pr_vec, lnz_vec, ln_kT_vec, lnx_erf_vec
      integer :: i
    
      if (allocated(MC_Compton_grids%x_pr_vec)) deallocate(MC_Compton_grids%x_pr_vec)
      if (allocated(MC_Compton_grids%z_vec)) deallocate(MC_Compton_grids%z_vec)
      if (allocated(MC_Compton_grids%mu_rel_pr_vec)) deallocate(MC_Compton_grids%mu_rel_pr_vec)
      if (allocated(MC_Compton_grids%kT_vec)) deallocate(MC_Compton_grids%kT_vec)
      if (allocated(MC_Compton_grids%x_erf_vec)) deallocate(MC_Compton_grids%x_erf_vec)
      if (allocated(MC_Compton_grids%mu_sc_erf_vec)) deallocate(MC_Compton_grids%mu_sc_erf_vec)
      if (allocated(lnx_pr_vec)) deallocate(lnx_pr_vec)
      if (allocated(lnz_vec)) deallocate(lnz_vec)
      if (allocated(ln_kT_vec)) deallocate(ln_kT_vec)
      if (allocated(lnx_erf_vec)) deallocate(lnx_erf_vec)
      allocate(MC_Compton_grids%x_pr_vec(n_x),MC_Compton_grids%z_vec(n_z),MC_Compton_grids%mu_rel_pr_vec(n_mu),MC_Compton_grids%kT_vec(n_kT),MC_Compton_grids%x_erf_vec(n_x_erf),MC_Compton_grids%mu_sc_erf_vec(n_mu_sc_erf))
      allocate(lnx_pr_vec(n_x),lnz_vec(n_z),ln_kT_vec(n_kT),lnx_erf_vec(n_x_erf))
    
      lnx_pr_min = log(x_pr_min)
      lnx_pr_max = log(x_pr_max)
      call Setup_grid(lnx_pr_vec,dlnx_pr,lnx_pr_min,lnx_pr_max)
      do i = 1,n_x
        MC_Compton_grids%x_pr_vec(i) = exp(lnx_pr_vec(i))
      end do
      
      lnz_min = log(z_min)
      lnz_max = log(z_max)
      call Setup_grid(lnz_vec,dlnz_pr,lnz_min,lnz_max)
      do i = 1,n_z
        MC_Compton_grids%z_vec(i) = exp(lnz_vec(i))
      end do
	  
      ln_kT_min = log(kT_min)
      ln_kT_max = log(kT_max)
      call Setup_grid(ln_kT_vec,d_ln_kT,ln_kT_min,ln_kT_max)
      do i = 1,n_kT
        MC_Compton_grids%kT_vec(i) = exp(ln_kT_vec(i))
      end do
      
      call Setup_grid(MC_Compton_grids%mu_rel_pr_vec,d_mu_rel,mu_rel_min,mu_rel_max)
      
      !!!#### For scattering in the electron rest frame ####################
      lnx_erf_min = log(x_erf_min)
      lnx_erf_max = log(x_erf_max)
      call Setup_grid(lnx_erf_vec,d_ln_x_erf,lnx_erf_min,lnx_erf_max)
      do i = 1,n_x_erf
        MC_Compton_grids%x_erf_vec(i) = exp(lnx_erf_vec(i))
      end do  
      
      call Setup_grid(MC_Compton_grids%mu_sc_erf_vec,d_mu_sc_erf,mu_sc_erf_min,mu_sc_erf_max)
      !!!###################################################################
      
      deallocate(lnx_pr_vec,lnz_vec,ln_kT_vec,lnx_erf_vec)
    
    End Subroutine Setup_grids_Compton_MC
    
    
    
    Subroutine KN_cross_tot(x_erf,sigma_tot)
  
      !!! total erf KN cross-section, normalized to sigma_T
    
      implicit none
      
      real(kind=rk), intent(in) :: x_erf
      real(kind=rk), intent(out) :: sigma_tot
      real(kind=rk) :: i_real, eps_chk, term
    
        !!!######## Klein-nishina cross_section #####################
          if (x_erf.gt.1.0e-1_rk) then     
            sigma_tot = 3.0_rk/(8.0_rk*x_erf)*((1.0_rk - 2.0_rk/x_erf - 2.0_rk/x_erf**2.0_rk)*log(1.0_rk+2.0_rk*x_erf) + 0.5_rk + 4.0_rk/x_erf - 0.5_rk/(1.0_rk+2.0_rk*x_erf)**2.0_rk)
          else
            sigma_tot = 0.0_rk
            i_real = 0.0_rk
            eps_chk = 1.0_rk
            do while (eps_chk.gt.1.0e-14_rk)
              i_real = i_real + 1.0_rk
              term = (-2.0_rk*x_erf)**(i_real-1.0_rk)*(1.0_rk/i_real + 4.0_rk/(i_real+1.0_rk) - 8.0_rk/(i_real+2.0_rk) + (i_real+1.0_rk)/2.0_rk)
              sigma_tot = sigma_tot + term
              eps_chk = abs(term/sigma_tot)
            end do
            sigma_tot = sigma_tot*3.0_rk/4.0_rk
          end if
        !!!###########################################################
    
    End Subroutine KN_cross_tot


    
    Subroutine P_cumul_distributions_Compton_MC
    
!       use grids_MC_Compton          !!! uses MC_Compton_grids, n_z, n_x, n_mu, n_x_erf, n_mu_sc_erf, n_kT
      use Cumulative_distributions, ONLY: P_cumul_MC_Compton,       &                         !!! P_cumul_MC_Compton% [kappa_pr_norm, P_mu_erf_ph, P_mu_pr_e, P_z]
                                          MC_Compton_grids, n_z, n_x, n_mu, n_x_erf, n_mu_sc_erf, n_kT
      use JCMaxwell, ONLY: f_MW_norm
      use auxiliary, ONLY: Findvalue_1dim_v3
    
      implicit none
      
      real(kind=rk), dimension(:,:), allocatable :: f_e_norm
      real(kind=rk), dimension(:), allocatable :: sigma_tot_vec, x_erf_loc_vec
      real(kind=rk) :: d_lnz, kT, I_norm, z, x_erf_loc_min, x_erf_loc_max, lnx_erf_min, lnx_erf_max, x_erf, sigma_tot
      real(kind=rk) :: x_pr, Int_z_mu, Int_z, mu_rel_pr, d_mu_rel, gamma, beta, Doppler_pr_erf, Integrand
      real(kind=rk) :: Int_mu_erf, mu_sc_erf, one_pl_mu, sin_sq, d_mu_erf, x_sc_ov_x, sigma_mu
      integer :: i, i_z, i_x, i_mu, i_kT, n_x_erf_loc, i_mu_sc_erf, i_x_erf, i_dummy
      
        print *,'P_cumul_distributions_Compton_MC: START'
      
      !!! taken from module grids
  !     n_z = size(z_vec)  
  !     n_x = size(x_pr_vec)  
  !     n_mu = size(mu_rel_pr_vec)		!!! Electron direction RELATIVE TO THE PHOTON in the flow frame 
  !     n_kT = size(theta_vec)
      
        if (allocated(P_cumul_MC_Compton%P_z)) deallocate(P_cumul_MC_Compton%P_z)
        if (allocated(P_cumul_MC_Compton%P_mu_pr_e)) deallocate(P_cumul_MC_Compton%P_mu_pr_e)
        if (allocated(P_cumul_MC_Compton%kappa_pr_norm)) deallocate(P_cumul_MC_Compton%kappa_pr_norm)
        if (allocated(P_cumul_MC_Compton%P_mu_erf_ph)) deallocate(P_cumul_MC_Compton%P_mu_erf_ph)
!         allocate(P_cumul_MC_Compton%P_z(n_x,n_kT,n_mu,n_z))
        allocate(P_cumul_MC_Compton%P_z(n_z,n_mu,n_kT,n_x))   !! inverted (20.10.20)
        allocate(P_cumul_MC_Compton%P_mu_pr_e(n_x,n_kT,n_mu))
        allocate(P_cumul_MC_Compton%kappa_pr_norm(n_x,n_kT))
        allocate(P_cumul_MC_Compton%P_mu_erf_ph(n_x_erf,n_mu_sc_erf))
        
        d_lnz = (log(MC_Compton_grids%z_vec(n_z)) - log(MC_Compton_grids%z_vec(1)))/(dble(n_z) - 1.0_rk)
        
        if (allocated(f_e_norm)) deallocate(f_e_norm)
        allocate(f_e_norm(n_z,n_kT))
        do i_kT = 1,n_kT
          kT = MC_Compton_grids%kT_vec(i_kT)
          I_norm = 0.0_rk
          do i_z = 1,n_z
        ! 	z = exp(lnz_vec(i_z))
            z = MC_Compton_grids%z_vec(i_z)
            f_e_norm(i_z,i_kT) = f_MW_norm(z,kT)
            I_norm = I_norm + f_e_norm(i_z,i_kT)*d_lnz
          end do
          print *,'f_e norm should be 1', I_norm
        end do  
    !     pause
        
        
        x_erf_loc_min = 1.0e-6_rk	! 1.0e-5_rk in LAT problem
        x_erf_loc_max = 1.0e6 		! 1.0e10_rk in LAT problem
        lnx_erf_min = log(x_erf_loc_min)
        lnx_erf_max = log(x_erf_loc_max)
        n_x_erf_loc = 2000
        if (allocated(sigma_tot_vec)) deallocate(sigma_tot_vec)
        if (allocated(x_erf_loc_vec)) deallocate(x_erf_loc_vec)
        allocate(sigma_tot_vec(n_x_erf_loc),x_erf_loc_vec(n_x_erf_loc))
        do i = 1,n_x_erf_loc
          x_erf = exp(lnx_erf_min + (i-1.0_rk)*(lnx_erf_max-lnx_erf_min)/(n_x_erf_loc-1.0_rk))
          call KN_cross_tot(x_erf,sigma_tot)
          x_erf_loc_vec(i) = x_erf
          sigma_tot_vec(i) = sigma_tot
  ! 	     print *,'sigma_tot', sigma_tot
        end do
        
        print *,'Precalc. MC_Compton_grids'    
        do i_x = 1,n_x
          if (i_x/100*100.eq.i_x) print *,'Precalc. MC_Compton_grids: i_x', i_x
          x_pr = MC_Compton_grids%x_pr_vec(i_x)
          do i_kT = 1,n_kT
            kT = MC_Compton_grids%kT_vec(i_kT)
            Int_z_mu = 0.0_rk
            P_cumul_MC_Compton%P_mu_pr_e(i_x,i_kT,1) = 0.0_rk
            do i_mu = 1,n_mu
              mu_rel_pr = MC_Compton_grids%mu_rel_pr_vec(i_mu)
              if (i_mu.eq.1) then
                d_mu_rel = (MC_Compton_grids%mu_rel_pr_vec(i_mu+1) - MC_Compton_grids%mu_rel_pr_vec(i_mu))/2.0_rk
              else if (i_mu.eq.n_mu) then
                d_mu_rel = (MC_Compton_grids%mu_rel_pr_vec(i_mu) - MC_Compton_grids%mu_rel_pr_vec(i_mu-1))/2.0_rk
              else
                d_mu_rel = (MC_Compton_grids%mu_rel_pr_vec(i_mu+1) - MC_Compton_grids%mu_rel_pr_vec(i_mu-1))/2.0_rk		! not sure if best for non-uniform
              end if
              Int_z = 0.0_rk
!               P_cumul_MC_Compton%P_z(i_x,i_kT,i_mu,1) = 0.0_rk
              P_cumul_MC_Compton%P_z(1,i_mu,i_kT,i_x) = 0.0_rk   !! inverted (20.10.20)
              
              do i_z = 1,n_z
                z = MC_Compton_grids%z_vec(i_z)				!! Shouldnt z be at halfstep???
                gamma = (z*z + 1.0_rk)**0.5_rk
                beta = z/gamma
                Doppler_pr_erf = 1.0_rk/(gamma*(1.0_rk - beta*mu_rel_pr))
                x_erf = x_pr/Doppler_pr_erf
                call Findvalue_1dim_v3(x_erf,x_erf_loc_vec,sigma_tot_vec,sigma_tot,i_dummy)
                Int_z = Int_z + f_e_norm(i_z,i_kT)*sigma_tot*(1.0_rk-beta*mu_rel_pr)
        ! 	    if (i_z.gt.1) P_cumul%P_z(i_x,i_kT,i_mu,i_z) = P_cumul%P_z(i_x,i_kT,i_mu,i_z-1) + f_e_norm(i_z,i_kT)*sigma_tot	!! WHY NO (1.0_rk-beta*mu_rel_pr)??
!                 if (i_z.gt.1) P_cumul_MC_Compton%P_z(i_x,i_kT,i_mu,i_z) = P_cumul_MC_Compton%P_z(i_x,i_kT,i_mu,i_z-1) + f_e_norm(i_z,i_kT)*sigma_tot	*(1.0_rk-beta*mu_rel_pr) !! 08.09.16
                  if (i_z.gt.1) P_cumul_MC_Compton%P_z(i_z,i_mu,i_kT,i_x) = P_cumul_MC_Compton%P_z(i_z-1,i_mu,i_kT,i_x) + f_e_norm(i_z,i_kT)*sigma_tot	*(1.0_rk-beta*mu_rel_pr) !! 08.09.16   !! inverted (20.10.20)
              end do
              
!               P_cumul_MC_Compton%P_z(i_x,i_kT,i_mu,:) = P_cumul_MC_Compton%P_z(i_x,i_kT,i_mu,:)/P_cumul_MC_Compton%P_z(i_x,i_kT,i_mu,n_z)
              P_cumul_MC_Compton%P_z(:,i_mu,i_kT,i_x) = P_cumul_MC_Compton%P_z(:,i_mu,i_kT,i_x)/P_cumul_MC_Compton%P_z(n_z,i_mu,i_kT,i_x)   !! inverted (20.10.20)
              Integrand = Int_z*d_lnz*d_mu_rel/2.0_rk   !!! /2.0 because averaging over mu_rel (since f_e_norm has already effectively been integrated over solid angle)
              Int_z_mu = Int_z_mu + Integrand
              if (i_mu.eq.n_mu) then
                P_cumul_MC_Compton%P_mu_pr_e(i_x,i_kT,i_mu) = Int_z_mu
              else if (i_mu.gt.1) then
                P_cumul_MC_Compton%P_mu_pr_e(i_x,i_kT,i_mu) = Int_z_mu - Integrand/2.0_rk
              end if
            end do
            
            P_cumul_MC_Compton%kappa_pr_norm(i_x,i_kT) = Int_z_mu
            P_cumul_MC_Compton%P_mu_pr_e(i_x,i_kT,:) = P_cumul_MC_Compton%P_mu_pr_e(i_x,i_kT,:)/P_cumul_MC_Compton%P_mu_pr_e(i_x,i_kT,n_mu)
          end do   
        end do
            
        do i_x_erf = 1,n_x_erf
          x_erf = MC_Compton_grids%x_erf_vec(i_x_erf)
          Int_mu_erf = 0.0_rk
          P_cumul_MC_Compton%P_mu_erf_ph(i_x_erf,1) = 0.0_rk
          do i_mu_sc_erf = 1,n_mu_sc_erf
            mu_sc_erf = MC_Compton_grids%mu_sc_erf_vec(i_mu_sc_erf)
            one_pl_mu = 1.0_rk + mu_sc_erf
            sin_sq = one_pl_mu*(1.0_rk - mu_sc_erf) 
            
            if (i_mu_sc_erf.eq.1) then
              d_mu_erf = (MC_Compton_grids%mu_sc_erf_vec(i_mu_sc_erf+1) - MC_Compton_grids%mu_sc_erf_vec(i_mu_sc_erf))/2.0_rk
            else if (i_mu_sc_erf.eq.n_mu_sc_erf) then
              d_mu_erf = (MC_Compton_grids%mu_sc_erf_vec(i_mu_sc_erf) - MC_Compton_grids%mu_sc_erf_vec(i_mu_sc_erf-1))/2.0_rk
            else
              d_mu_erf = (MC_Compton_grids%mu_sc_erf_vec(i_mu_sc_erf+1) - MC_Compton_grids%mu_sc_erf_vec(i_mu_sc_erf-1))/2.0_rk		! not sure if best for non-uniform
            end if
            
            x_sc_ov_x = 1.0_rk/(1.0_rk + x_erf*(1.0_rk - mu_sc_erf))
            sigma_mu = 3.0_rk/8.0_rk*x_sc_ov_x**2.0_rk*(x_sc_ov_x + 1.0_rk/x_sc_ov_x - sin_sq) 
                
            Int_mu_erf = Int_mu_erf + sigma_mu*d_mu_erf
            
            if (i_mu_sc_erf.eq.n_mu_sc_erf) then
              P_cumul_MC_Compton%P_mu_erf_ph(i_x_erf,i_mu_sc_erf) = Int_mu_erf
            else if (i_mu_sc_erf.gt.1) then 
              P_cumul_MC_Compton%P_mu_erf_ph(i_x_erf,i_mu_sc_erf) = Int_mu_erf - sigma_mu*d_mu_erf/2.0_rk
            end if
            
          end do
          P_cumul_MC_Compton%P_mu_erf_ph(i_x_erf,:) = P_cumul_MC_Compton%P_mu_erf_ph(i_x_erf,:)/P_cumul_MC_Compton%P_mu_erf_ph(i_x_erf,n_mu_sc_erf)
        end do
        
        deallocate(sigma_tot_vec,x_erf_loc_vec,f_e_norm)
        
        print *,'P_cumul_distributions_Compton_MC: DONE'
    
    End Subroutine P_cumul_distributions_Compton_MC

       
!!!################################### End: Compton scattering #########################################################################################


    Subroutine P_cumul_Thermal_brems_MC
    
      use Expintegral
      use Cumulative_distributions, ONLY: P_cumul_Therm_brems, y_min_brems, y_max_brems, n_y_brems
      
      use auxiliary, ONLY: Setup_grid
    
      implicit none
      
      real(kind=rk) :: d_lny, lny, y
      real(kind=rk) :: E_1_ymin, E_1_ymax, E_1
      integer :: i
        
        allocate(P_cumul_Therm_brems%lny(n_y_brems),P_cumul_Therm_brems%P_y(n_y_brems),P_cumul_Therm_brems%P_y_unnorm(n_y_brems))
        P_cumul_Therm_brems%P_y(:) = 0.0_rk
        P_cumul_Therm_brems%P_y_unnorm(:) = 0.0_rk
        call Setup_grid(P_cumul_Therm_brems%lny,d_lny,log(y_min_brems),log(y_max_brems))
        
        E_1_ymin = expint(1,y_min_brems)
        E_1_ymax = expint(1,y_max_brems)
        
        P_cumul_Therm_brems%Ndot_norm = E_1_ymin - E_1_ymax   !! Normalized total photon emission rate in the range [y_min_brems,y_max_brems]
        
        do i = 2,n_y_brems
          lny = P_cumul_Therm_brems%lny(i)
          y = exp(lny)
          E_1 = expint(1,y)
!           P_cumul_Therm_brems%P_y(i) = (E_1_ymin - E_1)/(E_1_ymin - E_1_ymax)
          P_cumul_Therm_brems%P_y_unnorm(i) = (E_1_ymin - E_1)
        end do  !!! do i = 1,n_y_brems
        
        P_cumul_Therm_brems%P_y(:) = P_cumul_Therm_brems%P_y_unnorm(:)/P_cumul_Therm_brems%Ndot_norm
        
        print *,'P_cumul_Therm_brems%P_y(:)', P_cumul_Therm_brems%P_y(:)
        print *,'Brems: 1 ==', P_cumul_Therm_brems%P_y(n_y_brems)
        
    End Subroutine P_cumul_Thermal_brems_MC








!!!################################################## End: Precalculation: Photon interactions #########################################################################


















End Module Initialize_physics
