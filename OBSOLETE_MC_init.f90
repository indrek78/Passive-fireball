! Module MC_init
! 
!   use realkind
!   use constants
!   
!   implicit none
! 
!     Contains
! 
! 
! !!!############################################# Initialize MC particle distributions #####################################################
! 
! 
!     Subroutine Initialize_Radiation()
!     
!       use Ejecta_choice_param, ONLY: i_choose_ejecta
!     
!       implicit none
!     
!         if (i_choose_ejecta.eq.1) then
!           call Initialize_Radiation_Passive_fireball_Planck()
!         else
!           stop 'Initialize_Radiation: cannot handle i_choose_ejecta'
!         end if
!             
!     End Subroutine Initialize_Radiation
!     
!     
!     
!     Subroutine Initialize_Radiation_Passive_fireball_Planck()
!       !!! Initialization specific to passife fireball, using type(Sphere_homol) :: Sphere_object
!     
!       use Photon_distribution, ONLY: Photons_vec
!       use Photon_distribution_escape, ONLY: Photons_escape, i_exist_max_escape
!       
! ! ! !       use Passive_fireball_params, ONLY: R_sph_0
!       use Ejecta_generic, ONLY: Ejecta_obj      !!! Only used here to obtain Ejecta_obj%i_choose_ejecta for a consistency check
!       use Ejecta_objects, ONLY: Sphere_object   !! Use specific object (passive fireball), which requires i_choose_ejecta = 1
!       
!       use Ejecta_choice_param, ONLY: i_choose_ejecta  !! Only used for double checking Ejecta_obj%i_choose_ejecta
!       
!       use MC_photon_params, ONLY: n_MC_ph
!       
!       use Ph_thermal, ONLY: f_Planck, Planck_moments
!       use auxiliary, ONLY: Setup_grid, Findvalue_1dim_v3
!       use Geometry, ONLY: get_v_vec_from_angles_in_r_system
! 
!       implicit none
!       
!         real(kind=rk) :: v_ej, xi_min, ct, d_Vol
!         real(kind=rk), dimension(:), allocatable :: kT_vec, xi_cell_bnd_vec
!         real(kind=rk) :: theta, n_ph, u_ph, N_ph_tot
!         real(kind=rk) :: x_min, x_max, lnx_min, lnx_max, d_lnx, x, lnx_val, mu, phi, mu_loc, phi_loc, r
!         real(kind=rk) :: r_cell_min, r_cell_max
!         real(kind=rk), dimension(:,:), allocatable :: P_cumul_r_x, lnx
!         real(kind=rk) :: w, rnd1, rnd2, rnd3, rnd4, rnd5, rnd6
!         real(kind=rk) :: n_MC_dV_real
!         real(kind=rk) :: R_sph_0  !!! NB: locally-defined (no longer rea from Module Passive_fireball_params)
!         integer :: n_r, i_r, n_x, i_x, i_MC, i_MC_dV, n_MC_dV, i_dummy
!       
!         if ((Ejecta_obj%i_choose_ejecta.ne.1).or.(i_choose_ejecta.ne.1)) then
!           print *,'Initialize_Radiation_Passive_fireball: Ejecta_obj%i_choose_ejecta.ne.1 or i_choose_ejecta.ne.1', Ejecta_obj%i_choose_ejecta, i_choose_ejecta
!           print *,'This routine is specific to passive fireball, which should have i_choose_ejecta = 1' 
!           stop
!         end if
!             
!         if (allocated(Photons_vec)) deallocate(Photons_vec)
!         
!         if (allocated(Photons_escape)) deallocate(Photons_escape)
!         i_exist_max_escape = 0
!         
!         n_r = Sphere_object%n_r
!         v_ej = Sphere_object%v_ej
!         xi_min = Sphere_object%xi_min
!         R_sph_0 = Sphere_object%R_sph_0
!         
!         allocate(kT_vec(n_r),xi_cell_bnd_vec(n_r+1))
!         kT_vec = Sphere_object%kT_vec
!         xi_cell_bnd_vec = Sphere_object%xi_cell_bnd_vec
!         
!         ct = R_sph_0/v_ej
!         
!         n_x = 1000
!         allocate(lnx(n_r,n_x))   !! Different grid for each location (allowing for different temperatures)
!         allocate(P_cumul_r_x(n_r,n_x))
!         P_cumul_r_x(:,:) = 0.0_rk
!         
! !         Vol_tot = 4.0_rk*pi/3.0_rk*(1.0_rk - xi_min**3.0_rk)*R_sph_0**3.0_rk   !!! Total volume of the ejecta 
!         
!         N_ph_tot = 0.0_rk  !!! Total number of physical photons
!         do i_r = 1,n_r
! !           print *,'i_r', i_r
! !           d_r = (xi_cell_bnd_vec(i_r+1) - xi_cell_bnd_vec(i_r))*R_sph_0
!           d_Vol = 4.0_rk*pi/3.0_rk*(xi_cell_bnd_vec(i_r+1)**3.0_rk - xi_cell_bnd_vec(i_r)**3.0_rk)*R_sph_0**3.0_rk
!           
!           theta = kT_vec(i_r)
!           call Planck_moments(theta,n_ph,u_ph)
!           
!           !!!############## Grids for cumulative distr. in ph. energy #################
!             x_min = theta*1.0e-7_rk
!             x_max = theta*1.0e2_rk
!             lnx_min = log(x_min)
!             lnx_max = log(x_max)
!             call Setup_grid(lnx(i_r,:),d_lnx,lnx_min,lnx_max)
!           !!!#########################################################################
!           
!           
!           N_ph_tot = N_ph_tot + n_ph*d_Vol
!           
!           P_cumul_r_x(i_r,:) = 0.0_rk   ! redundant
!           do i_x = 2,n_x
!             x = exp((lnx(i_r,i_x) + lnx(i_r,i_x-1))/2.0_rk)
!             P_cumul_r_x(i_r,i_x) = P_cumul_r_x(i_r,i_x-1) + f_Planck(x,theta)*d_lnx
!             
! !             print *,'x, d_lnx, theta, f_Planck(x,theta)', x, d_lnx, theta, f_Planck(x,theta)
! !             print *, 'P_cumul_r_x(i_r,i_x)', P_cumul_r_x(i_r,i_x)
! !             stop
!             
!           end do  !!!### do i_x = 1,n_x
!           
!           P_cumul_r_x(i_r,:) = P_cumul_r_x(i_r,:)/P_cumul_r_x(i_r,n_x) 
!                               
!         end do  !!!### do i_r = 1,n_r
!         
!         print *,'N_ph_tot', N_ph_tot
!         
! !         print *,'P_cumul_r_x(1,:)', P_cumul_r_x(1,:)
! !         stop
!         
!         !!!############# Draw photons ##################
!         !!########### Strategy 2 (in notebook) ####
!         
!         allocate(Photons_vec(n_MC_ph))
!         Photons_vec(:)%exist = .FALSE.
!         
!         w = N_ph_tot/n_MC_ph
!         
!         i_MC = 0
!         do i_r = 1,n_r
! !           print *,'i_r', i_r
!           d_Vol = 4.0_rk*pi/3.0_rk*(xi_cell_bnd_vec(i_r+1)**3.0_rk - xi_cell_bnd_vec(i_r)**3.0_rk)*R_sph_0**3.0_rk
!           theta = kT_vec(i_r)
!           call Planck_moments(theta,n_ph,u_ph)
!           
!           w = N_ph_tot/n_MC_ph
!           
!           n_MC_dV_real = n_ph*d_Vol/w
!           n_MC_dV = n_MC_dV_real
!           w = w*n_MC_dV_real/n_MC_dV
!           
! !           print *,'n_MC_dV', n_MC_dV
!  
!           do i_MC_dV = 1,n_MC_dV
!             i_MC = i_MC + 1
!           
!             call RANDOM_NUMBER(rnd1)
!             call RANDOM_NUMBER(rnd2)
!             call RANDOM_NUMBER(rnd3)
!             
!             call Findvalue_1dim_v3(rnd1,P_cumul_r_x(i_r,:),lnx(i_r,:),lnx_val,i_dummy)
!             
! !             print *,'P_cumul_r_x(i_r,:)', P_cumul_r_x(i_r,:)
! !             print *,'lnx(i_r,:)', lnx(i_r,:)
! !             print *,'lnx_val', lnx_val
! !             print *,'rnd1', rnd1
! !             stop
!             
!             mu = 2.0_rk*rnd2 - 1.0_rk
!             phi = 2.0_rk*pi*rnd3
!             
!             call RANDOM_NUMBER(rnd6)
!             r_cell_min = xi_cell_bnd_vec(i_r)*R_sph_0
!             r_cell_max = xi_cell_bnd_vec(i_r+1)*R_sph_0
!             r = (rnd6*r_cell_max**3.0_rk + (1.0_rk - rnd6)*r_cell_min**3.0_rk)**(1.0_rk/3.0_rk)      !! uniform distribution in shperical geometry
! ! !             r = (rnd6*r_cell_max + (1.0_rk - rnd6)*r_cell_min)  !! WRONG, for testing
!             
!             Photons_vec(i_MC)%ct = ct
!             Photons_vec(i_MC)%r = r
!             Photons_vec(i_MC)%x = exp(lnx_val)
!             Photons_vec(i_MC)%theta = acos(mu)
!             Photons_vec(i_MC)%phi = phi
!             Photons_vec(i_MC)%weight = w
!             
!             Photons_vec(i_MC)%ct_aux = 0.0_rk !!! unused
!             
!             Photons_vec(i_MC)%id = i_MC   !!! Not unique
!             Photons_vec(i_MC)%exist = .TRUE.
!             
!             !!!#### Draw location ####
!               call RANDOM_NUMBER(rnd4)
!               call RANDOM_NUMBER(rnd5)
!               mu_loc = 2.0_rk*rnd4 - 1.0_rk
!               phi_loc = 2.0_rk*pi*rnd5
!                       
!               Photons_vec(i_MC)%r_loc_vec(1) = r*(1.0_rk - mu_loc**2.0_rk)**0.5_rk*cos(phi_loc)
!               Photons_vec(i_MC)%r_loc_vec(2) = r*(1.0_rk - mu_loc**2.0_rk)**0.5_rk*sin(phi_loc)
!               Photons_vec(i_MC)%r_loc_vec(3) = r*mu_loc
!             !!!#######################
!             
!             call get_v_vec_from_angles_in_r_system(Photons_vec(i_MC)%theta,Photons_vec(i_MC)%phi,Photons_vec(i_MC)%r_loc_vec(:),Photons_vec(i_MC)%v_vec(:))   !! Calculates Photons_vec(i_MC)%v_vec
!                         
!             Photons_vec(i_MC)%ct_sourced = Photons_vec(i_MC)%ct
!             Photons_vec(i_MC)%r_sourced = Photons_vec(i_MC)%r
!             
!             Photons_vec(i_MC)%count_scatter = 0
!           
!           end do
!         
!         end do
!         
!         if (i_MC.ne.n_MC_ph) then
!           print *,'Initialize_Radiation_Passive_fireball: i_MC.ne.n_MC_ph', i_MC, n_MC_ph
! !           stop
!         end if
!         
!         deallocate(kT_vec,xi_cell_bnd_vec)
!         deallocate(P_cumul_r_x)
!         deallocate(lnx)
! 
!     End Subroutine Initialize_Radiation_Passive_fireball_Planck
!         
! !!!########################################################################################################################################
! 
! End Module MC_init








