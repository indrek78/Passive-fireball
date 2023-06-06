Module Printouts

  use realkind
  use constants
  
  implicit none
  
  Contains
  
  
      Subroutine Save_output(ct,sw_reset_Photons_escape)
      
        use Ejecta_generic, ONLY: Ejecta_obj
      
        implicit none
        
          real(kind=rk), intent(in) :: ct
          logical, intent(in) :: sw_reset_Photons_escape
        
          call Ejecta_obj%Print_spectrum(ct,sw_reset_Photons_escape)
          call Ejecta_obj%Print_ejecta(ct)
          call Ejecta_obj%Printout_Energy_conserv(ct,sw_reset_Photons_escape)

      End Subroutine Save_output
! 
!     Subroutine Printout_timedep_spectrum(ct,sw_reset_Photons_escape)
!           
!       use Photon_distribution_escape, ONLY: Photons_escape, i_exist_max_escape
!       use Photon_distribution, ONLY: Photons_vec
!       use auxiliary, ONLY: Setup_grid, bin
!       use Init_printout, ONLY: FF
!       
!       implicit none
!       
!       real(kind=rk), intent(in) :: ct
!       real(kind=rk), dimension(:), allocatable :: x_vec, lnx_vec, Nx, Nx_MC, Nx_esc
!       real(kind=rk) :: lnx_min, lnx_max, d_lnx, x
!       integer :: n_x_bin, i_ph, n_ph, n_ph_esc, i
!       integer, save :: i_save = 0
!       integer, save :: i_exist_max_escape_chk = 0
!       logical, intent(in) :: sw_reset_Photons_escape
!       
!       3 format(e14.7,' ',e14.7,' ',e14.7)
!       4 format(e14.7,' ',e14.7,' ',e14.7,' ',e14.7)
!       
!         n_x_bin = 20
!         lnx_min = log(1.0e-6_rk)
!         lnx_max = log(1.0e-2_rk)
!         allocate(x_vec(n_x_bin),lnx_vec(n_x_bin))        
!         call Setup_grid(lnx_vec,d_lnx,lnx_min,lnx_max)
!         x_vec(:) = exp(lnx_vec)
! 
!         if (.NOT.(allocated(Nx))) then  !!! 01.10.20
!           allocate(Nx(n_x_bin), Nx_MC(n_x_bin))
!           Nx(:) = 0.0_rk
!           Nx_MC(:) = 0.0_rk
!         end if
!         
!         if (allocated(Photons_vec)) then          
!           n_ph = size(Photons_vec)
!           do i_ph = 1,n_ph
!             x = Photons_vec(i_ph)%x
!             
!             if (Photons_vec(i_ph)%exist) then
!               call Bin(x,x_vec,Photons_vec(i_ph)%weight,Nx)   
!               call Bin(x,x_vec,1.0_rk,Nx_MC)    !! wastes time, should be merged with previous
!             end if
!             
!           end do
!         end if
!         
!         if (i_save.eq.0) then
!           open(100,FILE='results/Spectra/Spec_photons_' // trim(adjustl(FF)) // '.dat')
! !           open(100,FILE='results/Spectra/Spec_photons.dat')
!         else
!           open(100,POSITION='APPEND',FILE='results/Spectra/Spec_photons_' // trim(adjustl(FF)) // '.dat')
! !           open(100,POSITION='APPEND',FILE='results/Spectra/Spec_photons.dat')
!         end if
! 
!         do i = 1,n_x_bin
!           x = exp(log(x_vec(i)) + d_lnx/2.0_rk)
!           write(100,4) ct, x, Nx(i)/d_lnx, Nx_MC(i)
!         end do
!         close(100)
!         deallocate(Nx,Nx_MC)
!         
!         
! 
!         !!!############ Escaping photons #################
!           
!           if (.NOT.(allocated(Nx_esc))) then  !!! 01.10.20
!             allocate(Nx_esc(n_x_bin))
!             Nx_esc(:) = 0.0_rk
!           end if
!           
!           if (sw_reset_Photons_escape) then  !! Is Photons_escape reset at each timestep? (01.10.20)
!             Continue  !! DO NOT RESET Nx_esc, i_exist_max_escape_chk
!           else
!             i_exist_max_escape_chk = 0
!             Nx_esc(:) = 0.0_rk
!           end if
!                     
!           if (allocated(Photons_escape)) then          
!             n_ph_esc = size(Photons_escape)
!             do i = 1,n_ph_esc
!               x = Photons_escape(i)%x
!               if (Photons_escape(i)%exist) then
!               
!                 if (i.le.i_exist_max_escape) then   !!! SHOULD (!) be redundant
!                   call Bin(x,x_vec,Photons_escape(i)%weight,Nx_esc)              
!                 end if
!                 i_exist_max_escape_chk = i_exist_max_escape_chk + 1
!   
!               end if
!             end do
!           end if
!           
!           if (i_save.eq.0) then
!             open(110,FILE='results/Spectra/Spec_esc_cumul_' // trim(adjustl(FF)) // '.dat')
! !             open(110,FILE='results/Spectra/Spec_esc_cumul.dat')
!           else
!             open(110,POSITION='APPEND',FILE='results/Spectra/Spec_esc_cumul_' // trim(adjustl(FF)) // '.dat')
! !             open(110,POSITION='APPEND',FILE='results/Spectra/Spec_esc_cumul.dat')
!           end if
! 
!           do i = 1,n_x_bin
!             x = exp(log(x_vec(i)) + d_lnx/2.0_rk)
!             write(110,3) ct, x, Nx_esc(i)/d_lnx 
!           end do
!           close(110)
!           
!           
!           if (.NOT.(sw_reset_Photons_escape)) deallocate(Nx_esc)
!         
!         !!!################################################
!         
!         i_save = i_save + 1
!         
!     End Subroutine Printout_timedep_spectrum
!     
!     
!     
!     
!     Subroutine Printout_Sphere_qties(ct)
!     
!       use Ejecta_objects, ONLY: Sphere_object
!       use Ejecta_generic, ONLY: Ejecta_obj
!       use Init_printout, ONLY: FF
!       
!       implicit none
!       
!       real(kind=rk), intent(in) :: ct
!       integer :: n_r, i_r
!       integer, save :: i_save = 0
!       
!       4 format(e14.7,' ',e14.7,' ',e14.7,' ',e14.7)
!       
!       if (Ejecta_obj%i_choose_ejecta.eq.1) then   !!! Printouts specific to a passive homologous ejecta
!         if (i_save.eq.0) then
!           open(200,FILE='results/Ejecta/Ejecta_qties_' // trim(adjustl(FF)) // '.dat')
! !           open(200,FILE='results/Ejecta/Ejecta_qties.dat')
!         else
!           open(200,POSITION='APPEND',FILE='results/Ejecta/Ejecta_qties_' // trim(adjustl(FF)) // '.dat')
! !           open(200,POSITION='APPEND',FILE='results/Ejecta/Ejecta_qties.dat')
!         end if
!         
!         n_r = size(Sphere_object%kT_vec)
!         do i_r = 1,n_r
!         
! !         if (alpha_r.eq.1.0_rk) then
! !           tau_T = 2.0_rk*M_ej/(4.0_rk*pi*R_sph_0**2.0_rk)*Zav*sigma_T/m_p*(-log(xi_min))/(1.0_rk - xi_min**2.0_rk)
! !         else
! !           tau_T = (3.0_rk - alpha_r)/(1.0_rk - alpha_r)*M_ej/(4.0_rk*pi*R_sph_0**2.0_rk)*Zav*sigma_T/m_p*(1.0_rk - xi_min**(1.0_rk - alpha_r))/(1.0_rk - xi_min**(3.0_rk - alpha_r))
! !         end if
!         
!           write(200,4) ct, Sphere_object%xi_cell_vec(i_r), Sphere_object%kT_vec(i_r), Sphere_object%M_cell_vec(i_r)
!         end do
!         
!         close(200)
!         
!       else
!         !!! do nothing
!       end if
!     
!     
!       i_save = i_save + 1 
!     
!     End Subroutine Printout_Sphere_qties
!     
!     
    
    
    
    
End Module Printouts
        
