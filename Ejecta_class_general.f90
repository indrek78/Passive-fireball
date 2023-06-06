
!!!######################################################################################################################################################################
!!!##################### Generic ejecta class and instance definition (basically a wrapper for different possible ejecta configurations) ################################
!!!######################################################################################################################################################################

!!! NBNBNB: t is actually ct, v_rad is actually v_rad/c

Module class_Ejecta

  !!!! NBNBNB: ejecta object could be generalized further if the type-bound routines accepted r_loc_vec instead of r

  use RealKind
  use auxiliary, ONLY: Findvalue_1dim_v4
!   use class_Sphere_homol ! not needed here explicitly
  use Ejecta_objects, ONLY: Sphere_object, Some_yet_undefined_ejecta_quantity  !!! All implemented ejecta objects, only one at a time is used
  implicit none
  private
  
  type, public :: Ejecta
    integer :: i_choose_ejecta = 0
    
    real(kind=rk) :: ct_0
!     real(kind=rk) :: kT => 
    contains
      procedure, PASS :: rho => get_density
      procedure, PASS :: v_rad => get_velocity
!       procedure, PASS :: R_sph => get_Rsph     !!! Not so good, sphere-specific quantity (loss of generality)
      procedure, PASS :: domain_size => get_domain_size
      procedure, PASS :: chk_esc
      
      procedure, PASS :: M_cell => get_M_cell   !!! No corresponding procedure in the specific object (i.e. Sphere_object)
      procedure, PASS :: kT => get_kT           !!! No corresponding procedure in the specific object (i.e. Sphere_object)
      procedure, PASS :: set_kT                 !!! No corresponding procedure in the specific object (i.e. Sphere_object)
      
      procedure, PASS :: Zav => get_Zav         !!! No corresponding procedure in the specific object (i.e. Sphere_object)
      procedure, PASS :: Zz => get_Zz         !!! No corresponding procedure in the specific object (i.e. Sphere_object)
      procedure, PASS :: mu_mol => get_mu_mol         !!! No corresponding procedure in the specific object (i.e. Sphere_object)
      !!!procedure :: R_sph => Sphere_object%R_sph
      
      procedure, PASS :: initialize
      procedure, PASS :: initialize_radiation
      
      procedure, PASS :: Emit_thermal
      procedure, PASS :: Min_cool_time
      
      procedure, PASS :: Print_spectrum
      procedure, PASS :: Print_ejecta => Printout_Ejecta_qties
      procedure, PASS :: Printout_Energy_conserv
      
      !!!###########
      procedure, PASS :: c_Tambov => get_Tambov
      !!!###########
      
      procedure, PASS :: chk_within
      
  end type Ejecta
  
  Contains
  
  
    Subroutine initialize(self)
      implicit none
      class(Ejecta) :: self
        if (self%i_choose_ejecta.eq.1) then
          call Sphere_object%initialize()
          self%ct_0 = Sphere_object%R_sph_0/Sphere_object%v_ej  !! initial time
        else
          !print *,'get_density: cannot handle i_choose_ejecta', self%i_choose_ejecta
!           write(get_density,*) self%i_choose_ejecta
          stop 'problem'
        end if
    End Subroutine initialize
    
    Subroutine initialize_radiation(self)
      implicit none
      class(Ejecta), intent(in) :: self
        if (self%i_choose_ejecta.eq.1) then
          call Sphere_object%initialize_radiation()
        else
          stop 'problem'
        end if
    End Subroutine initialize_radiation
    
!     initialize_radiation
  
    Function get_density(self,t,r) result(rho)
      !!! t is actually ct
      implicit none
      class(Ejecta), intent(in) :: self
!       type(Sphere_homol) :: Sphere_object
      real(kind=rk), intent(in) :: t, r
      real(kind=rk) rho
        if (self%i_choose_ejecta.eq.1) then
          rho = Sphere_object%rho(t,r)
        else
          !print *,'get_density: cannot handle i_choose_ejecta', self%i_choose_ejecta
!           write(get_density,*) self%i_choose_ejecta
          stop 'class_Ejecta/get_density PROBLEM: self%i_choose_ejecta'
        end if
    End Function get_density
    
    Function get_velocity(self,t,r) result(v_rad)
      !!! t is actually ct
      implicit none
      class(Ejecta), intent(in) :: self
!       type(Sphere_homol) :: Sphere_object
      real(kind=rk), intent(in) :: t, r
      real(kind=rk) v_rad
        if (self%i_choose_ejecta.eq.1) then
          v_rad = Sphere_object%v_rad(t,r)
        else
          !print *,'get_density: cannot handle i_choose_ejecta', self%i_choose_ejecta
!           write(get_density,*) self%i_choose_ejecta
          stop
        end if
    End Function get_velocity
    
!     Function get_Rsph(self,t) result(R_sph)
    Function get_domain_size(self,t) result(domain_size)
      !!! t is actually ct
      implicit none
      class(Ejecta), intent(in) :: self
      real(kind=rk), intent(in) :: t
      real(kind=rk) :: R_sph
      real(kind=rk) :: domain_size
        if (self%i_choose_ejecta.eq.1) then
          R_sph = Sphere_object%v_ej*t
          domain_size = R_sph
        else
!           print *,'get_Rsph: cannot handle i_choose_ejecta', self%i_choose_ejecta
          stop
        end if
    End Function get_domain_size
    
    
    Logical Function chk_esc(self,t,r)   !! Why is chk_esc duplicated with Sphere_object but the latter not used here? ... now it is used
    !!! t is actually ct
      implicit none
      class(Ejecta), intent(in) :: self
      real(kind=rk), intent(in) :: t, r
      real(kind=rk) :: R_sph
        if (self%i_choose_ejecta.eq.1) then
        
          chk_esc = Sphere_object%chk_esc(t,r)
        
!           R_sph = Sphere_object%v_ej*t
!           if (r.gt.R_sph) then
!             chk_esc = .TRUE.
!           else
!             chk_esc = .FALSE.
!           end if
          
        else
          stop
        end if
        
    End Function chk_esc
    
    
    Logical Function chk_within(self,t,r)   !! Why is chk_esc duplicated with Sphere_object but the latter not used here? ... now it is used
    !!! t is actually ct
      implicit none
      class(Ejecta), intent(in) :: self
      real(kind=rk), intent(in) :: t, r
      real(kind=rk) :: R_sph
        if (self%i_choose_ejecta.eq.1) then
          chk_within = Sphere_object%chk_within(t,r)
        else
          stop
        end if
        
    End Function chk_within
    
    
    
    Function get_kT(self,t,r) result(kT)
      !!! t is actually ct
      implicit none
      class(Ejecta), intent(in) :: self
      real(kind=rk), intent(in) :: t, r
      real(kind=rk) :: kT, xi, R_sph, dummy
      integer :: i_cell
        if (self%i_choose_ejecta.eq.1) then
          R_sph = Sphere_object%v_ej*t
          xi = r/R_sph
          call Findvalue_1dim_v4(xi,Sphere_object%xi_cell_bnd_vec(:),Sphere_object%xi_cell_bnd_vec(:),dummy,i_cell)
          if (i_cell.ne.0) then          
            kT = Sphere_object%kT_vec(i_cell)
          else
            kT = 0.0_rk
          end if
        else
          stop
        end if
    End Function get_kT
    
    
    Subroutine set_kT(self,t,r,kT)
      !!! t is actually ct
      implicit none
      class(Ejecta), intent(in) :: self
      real(kind=rk), intent(in) :: t, r, kT
      real(kind=rk) :: xi, R_sph, dummy
      integer :: i_cell
      
      if (self%i_choose_ejecta.eq.1) then
        R_sph = Sphere_object%v_ej*t
        xi = r/R_sph
        call Findvalue_1dim_v4(xi,Sphere_object%xi_cell_bnd_vec(:),Sphere_object%xi_cell_bnd_vec(:),dummy,i_cell)
        if (i_cell.ne.0) then          
          Sphere_object%kT_vec(i_cell) = kT
        else
          !!! A bit dangerous to do nothing, transferred energy must go somewhere...
          stop 'set_kT: code shouldnt get here'
        end if
      else
        stop
      end if
    End Subroutine set_kT
    

    Function get_M_cell(self,t,r) result(M_cell)
      !!! t is actually ct
      implicit none
      class(Ejecta), intent(in) :: self
      real(kind=rk), intent(in) :: t, r
      real(kind=rk) :: M_cell, xi, R_sph, dummy
      integer :: i_cell
        if (self%i_choose_ejecta.eq.1) then
          R_sph = Sphere_object%v_ej*t
          xi = r/R_sph
          call Findvalue_1dim_v4(xi,Sphere_object%xi_cell_bnd_vec(:),Sphere_object%xi_cell_bnd_vec(:),dummy,i_cell)
          if (i_cell.ne.0) then          
            M_cell = Sphere_object%M_cell_vec(i_cell)
!             M_cell = Sphere_object%M_ej*(Sphere_object%xi_cell_bnd_vec(i_cell+1)**(3.0_rk - Sphere_object%alpha_r) - Sphere_object%xi_cell_bnd_vec(i_cell)**(3.0_rk - Sphere_object%alpha_r))/(1.0_rk - Sphere_object%xi_min**(3.0_rk - Sphere_object%alpha_r))
          else
            M_cell = 0.0_rk
          end if
        else
          stop
        end if
    End Function get_M_cell
    
    Function get_Zav(self) result(Zav)
      implicit none
      class(Ejecta), intent(in) :: self
!       real(kind=rk), intent(in) :: t, r
      real(kind=rk) :: Zav
        if (self%i_choose_ejecta.eq.1) then
          Zav = Sphere_object%Zav  !! single value for the netire region (could be generalized similar to get_kT)
        else
          !print *,'get_density: cannot handle i_choose_ejecta', self%i_choose_ejecta
!           write(get_density,*) self%i_choose_ejecta
          stop
        end if
    End Function get_Zav
    
    
    Function get_Zz(self) result(Zz)
      implicit none
      class(Ejecta), intent(in) :: self
      real(kind=rk) :: Zz
        if (self%i_choose_ejecta.eq.1) then
          Zz = Sphere_object%Zz  !! single value for the netire region (could be generalized similar to get_kT)
        else
          !print *,'get_density: cannot handle i_choose_ejecta', self%i_choose_ejecta
!           write(get_density,*) self%i_choose_ejecta
          stop
        end if
    End Function get_Zz
    
    
    Function get_mu_mol(self) result(mu_mol)
      implicit none
      class(Ejecta), intent(in) :: self
!       real(kind=rk), intent(in) :: t, r
      real(kind=rk) :: mu_mol
        if (self%i_choose_ejecta.eq.1) then
          mu_mol = Sphere_object%mu_mol  !! single value for the netire region (could be generalized similar to get_kT)
        else
          !print *,'get_density: cannot handle i_choose_ejecta', self%i_choose_ejecta
!           write(get_density,*) self%i_choose_ejecta
          stop
        end if
    End Function get_mu_mol
    
    
    Function get_Tambov(self) result(c_Tambov)
      implicit none
      class(Ejecta), intent(in) :: self
      real(kind=rk) :: c_Tambov
        if (self%i_choose_ejecta.eq.1) then
          c_Tambov = Sphere_object%c_Tambov  !! single value for the netire region (could be generalized similar to get_kT)
        else
          stop 'Ejecta_class_general/get_Tambov'
        end if
    End Function get_Tambov
    
    !!!######### Construction ##############
    
    Subroutine Emit_thermal(self,ct,d_ct)
      implicit none
      
      class(Ejecta), intent(in) :: self
      real(kind=rk), intent(in) :: ct, d_ct
      
      if (self%i_choose_ejecta.eq.1) then
        call Sphere_object%Emit_therm_brems(ct,d_ct)
      else
        stop
      end if
    
    End Subroutine Emit_thermal
    
    
    
    Function Min_cool_time(self,ct) result(ct_cool_min)
      implicit none
      
      class(Ejecta), intent(in) :: self
      real(kind=rk), intent(in) :: ct
      real(kind=rk) :: ct_cool_min
      
      if (self%i_choose_ejecta.eq.1) then
        ct_cool_min = Sphere_object%tcool_therm_brems(ct)
      else
        stop
      end if
    
    End Function Min_cool_time
    
    
    
    !!!################################### Printouts ############################################
      Subroutine Print_spectrum(self,ct,sw_reset_Photons_escape)
        implicit none
        class(Ejecta), intent(in) :: self
        real(kind=rk), intent(in) :: ct
        logical, intent(in) :: sw_reset_Photons_escape
        if (self%i_choose_ejecta.eq.1) then
          call Sphere_object%Print_spectrum(ct,sw_reset_Photons_escape)
        else
          stop
        end if      
      End Subroutine Print_spectrum  
    
    
      Subroutine Printout_Ejecta_qties(self,ct)
        implicit none
        class(Ejecta), intent(in) :: self
        real(kind=rk), intent(in) :: ct
        if (self%i_choose_ejecta.eq.1) then
          call Sphere_object%Print_Sphere(ct)
        else
          stop
        end if      
      End Subroutine Printout_Ejecta_qties
      
      
      Subroutine Printout_Energy_conserv(self,ct,sw_reset_Photons_escape)
        implicit none
        class(Ejecta), intent(in) :: self
        real(kind=rk), intent(in) :: ct
        logical, intent(in) :: sw_reset_Photons_escape
        if (self%i_choose_ejecta.eq.1) then
          call Sphere_object%Printout_Energy_conserv(ct,sw_reset_Photons_escape)
        else
          stop
        end if      
      End Subroutine Printout_Energy_conserv
      
      
      
                
    
    !!!##########################################################################################
        
End Module class_Ejecta


!!!################ Store Ejecta object to be used by the propagation routines ###############
Module Ejecta_generic
  use class_Ejecta
  implicit none
  type(Ejecta) :: Ejecta_obj
  ! class(Ejecta) :: Ejecta_obj ! CLASS variable ‘ejecta_obj’ at (1) must be dummy, allocatable or pointer
End Module Ejecta_generic

!!!###################################################################################################
!!!##################### End: Generic ejecta class and instance definition ###########################
!!!###################################################################################################


