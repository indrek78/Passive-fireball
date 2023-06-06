Module Geometry

  use Realkind

  implicit none
  
  Contains
  
    Subroutine Vector_norm(a,norm)
    
      implicit none
      
      real(kind=rk), dimension(3), intent(in) :: a
      real(kind=rk), intent(out) :: norm
      
      norm = (a(1)**2.0_rk + a(2)**2.0_rk + a(3)**2.0_rk)**0.5_rk
!       norm = sqrt(sum(a(:)**2.0_rk))
    
    End Subroutine Vector_norm
    
    
    
    Subroutine Scalar_product(a,b,ab)
    
      implicit none
      
      real(kind=rk), dimension(3), intent(in) :: a, b
      real(kind=rk), intent(out) :: ab
      
      ab = a(1)*b(1)+ a(2)*b(2) + a(3)*b(3)
      
    End Subroutine Scalar_product
    
    Subroutine Cross_product(a,b,a_cross_b)
    
      implicit none
      
      real(kind=rk), dimension(3), intent(in) :: a, b
      real(kind=rk), dimension(3), intent(out) :: a_cross_b
      
      a_cross_b(1) = a(2)*b(3) - a(3)*b(2)
      a_cross_b(2) = a(3)*b(1) - a(1)*b(3)
      a_cross_b(3) = a(1)*b(2) - a(2)*b(1)
      
    End Subroutine Cross_product
  
  
    Subroutine get_photon_angles_in_rz_system(r_loc_vec,v_ph_vec,theta_vr,phi_r_vz)
  
      !!! 15.08.18
      !!! Given the cartesian vectors r_loc_vec (radius vector), v_ph_vec (photon velocity), calculate angles theta_vr (angle betw. v and r),
      !!! phi_r_vz (angle between the projections of v and z onto the plane perpendicular to r; apparently, looking along the r direction, phi_r_vz is measured clockwise starting from the projection of z onto the plane perp. to r (CHECK!))
      !!! The angles theta_vr, phi_r_vz are the relevant ones for the scattering process if the bulk velocity of the medium is directed in the radial direction.
  
      use Constants, ONLY: pi
  
      implicit none
      
      real(kind=rk), dimension(3), intent(in) :: r_loc_vec, v_ph_vec
      real(kind=rk), intent(out) :: theta_vr, phi_r_vz
      
      real(kind=rk), dimension(3) :: z_vec, r_loc_vec_norm, v_ph_vec_norm
      real(kind=rk), dimension(3) :: eps_vec
      real(kind=rk), dimension(3) :: r_cross_z
      
      real(kind=rk) :: r_norm, v_norm
      real(kind=rk) :: cos_phi_r_vz, cos_theta_vr
      real(kind=rk) :: rv, rz, vz, rv2, rz2
      real(kind=rk) :: rzv
      real(kind=rk) :: eps_scal, eps_scal_2, veps, cos_phi_r_vz_backup, Num, Denom
      
      real(kind=rk), parameter :: tol = 1.0e-6_rk ! 1.0e-12_rk
      real(kind=rk), parameter :: tol2 = 1.0e-8_rk ! 1.0e-10_rk !! 1.0e-12_rk !!! set to 1e-10 on 07.08.19 (was 1e-12), to 1e-8 on 09.03.21 (was 1e-10)
      
      real(kind=rk) :: test
      
      z_vec = (/0.0_rk, 0.0_rk, 1.0_rk/)
      
!       call Vector_norm(r_loc_vec,r_norm)
!       call Vector_norm(v_ph_vec,v_norm)
      r_norm = norm2(r_loc_vec)
      v_norm = norm2(v_ph_vec)
      
      r_loc_vec_norm = r_loc_vec/r_norm
      v_ph_vec_norm = v_ph_vec/v_norm
      
      !!### Scalar products rv, rz of normalized z, v, r #############
        call Scalar_product(r_loc_vec_norm,v_ph_vec_norm,rv)
        call Scalar_product(r_loc_vec_norm,z_vec,rz)
        call Scalar_product(v_ph_vec_norm,z_vec,vz)
      !!!#############################################################
      
      rv2 = rv*rv
      rz2 = rz*rz
      
!       if ((abs(rz).lt.(1.0_rk + tol2)).and.(abs(rz).gt.(1.0_rk - tol2))) then
!         print *,'get_photon_angles_in_rz_system PROBLEM; cannot currently handle rz = +-1', rz
!         print *,'r_loc_vec_norm', r_loc_vec_norm
!         print *,'v_ph_vec_norm', v_ph_vec_norm
!         stop
!       end if
      
      if ((abs(rv).lt.(1.0_rk + tol)).and.(abs(rv).gt.(1.0_rk - tol))) then
        cos_phi_r_vz = 1.0_rk   !!! ARBITRARY!!!
      else
        cos_phi_r_vz = (vz - rz*rv)/((1.0_rk - rv**2.0_rk)*(1.0_rk - rz**2.0_rk))**0.5_rk
        
        
        !!!##### Construction ############
          eps_vec(:) = r_loc_vec_norm - z_vec
          eps_scal = norm2(eps_vec)
          eps_scal_2 = eps_scal*eps_scal
          
          call Scalar_product(v_ph_vec_norm,eps_vec,veps)
          
          Num = 0.5_rk*eps_scal_2*(vz + veps) - veps
          Denom = eps_scal*((1.0_rk + rz)/2.0_rk*(1.0_rk - rv2))**0.5_rk
        
          cos_phi_r_vz_backup = Num/Denom
          
          if ((abs(rz).lt.(1.0_rk + tol2)).and.(abs(rz).gt.(1.0_rk - tol2))) then
            cos_phi_r_vz = cos_phi_r_vz_backup
          else
            if (abs(cos_phi_r_vz - cos_phi_r_vz_backup).gt.tol) then
              print *,'get_photon_angles_in_rz_system PROBLEM; cos_phi_r_vz.ne.cos_phi_r_vz_backup', cos_phi_r_vz, cos_phi_r_vz_backup
              print *,'r_loc_vec_norm', r_loc_vec_norm
              print *,'v_ph_vec_norm', v_ph_vec_norm
              print *,'eps_scal, eps_scal_2, veps, vz, rz, rv2', eps_scal, eps_scal_2, veps, vz, rz, rv2
              print *,'eps_vec', eps_vec
              stop
            end if
          
          end if
        !!!###############################
        
        
        
      end if
      
      if ((cos_phi_r_vz.gt.1.0_rk)) then
        if (cos_phi_r_vz.lt.(1.0_rk + tol)) then
          cos_phi_r_vz = 1.0_rk
        else
          print *,'get_photon_angles_in_rz_system PROBLEM: cos_phi_r_vz > 1', cos_phi_r_vz
          print *,'vz, rz, rv', vz, rz, rv
          print *,'r_loc_vec_norm', r_loc_vec_norm
          print *,'v_ph_vec_norm', v_ph_vec_norm
          stop
        end if
      end if
      
      if ((cos_phi_r_vz.lt.-1.0_rk)) then
        if (cos_phi_r_vz.gt.(-1.0_rk - tol)) then
          cos_phi_r_vz = -1.0_rk
        else
          print *,'get_photon_angles_in_rz_system PROBLEM: cos_phi_r_vz < -1', cos_phi_r_vz
          print *,'vz, rz, rv', vz, rz, rv
          print *,'r_loc_vec_norm', r_loc_vec_norm
          print *,'v_ph_vec_norm', v_ph_vec_norm
          stop
        end if
      end if
      
      cos_theta_vr = rv
      
      if (cos_theta_vr.gt.1.0_rk) then
        if (cos_theta_vr.lt.(1.0_rk + tol)) then
          cos_theta_vr = 1.0_rk
          theta_vr = 0.0_rk
        else
          print *,'get_photon_angles_in_rz_system PROBLEM: cos_theta_vr', cos_theta_vr
          print *,'vz, rz, rv', vz, rz, rv
          print *,'r_loc_vec_norm', r_loc_vec_norm
          print *,'v_ph_vec_norm', v_ph_vec_norm
          stop
        end if
      else if (cos_theta_vr.lt.-1.0_rk) then  
        if (cos_theta_vr.ge.(-1.0_rk - tol)) then
          cos_theta_vr = -1.0_rk
          theta_vr = pi
        else
          print *,'get_photon_angles_in_rz_system PROBLEM: cos_theta_vr', cos_theta_vr
          print *,'vz, rz, rv', vz, rz, rv
          print *,'r_loc_vec_norm', r_loc_vec_norm
          print *,'v_ph_vec_norm', v_ph_vec_norm
          stop
        end if
      else
        theta_vr = acos(cos_theta_vr)
      end if
      
      !!!##### TEST, can be removed #############
        if ((theta_vr.gt.0.0_rk).or.(theta_vr.le.0.0_rk)) then
          Continue
        else
          print *,'theta_vr, cos_theta_vr', theta_vr, cos_theta_vr
          print *,'r_loc_vec', r_loc_vec
          print *,'v_ph_vec', v_ph_vec
          print *,'r_loc_vec_norm', r_loc_vec_norm
          print *,'v_ph_vec_norm', v_ph_vec_norm
          print *,'r_norm, v_norm', r_norm, v_norm
          call Scalar_product(r_loc_vec_norm,v_ph_vec_norm,test)
          print *,'test', test
        end if 
      !!!########################################            
      
      call Cross_product(r_loc_vec_norm,z_vec,r_cross_z)
      call Scalar_product(r_cross_z,v_ph_vec_norm,rzv)
      
      if (rzv.ge.0.0_rk) then
        phi_r_vz = acos(cos_phi_r_vz)
        if (phi_r_vz.lt.(0.0_rk-tol)) then
          print *,'get_photon_angles_in_rz_system PROBLEM; phi_r_vz, cos_phi_r_vz, rzv', phi_r_vz, cos_phi_r_vz, rzv
          stop
        end if
      else if (rzv.lt.0.0_rk) then
        phi_r_vz = -acos(cos_phi_r_vz)
        if (phi_r_vz.ge.(0.0_rk+tol)) then
          print *,'get_photon_angles_in_rz_system PROBLEM; phi_r_vz, cos_phi_r_vz, rzv', phi_r_vz, cos_phi_r_vz, rzv
          stop
        end if
      end if
      
      if ((phi_r_vz.gt.(pi + tol)).or.(phi_r_vz.lt.(-pi - tol))) then
        print *,'get_photon_angles_in_rz_system PROBLEM: phi_r_vz should be between -pi an pi; phi_r_vz', phi_r_vz
        stop
      end if
      
      if ((theta_vr.gt.(pi + tol)).or.(theta_vr.lt.(-tol))) then
        print *,'get_photon_angles_in_rz_system PROBLEM: theta_vr should be between 0 an pi; theta_vr', theta_vr
        stop
      end if
  
    End Subroutine get_photon_angles_in_rz_system
    
       
    
    
    Subroutine get_v_vec_from_angles_in_r_system(theta_vr,phi_r_vz,r_loc_vec,v_ph_vec_norm) 
    
      !!! same as v0, but less intermediate normalizations to save time
    
!       use timing, ONLY: cpu_t_1, cpu_t_2, cpu_t_3, cpu_t_4, cpu_t_5, cpu_t_6, cpu_t_7, cpu_t_8 
    
      implicit none
      
      real(kind=rk), dimension(3), intent(in) :: r_loc_vec
      real(kind=rk), intent(in) :: theta_vr, phi_r_vz
      real(kind=rk), dimension(3), intent(out) ::  v_ph_vec_norm
      
      real(kind=rk), dimension(3) :: z_vec, r_loc_vec_norm
      real(kind=rk), dimension(3) :: r_cross_z
      real(kind=rk), dimension(3) :: u_rz_para_norm, u_rz_perp_norm
      real(kind=rk) :: cos_theta_vr, sin_theta_vr
      
      real(kind=rk) :: r_norm, normaalne
      real(kind=rk) :: rz, sin_rz
      
!       real(kind=rk) :: test_norm
      real(kind=rk), parameter :: tol = 1.0e-4_rk !!! 1.0e-6_rk ! 1.0e-8 !! relaxed to 1e-6 on 21.01.19  !! relaxed to 1e-4 (1e-5) on 14.18.19
      
      real(kind=rk) :: cpu_t_0_add, cpu_t_1_add, cpu_t_2_add, cpu_t_3_add, cpu_t_4_add, cpu_t_5_add, cpu_t_6_add, cpu_t_7_add, cpu_t_8_add 
      
      z_vec = (/0.0_rk, 0.0_rk, 1.0_rk/)
     
!   call cpu_time(cpu_t_0_add)   
!   call cpu_time(cpu_t_1_add)
!   cpu_t_1 = cpu_t_1 + (cpu_t_1_add - cpu_t_0_add)
      
!       call Vector_norm(r_loc_vec,r_norm)
      r_norm = norm2(r_loc_vec)
      
      r_loc_vec_norm = r_loc_vec/r_norm
      
!   call cpu_time(cpu_t_2_add)
!   cpu_t_2 = cpu_t_2 + (cpu_t_2_add - cpu_t_1_add)
      
      call Scalar_product(r_loc_vec_norm,z_vec,rz)

!   call cpu_time(cpu_t_3_add)
!   cpu_t_3 = cpu_t_3 + (cpu_t_3_add - cpu_t_2_add)
      
      
      sin_rz = (1.0_rk - rz**2.0_rk)**0.5_rk
      
      !!!!!### u_rz_para_norm - unit vector in the direction of the projection of z onto the plane perp to r ######
      u_rz_para_norm = (z_vec - r_loc_vec_norm*rz)/sin_rz
      
      call Cross_product(r_loc_vec_norm,z_vec,r_cross_z)
      u_rz_perp_norm = r_cross_z/sin_rz

!   call cpu_time(cpu_t_4_add) 
!   cpu_t_4 = cpu_t_4 + (cpu_t_4_add - cpu_t_3_add)
  
       
      
      cos_theta_vr = cos(theta_vr)
      sin_theta_vr = sin(theta_vr)
     
!   call cpu_time(cpu_t_7_add) 
!   cpu_t_7 = cpu_t_7 + (cpu_t_7_add - cpu_t_4_add)
      
!       v_ph_vec_norm = r_loc_vec_norm*cos(theta_vr) + u_rz_para_norm*sin(theta_vr)*cos(phi_r_vz) + u_rz_perp_norm*sin(theta_vr)*sin(phi_r_vz)
      v_ph_vec_norm = r_loc_vec_norm*cos_theta_vr + u_rz_para_norm*sin_theta_vr*cos(phi_r_vz) + u_rz_perp_norm*sin_theta_vr*sin(phi_r_vz)
      
!       call cpu_time(cpu_t_5_add)
!       cpu_t_5 = cpu_t_5 + (cpu_t_5_add - cpu_t_7_add)
      
!        print *,'phi_r_vz', phi_r_vz
      
if (1.eq.1) then    !!! do not renormalize      
      if (1.eq.0) then
        call Vector_norm(v_ph_vec_norm,normaalne)
      else
!         normaalne = (v_ph_vec_norm(1)**2.0_rk + v_ph_vec_norm(2)**2.0_rk + v_ph_vec_norm(3)**2.0_rk)**0.5_rk
        normaalne = norm2(v_ph_vec_norm)
      end if
        
        
!   call cpu_time(cpu_t_8_add) 
!   cpu_t_8 = cpu_t_8 + (cpu_t_8_add - cpu_t_7_add)
      
      !!!### Check #####
        if ((normaalne.gt.(1.0_rk + tol)).or.(normaalne.lt.(1.0_rk - tol))) then
          print *,'get_v_vec_from_angles_in_r_system PROBLEM; normalization of v_ph_vec_norm', normaalne
          print *,'z_vec, r_loc_vec_norm, rz', z_vec, r_loc_vec_norm, rz
          print *,'v_ph_vec_norm', v_ph_vec_norm
          print *,'u_rz_para_norm', u_rz_para_norm
          print *,'u_rz_perp_norm', u_rz_perp_norm
          print *,'sin_theta_vr', sin_theta_vr
          print *,'cos_theta_vr', cos_theta_vr
          print *,'phi_r_vz', phi_r_vz
          stop
        end if
      !!!################
      
      v_ph_vec_norm(:) = v_ph_vec_norm(:)/normaalne

end if      
      
    End Subroutine get_v_vec_from_angles_in_r_system
    
    
  
End Module Geometry
