!!! Renamed from Maxw_routines.f90

 
Module JCMaxwell

  use Realkind
  
  implicit none
  
  contains
  
  function f_MW(z,theta)
  
    !!! maxwellian (arbitrary norm)
  
    implicit none
    
    real(kind=rk), intent(in) :: theta,z
    real(kind=rk) :: f_MW
    real(kind=rk) :: gamma, gamma_min_one

    gamma = (z*z + 1.0_rk)**0.5_rk
!     gamma_min_one = gamma - 1.0_rk
    gamma_min_one = z*z/(gamma + 1.0_rk)
!     f_MW = z**3.0_rk*exp(-gamma/theta)
    f_MW = z**3.0_rk*exp(-gamma_min_one/theta)
    
  end function f_MW  
  
  function f_MW_norm(z,theta)
  
    !!! maxwellian (normalized to a single electron)
  
    implicit none
    
    real(kind=rk), intent(in) :: theta, z
    real(kind=rk) :: f_MW_norm
    real(kind=rk) :: gamma, gamma_min_one, gminone_0, gminone_1, gminone_2

    gamma = (z*z + 1.0_rk)**0.5_rk
!     gamma_min_one = gamma - 1.0_rk
    gamma_min_one = z*z/(gamma + 1.0_rk)
!     f_MW = z**3.0_rk*exp(-gamma/theta)
    f_MW_norm = z**3.0_rk*exp(-gamma_min_one/theta)
    
    call MW_moments(theta,gminone_0,gminone_1,gminone_2)
    f_MW_norm = f_MW_norm/gminone_0
    
  end function f_MW_norm  
  
  
  subroutine MW_moments(theta,gminone_0,gminone_1,gminone_2)
  
    !! moments of (\gamma-1), norm to 0-th moment (exc. the 0-th moment itself)
  
    implicit none
    
    real(kind=rk), intent(in) :: theta
    real(kind=rk), intent(out) :: gminone_1, gminone_2
    real(kind=rk) :: g_min_one_min, g_min_one_max, z_min, z_max, lnz_min, lnz_max, d_lnz
    real(kind=rk) :: gminone_0, gamma_min_one, z, gamma, f
    real(kind=rk), dimension(:), allocatable :: lnz
    integer :: i, n_z
  
    n_z = 10000
    g_min_one_min = theta*1.0e-7_rk
    g_min_one_max = theta*1.0e2_rk
    
    z_min = ((g_min_one_min+1.0_rk)**2.0_rk - 1.0_rk)**0.5_rk
    z_max = ((g_min_one_max+1.0_rk)**2.0_rk - 1.0_rk)**0.5_rk
    lnz_min = log(z_min)
    lnz_max = log(z_max)

    allocate(lnz(n_z))
    d_lnz = (lnz_max - lnz_min)/(dble(n_z) - 1.0_rk)
    do i = 1,n_z
      lnz(i) = lnz_min + (dble(i) - 1.0_rk)*d_lnz
    end do
    
    
    gminone_0 = 0.0_rk
    gminone_1 = 0.0_rk
    gminone_2 = 0.0_rk
    do i = 1,n_z
      z = exp(lnz(i))
      gamma = (z*z + 1.0_rk)**0.5_rk
      gamma_min_one = z*z/(gamma + 1.0_rk)
      f = f_MW(z,theta)
      gminone_2 = gminone_2 + f*gamma_min_one**2.0_rk ! *gamma**2.0_rk
      gminone_1 = gminone_1 + f*gamma_min_one	 ! *gamma ! *(gamma-1.0_rk)
      gminone_0 = gminone_0 + f
    end do  
    gminone_2 = gminone_2/gminone_0
    gminone_1 = gminone_1/gminone_0
    gminone_0 = gminone_0*d_lnz
  
  end subroutine MW_moments
  
End Module JCMaxwell




Module Ph_thermal

  use Realkind
  use constants

  implicit none
  
  Contains
  
    Function f_Planck(x,theta)
      !!! dN/(dlnx*dV) (angle-integrated)
      implicit none
      real(kind=rk), intent(in) :: x, theta
      real(kind=rk) :: f_Planck
        f_Planck = 8.0_rk*pi/(lambda_C**3.0_rk)*x**3.0_rk/(exp(x/theta) - 1.0_rk)
    End Function f_Planck
    
    
    Subroutine Planck_moments(theta,n_ph,u_ph)
      !!! Dimensionless Planck radiation energy density (u_ph) and number density (n_ph)
    
      implicit none
      
      real(kind=rk), intent(in) :: theta
      real(kind=rk), intent(out) :: n_ph, u_ph
      real(kind=rk) :: T
      
        T = m_e*c0**2.0_rk*theta/k_B
        u_ph = a_rad*T**4.0_rk/(m_e*c0**2.0) !! energy density in m_e*c^2 units
        n_ph = a_rad*T**3.0_rk/(2.701*k_B)
    
    End Subroutine Planck_moments
    
    
! ! !     Function f_Wien(x,theta,n_ph)
    Function f_Wien(x,theta,u_ph)
      !!! dN/(dlnx*dV) (angle-integrated)
      !!! u_ph is without energy dimension (i.e. units of m_e*c^2)
      implicit none
      real(kind=rk), intent(in) :: x, theta, u_ph  !!! n_ph
      real(kind=rk) :: Cc
      real(kind=rk) :: f_Wien
      
!         Cc = n_ph/(2.0_rk*theta**3.0_rk)
        Cc = u_ph/(6.0_rk*theta**4.0_rk)
        
        f_Wien = Cc*x**3.0_rk*exp(-x/theta)
    End Function f_Wien
    
End Module Ph_thermal










  
