Module RealKind
  implicit none
 integer, parameter :: rk = 8 ! 16 for precalculating pair production
end module RealKind


Module Constants
use RealKind
real(kind=rk), parameter :: m_e = 9.1093826E-28_rk ! g
real(kind=rk), parameter :: e = 4.803250E-10_rk
real(kind=rk), parameter :: c0 = 2.99792458E10_rk ! cm/s
real(kind=rk), parameter :: pi = 3.14159265358979323846264_rk
real(kind=rk), parameter :: r_0 = 2.817940325E-13_rk ! cm
real(kind=rk), parameter :: h_pc = 6.6260693E-27_rk ! erg*s (Plank constant)
real(kind=rk), parameter :: sigma_t = 6.65245873E-25_rk ! cm**2
real(kind=rk), parameter :: sigma_SB = 5.670400e-5_rk ! erg/(s*cm**2*K**4)		Check; CHANGED sigma -> sigma_SB
real(kind=rk), parameter :: k_B = 1.3806505e-16_rk ! erg/K
real(kind=rk), parameter :: sigma_n = 3.0e-26_rk ! cm**2. Approximate nuclear cross-section
real(kind=rk), parameter :: lambda_C = 2.42631021e-10_rk			! Compton wavelength

real(kind=rk), parameter :: m_p = 1.672621637E-24_rk ! g
real(kind=rk), parameter :: a_rad = 7.565767E-15_rk ! erg/(cm**3*K**4)
real(kind=rk), parameter :: afs = 7.2973525698e-3_rk

real(kind=rk), parameter :: Msun = 1.9891e33_rk	! (g)
real(kind=rk), parameter :: eV = 1.602176565e-12_rk	! erg

real(kind=rk), parameter :: G_const = 6.67384e-8_rk
real(kind=rk), parameter :: pc = 3.08568e18_rk

End Module constants 
