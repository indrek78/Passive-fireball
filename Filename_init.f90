Module Init_printout

  use Realkind
  
  implicit none

  character(LEN=100) :: FF

  Contains
  
    Subroutine Read_Filename
      implicit none
      
        FF = 'test_1'  !! 1st test with initial powerlaw
        FF = 'test_2'  !!! steeper PL: alpha_x = 1
        FF = 'test_3'  !! longer sim: ct_max = ct_in*10.0 (3.0)
        
        FF = 'test_4'
        FF = 'test_5'  !!! !st working version with Compt and brems
        FF = 'test_5a'
!         FF = 'test_1a'

        FF = 'test_10'   !! 5a, but with "more" relativistically correct (except in the initialization phase) (08.05.23)
        
        FF = 'test_10a'  !! Introduced MC_photon_cleanup, MC_photon_cull

        FF = 'test_11'   !! 11, but R_sph_0 = 1.0e11 (1.0e12), n_MC_ph = 300 (3000), n_ph_cull = 100*n_MC_ph; n_t_nominal = 300 (30); [ct_max = ct_in*10.0_rk];  PROBLEM, energy non-conser at the beginning!  INTERRUPTED
        
        !!! Discovered ERROR: brems photon angles were not drawn properly
        
        FF = 'test_11a'  !! 11, but in Evolve_photon_diffusion: eps = 0.1 (0.5), and l_rad = 1.0e4 (1.0e1); ... 11a somewhat of a misnomer
        
        !!!### Introduced c_Tambov to artificially increase heat capacity ########
        
        FF = 'test_11b' !! 11a, but c_Tambov = 10
        
        !!!######## more printouts ### 
        
        FF = 'test_12'  !! 11b, but l_rad = 1.0e5, c_Tambov = 100; [n_t_nominal = 300, ct_max = ct_in*10, R_sph_0 = 1e11 cm]
                        !! Don't understand why adiabatic cooling scaling is not attained
        FF = 'test_13'  !! 12, but xi_min = 0.1 (0.5), l_rad = 5e5 (1e5)
        
        !!!## 10.05.23
        !!!## Added option for the case when formally less than one photon would be emitted; probably renders previous runs incorrect, because formally <1 photons emitted quite frequently, which previously led to no emission
        !!!## introduced chk_within
        !!!## Adjusted brems cooling for the range within which photons are actually emitted in Emit_thermal_brems_Sphere_homol
                
        FF = 'test_13a'
        
        !!!###### 11.05.23 #############
        !!!###### Longer runs: ct_max = ct_in*300 (ct_in*10) ###########
        !!!###### For ejecta: uniform grid in ejecta mass, rather than radius ######
        
        FF = 'test_14'
        !!FF = 'test_14a'  !! 14, but R_sph_0 = 3e11 (1e11) and c_Tambov = 40 (lower crashes for kT<0; somehow 14 with smaller starting radius did not.) !! CRASHED, wrong params, OVERWRITTEN
        FF = 'test_14a'   !! 14, but R_sph_0 = 3e11 (1e11), kT_in = 1.0e-3 --> 1.0e-3/3.0 = 3.33e-4 (accounting for adiab cooling betw 1e11 and 3e11), l_rad = 5e5 --> 5e5/3^3 = 1.85e4, so that total energy density would remain the same; [c_Tambov = 10]
        
        FF = 'test_14b'   !! 14a, but with en. conserv prinouts, varying target number of photons, and ct_max = ct_in*100 (300), n_t_nominal = 100 (300)
                          !! Result: looks good, energy conserved, but crashed the terminal (memory? doubtful since only 2e6 photons at the end)
                          
        FF = 'test_14c'  !! repeat of 14b, but with sw_reset_Photons_escape = .TRUE. (escaping photons reset each step, keeping only binned qties); also limit n_MC_ph (maximally) to 30*n_MC_ph_0   
        
                        !!! NBNBNBNB: setting sw_reset_Photons_escape = .TRUE. causes energy conservation to appear violated in Printout_Energy_conserv and python notebook (since Photons_escape is reset each step)
        
        
        !!!### Looks like initial n_ph/n_e is not terribly different from unity, i.e. initial temperature could be much higher than used here ###
        
        !!! Increase temperature
        
        FF = 'test_15'  !! 14c, but kT_in = 3.33e-3 (3.33e-4), and c_Tambov = 1 (10)
        
        FF = 'test_16'  !! 14, but kT_in = 3.33e-2 (3.33e-3); also, in Printout_timedep_spectrum,  broadened range of binning up to 1e2 (1e-2) and increased the number of bins to n_x_bin = 40 (20); [c_Tambov = 1]
                        !! increased n_kT = 60 (30) in Cumulative_distributions  
                        
        !!!### Setting alpha_r = 0 (1.1)
        
        FF = 'test_20' !!! 16, but alpha_r = 0 (1.1) 
        
        FF = 'test_30' !!! 20, but M_ej_norm = 3e-6 (1e-6), and n_MC_ph_0 = 10000 (30000); very slow
        
        FF = 'test_30a' !!! 30, but n_MC_ph_0 = 1000 (10000), c_Tambov = 3 (1); Ignoring ERROR 'Initialize_Radiation_Passive_fireball_Wien: self%c_Tambov.ge.(0.1_rk*n_ph/n_tot); self%c_Tambov, n_ph/n_tot: 3, 16' for this simul only
                        !!! Derived params.: E_rad[0] = 9.5e+45, E_kin = M_ej*v_ej^2/2 = 2.7e46 (not actually correct, v changes within ejecta; 1.6e46 if alpha_r = 0, xi_min -> 0)
                        
        FF = 'test_40'  !!! 30a, but v_ej = 1.5e-1 (1.0e-1); [n_MC_ph_0 = 1000, c_Tambov = 3]; Ignoring ERROR 'Initialize_Radiation_Passive_fireball_Wien: self%c_Tambov.ge.(0.1_rk*n_ph/n_tot); self%c_Tambov, n_ph/n_tot: 3, 16' for this simul only
        
        !!!## Do we need smaller initial radii (Itai, Metzger: X-ray spectrum looks thermal, with radius comparable to solar)?
        
        !!!### Must one reintroduce a sizable inner radius? since the driver (star) is of comparable size?
        
!         FF = 'del'
        
    End Subroutine Read_Filename
  
  
  
End Module Init_printout
