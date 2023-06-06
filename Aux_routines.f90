Module auxiliary  

  use Realkind

  implicit none
  
    Contains

    
  
    Subroutine Findvalue_1dim_v3_nearestindex(x,xvec,f,value,i_x)

! ! 	v2: based on Findvalue_1dim: now returns the index kx1 for which xvec(kx1) is right below x.
!!	v3 : if x out of bounds, value is set to the value at the boundary instead of 0
!!	Findvalue_1dim_v3_nearestindex: based on v3: returns kx1 as the index of the nearest gridpoint

      implicit none

      real(kind=rk), dimension(:), intent(in) :: f 
      real(kind=rk), dimension(:), intent(in) :: xvec
      real(kind=rk), intent(out) :: value
      real(kind=rk) :: dx1, dx2, cfx, x
      integer :: i, imin, imax, kx1, kx2, n_x
      integer, intent(out) :: i_x
      logical :: cycl

    !       n_x = size(f)
        n_x = size(xvec)

        if (n_x.ne.(size(f))) then
          print *,'Dimesnions do not match'
          print *,'n_x, size(x)', n_x, size(xvec)
          print *,'size(f)', size(f)
          stop
        end if

        i = n_x/2
        imin = 1
        imax = n_x

      ! 	if ((x.lt.xvec(1)).or.(x.gt.xvec(n_x))) then
        if (x.lt.xvec(1)) then
        
      ! 		     print *,'Findvalue_1dim_v3: out of bounds;x,xvec(1),xvec(n_x)', x,xvec(1),xvec(n_x)	!!! NBNBNBNBNBNBNB: only for test, temporary
      ! 		     stop											!!! NBNBNBNBNBNBNB: only for test, temporary
      ! 	value = 0.0_rk
          value = f(1)
          i_x = 1
          return
        end if
        if (x.gt.xvec(n_x)) then
      ! 		      print *,'Findvalue_1dim_v3: out of bounds;x,xvec(1),xvec(n_x)', x,xvec(1),xvec(n_x)	!!! NBNBNBNBNBNBNB: only for test, temporary
      ! 		     stop											!!! NBNBNBNBNBNBNB: only for test, temporary
        
      ! 	value = 0.0_rk
          value = f(n_x)
          i_x = n_x
          return
        end if

        cycl = .TRUE.
        do while (cycl)
          if (x.gt.xvec(i)) then
            imin = i
            i = i + (imax - i)/2
          else if (x.lt.xvec(i)) then
            imax = i
            i = i - (i - imin)/2 
          else if (x.eq.xvec(i)) then
            imin = i
            imax = i+1
          else
            print *,'Findvalue_1dim_v3_nearestindex ERROR; x, xvec(i)', x, xvec(i)
            print *,'i, imin, imax', i, imin, imax
            print *,'xvec(:)', xvec(:)
            stop
          end if
          if ((imax-imin).eq.1) then
            kx1 = imin
            kx2 = imax
            cycl = .FALSE.
          end if
        end do

        cfx = (x - xvec(kx1))/(xvec(kx2) - xvec(kx1))

        ! value = M(kt1,kx1) + cfx*(M(kt1,kx2) - M(kt1,kx1)) + cft*(M(kt2,kx1) - M(kt1,kx1)) + cft*cfx*(M(kt2,kx2) + M(kt1,kx1) - M(kt2,kx1) - M(kt1,kx2))
        value = f(kx1) + cfx*(f(kx2) - f(kx1))
        
        if (abs(cfx).lt.0.5_rk) then
          i_x = kx1
        else
          i_x = kx2
        end if

        if((x - xvec(kx1))/(xvec(kx2) - xvec(kx1)) .gt. 1.0_rk) then
          print *,'Error in findvalue_1dim'
          stop
        end if

        if((x - xvec(kx1))/(xvec(kx2) - xvec(kx1)) .lt. 0.0_rk) then
          print *,'Error in findvalue_1dim'
          stop
        end if

        !!!***********************************************************************************************************

    End Subroutine Findvalue_1dim_v3_nearestindex
    
  

    Subroutine Findvalue_1dim_v3(x,xvec,f,value,kx1)
        ! (lnz,lnz_vec,f,w)

        ! ! 	v2: based on Findvalue_1dim: now returns the index kx1 for which xvec(kx1) is right below x.
        !!	v3 : if x out of bounds, value is set to the value at the boundary instead of 0

      implicit none

      real(kind=rk), dimension(:), intent(in) :: f 
      real(kind=rk), dimension(:), intent(in) :: xvec
      real(kind=rk), intent(out) :: value
      real(kind=rk) :: dx1, dx2, cfx, x
      integer :: i, imin, imax, kx2, n_x
      integer, intent(out) :: kx1
      logical :: cycl
      
        n_x = size(xvec)

        if (n_x.ne.(size(f))) then
          print *,'Dimesnions do not match'
          print *,'n_x, size(x)', n_x, size(xvec)
          print *,'size(f)', size(f)
          stop
        end if

        i = n_x/2
        imin = 1
        imax = n_x

        if (x.lt.xvec(1)) then
      ! 	value = 0.0_rk
          value = f(1)
          kx1 = 1
          return
        end if
        if (x.gt.xvec(n_x)) then      
      ! 	value = 0.0_rk
          value = f(n_x)
          kx1 = n_x
          return
        end if

        cycl = .TRUE.
        do while (cycl)
          if (x.gt.xvec(i)) then
            imin = i
            i = i + (imax - i)/2
          else if (x.lt.xvec(i)) then
            imax = i
            i = i - (i - imin)/2 
          else if (x.eq.xvec(i)) then
            imin = i
            imax = i+1
          end if
          if ((imax-imin).eq.1) then
            kx1 = imin
            kx2 = imax
            cycl = .FALSE.
          end if
        end do

        cfx = (x - xvec(kx1))/(xvec(kx2) - xvec(kx1))

        ! value = M(kt1,kx1) + cfx*(M(kt1,kx2) - M(kt1,kx1)) + cft*(M(kt2,kx1) - M(kt1,kx1)) + cft*cfx*(M(kt2,kx2) + M(kt1,kx1) - M(kt2,kx1) - M(kt1,kx2))
        value = f(kx1) + cfx*(f(kx2) - f(kx1))

        if((x - xvec(kx1))/(xvec(kx2) - xvec(kx1)) .gt. 1.0_rk) then
          print *,'Error in findvalue_1dim'
          stop
        end if

        if((x - xvec(kx1))/(xvec(kx2) - xvec(kx1)) .lt. 0.0_rk) then
          print *,'Error in findvalue_1dim'
          stop
        end if
        
        !!!***********************************************************************************************************

    End Subroutine Findvalue_1dim_v3
    
        
    
       
    
    Subroutine Findvalue_1dim_v4(x,xvec,f,value,kx1)
            !!! (lnz,lnz_vec,f,w)

! ! 	v2: based on Findvalue_1dim: now returns the index kx1 for which xvec(kx1) is right below x.
!!	v3 : if x out of bounds, value is set to the value at the boundary instead of 0
!!	v4 : if x out of bounds, value and index set to 0

      implicit none

      real(kind=rk), dimension(:), intent(in) :: f 
      real(kind=rk), dimension(:), intent(in) :: xvec
      real(kind=rk), intent(out) :: value
      real(kind=rk) :: dx1, dx2, cfx, x
      integer :: i, imin, imax, kx2, n_x
      integer, intent(out) :: kx1
      logical :: cycl

        n_x = size(f)

        if (n_x.ne.(size(f))) then
          print *,'Dimesnions do not match'
          print *,'n_x, size(x)', n_x, size(xvec)
          print *,'size(f)', size(f)
          stop
        end if

        i = n_x/2
        imin = 1
        imax = n_x

    ! 	if ((x.lt.xvec(1)).or.(x.gt.xvec(n_x))) then
        if (x.lt.xvec(1)) then
          value = 0.0_rk
          kx1 = 0
    ! 	value = f(1)
    ! 	kx1 = 1
          return
        end if
        if (x.gt.xvec(n_x)) then
          value = 0.0_rk
          kx1 = 0
    ! 	value = f(n_x)
    ! 	kx1 = n_x
          return
        end if

        cycl = .TRUE.
        do while (cycl)
          if (x.gt.xvec(i)) then
            imin = i
            i = i + (imax - i)/2
          else if (x.lt.xvec(i)) then
            imax = i
            i = i - (i - imin)/2 
          else if (x.eq.xvec(i)) then
            imin = i
            imax = i+1
          end if
          if ((imax-imin).eq.1) then
            kx1 = imin
            kx2 = imax
            cycl = .FALSE.
          end if
        end do

        cfx = (x - xvec(kx1))/(xvec(kx2) - xvec(kx1))

        ! value = M(kt1,kx1) + cfx*(M(kt1,kx2) - M(kt1,kx1)) + cft*(M(kt2,kx1) - M(kt1,kx1)) + cft*cfx*(M(kt2,kx2) + M(kt1,kx1) - M(kt2,kx1) - M(kt1,kx2))
        value = f(kx1) + cfx*(f(kx2) - f(kx1))

        if((x - xvec(kx1))/(xvec(kx2) - xvec(kx1)) .gt. 1.0_rk) then
          print *,'Error in findvalue_1dim'
          stop
        end if

        if((x - xvec(kx1))/(xvec(kx2) - xvec(kx1)) .lt. 0.0_rk) then
          print *,'Error in findvalue_1dim'
          stop
        end if

        !!!***********************************************************************************************************

    End Subroutine Findvalue_1dim_v4
  
    
    
    Subroutine Bin(x,x_vec,weight,Nx)
    
      use realkind
      
      implicit none
      
      real(kind=rk), intent(in) :: x, weight
      real(kind=rk), dimension(:), intent(in) :: x_vec
      real(kind=rk), dimension(:) :: Nx
      real(kind=rk) :: x_min, x_max
      integer :: i, n

      n = size(x_vec)
      
      do i=1,n-1
        x_min = x_vec(i)
        x_max = x_vec(i+1)
        if ((x.ge.x_min).and.(x.lt.x_max)) then
          Nx(i) = Nx(i) + weight ! 1.0_rk
        return
        end if
      end do
    
    End Subroutine Bin
        
    
    

    Subroutine Setup_grid(x_vec,dx,x_min,x_max)
      
        use realkind
      
        implicit none
        
        real(kind=rk), dimension(:), intent(out) :: x_vec
        real(kind=rk), intent(out) :: dx
        real(kind=rk), intent(in) :: x_min, x_max
        integer :: i, n
        
        n = size(x_vec)
        
        dx = (x_max - x_min)/(dble(n) - 1.0_rk)
        do i=1,n
          x_vec(i) = x_min + (dble(i) - 1.0_rk)*dx
        end do
        
    End Subroutine Setup_grid  
    
    
    
    
    
    
    
    
    
     
  
    Subroutine MC_append_escape(Photon)
    
      use data_type_def, ONLY: Photon_param
      use Photon_distribution_escape, ONLY: Photons_escape, i_exist_max_escape
    
      implicit none
      
      type(Photon_param), intent(in) :: Photon
      type(Photon_param), dimension(:), allocatable :: Photons_tmp
!       real(kind=rk) :: 
      integer :: n_size, n_size_neu
      
        if (allocated(Photons_escape)) then
          n_size = size(Photons_escape)
        else
          n_size = 1000  !! some arbitrary initial size
          allocate(Photons_escape(n_size))
          i_exist_max_escape = 0
          
          Photons_escape(:)%x = 0.0_rk
          Photons_escape(:)%exist = .FALSE.
          Photons_escape(:)%count_scatter = 0   !!! diagnostics only
        end if
        
        if ((i_exist_max_escape+1).gt.n_size) then
          allocate(Photons_tmp(n_size))
          Photons_tmp(1:n_size) = Photons_escape(1:n_size)
          deallocate(Photons_escape)
          n_size_neu = n_size*2
          allocate(Photons_escape(n_size_neu))
          Photons_escape(:)%exist = .FALSE.
          Photons_escape(:)%count_scatter = 0   !!! diagnostics only
          Photons_escape(1:n_size) = Photons_tmp(1:n_size)
          deallocate(Photons_tmp)
        end if
        
        Photons_escape(i_exist_max_escape+1) = Photon
        i_exist_max_escape = i_exist_max_escape + 1
        
      
    End Subroutine MC_append_escape
    
    

    Subroutine MC_photon_cleanup
    
      use Photon_distribution, ONLY: Photons_vec
      use data_type_def, ONLY: Photon_param
    
      implicit none
      
      type(Photon_param), dimension(:), allocatable :: Photons_tmp
      integer :: n_ph_bin, i, i_exist
    
        !!!#### Update photon array, remove dead souls #####
        if (allocated(Photons_vec)) then ! 16.09.19      
          n_ph_bin = size(Photons_vec)
          allocate(Photons_tmp(n_ph_bin))
          i_exist = 0
          do i = 1,n_ph_bin
            if (Photons_vec(i)%exist) then
              i_exist = i_exist + 1
              Photons_tmp(i_exist) = Photons_vec(i)
            end if
          end do
          deallocate(Photons_vec)
          allocate(Photons_vec(1:i_exist))
          Photons_vec(1:i_exist) = Photons_tmp(1:i_exist)
          deallocate(Photons_tmp)
        end if
        !!!##################################################
    
    End Subroutine MC_photon_cleanup
    
    
    Subroutine MC_photon_cull
    
      use MC_photon_params, ONLY: n_ph_cull
      use Photon_distribution, ONLY: Photons_vec
      use data_type_def, ONLY: Photon_param
    
      implicit none
      
        type(Photon_param), dimension(:), allocatable :: Photons_tmp
        real(kind=rk) :: frac_keep, rnd1
        integer :: n_MC, i, i_exist
        logical :: sw_keep
        
      
        if (allocated(Photons_vec)) then
          n_MC = size(Photons_vec)
          
          frac_keep = 1.0_rk
          if ((n_ph_cull.gt.0).and.(n_MC.gt.(n_ph_cull*2))) then
            frac_keep = min(1.0_rk,dble(n_ph_cull)/dble(n_MC))
            print *,'Reducing Photons_therm size to n_ph_cull; n_MC, n_ph_cull', n_MC, n_ph_cull
!           end if  !!! MOVED 'end if' below, to do nothing if conditions for culling are not satisfied
          
            allocate(Photons_tmp(n_MC))
  !           Photons_tmp(:) = .FALSE.
            i_exist = 0
            do i = 1,n_MC
              if (Photons_vec(i)%exist) then
              
                sw_keep = .TRUE.
                if (frac_keep.lt.1.0_rk) then
                  call RANDOM_NUMBER(rnd1)
                  if (rnd1.lt.frac_keep) then
                    sw_keep = .TRUE.
                  else
                    sw_keep = .FALSE.
                  end if
                end if
                
                if (sw_keep) then
                  i_exist = i_exist + 1
                  Photons_tmp(i_exist) = Photons_vec(i)
                end if
            
              end if    !!! if (Photons_vec(i)%exist) then
            end do  !!! do i = 1,n_MC
            
            deallocate(Photons_vec)
            allocate(Photons_vec(1:i_exist))
            Photons_vec(1:i_exist) = Photons_tmp(1:i_exist)
            deallocate(Photons_tmp)
            
            Photons_vec(1:i_exist)%weight = Photons_vec(1:i_exist)%weight/frac_keep
            
          
          end if
          
        end if  !! if (allocated(Photons_vec)) then
    
    End Subroutine MC_photon_cull
    
    
    Subroutine MC_photon_split(ct,sw_force,n_cf)
    
      use Photon_distribution, ONLY: Photons_vec
      use data_type_def, ONLY: Photon_param
      use Init_printout, ONLY: FF
      
      use MC_photon_params, ONLY: n_MC_ph, w_fiducial      
    
      implicit none
    
      real(kind=rk), intent(in) :: ct
      type(Photon_param), dimension(:), allocatable :: Photons_tmp
      logical, intent(in) :: sw_force   !! Force a split?
      integer, intent(in) :: n_cf
      
      integer :: n_MC, i_cf, i_min, i_max  !! , n_cf
      integer, save :: i_save = 0
      
      logical :: sw_split
      
        4 format(e14.7,' ',e14.7,' ',e14.7,' ',e14.7)
      
        call MC_photon_cleanup  !! remove dead souls
        
! !         n_cf = 2
      
        !!!## Split if number of MC photons drops too low
        sw_split = .FALSE.
        if (allocated(Photons_vec)) then
          n_MC = size(Photons_vec)
          if (n_MC.lt.(n_MC_ph/n_cf)) then   !!! Split photons if decrease below n_cf times the fiducial number
            sw_split = .TRUE.
          end if
        end if
      
      
        sw_split = sw_force.OR.sw_split    !!! is sw_force = .TRUE., then always split
      
        if (sw_split) then  !!! Do the split
          allocate(Photons_tmp(n_cf*n_MC))
          Photons_tmp(:)%exist = .FALSE.
          
          !!!### Create n_cf times more photons than before, adjust weights accordingly ####
          do i_cf = 1,n_cf
            i_min = (i_cf-1)*n_MC + 1
            i_max = i_cf*n_MC
            Photons_tmp(i_min:i_max) = Photons_vec(1:n_MC)
          end do
          Photons_tmp(:)%weight = Photons_tmp(:)%weight/dble(n_cf)
          
          w_fiducial = w_fiducial/dble(n_cf)   !!! Also reduce the fiducial weight so that more photons are emitted
          !!!###############################################################################
      
          deallocate(Photons_vec)
          allocate(Photons_vec(n_cf*n_MC))
          Photons_vec = Photons_tmp
          deallocate(Photons_tmp)
          
          if (i_save.eq.0) then
            open(500,FILE='results/Misc/Photon_splitting_' // trim(adjustl(FF)) // '.dat')
          else
            open(500,POSITION='APPEND',FILE='results/Misc/Photon_splitting_' // trim(adjustl(FF)) // '.dat')
          end if
          write(500,*) '############ Splitting photons ###############'
          write(500,*) 'Count'
          write(500,*) i_save + 1
          write(500,*) 'ct'
          write(500,*) ct
          write(500,*) 'n_MC'
          write(500,*) n_MC
          write(500,*) 'New n_MC'
          write(500,*) n_cf*n_MC
          write(500,*) 'n_MC_ph'
          write(500,*) n_MC_ph
          close(500)
          
          i_save = i_save + 1   
        end if  !! if (sw_split) then
      
    End Subroutine MC_photon_split
    
    
    
    Subroutine MC_target_number(domain_size,domain_size_in)
    
      use MC_photon_params, ONLY: n_MC_ph, n_MC_ph_0, n_ph_cull
    
      implicit none
    
      real(kind=rk), intent(in) :: domain_size, domain_size_in
    
        n_MC_ph = n_MC_ph_0*(domain_size/domain_size_in)**1.5_rk  !! Adjust target number of MC photons (power somewhat arbitrary)
        n_MC_ph = min(n_MC_ph,30*n_MC_ph_0)   !!! Limit how much n_MC_ph can grow
        
        
        n_ph_cull = 100*n_MC_ph   !! 
    
    End Subroutine MC_target_number
    
    
End Module auxiliary 
