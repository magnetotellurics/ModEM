    module sg_vector_mg !
      ! created by Cherevatova (Jan, 2012)
      ! This module creates(allocates) cvectors/rvectors on staggered multi-grid
      ! Deallocates  and does some basic algebraic operations.
      ! Old routines, from sg_cvector are used

      ! Derived cvector_mg and rvector_mg data types are defined in sg_vector
      ! module (after cvector, rvector respectively)

    use math_constants        ! math/ physics constants
    use utilities
    use griddef
    use sg_vector
    implicit none

      ! Important - overloading the '=' assignment
      interface assignment (=)
        module procedure copy_cvector_mg
        module procedure c2mg ! copies cvector to cvector_mg
        module procedure mg2c ! cvector_mg to cvector
      end interface

      ! Generic interfaces are done through subroutines
      ! creates edge/ face nodes
      interface create
        module procedure create_rvector_mg
        module procedure create_cvector_mg
      end interface

      ! deallocates the edge/ face nodes
      interface deall
        module procedure deall_rvector_mg
        module procedure deall_cvector_mg
      end interface

      ! set the values to zero
      interface zero
        module procedure zero_cvector_mg
      end interface

      ! scalar value multiplies the edge/ face nodes
      interface scMult
        module procedure scMult_cvector_mg
        module procedure scMultReal_cvector_mg
      end interface

      INTERFACE scMultAdd
         module procedure scMultAdd_cvector_mg
      END INTERFACE

      interface linComb
        module procedure linComb_cvector_mg
      end interface

        ! adds the edge/ face nodes
        interface add
          module procedure add_cvector_mg
          module procedure add_rvector_mg
        end interface
        ! pointwise vector (two vector data types) multiplication of edge/ face
        ! nodes
        ! and pointwise real-complex (mixed) multiplication of edge/ face nodes
        ! Both are vector data types

        interface diagMult
          module procedure diagMult_cvector_mg
          module procedure diagMult_rvector_mg
          module procedure diagMult_crvector_mg
          module procedure diagMult_rcvector_mg
        end interface

        ! subtracts the edge/ face nodes
        interface subtract
          module procedure subtract_rvector_mg
          module procedure subtract_cvector_mg
        end interface

        interface dotProd
          module procedure dotProd_cvector_mg_f
        end interface

        interface dotProd_noConj
          module procedure dotProd_noConj_cvector_mg_f
        end interface

        interface  diagDiv
          module procedure diagDiv_rcvector_mg
          module procedure diagDiv_crvector_mg
        end interface

          ! overload some intrinsic functions for complex numbers
          INTERFACE conjg
             MODULE PROCEDURE conjg_cvector_mg_f
          END INTERFACE

          INTERFACE cmplx
             MODULE PROCEDURE cmplx_rvector_mg_f
          END INTERFACE

          INTERFACE real
             MODULE PROCEDURE real_cvector_mg_f
          END INTERFACE

          INTERFACE imag
             MODULE PROCEDURE imag_cvector_mg_f
          END INTERFACE


        public  :: create_rvector_mg, create_cvector_mg, &
                   deall_rvector_mg, deall_cvector_mg, &
                   copy_cvector_mg, c2mg, mg2c, &
                   zero_cvector_mg, &
                   subtract_rvector_mg, subtract_cvector_mg, &
                   scMult_cvector_mg, &
                   add_cvector_mg,add_rvector_mg, &
                   diagMult_cvector_mg, diagMult_crvector_mg, &
                   diagMult_rvector_mg, diagMult_rcvector_mg, &
                   linComb_cvector_mg, &
                   dotProd_cvector_mg_f,  dotProd_noConj_cvector_mg_f, &
                   diagDiv_rcvector_mg, diagDiv_crvector_mg, &
                   conjg_cvector_mg_f, cmplx_rvector_mg_f, real_cvector_mg_f, imag_cvector_mg_f, &
                   plotCvector_mg

    ! *******************************************************************************
    ! this type defines cvectors on multi-grid

    type :: cvector_mg
       ! number of subgrids
       integer  :: mgridSize
       ! cvectors for subgrids
       type(cvector), pointer :: cvArray(:)
       ! coarseness
       integer, allocatable  :: coarseness(:)
       ! allocated:  .true.  x, y, z arrays have been allocated
       logical  :: allocated = .false.
       ! store the intention of the use in a character string defined
       ! in GridDef as a parameter: EDGE or FACE are two possibilities

       ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
       ! (probably will not be needed in the future when compilers will support
       ! ISO/IEC 15581 - the "allocatable array extension")
       logical                                        :: temporary = .false.

       character (len=80)   :: gridType

       ! pointer to parent grid
       type (grid_t), pointer  :: grid
    end type cvector_mg

    ! ******************************************************************************
    ! this type defines multigrids rvector
    type rvector_mg
      ! number of subgrids
      integer  :: mgridSize
      ! store the intention of the use in a character string defined
      ! in GridDef as a parameter: EDGE or FACE are two possibilities
      character (len=80)  :: gridType
      ! rvectors for subgrids
      type(rvector), pointer :: rvArray(:)
      ! coarseness
      integer, allocatable  :: coarseness(:)
      ! allocated:  .true.  x, y, z arrays have been allocated
      logical  :: allocated = .false.
      ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
      ! (probably will not be needed in the future when compilers will support
      ! ISO/IEC 15581 - the "allocatable array extension")
      logical                                        :: temporary = .false.

       ! pointer to parent grid
       type (grid_t), pointer  :: grid

    end type rvector_mg


    Contains

      ! *******************************************************************************
      ! allocates multigrids rvector
      subroutine create_rvector_mg (mgrid, e, gridType)

      implicit none
        type(grid_t), target, intent(in)  :: mgrid
        ! the grid for which an edge/ face node field is being initialized
        type (rvector_mg), intent(inout)  :: e
        character (len=80), intent(in)    :: gridType
        !local
        integer  :: status, imgrid


        ! First deallocate anything, that's allocated
        call deall(e) !deall_rvector_mg
        ! Set pointer
        e%grid => mgrid
       ! allocate memory for rvector arrays
        allocate(e%rvArray(mgrid%mgridSize), STAT = status)
        ! allocate rvectors
        do imgrid = 1, mgrid%mgridSize
          call create(mgrid%gridArray(imgrid),e%rvArray(imgrid),gridType) !create_rvector
        enddo
        ! allocate memory for coarseness array
        allocate(e%coarseness(mgrid%mgridSize), STAT = status)
        ! number of subgrids
        e%mgridSize = mgrid%mgridSize
        ! gridType
        e%gridType = gridType
        ! e%allocated will be true if all allocations succeed
        e%allocated = .true.

        ! copy coarseness from multigrid
        e%coarseness = mgrid%coarseness

      end subroutine create_rvector_mg

      ! ******************************************************************************
      ! allocates cvector for multigrid
      subroutine create_cvector_mg (mgrid, e, gridType)

      implicit none
        type(grid_t), target, intent(in)  :: mgrid
        ! the grid for which an edge/ face node field is being initialized
        type (cvector_mg), intent(inout) :: e
        character (len=80), intent(in) :: gridType
        !local
        integer  :: status, imgrid

        ! First deallocate anything, that's allocated
        call deall(e)!deall_cvector_mg
        ! Set pointer
        e%grid => mgrid
        ! allocate memory for cvector arrays
        allocate(e%cvArray(mgrid%mgridSize), STAT = status)
        ! allocate cvectors
        do imgrid = 1, mgrid%mgridSize
          call create(mgrid%gridArray(imgrid),e%cvArray(imgrid),gridType) !create_cvector
        enddo
        ! allocate memory for coarseness array
        allocate(e%coarseness(mgrid%mgridSize), STAT = status)
        ! number of subgrids
        e%mgridSize = mgrid%mgridSize
        ! gridType
        e%gridType = gridType
        ! e%allocated will be true if all allocations succeed
        e%allocated = .true.

        ! copy coarseness from multigrid
        e%coarseness = mgrid%coarseness

      end subroutine create_cvector_mg

    ! *****************************************************************************
      subroutine deall_rvector_mg(e)

      implicit none
        type (rvector_mg)  :: e
        integer  :: status, imgrid


          if (e%allocated) then
        ! deallocate cvectors of subgrids
          do imgrid = 1, e%mgridSize
            call deall(e%rvArray(imgrid))!deall_rvector
          enddo
          deallocate(e%coarseness, STAT=status)
          deallocate(e%rvArray, STAT = status)
          end if
          if(associated(e%grid)) nullify(e%grid)
          e%mgridSize = 0
          E%gridType = ''
          E%allocated = .false.

      end  subroutine deall_rvector_mg

      ! *****************************************************************************
      subroutine deall_cvector_mg(e)

      implicit none
        type (cvector_mg)  :: e
        integer  :: status, imgrid


          if (e%allocated) then
        ! deallocate cvectors of subgrids
          do imgrid = 1, e%mgridSize
            call deall(e%cvArray(imgrid))!deall_cvector
          enddo
          deallocate(e%coarseness, STAT=status)
          deallocate(e%cvArray, STAT = status)
          end if

          if(associated(e%grid)) nullify(e%grid)


          e%mgridSize = 0
          E%gridType = ''
          E%allocated = .false.

      end  subroutine deall_cvector_mg

      ! *****************************************************************************
      ! copy cvector_mg to cvector_mg
      subroutine copy_cvector_mg(e2,e1)

      ! first argument is output
      implicit none
        type (cvector_mg), intent(in)  :: e1
        type (cvector_mg), intent(inout)  :: e2
        !local
        integer  :: imgrid

        if(.not.e1%allocated) then
          print *, 'RHS not allocated yet for copy_cvector_mg'
        endif

        if(.not.e2%allocated) then
          call create(e1%grid, e2, e1%gridType)
        endif

         if (e1%gridType == e2%gridType.and.e1%mgridSize == e2%mgridSize) then
            do imgrid = 1, e1%mgridSize
              if(e2%cvArray(imgrid)%nx ==  e1%cvArray(imgrid)%nx.and. &
                 e2%cvArray(imgrid)%ny ==  e1%cvArray(imgrid)%ny.and. &
                 e2%cvArray(imgrid)%nz ==  e1%cvArray(imgrid)%nz) then

                      e2%cvArray(imgrid)%x = e1%cvArray(imgrid)%x
                      e2%cvArray(imgrid)%y = e1%cvArray(imgrid)%y
                      e2%cvArray(imgrid)%z = e1%cvArray(imgrid)%z
                 e2%cvArray(imgrid)%gridType = e1%cvArray(imgrid)%gridType
                 e2%cvArray(imgrid)%grid => e1%cvArray(imgrid)%grid
              else
                print *, 'e1 and e2 are not the same size; copy_cvector_mg'
              endif
            enddo
            e2%gridType = e1%gridType
            e2%grid => e1%grid
            e2%coarseness = e1%coarseness
          else
            print *, 'not compatible usage for copy_cvector'
          endif

      end subroutine copy_cvector_mg
      ! *****************************************************************************
      ! copy fields from cvector to cvector_mg
      subroutine c2mg(e2,e1)

      implicit none

        type(cvector), intent(in)  :: e1       ! input cvector
        type(cvector_mg), intent(inout)  ::e2  ! output cvector

        ! local
        integer  :: imgrid,ifine, ix,iy,iz,izv,ic
        integer  :: nx,ny,nz,nzCum,ccoeff_current
        integer  :: errAll

        if (.not.e1%allocated)then
          print *, 'Error c2mg; e1 (cvector) is not allocated'
        endif

        if (.not.e2%allocated)then
          print *, 'Error c2mg; e2 (cvector_mg) is not allocated'
        endif

        ! which gridType and check
        if (e1%gridType == EDGE .and. e1%gridType == e2%gridType) then
          nzCum = 0
          do imgrid = 1, e2%mgridSize  ! Global loop on sub-grids
            nx = e2%cvArray(imgrid)%nx
            ny = e2%cvArray(imgrid)%ny
            nz = e2%cvArray(imgrid)%nz
            ccoeff_current = 2**e2%coarseness(imgrid)
            ! re-count x component
            do iz = 1, nz+1
               izv = iz + nzCum
              do iy = 1, ny+1
                do ix =1, nx
                  ! do ic = 1, ccoeff_current
                  ! average fields
                  ! does not make any difference in the final solution, compare with copy, but
                  ! leads to numerical unstability, therefore not used
                  !  e2%cvArray(imgrid)%x(ix,iy,iz) = e2%cvArray(imgrid)%x(ix,iy,iz) + &
                  !                e1%x(ccoeff_current*(ix-1)+ic,ccoeff_current*(iy-1)+1,izv)*e1%grid%dx(ccoeff_current*(ix-1)+ic)
                  ! enddo ! ic
                  !  e2%cvArray(imgrid)%x(ix,iy,iz) = e2%cvArray(imgrid)%x(ix,iy,iz)/e2%grid%gridArray(imgrid)%dx(ix)
                  ! copy fields
                    e2%cvArray(imgrid)%x(ix,iy,iz) = e1%x(ccoeff_current*(ix-1)+1,ccoeff_current*(iy-1)+1,izv)
                enddo  !ix
              enddo !iy
            ! re-count y component
              do ix = 1, nx+1
                do iy = 1, ny
!                   do ic = 1, ccoeff_current
!                    e2%cvArray(imgrid)%y(ix,iy,iz) = e2%cvArray(imgrid)%y(ix,iy,iz) +  &
!                                     e1%y(ccoeff_current*(ix-1)+1,ccoeff_current*(iy-1)+ic,izv)*e1%grid%dy(ccoeff_current*(iy-1)+ic)
!                   enddo !ic
!                    e2%cvArray(imgrid)%y(ix,iy,iz) = e2%cvArray(imgrid)%y(ix,iy,iz)/e2%grid%gridArray(imgrid)%dy(iy)
                   e2%cvArray(imgrid)%y(ix,iy,iz) = e1%y(ccoeff_current*(ix-1)+1,ccoeff_current*(iy-1)+1,izv)
                enddo !iy
              enddo !ix
            enddo ! iz
            ! re-count z component
            do iy =1, ny+1
              do ix= 1, nx+1
                do iz = 1, nz
                   izv = iz + nzCum
                    e2%cvArray(imgrid)%z(ix,iy,iz) = e1%z(ccoeff_current*(ix-1)+1,ccoeff_current*(ix-1)+1,izv)
               enddo !iz
              enddo !ix
            enddo !iy
            nzCum =  nzCum + nz
          enddo  ! Global loop over sub-grids

        else if(e1%gridType == FACE .and. e1%gridType == e2%gridType) then
          nzCum = 0
          do imgrid = 1, e2%mgridSize  ! Global loop on sub-grids
            nx = e2%cvArray(imgrid)%nx
            ny = e2%cvArray(imgrid)%ny
            nz = e2%cvArray(imgrid)%nz
            ccoeff_current = 2**e2%coarseness(imgrid)
            ! re-count x component
            do iz = 1, nz
               izv = iz + nzCum
              do iy = 1, ny
                do ix =1, nx+1
                  !do ic = 1, ccoeff_current
                  ! average fields
                  !  e2%cvArray(imgrid)%x(ix,iy,iz) = e2%cvArray(imgrid)%x(ix,iy,iz) + &
                  !            e1%x(ccoeff_current*(ix-1)+ic,ccoeff_current*(iy-1)+1,izv)*e1%grid%dx(ccoeff_current*(ix-1)+ic)
                  !enddo ! ic
                  !  e2%cvArray(imgrid)%x(ix,iy,iz) = e2%cvArray(imgrid)%x(ix,iy,iz)/e2%grid%gridArray(imgrid)%dx(ix)
                    e2%cvArray(imgrid)%x(ix,iy,iz) = e1%x(ccoeff_current*(ix-1)+ic,ccoeff_current*(iy-1)+1,izv)
                enddo  !ix
              enddo !iy
            ! re-count y component
              do ix = 1, nx
                do iy = 1, ny+1
                  ! do ic = 1, ccoeff_current
                  !   e2%cvArray(imgrid)%y(ix,iy,iz) = e2%cvArray(imgrid)%y(ix,iy,iz)  + &
                  !                   e1%y(ccoeff_current*(ix-1)+1,ccoeff_current*(iy-1)+ic,izv)*e1%grid%dy(ccoeff_current*(iy-1)+ic)
                  !enddo !ic
                  !  e2%cvArray(imgrid)%y(ix,iy,iz) = e2%cvArray(imgrid)%y(ix,iy,iz)/e2%grid%gridArray(imgrid)%dy(iy)
                    e2%cvArray(imgrid)%y(ix,iy,iz) = e1%y(ccoeff_current*(ix-1)+1,ccoeff_current*(iy-1)+ic,izv)
                enddo !iy
              enddo !ix
            enddo ! iz
            ! re-count z component
            do iy =1, ny
              do ix= 1, nx
                do iz = 1, nz+1
                   izv = iz + nzCum
                    e2%cvArray(imgrid)%z(ix,iy,iz) = e1%z(ccoeff_current*(ix-1)+1,ccoeff_current*(ix-1)+1,izv)
               enddo !iz
              enddo !ix
            enddo !iy
            nzCum =  nzCum + nz
          enddo  ! Global loop over sub-grids
        else
          print *, 'Error c2mg; cvector and cvector_mg not are the same gridType'
        endif

      end subroutine c2mg

      ! *****************************************************************************

      ! convert cvector_mg to cvector. Copy fields
      ! used only to plot fields

      subroutine mg2c(e2, e1)

      implicit none
        type(cvector_mg), intent(in)  :: e1  ! cvector_mg in
        type(cvector), intent(inout)  :: e2  ! cvector out

        ! local
        integer :: nzCum, nc, nx,ny,nz
        integer  :: imgrid,ix,iy,iz,izv, ic,icx,icy
        integer  :: status

        if (.not.e1%allocated)then
          print *, 'Error mg2c; e1 (cvector_mg) is not allocated'
        endif

        if (.not.e2%allocated)then
          print *, 'Error mg2c; e2 (cvector) is not allocated'
        endif

        ! check gridType
        if (e1%gridType == e2%gridType) then

            nzCum = 0
            do imgrid = 1, e1%mgridSize  ! Global loop over sub-grids
               nx = e1%cvArray(imgrid)%nx
               ny = e1%cvArray(imgrid)%ny
               nz = e1%cvArray(imgrid)%nz
               nc= 2**e1%coarseness(imgrid)
               do iz = 1, nz+1
                  izv = iz + nzCum
                  ! Ex components
                 do iy = 1, ny
                     do ix =1, nx
!                           e2%x(nc*(ix-1)+1,nc*(iy-1)+1,izv) = e1%cvArray(imgrid)%x(ix,iy,iz)
                       do icx = 1, nc
!                         do icy = 1, nc
                           e2%x(nc*(ix-1)+icx,nc*(iy-1)+1,izv) = e1%cvArray(imgrid)%x(ix,iy,iz)
!                           e2%x(nc*(ix-1)+icx,nc*(iy-1)+icy,izv) = e1%cvArray(imgrid)%x(ix,iy,iz)
!                         enddo ! icy
                       enddo  ! icx
                     enddo ! ix
                  enddo ! iy
                  ! Ey components
                 do ix = 1, nx
                     do iy =1, ny
                      !      e2%y(nc*(ix-1)+1,nc*(iy-1)+1,izv) = e1%cvArray(imgrid)%y(ix,iy,iz)
                       do icx =1, nc
!                         do icy = 1, nc
                           e2%y(nc*(ix-1)+icx,nc*(iy-1)+1,izv) = e1%cvArray(imgrid)%y(ix,iy,iz)
!                         enddo  ! icx
                       enddo  ! icy
                     enddo ! iy
                  enddo ! ix
               enddo !iz

              ! Ez components
              do iz = 1, nz
                  izv = iz + nzCum
                  do iy = 1, ny
                     do ix =1, nx
                     !     e2%z(nc*(ix-1)+1,nc*(iy-1)+1,izv)= e1%cvArray(imgrid)%z(ix,iy,iz)
                       do icx = 1, nc
!                         do icy = 1, nc
                           e2%z(nc*(ix-1)+icx,nc*(iy-1)+1,izv)= e1%cvArray(imgrid)%z(ix,iy,iz)
!                        enddo ! icx
                       enddo ! icy
                     enddo ! ix
                 enddo !iy
              enddo !iz
            nzCum = nzCum + nz
            enddo   ! Global loop over sub-grids

        else
          print *, 'Error mg2c; cvector and cvector_mg are not the same gridType'
        endif

      end subroutine mg2c

      ! *****************************************************************************
      subroutine random_rvector_mg(e,eps)

        implicit none
        type(rvector_mg), intent(inout)  :: e
        real(8), intent(in), optional    :: eps
        ! local
        integer  ::imgrid

        do imgrid = 1, e%mgridSize
            if (present(eps)) then
                call random_rvector(e%rvArray(imgrid),eps)
            else
                call random_rvector(e%rvArray(imgrid))
            endif
        enddo
     end subroutine random_rvector_mg

      ! *****************************************************************************
      subroutine random_cvector_mg(e,eps)

        implicit none
        type(cvector_mg), intent(inout)  :: e
        real(8), intent(in), optional    :: eps
        ! local
        integer  ::imgrid

        do imgrid = 1, e%mgridSize
            if (present(eps)) then
                call random_cvector(e%cvArray(imgrid),eps)
            else
                call random_cvector(e%cvArray(imgrid))
            endif
        enddo
     end subroutine random_cvector_mg


      ! *****************************************************************************
      subroutine zero_rvector_mg(e)

        implicit none
        type(rvector_mg), intent(inout)  :: e
        ! local
        integer  ::imgrid

        do imgrid = 1, e%mgridSize
            call zero(e%rvArray(imgrid)) !zero_rvector
       enddo
     end subroutine zero_rvector_mg

      ! *****************************************************************************
      subroutine zero_cvector_mg(e)

        implicit none
        type(cvector_mg), intent(inout)  :: e
        ! local
        integer  ::imgrid

        do imgrid = 1, e%mgridSize
            call zero(e%cvArray(imgrid)) !zero_cvector
       enddo
     end subroutine zero_cvector_mg

     ! ***************************************************************************
     subroutine scMult_cvector_mg(c, e1, e2)

       implicit none
       ! a complex scalar to be multiplied with
       complex(kind=prec), intent(in)  :: c
       type (cvector_mg), intent(in)   :: e1
       type (cvector_mg), intent(inout)  :: e2
        !local
        integer  :: imgrid

       ! check whether e1 and e2 have the same number of subgrids
         if (e1%mgridSize == e2%mgridSize) then
            do imgrid = 1, e1%mgridSize
              call scMult(c, e1%cvArray(imgrid),e2%cvArray(imgrid)) !scMult_cvector
            enddo

         else
           write (0, *) 'Error:scMult_cvector_mg: vectors not same subgrids size'
         endif
     end subroutine scMult_cvector_mg
     ! *****************************************************************************
     subroutine scMultReal_cvector_mg(c, e1, e2)

       implicit none
        real (kind=prec), intent(in)   :: c
        ! a real scalar to be multiplied with
        type (cvector_mg), intent(in)  :: e1
        type (cvector_mg), intent(inout)  :: e2
        ! local
        integer  :: imgrid

          if(.not.e1%allocated) then
            print *, 'RHS not allocated yet for scMultReal_cvector_mg'
            return
          endif

         ! check to see if LHS (E2) is active (allocated)
         if(.not.e2%allocated) then
           print *,'LHS was not allocated yet for scMultReal_cvector_mg'
         else

           if ((e1%gridType == e2%gridType).and.(e1%mgridSize == e2%mgridSize)) then

             do imgrid = 1, e1%mgridSize
                call scMult(c, e1%cvArray(imgrid), e2%cvArray(imgrid)) ! scMultReal_cvector
             enddo

           else
             print *, 'not compatible usage for scMultReal_cvector_mg'
           end if
         endif

      end subroutine scMultReal_cvector_mg ! scMultReal_cvector_mg
   ! *****************************************************************************
     subroutine linComb_cvector_mg(inc1, e1, inc2, e2, e3)

      implicit none
      !   input vectors
      type (cvector_mg), intent(in)  :: e1, e2
      !  input complex scalars
      complex (kind=prec), intent(in)  :: inc1, inc2
      ! lin comp cvector
      type (cvector_mg), intent(inout)          :: e3
      !local
      integer  :: imgrid

      ! check whether e1 and e2 have the same number of subgrids
        if (e1%mgridSize == e2%mgridSize) then
          e3%mgridSize = e1%mgridSize
          do imgrid = 1, e1%mgridSize
            call linComb(inc1, e1%cvArray(imgrid), inc2, &
                              e2%cvArray(imgrid), e3%cvArray(imgrid))!linComb_cvector
          enddo
        else
          write (0, *) 'Error:linComb_cvector_mg: vectors not same subgrids size'
        endif
     end subroutine linComb_cvector_mg
   ! *****************************************************************************************
     subroutine add_cvector_mg(e1, e2, e3)

      implicit none
        type (cvector_mg), intent(in)  :: e1, e2
        type (cvector_mg), intent(inout)  :: e3
        !local
        integer  :: imgrid

        ! check whether e1 and e2 have the same number of subgrids
          if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                           (e2%mgridSize == e3%mgridSize)) then

            do imgrid = 1, e1%mgridSize
              call add(e1%cvArray(imgrid), e2%cvArray(imgrid), e3%cvArray(imgrid)) !add_cvector
            enddo

          else
             write(0, *) 'Error:add_cvector_mg: vectors not same subgrids size'
          endif
     end subroutine add_cvector_mg
     ! *****************************************************************************************
     subroutine add_rvector_mg(e1, e2, e3)

        implicit none
        type (rvector_mg), intent(in)  :: e1, e2
        type (rvector_mg), intent(inout)  :: e3
        !local
        integer  :: imgrid

        ! check whether e1 and e2 have the same number of subgrids
          if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                           (e2%mgridSize == e3%mgridSize)) then

            do imgrid = 1, e1%mgridSize
              call add(e1%rvArray(imgrid), e2%rvArray(imgrid), e3%rvArray(imgrid))!add_rvector
            enddo

          else
             write(0, *) 'Error:add_rvector_mg: vectors not same subgrids size'
          endif
     end subroutine add_rvector_mg
     ! **********************************************************************************************
     subroutine diagMult_cvector_mg(e1, e2, e3)

        implicit none
        type (cvector_mg), intent(in)  :: e1, e2
        type (cvector_mg), intent(inout)  :: e3
        !local
        integer  :: imgrid

        ! check whether e1 and e2 have the same number of subgrids
          if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                           (e2%mgridSize == e3%mgridSize)) then

            do imgrid = 1, e1%mgridSize
              call diagMult(e1%cvArray(imgrid), e2%cvArray(imgrid), e3%cvArray(imgrid))!diagMult_cvector
            enddo

          else
             write(0, *) 'Error::diagMult_cvector_mg vectors not same subgrids size'
          endif
     end subroutine diagMult_cvector_mg
     ! *****************************************************************************************
     subroutine diagMult_rvector_mg(e1, e2, e3)

        implicit none
        type (rvector_mg), intent(in)  :: e1, e2
        type (rvector_mg), intent(inout)  :: e3
        !local
        integer  :: imgrid

        ! check whether e1 and e2 have the same number of subgrids
          if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                           (e2%mgridSize == e3%mgridSize)) then

            do imgrid = 1, e1%mgridSize
              call diagMult(e1%rvArray(imgrid), e2%rvArray(imgrid), e3%rvArray(imgrid)) !diagMult_rvector
            enddo

          else
             write(0, *) 'Error::diagMult_rvector_mg vectors not same subgrids size'
          endif
     end subroutine diagMult_rvector_mg
     ! *****************************************************************************************
     subroutine diagMult_crvector_mg(e1, e2, e3)

        implicit none
        type (cvector_mg), intent(in)  :: e1
        type (rvector_mg), intent(in)  :: e2
        type (cvector_mg), intent(inout)  :: e3
        !local
        integer  :: imgrid

        ! check whether e1 and e2 have the same number of subgrids
          if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                           (e2%mgridSize == e3%mgridSize)) then

            do imgrid = 1, e1%mgridSize
              call diagMult_crvector(e1%cvArray(imgrid), e2%rvArray(imgrid), e3%cvArray(imgrid))!diagMult_crvector
            enddo
          else
             write(0, *) 'Error::diagMult_crvector_mg vectors not same subgrids size'
          endif

      end subroutine diagMult_crvector_mg
      ! **********************************************************************************************
      subroutine diagMult_rcvector_mg(e1, e2, e3)

            implicit none
            type (rvector_mg), intent(in)  :: e1
            type (cvector_mg), intent(in)  :: e2
            type (cvector_mg), intent(inout)  :: e3
            ! local
            integer  :: imgrid

            ! check whether e1 and e2 have the same number of subgrids
              if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                               (e2%mgridSize == e3%mgridSize)) then

                do imgrid = 1, e1%mgridSize
                  call diagMult(e1%rvArray(imgrid), e2%cvArray(imgrid), e3%cvArray(imgrid))! diagMult_rcvector
                enddo

              else
                 write(0, *) 'Error::diagMult_rcvector vectors not same subgrids size'
              endif
      end subroutine diagMult_rcvector_mg
      ! ********************************************************************************************88
      subroutine subtract_rvector_mg(e1, e2, e3)

        implicit none
        type (rvector_mg), intent(in)  :: e1, e2
        type (rvector_mg), intent(inout)  :: e3
        ! local variables
        integer  :: imgrid

        ! check whether e1 and e2 have the same number of subgrids
          if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                           (e2%mgridSize == e3%mgridSize)) then
            do imgrid = 1, e1%mgridSize
               call subtract(e1%rvArray(imgrid), e2%rvArray(imgrid), e3%rvArray(imgrid)) !subtract_rvector
             enddo

          else
             write(0, *) 'Error::subtract_rvector_mg vectors not same subgrids size'
          endif
     end subroutine subtract_rvector_mg
     ! ********************************************************************************************88
     subroutine subtract_cvector_mg(e1, e2, e3)

        implicit none
        type (cvector_mg), intent(in)  :: e1, e2
        type (cvector_mg), intent(inout)  :: e3
        ! local variables
        integer  :: imgrid

        ! check whether e1 and e2 have the same number of subgrids
          if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                           (e2%mgridSize == e3%mgridSize)) then
            do imgrid = 1, e1%mgridSize
               call subtract(e1%cvArray(imgrid), e2%cvArray(imgrid), e3%cvArray(imgrid))!subtract_cvector
             enddo

          else
             write(0, *) 'Error::subtract_cvector_mg vectors not same subgrids size'
          endif
     end subroutine subtract_cvector_mg

     ! ***************************************************************************
     ! dotProd_cvector_mg computes dot product of two vectors stored
     ! as derived data type cvector_mg, returning a complex number
     function dotProd_cvector_mg_f(e1, e2) result(c)

          implicit none

            type (cvector_mg), intent(in)  :: e1, e2
            complex(kind=prec)  :: c
            ! local
            type (cvector_mg) :: e3
            integer  :: imgrid
            complex (kind=prec),allocatable  :: ctemp(:)

            c = R_ZERO
            if((.not.e1%allocated).or.(.not.e2%allocated)) then
              write(0,*) 'RHS not allocated yet for dotProdCC multi-grid'
              return
            endif
            e3 = e1
           ! check whether e1 and e2 have the same number of subgrids
           if (e1%mgridSize == e2%mgridSize) then
             allocate(ctemp(e1%mgridSize))

             do imgrid = 1, e1%mgridSize
               ! delete duplicate edges
               e3%cvarray(imgrid)%x(:,:,e1%cvarray(imgrid)%nz+1) = C_ZERO
               e3%cvarray(imgrid)%y(:,:,e1%cvarray(imgrid)%nz+1) = C_ZERO
               ctemp(imgrid) = dotProd(e3%cvArray(imgrid), e2%cvArray(imgrid)) !dotProd_cvector_f
             enddo
          else
            write(0,*) 'Error :: dotProd_cvector_mg_f: e1 and e2 not the same subgrids'
          end if
          c = sum(ctemp)
          deallocate(ctemp)
          call deall(e3)!deall_cvector_mg
      end function dotProd_cvector_mg_f

     ! **************************************************************************8
     ! dotProd_noConj_cvector_mg computes dot product of two vectors stored
     ! as derived data type cvector_mg, returning a complex number
     function dotProd_noConj_cvector_mg_f(e1, e2) result(c)

        implicit none
        type (cvector_mg), intent(in)  :: e1, e2
        complex(kind=prec)  :: c
        ! local
        integer  :: imgrid
        complex (kind=prec),allocatable  :: ctemp(:)

        ! check whether e1 and e2 have the same number of subgrids
         if (e1%mgridSize == e2%mgridSize) then
           allocate(ctemp(e1%mgridSize))
           do imgrid = 1, e1%mgridSize
             ctemp(imgrid) = dotProd_noConj(e1%cvArray(imgrid), e2%cvArray(imgrid))!dotProd_noConj_cvector_f
           enddo
         else
            write(0,*) 'Multigrid cvector 1 # of subgrids: ',e1%mgridSize
            write(0,*) 'Multigrid cvector 2 # of subgrids: ',e2%mgridSize
            write(0,*) 'Error :: dotProd_noConj_cvector_mg_f: e1 and e2 not the same subgrids'
         end if
         c = sum(ctemp)

          deallocate(ctemp)
     end function dotProd_noConj_cvector_mg_f
     ! ***************************************************************************8
     subroutine diagDiv_crvector_mg(e1,e2,e3)

         implicit none
         type (cvector_mg), intent(in)  :: e1
         type (rvector_mg), intent(in)  :: e2
         type (cvector_mg), intent(inout)  :: e3

         integer   :: imgrid

         if((.not.e1%allocated).or.(.not.e2%allocated)) then
           print*, 'RHS not allocated yet for diagDiv_crvector_mg'
           return
         endif

         ! check to see if LHS (e3) is active (allocated)
         if(.not.e3%allocated) then
           print *, 'LHS was not allocated for diagDiv_crvector_mg'
         else
           do imgrid = 1, e1%mgridSize ! loop on subgrids
             Call diagDiv(e1%cvArray(imgrid), e2%rvArray(imgrid), e3%cvArray(imgrid)) !diagDiv_crvector

           enddo
         endif

     end subroutine diagDiv_crvector_mg
     ! ******************************************************************************
     subroutine diagDiv_rcvector_mg(e1, e2, e3)

         implicit none
         type (rvector_mg), intent(in)  :: e1
         type (cvector_mg), intent(in)  :: e2
         type (cvector_mg), intent(inout)  :: e3

         integer   :: imgrid

         if((.not.e1%allocated).or.(.not.e2%allocated)) then
           print*, 'RHS not allocated yet for diagDiv_rcvector_mg'
           return
         endif

         ! check to see if LHS (e3) is active (allocated)
         if(.not.e3%allocated) then
           print *, 'LHS was not allocated for diagDiv_rcvector_mg'
         else
           do imgrid = 1, e1%mgridSize ! loop on subgrids
             Call diagDiv(e1%rvArray(imgrid), e2%cvArray(imgrid), e3%cvArray(imgrid))!diagDiv_rcvector
           enddo
         endif
     end subroutine diagDiv_rcvector_mg
     ! *****************************************************************************
      subroutine scMultAdd_cvector_mg(c, e1, e2)

       implicit none
            complex(kind=prec), intent(in)  :: c
            ! a complex scalar to be multiplied with
            type (cvector_mg), intent(in)   :: e1
            type (cvector_mg)  :: e2

        ! local variables
           integer  :: imgrid

            if(.not.e1%allocated) then
               print *, 'RHS not allocated yet for scMultAdd_cvector_mg'
               return
            endif

            ! check to see if LHS (e2) is active (allocated)
            if(.not.e2%allocated) then
               print *, 'LHS was not allocated for scMultAdd_cvector_mg'
            else

              if ((e1%gridType == e2%gridType)) then
                do imgrid = 1, e1%mgridSize  ! loop on subgrids
                  ! Check whether both vectors are of the same size
                  if((e1%cvArray(imgrid)%nx == e2%cvArray(imgrid)%nx).and.(e1%cvArray(imgrid)%ny == e2%cvArray(imgrid)%ny) &
                                                          .and.(e1%cvArray(imgrid)%nz == e2%cvArray(imgrid)%nz)) then

                     ! complex scalar multiplication for x,y,z-components

                     e2%cvArray(imgrid)%x = e2%cvArray(imgrid)%x + e1%cvArray(imgrid)%x * c
                     e2%cvArray(imgrid)%y = e2%cvArray(imgrid)%y + e1%cvArray(imgrid)%y * c
                     e2%cvArray(imgrid)%z = e2%cvArray(imgrid)%z + e1%cvArray(imgrid)%z * c

               else
                 print *, 'Error:scMultAdd_cvector_mg: vectors not same size'

               end if
                enddo ! loop on subgrids
              else
                print *, 'not compatible usage for scMultAdd_cvector'
              end if
            end if

     end subroutine scMultAdd_cvector_mg

! ****************************************************************************************************************
    ! Created by Cherevatova (Jan, 2013)
    ! This routine computes Cvector_mg on the coarse side
    ! taking Cvector_mg on the fine side.
    ! Used in some routines like M1Interface
    ! to fill in nz+1 layer on the coarse side from 1 layer of fine side
    ! or vice-versa
    subroutine average_cvector_mg(E)

    implicit none
     type(cvector_mg), intent(inout)    :: E

     !local
     integer                            :: imgrid
     integer                            :: ix, iy, iz
     integer                            :: nz

     if(.not.E%allocated)then
       print*, 'average_cvector_mg: E not allocated'
       stop
     endif

     do imgrid = 1, E%mGridSize-1
       nz = E%cvArray(imgrid)%nz
       select case (E%grid%interfaceType(imgrid))

       case(f2c) ! interface: fine to coarse
         ! Ex component
         do iy = 2, E%cvArray(imgrid+1)%ny
           do ix = 1, E%cvArray(imgrid+1)%nx
             E%cvArray(imgrid+1)%x(ix,iy,1) =  (E%cvArray(imgrid)%x(2*ix-1,2*iy-1,nz+1)*E%grid%gridArray(imgrid)%dx(2*ix-1)&
                                                 +E%cvArray(imgrid)%x(2*ix,2*iy-1,nz+1)*E%grid%gridArray(imgrid)%dx(2*ix))/&
                                                 (E%grid%gridArray(imgrid)%dx(2*ix-1)+E%grid%gridArray(imgrid)%dx(2*ix))
           enddo ! ix
        enddo ! iy
        ! Ey component
        do iy = 1, E%cvArray(imgrid+1)%ny
          do ix = 2, E%cvArray(imgrid+1)%nx
            E%cvArray(imgrid+1)%y(ix,iy,1) =  (E%cvArray(imgrid)%y(2*ix-1,2*iy-1,nz+1)*E%grid%gridArray(imgrid)%dy(2*iy-1)&
                                                  +E%cvArray(imgrid)%y(2*ix-1,2*iy,nz+1)*E%grid%gridArray(imgrid)%dy(2*iy))/&
                                                  (E%grid%gridArray(imgrid)%dy(2*iy-1)+E%grid%gridArray(imgrid)%dy(2*iy))
          enddo ! ix
        enddo ! iy

       case(c2f) ! interface : coarse to fine
       ! Ex component
       do iy = 2, E%cvArray(imgrid)%ny
         do ix = 1, E%cvArray(imgrid)%nx
            E%cvArray(imgrid)%x(ix,iy,nz+1) =  (E%cvArray(imgrid+1)%x(2*ix-1,2*iy-1,1)*E%grid%gridArray(imgrid+1)%dx(2*ix-1)&
                                                  +E%cvArray(imgrid+1)%x(2*ix,2*iy-1,1)*E%grid%gridArray(imgrid+1)%dx(2*ix))/&
                                                  (E%grid%gridArray(imgrid+1)%dx(2*ix-1)+E%grid%gridArray(imgrid+1)%dx(2*ix))
         enddo ! ix
       enddo  ! iy
       ! Ey component
       do iy = 1, E%cvArray(imgrid)%ny
         do ix = 2, E%cvArray(imgrid)%nx
            E%cvArray(imgrid)%y(ix,iy,nz+1) =  (E%cvArray(imgrid+1)%y(2*ix-1,2*iy-1,1)*E%grid%gridArray(imgrid+1)%dy(2*iy-1)&
                                                  +E%cvArray(imgrid+1)%y(2*ix-1,2*iy,1)*E%grid%gridArray(imgrid+1)%dy(2*iy))/&
                                                  (E%grid%gridArray(imgrid+1)%dy(2*iy-1)+E%grid%gridArray(imgrid+1)%dy(2*iy))
         enddo ! ix
       enddo  ! iy
       case(f2f) ! interface fine to fine/ coarse to coarse grid
         E%cvArray(imgrid)%x(:, :, nz+1) = E%cvArray(imgrid+1)%x(:, :, 1)
         E%cvArray(imgrid)%y(:, :, nz+1) = E%cvArray(imgrid+1)%y(:, :, 1)
       case (orig)
         ! this is original ModEM, no sub-grids
         ! nothing to do with the interfaces
         return
       case default
         print*, 'average_cvector_mg: select case statement error'
       end select


     enddo ! imgrid

    end subroutine average_cvector_mg
! **************************************************************************************************************
    ! Created by Cherevatova (Jan, 2013)
    ! This routine computes Rvector_mg on the coarse side
    ! taking Rvector_mg on the fine side.
    ! Used in some routines like M1Interface
    ! to fill in nz+1 layer on the coarse side from 1 layer of fine side
    ! or vice-versa
    subroutine average_rvector_mg(E)

    implicit none
     type(rvector_mg), intent(inout)    :: E

     !local
     integer                            :: imgrid
     integer                            :: ix, iy, iz
     integer                            :: nz

     if(.not.E%allocated)then
       print*, 'average_rvector_mg: E not allocated'
       stop
     endif

     do imgrid = 1, E%mGridSize-1
       nz = E%rvArray(imgrid)%nz
       select case (E%grid%interfaceType(imgrid))

       case(f2c) ! interface: fine to coarse
         ! Ex component
         do iy = 2, E%rvArray(imgrid+1)%ny
           do ix = 1, E%rvArray(imgrid+1)%nx
             E%rvArray(imgrid+1)%x(ix,iy,1) =  (E%rvArray(imgrid)%x(2*ix-1,2*iy-1,nz+1)*E%grid%gridArray(imgrid)%dx(2*ix-1)&
                                                 +E%rvArray(imgrid)%x(2*ix,2*iy-1,nz+1)*E%grid%gridArray(imgrid)%dx(2*ix))/&
                                                 (E%grid%gridArray(imgrid)%dx(2*ix-1)+E%grid%gridArray(imgrid)%dx(2*ix))
           enddo ! ix
        enddo ! iy
        ! Ey component
        do iy = 1, E%rvArray(imgrid+1)%ny
          do ix = 2, E%rvArray(imgrid+1)%nx
            E%rvArray(imgrid+1)%y(ix,iy,1) =  (E%rvArray(imgrid)%y(2*ix-1,2*iy-1,nz+1)*E%grid%gridArray(imgrid)%dy(2*iy-1)&
                                                  +E%rvArray(imgrid)%y(2*ix-1,2*iy,nz+1)*E%grid%gridArray(imgrid)%dy(2*iy))/&
                                                  (E%grid%gridArray(imgrid)%dy(2*iy-1)+E%grid%gridArray(imgrid)%dy(2*iy))
          enddo ! ix
        enddo ! iy

       case(c2f) ! interface : coarse to fine
       ! Ex component
       do iy = 2, E%rvArray(imgrid)%ny
         do ix = 1, E%rvArray(imgrid)%nx
            E%rvArray(imgrid)%x(ix,iy,nz+1) =  (E%rvArray(imgrid+1)%x(2*ix-1,2*iy-1,1)*E%grid%gridArray(imgrid+1)%dx(2*ix-1)&
                                                  +E%rvArray(imgrid+1)%x(2*ix,2*iy-1,1)*E%grid%gridArray(imgrid+1)%dx(2*ix))/&
                                                  (E%grid%gridArray(imgrid+1)%dx(2*ix-1)+E%grid%gridArray(imgrid+1)%dx(2*ix))
         enddo ! ix
       enddo  ! iy
       ! Ey component
       do iy = 1, E%rvArray(imgrid)%ny
         do ix = 2, E%rvArray(imgrid)%nx
            E%rvArray(imgrid)%y(ix,iy,nz+1) =  (E%rvArray(imgrid+1)%y(2*ix-1,2*iy-1,1)*E%grid%gridArray(imgrid+1)%dy(2*iy-1)&
                                                  +E%rvArray(imgrid+1)%y(2*ix-1,2*iy,1)*E%grid%gridArray(imgrid+1)%dy(2*iy))/&
                                                  (E%grid%gridArray(imgrid+1)%dy(2*iy-1)+E%grid%gridArray(imgrid+1)%dy(2*iy))
         enddo ! ix
       enddo  ! iy
       case(f2f) ! interface fine to fine/ coarse to coarse grid
         E%rvArray(imgrid)%x(:, :, nz+1) = E%rvArray(imgrid+1)%x(:, :, 1)
         E%rvArray(imgrid)%y(:, :, nz+1) = E%rvArray(imgrid+1)%y(:, :, 1)
       case (orig)
         ! this is original ModEM, no sub-grids
         ! nothing to do with the interfaces
         return
       case default
         print*, 'average_rvector_mg: select case statement error'
       end select


     enddo ! imgrid

    end subroutine average_rvector_mg


     ! *********************************************************************************************************

      function conjg_cvector_mg_f(e1) result (e2)
        ! conjg_cvector_mg_f computes a conjugate of a derived data type cvector_mg
        ! A.K.
        implicit none
        type (cvector_mg), intent(in)            :: e1
        type (cvector_mg)                        :: e2

        integer                                  :: status
        integer                                  :: imgrid

        ! check to see if RHS (e1) is active (allocated)
        if(.not.e1%allocated) then
           print*, 'input not allocated yet for conjg_cvector_mg_f'
        else

          if  (e1%gridType == e2%gridType) then
            if(e2%mgridSize /= e1%mgridSize) then
              print*, 'Warning: e2 and e1 are not the same size in function conjg_cvector_mg_f. Reallocated! '
              call deall(e2)
              call create(e1%grid, e2,e1%gridType)
            else
              do imgrid = 1, e1%mgridSize ! Global loop over sub-grid
                if((e2%cvArray(imgrid)%nx == e1%cvArray(imgrid)%nx).and.(e2%cvArray(imgrid)%ny == e1%cvArray(imgrid)%ny) &
                                                                 .and.(e2%cvArray(imgrid)%nz == e1%cvArray(imgrid)%nz)) then

                  ! just conjugate components
                  e2%cvArray(imgrid)%x = conjg(e1%cvArray(imgrid)%x)
                  e2%cvArray(imgrid)%y = conjg(e1%cvArray(imgrid)%y)
                  e2%cvArray(imgrid)%z = conjg(e1%cvArray(imgrid)%z)
                  e2%cvArray(imgrid)%gridType = e1%cvArray(imgrid)%gridType
                else
                  if(e2%allocated)then
                    ! first deallocate memory for x,y,z
                    deallocate(e2%cvArray(imgrid)%x,e2%cvArray(imgrid)%y,e2%cvArray(imgrid)%z, STAT=status)
                  endif

                  !  then allocate E3 as correct size ...
                  call create(e1%cvarray(imgrid)%grid, e2%cvarray(imgrid),e1%gridType)

                  !   .... and conjugate E1
                  e2%cvArray(imgrid)%x = conjg(e1%cvArray(imgrid)%x)
                  e2%cvArray(imgrid)%y = conjg(e1%cvArray(imgrid)%y)
                  e2%cvArray(imgrid)%z = conjg(e1%cvArray(imgrid)%z)
                  e2%cvArray(imgrid)%gridType = e1%cvArray(imgrid)%gridType
                endif ! check nx,ny,nz
              enddo  ! Global loop over sub-grid
            endif  ! check mgridSize
          else
           print*, 'not compatible usage for conjg_cvector_mg_f'
         end if ! check gridType
        endif ! check e1 allocation
        e2%temporary = .true.

      end function conjg_cvector_mg_f  ! conjg_cvector_mg_f

      ! ***************************************************************************
      function cmplx_rvector_mg_f(e1, e2) result (e3)
      ! inputs two real vectors defined on the multi-grid, merges them as real1 + imag(real2), a complex
      ! vector on mg
        implicit none
        type (rvector_mg), intent(in)            :: e1
        type (rvector_mg), intent(in)            :: e2
        type (cvector_mg)                        :: e3

        integer                               :: status
        integer                               :: imgrid

        ! check to see if RHS (E1 and E2) are active (allocated)
        if((.not.e1%allocated).or.(.not.e2%allocated)) then
           write(0,*) 'RHS not allocated yet for cmplx_rvector_mg_f'
        else
          if ((e1%gridType == e2%gridType).and.(e1%gridType == e3%gridType)) then
            if (e1%mgridSize == e2%mgridSize .and. e3%mgridSize == e2%mgridSize .and. e3%mgridSize == e1%mgridSize)then
              do imgrid = 1, e1%mgridSize ! Global loop over sub-grids
                if((e3%cvArray(imgrid)%nx == e1%rvArray(imgrid)%nx).and.(e3%cvArray(imgrid)%ny == e1%rvArray(imgrid)%ny).and.(e3%cvArray(imgrid)%nz == e1%rvArray(imgrid)%nz).and.&
                  (e3%cvArray(imgrid)%nx == e2%rvArray(imgrid)%nx).and.(e3%cvArray(imgrid)%ny == e2%rvArray(imgrid)%ny).and.(e3%cvArray(imgrid)%nz == e2%rvArray(imgrid)%nz))  then

                   ! create a complex pair
                  e3%cvArray(imgrid)%x = cmplx(e1%rvArray(imgrid)%x, e2%rvArray(imgrid)%x, prec)
                  e3%cvArray(imgrid)%y = cmplx(e1%rvArray(imgrid)%y, e2%rvArray(imgrid)%y, prec)
                  e3%cvArray(imgrid)%z = cmplx(e1%rvArray(imgrid)%z, e2%rvArray(imgrid)%z, prec)
                  e3%cvArray(imgrid)%gridType = e1%rvArray(imgrid)%gridType
                else
                  if(e3%allocated) then
                    ! first deallocate memory for x,y,z
                    deallocate(e3%cvArray(imgrid)%x,e3%cvArray(imgrid)%y,e3%cvArray(imgrid)%z, stat=status)
                  endif
                  !  then allocate e3 as correct size ...
                  call create(e1%rvArray(imgrid)%grid, e3%cvArray(imgrid), e1%gridType)
                  !   .... and create a complex pair
                  e3%cvArray(imgrid)%x = cmplx(e1%rvArray(imgrid)%x, e2%rvArray(imgrid)%x, prec)
                  e3%cvArray(imgrid)%y = cmplx(e1%rvArray(imgrid)%y, e2%rvArray(imgrid)%y, prec)
                  e3%cvArray(imgrid)%z = cmplx(e1%rvArray(imgrid)%z, e2%rvArray(imgrid)%z, prec)
                  e3%cvArray(imgrid)%gridType = e1%rvArray(imgrid)%gridType
                end if ! check nx, ny, nz
              enddo ! Global loop over sub-grids
            else
              write (0, *) 'not compatible usage for cmplx_rvector_mg_f'
            end if  ! check mgridSize
          else
           write (0, *) 'not compatible usage for cmplx_rvector_mg_f'
          end if  ! check gridType
        endif ! check allocate
        e3%temporary = .true.

      end function cmplx_rvector_mg_f  ! cmplx_rvector_f


      ! ***************************************************************************
      ! real_cvector_mg_f copies the real part of the derived data type cvector_mg variable;
      ! to produce a derived data type rvector_mg
      function real_cvector_mg_f(e1) result (e2)

        implicit none
        type (cvector_mg), intent(in)            :: e1
        type (rvector_mg)                        :: e2

        integer                               :: status
        integer                               :: imgrid

        ! check to see if RHS (e1) is active (allocated)
        if(.not.e1%allocated) then
           write(0,*) 'input not allocated yet for real_cvector_mg_f'
        else

          ! we know nothing about e2 ... deallocate just in case
          Call deall(e2) !deall_rvector_mg
          !  then allocate e2 as correct size ...
          Call create(e1%grid, e2, e1%gridType) !create_rvector_mg
          !   .... and copy E1
          do imgrid = 1, e1%mgridSize
            e2%rvArray(imgrid)%x = real(e1%cvArray(imgrid)%x)
            e2%rvArray(imgrid)%y = real(e1%cvArray(imgrid)%y)
            e2%rvArray(imgrid)%z = real(e1%cvArray(imgrid)%z)
            e2%rvArray(imgrid)%gridType = e1%cvArray(imgrid)%gridType
          enddo
        end if
        e2%temporary = .true.

      end function real_cvector_mg_f  ! real_cvector_mg_f


      ! ***************************************************************************
      ! imag_cvector_mg_f copies the imag part of the derived data type cvector_mg variable;
      ! to produce a derived data type rvector_mg
      function imag_cvector_mg_f(e1) result (e2)

        implicit none
        type (cvector_mg), intent(in)            :: e1
        type (rvector_mg)                        :: e2

        integer                               :: status
        integer                               :: imgrid

        ! check to see if RHS (e1) is active (allocated)
        if(.not.e1%allocated) then
           write(0,*) 'input not allocated yet for imag_cvector_mg_f'
        else

          ! we know nothing about e2 ... deallocate just in case
          Call deall(e2) !deall_rvector_mg
          !  then allocate e2 as correct size ...
          Call create(e1%grid, e2, e1%gridType) !create_rvector_mg
          !   .... and copy e1
          do imgrid = 1, e1%mgridSize
            e2%rvArray(imgrid)%x = aimag(e1%cvArray(imgrid)%x)
            e2%rvArray(imgrid)%y = aimag(e1%cvArray(imgrid)%y)
            e2%rvArray(imgrid)%z = aimag(e1%cvArray(imgrid)%z)
            e2%rvArray(imgrid)%gridType = e1%cvArray(imgrid)%gridType
          enddo
        end if

        e2%temporary = .true.

      end function imag_cvector_mg_f  ! imag_cvector_mg_f





      ! ******************************************************************************************************
      subroutine plotCvector_mg(inE) ! need for testing
        ! created by Cherevatova (Oct, 2012)

        implicit none
          type (cvector_mg), intent(in)         :: inE
          ! local
          character(len=10)                     :: XZYZ,XY
          real (kind=prec)                      :: tempD
          type(cvector)                         :: Etemp
          integer                               :: Io, ix, iy, iz, nx, ny, nz, imgrid

          nx = inE%grid%nx
          ny = inE%grid%ny
          nz = inE%grid%nz

          tempD= 0
          do ix = 1,nx-1
            tempD= tempD+inE%grid%dx(ix)
            Io = ix
            if (tempD == -inE%grid%ox) exit
          enddo

          ! plot vertical cut
          call create(inE%grid,Etemp,inE%gridType)
          open(001, file = 'XZYZ')
            write(001,*)'X, Y, Z, Ex, Ey, Ez'
                Etemp = inE !mg2c
                write(001,'(i3,2x,i3,2x,i3,2x,es13.5,2x,es13.5,2x,es13.5)') ((Io, iy,iz, &
                                     abs(Etemp%x(Io,iy,iz)), &
                                     abs(Etemp%y(Io,iy,iz)), &
                                     abs(Etemp%z(Io,iy,iz)),iy=1,ny),iz=1,nz)
                call zero(Etemp)
            close(001)

          ! plot horizontal cut
          open(002, file = 'XY')
            write(002,*)'X, Y, Z,Ex, Ey, Ez'
              ! specify number of iz manually
              imgrid = 1
              iz = 20
               write(002,'(i3,2x,i3,2x,i3,2x,es13.5,2x,es13.5,2x,es13.5)') ((ix, iy, iz, &
                                     abs(inE%cvArray(imgrid)%x(ix,iy,iz)), &
                                     abs(inE%cvArray(imgrid)%y(ix,iy,iz)), &
                                     abs(inE%cvArray(imgrid)%z(ix,iy,iz)),ix=1,inE%cvArray(imgrid)%nx),iy=1,inE%cvArray(imgrid)%ny)
            close(002)

        call deall(Etemp)

      end subroutine plotCvector_mg


    end module sg_vector_mg
