! *****************************************************************************
module output_orig
  ! Module containing the subroutines for a-posteriori analysis of the output
  ! vector h=lH from the routine SolveMaxwells

  use math_constants
  use elements
  use interp_orig
  use iospec
  implicit none


Contains

  ! ***************************************************************************
  ! * The following subroutines are all inherited from the older version of the
  ! * forward solver. They compute the data functionals and output them. However,
  ! * we aim to separate these two actions. The data functionals will be computed
  ! * in module dataFunc.f90, while this module shall contain output subroutines
  ! * only.

      subroutine avg_cresp_out(freq,num,fn_avg_cresp,cresp)

! ---- output calculated model response
	  integer, intent(in)						:: num
	  real(8), dimension(num), intent(in)		:: freq
      complex(8), dimension(num), intent(in)	:: cresp
	  character(80), intent(in)					:: fn_avg_cresp
	  integer									:: i

      open(ioAvgC,file=fn_avg_cresp,status='unknown',form='formatted')

      write(ioAvgC,'(a)') 'Global average C response [km]'
      do i=1,num
         write(ioAvgC,'(i5,1x,3g15.7)') i,freq(i),cresp(i)
      end do
 
      close(ioAvgC)

	  end subroutine avg_cresp_out


      subroutine allh_out(l,m,n,hx,hy,hz,x,y,z,nair)

!-------------------------------------------------------------
!     to output all the magnetic fields with position
!-------------------------------------------------------------

	  implicit none

	  integer						   	 :: l,m,n
      real(8),dimension(l)	             :: x
      real(8),dimension(m+1)             :: y
      real(8),dimension(n+1)             :: z
      complex(8),dimension(:)			 :: hx,hy,hz
      integer				             :: nair
      real(8)                            :: xx0,xx1,xx,yy,zz
      integer                            :: k,j,i,ii

	  inquire(ioHx,OPENED=opened)
	  if (.not.opened) then
		write (0,*) 'Error: (allh_out) output file for Hx not opened'
		return
	  end if
	  inquire(ioHy,OPENED=opened)
	  if (.not.opened) then
		write (0,*) 'Error: (allh_out) output file for Hy not opened'
		return
	  end if
	  inquire(ioHz,OPENED=opened)
	  if (.not.opened) then
		write (0,*) 'Error: (allh_out) output file for Hz not opened'
		return
	  end if

!
!   First Hx....
!
      do k=nair+1,n+1
!      do k=nair+1,nair+1 !temporarily only output surface fields
        zz=z(k)/1000.d0
        do j=2,m
          yy=y(j)*r2d
          xx0=0.d0
          do i=1,l
            xx1=(xx0+x(i)*r2d)
            xx=(xx0+xx1)/2.d0
            xx0=xx1
            call n_allhxijk(l,m,n,i,j,k,ii)
            write(ioHx,'(5e15.7)')xx,yy,zz,hx(ii)
          end do
        end do
      end do
!
!   Now Hy....
!
      do k=nair+1,n+1
!      do k=nair+1,nair+1 !temporarily only output surface fields
        zz=z(k)/1000.d0
        do j=1,m
          yy=( y(j)+y(j+1) )*r2d/2.d0
          xx=0.d0
          do i=1,l
             call n_allhyijk(l,m,n,i,j,k,ii)
             write(ioHy,'(5e15.7)')xx,yy,zz,hy(ii)
             xx=xx+x(i)*r2d
          end do
        end do
      end do

!
!   Now Hz....
!
      do k=nair+1,n
!      do k=nair+1,nair+1 !temporarily only output surface fields
        zz=( z(k)+z(k+1) )/1000.d0/2.d0
        j=1
          yy=y(j)*r2d
          i=1
            xx=0.d0
            call n_allhzijk(l,m,n,i,j,k,ii)
            write(ioHz,'(5e15.7)')xx,yy,zz,hz(ii)

        do j=2,m
          yy=y(j)*r2d
          xx=0.d0
          do i=1,l
            call n_allhzijk(l,m,n,i,j,k,ii)
            write(ioHz,'(5e15.7)')xx,yy,zz,hz(ii)
            xx=xx+x(i)*r2d
          end do
        end do

        j=m+1
          yy=y(j)*r2d
          i=1
            xx=0.d0
            call n_allhzijk(l,m,n,i,j,k,ii)
            write(ioHz,'(5e15.7)')xx,yy,zz,hz(ii)
      end do
      return
      end subroutine allh_out


! *****************************************************************************


      subroutine allh_out_geom(l,m,n,hx,hy,hz,x,y,z,nair)

!-------------------------------------------------------------
!     to output all the magnetic fields including air layers
!     last mod: A.S. 10.7.98 - fits format required by display
!     program natural3dsurf
!-------------------------------------------------------------

	  implicit none

	  integer								:: l,m,n
      real(8),dimension(l)		            :: x
      real(8),dimension(m+1)		        :: y
      real(8),dimension(n+1)                :: z
      real(8)                               :: x0,y0
      complex(8), dimension(:)				:: hx,hy,hz
      integer					            :: nair
      integer                               :: i,j,k,ii_hx,ii_hy,ii_hz
      real(8)                               :: xx0,xx1,xx,yy,zz

!      pi=datan(1.d0)*4.d0
!      r2d=180.d0/pi

!
!   output staggered grid of hx,hy,hz - but simplify output so all
!   values appear to be on a non-staggered grid defined at the lat,lon,r
!   nodes. In reality, Hx is not defined at the poles, and is defined at
!   longitudes intermediate between the nodes. Hy is not defined at
!   poles, and is defined at latitudes intermediate between the nodes.
!   Hz is defined at the poles, and at the lat,lon nodes - but at depths
!   intermediate between the shells.
!
!   Note - the indexing scheme is a bit odd here - specifically for purposes
!   of graphical output. There are, in reality, m+1 latitude nodes. At each of
!   these nodes, the Hz value is defined. Node 0 is the N pole, and node m+1
!   is the S pole. At a total of m nodes at latitudes intermediate between
!   the Hz nodes, the Hy field is defined. At m-1 nodes at the same latitudes
!   as the Hz nodes - but skipping the poles and at intermediate longitudes, the
!   Hz field is defined. To complicate matters further, the electrical
!   conductivity model is defined at m nodes at latitudes intermediate between
!   the Hz nodes. For graphical purposes, to these points we augment the
!   grid artificially with N and S pole points, otherwise the geometry of
!   the output object will no longer be convex and will appear odd when
!   visualised.
!
	  inquire(ioH,OPENED=opened)
	  if (.not.opened) then
		write (0,*) 'Error: (allh_out_geom) output file not opened'
		return
	  end if
	  write(ioH,*) n,l*m+2
      do k=1,n                                 ! loop over radial shells
        zz=z(k)/6371000.0d0
        write(ioH,*) zz                         ! radial distance of shell
        call n_allhxijk(l,m,n,1,2,k,ii_hx)     ! Hx N-pole value - not def - so one node down
        call n_allhyijk(l,m,n,1,1,k,ii_hy)     ! Hy N-pole value
        call n_allhzijk(l,m,n,1,1,k,ii_hz)     ! Hz N-pole value

        write(ioH,'(6g15.7)') hx(ii_hx),hy(ii_hy),hz(ii_hz)

        do i=1,l
          do j=1,m
            if (j == 1) then
               call n_allhxijk(l,m,n,i,2,k,ii_hx)
            else
               call n_allhxijk(l,m,n,i,j,k,ii_hx)
            end if
            call n_allhyijk(l,m,n,i,j,k,ii_hy)
            call n_allhzijk(l,m,n,i,j,k,ii_hz)

            write(ioH,'(6g15.7)') hx(ii_hx),hy(ii_hy),hz(ii_hz)

          end do
        end do
        call n_allhxijk(l,m,n,1,m,k,ii_hx)     ! Hx S-pole value
        call n_allhyijk(l,m,n,1,m,k,ii_hy)     ! Hy S-pole value
        call n_allhzijk(l,m,n,1,m+1,k,ii_hz)   ! Hz S-pole value

        write(ioH,'(6g15.7)') hx(ii_hx),hy(ii_hy),hz(ii_hz)

      end do

      return
      end subroutine allh_out_geom


! *****************************************************************************


      subroutine calc_response(omega,l,m,n,hx,hy,hz,x,y,z,nair,cdat,ddat)

!-------------------------------------------------------------
!   to output the c and d response vs position within the model domain
!
!   Note - the indexing scheme is a bit odd here - specifically for purposes
!   of graphical output. There are, in reality, m+1 latitude nodes. At each of
!   these nodes, the Hz value is defined. Node 1 is the N pole, and node m+1
!   is the S pole. At a total of m nodes at latitudes intermediate between
!   the Hz nodes, the Hy field is defined. At m-1 nodes at the same latitudes
!   as the Hz nodes - but skipping the poles and at intermediate longitudes, 
!   the Hx field is defined. To complicate matters further, the electrical
!   conductivity model is defined at m nodes at latitudes intermediate between
!   the Hz nodes. For graphical purposes, to these points we augment the
!   grid artificially with N and S pole points, otherwise the geometry of
!   the output object will no longer be convex and will appear odd when
!   visualised.
!-------------------------------------------------------------

	  implicit none

      real(8)                         :: omega
	  integer						  :: l,m,n
      real(8),    dimension(l)	      :: x,x_hx,x_hy,x_hz
      real(8),    dimension(m+1)      :: y,y_hx,y_hy,y_hz
      real(8),    dimension(n+1)      :: z,z_hx,z_hy,z_hz
      real(8),    dimension(3)        :: xa,ya,y2
      complex(8)                      :: cdat, ddat
      integer					      :: nair
      integer                         :: nc,nspline
      complex(8), dimension(:)		  :: hx,hy,hz
      complex(8)                      :: hz_up_n,hz_up_s,hz_i,hx_av
      complex(8)                      :: hx_i,hx_n,hx_s
      complex(8)                      :: c,d
! --  intel compiler produces a run-time stack overflow if this is used
!     complex(8)                      :: h(n+1,m+1,l,3)
      complex(8), allocatable         :: h(:,:,:,:)
      real(8)                         :: hz_up_r,hz_up_i,d_n,d_s
      real(8)                         :: d_ne,d_nw,d_se,d_sw,c_w
      real(8)                         :: eps
      real(8)                         :: xx0,xx1,xx,yy,zz,xz,dy
      real(8)                         :: colat,yp1,ypn
      integer                         :: i,ii,j,k,kk,im1

      data yp1/1.d99/,ypn/1.d99/

      allocate (h(1:n+1,1:m+1,1:l,3))

	  inquire(ioC,OPENED=opened)
	  if (.not.opened) then
		write (0,*) 'Error: (calc_response) output file for C responses not opened'
		return
	  end if
	  inquire(ioD,OPENED=opened)
	  if (.not.opened) then
		write (0,*) 'Error: (calc_response) output file for D responses not opened'
		return
	  end if


!      pi  = datan(1.d0)*4.d0
!      r2d = 180.d0/pi
      eps = 1.0D5*epsilon(1.0D0)
      nspline   = 3

!
!   First Hx....
!
      do k=nair+1,n+1
        z_hx(k) = z(k)/6371000.d0
        do j=2,m                                   ! Hx not actually defined at poles
          y_hx(j) = y(j)*r2d
          xx0     = 0.d0
          do i=1,l
            call n_allhxijk(l,m,n,i,j,k,ii)
            xx1        =  xx0 + x(i)*r2d
            x_hx(i)    = (xx0 + xx1)/2.d0
            xx0        = xx1
            h(k,j,i,1) = hx(ii)
          end do
        end do
      end do
!
!   Now Hy....
!
      do k=nair+1,n+1
        z_hy(k) = z(k)/6371000.d0
        do j=1,m                                  ! defined at intermediate latitudes
          y_hy(j) = (y(j) + y(j+1))*r2d/2.d0
          xx      = 0.0d0
          do i=1,l
             call n_allhyijk(l,m,n,i,j,k,ii)
             x_hy(i)    = xx
             xx         = xx + x(i)*r2d
             h(k,j,i,2) = hy(ii)
          end do
        end do
      end do
!
!   Now Hz....
!
!Schultz: For vertical interpolation of Hz, Hz should start from 1/2 above
!         the earth's surface. 24/FEB/99

      do k=(nair+1)-1,n
!Schultz
        z_hz(k)    = (z(k) + z(k+1))/6371000.d0/2.d0

        j          = 1                             ! N-pole
        y_hz(j)    = y(j)*r2d
        i          = 1
        xx         = 0.0d0
        call n_allhzijk(l,m,n,i,j,k,ii)
        h(k,j,i,3) = hz(ii)

!toh: Copy the N-pole Hz to other meridians. So, set the same values to h(k,1,i).

        do i = 2, l
           h(k,1,i,3) = hz(ii)
        end do

!toh: This seems redundant. However, this is essential for right C-response
!     evaluation at the N-most Hy latitude to which the N-pole Hz is
!     laterally interpolated. 24/FEB/99

        do j=2,m                                   ! non-polar
          y_hz(j) = y(j)*r2d
          xx      = 0.d0
          do i=1,l
            x_hz(i) = xx
            call n_allhzijk(l,m,n,i,j,k,ii)
            h(k,j,i,3) = hz(ii)
            xx=xx+x(i)*r2d
          end do
        end do

        j        = m+1                             ! S-pole
        y_hz(j)  = y(j)*r2d
        i        = 1
        x_hz(i)  = 0.0d0
        call n_allhzijk(l,m,n,i,j,k,ii)
        h(k,j,i,3) = hz(ii) 

!toh: Copy the S-pole Hz to other meridians. So, set the same values to h(k,m+1,i).

        do i = 2, l
           h(k,m+1,i,3) = hz(ii)
        end do

!toh: This seems redundant. However, this is essential for right C-response
!     evaluation at the S-most Hy latitude to which the S-pole Hz is
!     laterally interpolated. 24/FEB/99

      end do

!
!   Now calculate C and D responses. 
!   c = Re/2 * tan(theta) * Hz/Hy 
!   d = Re/2 * cos(theta) * Hx/Hy [theta is colatitude]
!
!   but, note that
!
!   while Hz is defined at the poles, Hy is defined at the
!   latitudes in between the Hz nodes and at the same longitudes
!   as Hz; and Hx is defined at the Hz nodes in latitude, but at staggered
!   longitudes and also not at the poles. The Hz nodes are at depths
!   intermediate between the Hx, Hy nodes for the current shell and
!   the shell below. The zeta response requires calculation of the theta and
!   phi derivatives of Hx, which are defined ambiguously at the poles.
!   We therefore apply the convention that none of the responses
!   will be defined at the poles;
!   thus, for reasons of graphical output, the response values will be set
!   identically to (0.0d0,0.0d0) at the poles.

! ----- generate C response by upward interpolating Hz values to lie on same shell plane
!       as Hx,Hy. The Hy nodes are defined at m points - which are at intermediate latitudes
!       between the m+1 Hz nodes; and are also coincident with the intermediate latitudes
!       at which the conductivity model is specified. Therefore, after upward interpolating the
!       Hz values, it is necessary to interplate along lines of constant longitude each pair
!       of Hz values - i.e. the one to the N and the one to the S of the Hy node location.
!       C is formed at this point from the local Hy and the upward and latitudinally averaged Hz
!       values.
!toh: # of C-response per frequency = l*m, ie no poles, as of 22/FEB/99

      !AK ---- outputting l and m for debug purposes
      ! (AK) write(ioC,*) "no. iter over r, phi, theta resp.",n-ksmap(1,1)+1,l,m
      ! (AK) write(ioD,*) "no. iter over r, phi, theta resp.",n-ksmap(1,1)+1,l,m

      cdat = dcmplx(0.0d0,0.0d0)
      ddat = dcmplx(0.0d0,0.0d0)
      nc   = 0

      write(ioC,*) "freq = ",omega/(2.0d0*pi)
      write(ioD,*) "freq = ",omega/(2.0d0*pi)
      do k=nair+1,n                              ! loop over non-air layers

         xa(1) = z_hz(k-1)                           ! depth of bottom air layer - for subsequent interpolation
         xa(2) = z_hz(k)                             ! first earth layer
         xa(3) = z_hz(k+1)                           ! shell below that
         xz    = z_hy(k)                             ! level where Hy is defined, and to which Hz will be interp.

!AK ---- simplifying the output
         !write(ioC,*) xz
         !write(ioD,*) xz

         kk = k
         if (kk > n-3) kk = n-3                       ! don't calculate response for shells > n-3th deepest.

         do i=1,l                                     ! all longitudes

            do j=1,m                                  ! intermediate latitudes where Hy and sigma is defined (non-polar)

               ya(1) = dreal(h(kk-1,j,i,3))           ! Real Hz N of Hy at level 1/2 below level above current Hx,Hy shell
               ya(2) = dreal(h(kk,j,i,3))             ! Real Hz N of Hy at level 1/2 below current Hx,Hy shell
               ya(3) = dreal(h(kk+1,j,i,3))           ! Real Hz N of Hy at level 1/2 below level below current

               call spline(xa,ya,nspline,yp1,ypn,y2)        ! spline interpolation of Real Hz upward to Hx,Hy node level
               call splint(xa,ya,y2,nspline,xz,hz_up_r)

               ya(1) = dimag(h(kk-1,j,i,3))           ! Imag Hz N of Hy at level 1/2 below level above current Hx,Hy shell
               ya(2) = dimag(h(kk,j,i,3))             ! Imag Hz N of Hy at level 1/2 below current Hx,Hy shell
               ya(3) = dimag(h(kk+1,j,i,3))           ! Imag Hz N of Hy at level 1/2 below level below current

               call spline(xa,ya,nspline,yp1,ypn,y2)        ! spline interpolation of Imag Hz upward to Hx,Hy node level
               call splint(xa,ya,y2,nspline,xz,hz_up_i)

               hz_up_n = dcmplx(hz_up_r,hz_up_i)      ! build new complex radially-upward interpolated Hz just N of Hy

               ya(1) = dreal(h(kk-1,j+1,i,3))         ! Real Hz S of Hy at level 1/2 below level above current Hx,Hy shell
               ya(2) = dreal(h(kk,j+1,i,3))           ! Real Hz S of Hy at level 1/2 below current Hx,Hy shell
               ya(3) = dreal(h(kk+1,j+1,i,3))         ! Real Hz S of Hy at level 1/2 below level below current

               call spline(xa,ya,nspline,yp1,ypn,y2)        ! spline interpolation of Real Hz upward to Hx,Hy node level
               call splint(xa,ya,y2,nspline,xz,hz_up_r)

               ya(1) = dimag(h(kk-1,j+1,i,3))         ! Imag Hz S of Hy at level 1/2 below level above current Hx,Hy shell
               ya(2) = dimag(h(kk,j+1,i,3))           ! Imag Hz S of Hy at level 1/2 below current Hx,Hy shell
               ya(3) = dimag(h(kk+1,j+1,i,3))         ! Imag Hz S of Hy at level 1/2 below level below current

               call spline(xa,ya,nspline,yp1,ypn,y2)        ! spline interpolation of Imag Hz upward to Hx,Hy node level
               call splint(xa,ya,y2,nspline,xz,hz_up_i)

               hz_up_s = dcmplx(hz_up_r,hz_up_i)      ! build new complex radially-upward interpolated Hz just S of Hy

               colat = y_hy(j)                        ! colatitude of Hy node

               d_n  = dabs(y_hy(j) - y_hz(j))         ! distance between Hy node & Hz node immediately N assume same lon.
               d_s  = dabs(y_hy(j) - y_hz(j+1))       ! distance between Hy node & Hz node immediately S assume same lon.
                                                      ! Hz interpolated laterally by L2 distance weighting

               hz_i = ( hz_up_n/(d_n**2) + hz_up_s/(d_s**2) )/( 1.0d0/(d_n**2) + 1.0d0/(d_s**2) )

               if (cdabs(h(kk,j,i,2)) > eps .and. colat /= 90.0d0) then

                  c = 6371.0d0*z_hy(kk)*(hz_i/dcos(colat))/(2.0d0*(h(kk,j,i,2)/dsin(colat)))
               else

                  c = dcmplx(0.0d0,0.0d0)

               end if

!toh: to interpolate Hx to where Hy resides

               if ( i == 1 ) then

                  im1 = l
                  c_w = 360.0d0

               else

                  im1 = i - 1
                  c_w = 0.0d0

               end if

! d_ne: distance along the latitude from Hx(i,j) to 1/2 north of Hy(i,j)
! d_nw: distance along the latitude from Hx(i-1,j) to 1/2 north of Hy(i,j)
! d_se: distance along the latitude from Hx(i,j+1) to 1/2 south of Hy(i,j)
! d_sw: distance along the latitude from Hx(i-1,j+1) to 1/2 south of Hy(i,j)

! hx_n: Hx interpolated along the same latitude 1/2 north of Hy (L2 weighting)
! hx_s: Hx interpolated along the same latitude 1/2 south of Hy (L2 weighting)

               if ( j == 1 ) then

                  d_se = dabs(dsin(y_hx(j+1))*(x_hx(i)-x_hz(i)))
                  d_sw = dabs(dsin(y_hx(j+1))*(x_hz(i)-x_hx(im1)+c_w))

                  hx_n = dcmplx(0.0d0,0.0d0)
                  hx_s = ( h(kk,j+1,i,1)/(d_se**2) + &
                    h(kk,j+1,im1,1)/(d_sw**2) )/( 1.0d0/(d_se**2) + 1.0d0/(d_sw**2) )

               else if ( j == m ) then

                  d_ne = dabs(dsin(y_hx(j))*(x_hx(i)-x_hz(i)))
                  d_nw = dabs(dsin(y_hx(j))*(x_hz(i)-x_hx(im1)+c_w))

                  hx_n = ( h(kk,j,i,1)/(d_ne**2) + &
                      h(kk,j,im1,1)/(d_nw**2) )/( 1.0d0/(d_ne**2) + 1.0d0/(d_nw**2) )
                  hx_s = dcmplx(0.0d0,0.0d0)

               else

                  d_ne = dabs(dsin(y_hx(j))*(x_hx(i)-x_hz(i)))
                  d_nw = dabs(dsin(y_hx(j))*(x_hz(i)-x_hx(im1)+c_w))
                  d_se = dabs(dsin(y_hx(j+1))*(x_hx(i)-x_hz(i)))
                  d_sw = dabs(dsin(y_hx(j+1))*(x_hz(i)-x_hx(im1)+c_w))

                  hx_n = ( h(kk,j,i,1)/(d_ne**2) + &
                      h(kk,j,im1,1)/(d_nw**2) )/( 1.0d0/(d_ne**2) + 1.0d0/(d_nw**2) )
                  hx_s = ( h(kk,j+1,i,1)/(d_se**2) + &
                    h(kk,j+1,im1,1)/(d_sw**2) )/( 1.0d0/(d_se**2) + 1.0d0/(d_sw**2) )

               end if

! Hx interpolated along the meridian by L2 distance weighting

               hx_i = ( hx_n/(d_n**2) + hx_s/(d_s**2) ) &
                     /( 1.0d0/(d_n**2) + 1.0d0/(d_s**2) )

!toh: to compute D-responses

               if (cdabs(h(kk,j,i,2)) > eps) then

                  d = 6371.0d0*z_hy(kk)*hx_i*dcos(colat)/(2.0d0*h(kk,j,i,2))
				  !print *, '!!!1 i,j,Hx,Hy,Hz=',i,j,hx_i,h(kk,j,i,2),hz_i
	!write(*,'(a40,2i6,2g15.7)') 'i,j,lon,colat = ',&
		!i,j,x_hx(i),colat

               else

                  d = dcmplx(0.0d0,0.0d0)

               end if
!toh
!AK: changes in the output
               ! (AK) write(ioC,'(3i5,2g15.7)') k,i,j,c

!toh: write out D-responses
!AK: changes in the output
               ! (AK) write(ioD,'(3i5,2g15.7)') k,i,j,d
!toh

!toh:
!toh: Determine the global mean C value at the Earth's surface using l times m
!toh: C values. 24/FEB/99
!toh:
               if ( k == nair+1 ) then	  !	k==ksmap(i,j)
                  write(ioC,'(2i5,2g15.7)') i,j,c
                  write(ioD,'(2i5,2g15.7)') i,j,d
                  cdat = cdat + c
                  ddat = ddat + d
                  nc = nc+1

               end if

            end do

         end do

      end do

      if ( nc == 0 ) then

         write(6,*) 'DEBUG response nc = 0 !'

      else if ( nc < l*m ) then

         write(6,*) 'DEBUG response nc =! l x m',nc

      else

         cdat = cdat/dfloat(nc)
         ddat = ddat/dfloat(nc)

      end if

!toh: added on 28/MAY/99 to save field values into files
!
!   First Hx....
!
      do k=1,n+1

        do j=2,m                             ! Hx not actually defined at poles

          do i=1,l

            call n_allhxijk(l,m,n,i,j,k,ii)
            write(ioHx,'(5g15.7)') x_hx(i),y_hx(j),z(k),hx(ii)

          end do

        end do

      end do
!
!   Now Hy....
!
      do k=1,n+1

        do j=1,m                            ! defined at intermediate latitudes

          do i=1,l

             call n_allhyijk(l,m,n,i,j,k,ii)
             write(ioHy,'(5g15.7)') x_hy(i),y_hy(j),z(k),hy(ii)

          end do

        end do

      end do
!
!   Now Hz....
!

      do k=1,n ! Hz's are defined at intermediate depths between two shells.

        do j=1,m+1 ! At poles, there're L same values per surface.

          do i=1,l

             if ( j == 1 .or. j == m+1 ) then

                call n_allhzijk(l,m,n,1,j,k,ii)

             else

                call n_allhzijk(l,m,n,i,j,k,ii)

             end if

             write(ioHz,'(5g15.7)') x_hz(i),y_hz(j),(z(k)+z(k+1))/2.d0,hz(ii)

          end do

        end do

      end do

!      write(ioH) hx
!      write(ioH) hy
!      write(ioH) hz

!toh: on 28/MAY/99

!toh: Deallocate h(n+1,m+1,l,3) explicitly. 25/FEB/99

      deallocate (h)
!toh
      return

      end subroutine calc_response



! *****************************************************************************



      subroutine allh_out_prin(l,m,n,hx,hy,hz,x,y,z)

!-------------------------------------------------------------
!     to output all the magnetic fields with position
!-------------------------------------------------------------

	  implicit none

	  integer						   :: l,m,n
      real(8),dimension(l)	           :: x
      real(8),dimension(m+1)           :: y
      real(8),dimension(n+1)           :: z
      complex(8),dimension(:)		   :: hx,hy,hz
      real(8)                          :: xx0,xx1,xx,yy,zz
      integer                          :: k,j,i,ii

	  inquire(ioHx,OPENED=opened)
	  if (.not.opened) then
		write (0,*) 'Error: (allh_out_prin) output file for Hx not opened'
		return
	  end if
	  inquire(ioHy,OPENED=opened)
	  if (.not.opened) then
		write (0,*) 'Error: (allh_out_prin) output file for Hy not opened'
		return
	  end if
	  inquire(ioHz,OPENED=opened)
	  if (.not.opened) then
		write (0,*) 'Error: (allh_out_prin) output file for Hz not opened'
		return
	  end if
      
!
!   First Hx....
!
      do k=1,n+1
        zz=z(k)/1000.d0
!       do j=2,m
        do j=m/2,m/2
          yy=y(j)*r2d
          xx0=0.d0
!         do i=1,l
          do i=1,l/2-1
            xx1=(xx0+x(i)*r2d)
            xx=(xx0+xx1)/2.d0
            xx0=xx1
          end do
          do i=l/2,l/2
            xx1=(xx0+x(i)*r2d)
            xx=(xx0+xx1)/2.d0
            xx0=xx1
            call n_allhxijk(l,m,n,i,j,k,ii)
            write(ioHx,'(5e12.4)')xx,yy,zz,hx(ii)
          end do
        end do
      end do
!
!   Now Hy....
!
      do k=1,n+1
        zz=z(k)/1000.d0
!       do j=1,m
        do j=m/2,m/2
          yy=( y(j)+y(j+1) )*r2d/2.d0
          xx=0.d0
!         do i=1,l
          do i=1,l/2-1
             xx=xx+x(i)*r2d
          end do
          do i=l/2,l/2
             call n_allhyijk(l,m,n,i,j,k,ii)
             write(ioHy,'(5e12.4)')xx,yy,zz,hy(ii)
             xx=xx+x(i)*r2d
          end do
        end do
      end do
!
!   Now Hz....
!
      do k=1,n
        zz=( z(k)+z(k+1) )/1000.d0/2.d0
        j=1
          yy=y(j)*r2d
          i=1
          xx=0.d0
          call n_allhzijk(l,m,n,i,j,k,ii)
          write(ioHz,'(5e12.4)')xx,yy,zz,hz(ii)
      end do
      return
      end subroutine allh_out_prin



! *****************************************************************************



      subroutine allj_out(iow,l,m,n,nair,omega,x,y,z,resist,hx,hy,hz)

!-------------------------------------------------------------
!
!     compute block surface centered j field
!     input h values are given on edges (parralel to edge) in 
!           real space units
!     output jx, jy and jz are given in real space units
!
!     by H Toh on 10/JUN/99 based on 'surfeh_out'
!
!-------------------------------------------------------------

	  implicit none

      integer                            :: l,m,n,iow
      complex(8),dimension(:)			 :: hx,hy,hz
      real(8),dimension(l)	             :: x
      real(8),dimension(m+1)             :: y
      real(8),dimension(n+1)             :: z
      real(8)                            :: omega
      real(8),dimension(l,m,n)	         :: resist

      complex(8)                         :: jx1,jx2,jy1,jy2
      complex(8)                         :: ciuw,jx,jy,jz
      real(8)                            :: uw
      real(8)                            :: yp,zp
      real(8)                            :: xijk,xiek,xijf,xief,xipp,xipk
      real(8)                            :: yjk,yjf,yjp
      real(8)                            :: zc,zk
      real(8)                            :: sijk,sjki,seki,skij,sjpi,sijp
      integer				             :: nair
      integer                            :: i,j,k,c,d,e,f
      integer                            :: ijkx,iekx,ijfx,iefx
      integer                            :: ijky,djky,ijfy,djfy
      integer                            :: ijcz,djcz,iecz,decz,ijkz
      integer                            :: djkz,iekz,dekz

      real(8)                            :: xxw0,xxw,yyw,zzw

	  inquire(ioJ,OPENED=opened)
	  if (.not.opened) then
		write (0,*) 'Error: (allj_out) output file for J not opened'
		return
	  end if

      uw=omega * MU_0
      ciuw=dcmplx(0.0d0,uw)

      !write(iow,*) n-ksmap(1,1)+1,l*m

      do k=nair+1,n

         write(iow,*) z(k)/6371000.d0

         c=k-1
         f=k+1

         call leng_zijk(c,z,zc)
         call leng_zijk(k,z,zk)

         xxw0=0.0d0
         zzw=z(k)/1000.0d0

         do i=1,l

            xxw=( xxw0+x(i)/2.0d0 )*r2d
            xxw0=xxw0+x(i)

            if (i.ne.l) then

               d=d+1

            else

               d=1

            end if

            do j=1,m

               e=j+1
               yyw=(y(j)+y(e))/2.0d0*r2d

               if (j.ne.1.and.j.ne.m) then

                  call leng_xijk(i,j,k,x,y,z,xijk)
                  call leng_xijk(i,e,k,x,y,z,xiek)

               else if (j == 1) then

                  call leng_xijk(i,e,k,x,y,z,xiek)

               else if (j == m) then

                  call leng_xijk(i,j,k,x,y,z,xijk)

               end if

!--------------------------------------------------------------
!
! compute current densities
!
! Vertical 3-point spline interpolation is definitely needed here.
! I just don't have time to do that, so jx & jy are always defined
! at half  shell thickness below jz in its present form.

! Also, horizontal interpolation should be of the form of L2
! weighing in the case of a variable mesh rather than the present
! mere arithmetic mean.
!
! by H TOH @ ITG, Dept Earth Sciences, Univ. Cambridge on 10/JUN/'99
!
!-----------------computation for Jx --------------------------

               call leng_yijk(j,k,y,z,yjk)
               call leng_yijk(j,f,y,z,yjf)

               call area_sijk(j,k,y,z,sijk)

               call n_allhyijk(l,m,n,i,j,k,ijky)
               call n_allhyijk(l,m,n,i,j,f,ijfy)
               call n_allhzijk(l,m,n,i,j,k,ijkz)
               call n_allhzijk(l,m,n,i,e,k,iekz)

               jx1=-hz(ijkz)*zk-hy(ijfy)*yjf+hz(iekz)*zk+hy(ijky)*yjk
               jx1=jx1/sijk

               call n_allhyijk(l,m,n,d,j,k,djky)
               call n_allhyijk(l,m,n,d,j,f,djfy)
               call n_allhzijk(l,m,n,d,j,k,djkz)
               call n_allhzijk(l,m,n,d,e,k,dekz)

               jx2=-hz(djkz)*zk-hy(djfy)*yjf+hz(dekz)*zk+hy(djky)*yjk
               jx2=jx2/sijk

               jx = ( jx1 + jx2 )/2.d0

!-----------------computation for Jy --------------------------

               if (j.ne.1) then

                  call n_allhxijk(l,m,n,i,j,k,ijkx)
                  call n_allhxijk(l,m,n,i,j,f,ijfx)
                  call n_allhzijk(l,m,n,i,j,k,ijkz)
                  call n_allhzijk(l,m,n,d,j,k,djkz)

                  call leng_xijk(i,j,k,x,y,z,xijk)
                  call leng_xijk(i,j,f,x,y,z,xijf)

                  jy1= hz(ijkz)*zk+hx(ijfx)*xijf-hz(djkz)*zk-hx(ijkx)*xijk

                  call area_sjki(i,j,k,x,y,z,sjki)

                  jy1=jy1/sjki

               end if

               if (j.ne.m) then

                  call n_allhxijk(l,m,n,i,e,k,iekx)
                  call n_allhxijk(l,m,n,i,e,f,iefx)
                  call n_allhzijk(l,m,n,i,e,k,iekz)
                  call n_allhzijk(l,m,n,d,e,k,dekz)

                  call leng_xijk(i,e,k,x,y,z,xiek)
                  call leng_xijk(i,e,f,x,y,z,xief)

                  jy2= hz(iekz)*zk+hx(iefx)*xief-hz(dekz)*zk-hx(iekx)*xiek

                  call area_sjki(i,e,k,x,y,z,seki)

                  jy2=jy2/seki

               end if

               jy = ( jy1 + jy2 )/2.d0

!-----------------computation for Jz --------------------------

               if (j.ne.1.and.j.ne.m) then

                  call n_allhxijk(l,m,n,i,j,k,ijkx)
                  call n_allhxijk(l,m,n,i,e,k,iekx)
                  call n_allhyijk(l,m,n,i,j,k,ijky)
                  call n_allhyijk(l,m,n,d,j,k,djky)

                  call leng_xijk(i,j,k,x,y,z,xijk)
                  call leng_xijk(i,e,k,x,y,z,xiek)
                  call leng_yijk(j,k,y,z,yjk)

                  jz=-hy(ijky)*yjk-hx(iekx)*xiek+hy(djky)*yjk+hx(ijkx)*xijk

                  call area_skij(i,j,k,x,y,z,skij)

                  jz=jz/skij

               else if (j == 1) then

                  call n_allhxijk(l,m,n,i,e,k,iekx)
                  call n_allhyijk(l,m,n,i,j,k,ijky)
                  call n_allhyijk(l,m,n,d,j,k,djky)

                  call leng_xijk(i,e,k,x,y,z,xiek)
                  call leng_yijk(j,k,y,z,yjk)

                  jz=-hy(ijky)*yjk-hx(iekx)*xiek+hy(djky)*yjk

                  call area_skij(i,j,k,x,y,z,skij)

                  jz=jz/skij

               else if (j == m) then

                  call n_allhxijk(l,m,n,i,j,k,ijkx)
                  call n_allhyijk(l,m,n,i,j,k,ijky)
                  call n_allhyijk(l,m,n,d,j,k,djky)

                  call leng_xijk(i,j,k,x,y,z,xijk)
                  call leng_yijk(j,k,y,z,yjk)

                  jz=-hy(ijky)*yjk+hy(djky)*yjk+hx(ijkx)*xijk

                  call area_skij(i,j,k,x,y,z,skij)

                  jz=jz/skij

               end if

!toh: note here that J is defined at the centre of a surface on a shell.
!toh: the point is shifted by dx/2 from that of C & D responses.

               ! (AK) write(iow,'(6g15.7)') jx,jy,jz

            end do

         end do

      end do

      return

      end subroutine allj_out
  

! *****************************************************************************



      subroutine surfeh_out(iow,l,m,n,nair,omega, &
                            x,y,z,resist,hx,hy,hz)


!-------------------------------------------------------------
!
!     compute block surface centered h and e field at top of ksurf layer
!     input h values are given on edges (parralel to edge) in real space
!           units
!     output hsx,hys,hzs,exs and eys are given in real space units
!     e field computation uses only fields in ksurf layer
!          extrapolation to top of layer is done using curle=iuw*h
!          and Ez=0
!-------------------------------------------------------------

      implicit none

      integer                            :: l,m,n,nair,iow
      complex(8),dimension(:)			 :: hx,hy,hz
      real(8),dimension(l)	             :: x
      real(8),dimension(m+1)             :: y
      real(8),dimension(n+1)             :: z
      real(8)                            :: omega
      real(8),dimension(l,m,n)	         :: resist

      complex(8)                         :: jx1,jx2,jy1,jy2,jz
      complex(8)                         :: exk,eyk,hzc,hzk,ciuw
      real(8)                            :: uw
      real(8)                            :: yp,zp
      real(8)                            :: xijk,xiek,xijf,xief,xipp,xipk
      real(8)                            :: yjk,yjf,yjp
      real(8)                            :: zc,zk
      real(8)                            :: sijk,sjki,seki,skij,sjpi,sijp
      integer,dimension(l,m+1)			 :: ksmap
      integer                            :: i,j,k,c,d,e,f
      integer                            :: ijkx,iekx,ijfx,iefx
      integer                            :: ijky,djky,ijfy,djfy
      integer                            :: ijcz,djcz,iecz,decz,ijkz
      integer                            :: djkz,iekz,dekz

      complex(8)                         :: hxs,hys,hzs,exs,eys,ezs
      real(8)                            :: xxw0,xxw,yyw,zzw

!      pi=4.0d0*datan(1.0d0)
!      r2d=180.0d0/pi
      uw=omega * MU_0
      ciuw=dcmplx(0.0d0,uw)

      do j=1,m
        e=j+1
        yyw=(y(j)+y(e))/2.0d0*r2d
        xxw0=0.0d0
        do i=1,l

          k=nair+1	  !k=ksmap(i,j)
          c=k-1
          f=k+1
          call leng_zijk(c,z,zc)
          call leng_zijk(k,z,zk)
          zzw=z(k)/1000.0d0

          xxw=( xxw0+x(i)/2.0d0 )*r2d
          xxw0=xxw0+x(i)
          if (i.ne.l) then
            d=d+1
          else
            d=1
          end if

          if (j.ne.1.and.j.ne.m) then
            call leng_xijk(i,j,k,x,y,z,xijk)
            call leng_xijk(i,e,k,x,y,z,xiek)
          else if (j == 1) then
            call leng_xijk(i,e,k,x,y,z,xiek)
          else if (j == m) then
            call leng_xijk(i,j,k,x,y,z,xijk)
          end if

!--------------------------------------------------------------
!
! compute surface magnetic fields
!
!-----------------computation for Hx --------------------------

          if (j.ne.1.and.j.ne.m) then
            call n_allhxijk(l,m,n,i,j,k,ijkx)
            call n_allhxijk(l,m,n,i,e,k,iekx)
            hxs=(hx(ijkx)*xijk+hx(iekx)*xiek)/(xijk+xiek)
!c          hsx(i,j)=(hx(ijkx)*xijk+hx(iekx)*xiek)/(xijk+xiek)
          else if (j == 1) then
            call n_allhxijk(l,m,n,i,e,k,iekx)
            hxs=hx(iekx)
!c          hsx(i,j)=hx(iekx)
          else if (j == m) then
            call n_allhxijk(l,m,n,i,j,k,ijkx)
            hxs=hx(ijkx)
!c          hsx(i,j)=hx(ijkx)
          end if

!-----------------computation for Hy --------------------------

          call n_allhyijk(l,m,n,i,j,k,ijky)
          call n_allhyijk(l,m,n,d,j,k,djky)
          hys=(hy(ijky)+hy(djky))/2.0d0
!c        hsy(i,j)=(hy(ijky)+hy(djky))/2.0

!-----------------computation for Hz --------------------------

          if (j.ne.1.and.j.ne.m) then
            call n_allhzijk(l,m,n,i,j,c,ijcz)
            call n_allhzijk(l,m,n,d,j,c,djcz)
            call n_allhzijk(l,m,n,i,e,c,iecz)
            call n_allhzijk(l,m,n,d,e,c,decz)
            hzc=(hz(ijcz)+hz(djcz)+hz(iecz)+hz(decz))/4.0d0
            call n_allhzijk(l,m,n,i,j,k,ijkz)
            call n_allhzijk(l,m,n,d,j,k,djkz)
            call n_allhzijk(l,m,n,i,e,k,iekz)
            call n_allhzijk(l,m,n,d,e,k,dekz)
            hzk=(hz(ijkz)+hz(djkz)+hz(iekz)+hz(dekz))/4.0d0
          else if (j == 1) then
            call n_allhzijk(l,m,n,i,j,c,ijcz)
            call n_allhzijk(l,m,n,i,e,c,iecz)
            call n_allhzijk(l,m,n,d,e,c,decz)
            hzc=(hz(ijcz)+hz(iecz)+hz(decz))/3.0d0
            call n_allhzijk(l,m,n,i,j,k,ijkz)
            call n_allhzijk(l,m,n,i,e,k,iekz)
            call n_allhzijk(l,m,n,d,e,k,dekz)
            hzk=(hz(ijkz)+hz(iekz)+hz(dekz))/3.0d0
          else if (j == m) then
            call n_allhzijk(l,m,n,i,j,c,ijcz)
            call n_allhzijk(l,m,n,d,j,c,djcz)
            call n_allhzijk(l,m,n,i,e,c,iecz)
            hzc=(hz(ijcz)+hz(djcz)+hz(iecz))/3.0d0
            call n_allhzijk(l,m,n,i,j,k,ijkz)
            call n_allhzijk(l,m,n,d,j,k,djkz)
            call n_allhzijk(l,m,n,i,e,k,iekz)
            hzk=(hz(ijkz)+hz(djkz)+hz(iekz))/3.0d0
          end if
          hzs=(hzc*zc+hzk*zk)/(zc+zk)
!c        hsz(i,j)=(hzc*zc+hzk*zk)/(zc+zk)

!--------------------------------------------------------------
!
! compute current densities
!
!-----------------computation for Jx --------------------------

          call leng_yijk(j,k,y,z,yjk)
          call leng_yijk(j,f,y,z,yjf)
          call area_sijk(j,k,y,z,sijk)

          call n_allhyijk(l,m,n,i,j,k,ijky)
          call n_allhyijk(l,m,n,i,j,f,ijfy)
          call n_allhzijk(l,m,n,i,j,k,ijkz)
          call n_allhzijk(l,m,n,i,e,k,iekz)
          jx1=-hz(ijkz)*zk-hy(ijfy)*yjf+hz(iekz)*zk+hy(ijky)*yjk
          jx1=jx1/sijk

          call n_allhyijk(l,m,n,d,j,k,djky)
          call n_allhyijk(l,m,n,d,j,f,djfy)
          call n_allhzijk(l,m,n,d,j,k,djkz)
          call n_allhzijk(l,m,n,d,e,k,dekz)
          jx2=-hz(djkz)*zk-hy(djfy)*yjf+hz(dekz)*zk+hy(djky)*yjk
          jx2=jx2/sijk

!-----------------computation for Jy --------------------------

          if (j.ne.1) then
            call n_allhxijk(l,m,n,i,j,k,ijkx)
            call n_allhxijk(l,m,n,i,j,f,ijfx)
            call n_allhzijk(l,m,n,i,j,k,ijkz)
            call n_allhzijk(l,m,n,d,j,k,djkz)
            call leng_xijk(i,j,k,x,y,z,xijk)
            call leng_xijk(i,j,f,x,y,z,xijf)
            jy1= hz(ijkz)*zk+hx(ijfx)*xijf-hz(djkz)*zk-hx(ijkx)*xijk
            call area_sjki(i,j,k,x,y,z,sjki)
            jy1=jy1/sjki
          end if

          if (j.ne.m) then
            call n_allhxijk(l,m,n,i,e,k,iekx)
            call n_allhxijk(l,m,n,i,e,f,iefx)
            call n_allhzijk(l,m,n,i,e,k,iekz)
            call n_allhzijk(l,m,n,d,e,k,dekz)
            call leng_xijk(i,e,k,x,y,z,xiek)
            call leng_xijk(i,e,f,x,y,z,xief)
            jy2= hz(iekz)*zk+hx(iefx)*xief-hz(dekz)*zk-hx(iekx)*xiek
            call area_sjki(i,e,k,x,y,z,seki)
            jy2=jy2/seki
          end if

!-----------------computation for Jz --------------------------

          if (j.ne.1.and.j.ne.m) then
            call n_allhxijk(l,m,n,i,j,k,ijkx)
            call n_allhxijk(l,m,n,i,e,k,iekx)
            call n_allhyijk(l,m,n,i,j,k,ijky)
            call n_allhyijk(l,m,n,d,j,k,djky)
            call leng_xijk(i,j,k,x,y,z,xijk)
            call leng_xijk(i,e,k,x,y,z,xiek)
            call leng_yijk(j,k,y,z,yjk)
            jz=-hy(ijky)*yjk-hx(iekx)*xiek+hy(djky)*yjk+hx(ijkx)*xijk
            call area_skij(i,j,k,x,y,z,skij)
            jz=jz/skij
          else if (j == 1) then
            call n_allhxijk(l,m,n,i,e,k,iekx)
            call n_allhyijk(l,m,n,i,j,k,ijky)
            call n_allhyijk(l,m,n,d,j,k,djky)
            call leng_xijk(i,e,k,x,y,z,xiek)
            call leng_yijk(j,k,y,z,yjk)
            jz=-hy(ijky)*yjk-hx(iekx)*xiek+hy(djky)*yjk
            call area_skij(i,j,k,x,y,z,skij)
            jz=jz/skij
          else if (j == m) then
            call n_allhxijk(l,m,n,i,j,k,ijkx)
            call n_allhyijk(l,m,n,i,j,k,ijky)
            call n_allhyijk(l,m,n,d,j,k,djky)
            call leng_xijk(i,j,k,x,y,z,xijk)
            call leng_yijk(j,k,y,z,yjk)
            jz=-hy(ijky)*yjk+hy(djky)*yjk+hx(ijkx)*xijk
            call area_skij(i,j,k,x,y,z,skij)
            jz=jz/skij
          end if

!--------------------------------------------------------------
!
! compute surface electric fields
!
!-----------------computation for Ex --------------------------

          zp=(z(k)+z(f))/2.0d0
          yp=(y(j)+y(e))/2.0d0
          
!         exk=(jx1+jx2)/2.0*resist(i,j,k)
          exk=(jx1+jx2)
          exk=exk*0.5d0
          exk=exk*dcmplx(resist(i,j,k))

          xipk=z(k)*dsin(yp)*x(i)
          xipp=zp  *dsin(yp)*x(i)
          sjpi=( z(k)**2 - zp**2 )*dsin(y(j))*x(i)/2.0d0
          exs=(xipp*exk+ciuw*hys*sjpi)/xipk
!c        esx(i,j)=(xipp*exk+ciuw*hsy(i,j)*sjpi)/xipk

!-----------------computation for Ey --------------------------

! NOTE - when using the -mp or -omp libraries, it appears that the
! compiler option -fpe1 (usually also with -fast) MUST be used
! or the underflow produces a fatal error. This does not happen
! under -wsf.

!         eyk=(jy1+jy2)/2.0*dcmplx(resist(i,j,k))
          eyk=jy1+jy2
          eyk=eyk*0.5d0
          eyk=eyk*dcmplx(resist(i,j,k))

          call leng_yijk(j,k,y,z,yjk)
          yjp=zp*(y(e)-y(j))
          sijp=( z(k)**2 - zp**2 )*( y(e)-y(j) )/2.0d0
          eys=(yjp*eyk-ciuw*hxs*sijp)/yjk
!c        esy(i,j)=(yjp*eyk-ciuw*hsx(i,j)*sijp)/yjk

!-----------------computation for Ez --------------------------

          ezs=jz*resist(i,j,k)
!c        esz(i,j)=jz*resist(i,j,k)
!--------------------------------------------------------------

        ! (AK) write(iow,'(15e15.7)')xxw,yyw,zzw,hxs,hys,hzs,exs,eys,ezs
       end do
      end do
      return
      end subroutine surfeh_out

end module output_orig
