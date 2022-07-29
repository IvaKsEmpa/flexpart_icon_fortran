! ----------------------------------------------------------------------
! Interpolation to point location byc 
! ----------------------------------------------------------------------
!
! o5/2022 - K.Ivanova, EMPA
! ----------------------------------------------------------------------
!
module precisions
  implicit none
  INTEGER, parameter :: dp=kind(1.0d0)
  logical     :: logik, ouvert
end module precisions

MODULE GlobalParam
  USE precisions
  IMPLICIT NONE
  real(kind=dp), parameter :: pi = 4._dp * atan(1._dp)
  REAL(kind=dp),PARAMETER :: plam_phi(2)=(/ 8.d0, 47.d0 /)
  ! local constants
  CHARACTER (LEN=*), PARAMETER :: grid_filename   = "/input/ICON/ivme/grid/ICON-1E_DOM01_tri.nc"
  CHARACTER (LEN=*), PARAMETER :: in_filename     = "/project/ivme/MCH-1/icon-art-BRM/icon_scripts_output_code/icon_art_output/ICON-ART-OEM_1504Unsrt_DOM01_00140000.grb"
  CHARACTER(LEN=10), PARAMETER :: open_mode='r'
  integer :: ncells, ncells_delaunay, nvertex, ii, i
  integer :: ncfileID, status, nn_tri, istat, ifield, arsize
  LOGICAL :: is_in, use_random
  INTEGER :: rfile, igrib
  CHARACTER(LEN=129) ::sn
  INTEGER :: j,varid, dimid
  CHARACTER (LEN=26):: dimname
END MODULE GlobalParam

PROGRAM interpol_delaunay
  USE precisions 
  USE GlobalParam
  use eccodes
  use netcdf
  IMPLICIT NONE
  REAL(kind=dp) :: point_loc(2), point_loc_cc(3), weights(3)
  ! local variables:
  REAL(kind=dp) :: r, v1(3), v2(3), v3(3), temp(1)
  real(kind=dp), allocatable:: lon(:), lat(:), cell_centers(:,:), vertex(:,:), vertex_of_cell(:,:)
  real(kind=dp), allocatable:: cell_cc(:,:), vertex_cc(:,:)
  real(kind=dp), allocatable   :: cell_of_triangle(:,:)
  v1=0.d0;v2=0.d0;v3=0.d0; is_in=.False.;use_random=.True.

  status = nf90_open(grid_filename, nf90_nowrite, ncfileID)
  CALL nf90_err(status)
  CALL READ_GRID()

  allocate(lon(ncells), &
        lat(ncells),    &
        cell_centers(ncells, 2),     &
        vertex(nvertex,2 ),          &
        vertex_of_cell(ncells,3),    &
        cell_of_triangle(ncells_delaunay,3))


  allocate(cell_cc(ncells,3), vertex_cc(nvertex,3))

  ! Get the varid of the data variable, based on its name.
  
  CALL GET_VAR_ID(lon, lat)
  
  ! cell=(lon,lat)
  cell_centers(:,1)=lon(1:ncells)
  cell_centers(:,2)=lat(1:ncells)

  deallocate(lon, lat)
  allocate(lon(nvertex), lat(nvertex))

  CALL GET_VLON_VLAT(LON,LAT)
  
  vertex(:,1)=lon(1:nvertex)
  vertex(:,2)=lat(1:nvertex)

  CALL GET_VERTEX_OF_CELL()
  CALL GET_CELL_OF_TRIANGLE(cell_of_triangle)

  nn_tri=ncells_delaunay  
  write(6,*) "nn_tri", nn_tri, "ncells_delaunay=", ncells_delaunay
  status = nf90_close(ncfileID)
  call nf90_err(status)

  ! convert to cartesian 
  vertex_cc = gc2cc(vertex, nvertex)
  cell_cc = gc2cc(cell_centers, ncells)

  write(*,*) "First vertex (lon/lat):", vertex(1,:) * 180/pi
  write(*,*) "First vertex (cart):   ", vertex_cc(1,:)
  
  write(*,*) "First cell (lon/lat):", cell_centers(1,:) * 180/pi
  write(*,*) "First cell (cart):   ", cell_cc(1,:)
  ! pick up random location from domain
if (use_random .eqv. .TRUE.) then
        call random_seed
        call random_number(r)
        write(*,*) "r", r
        
        temp=r*cell_centers(int(lbound(cell_centers(:,1))),1) + (1.d0-r)*cell_centers(int(ubound(cell_centers(:,1))),1)
        point_loc(1)=temp(1)
        
        temp=r*cell_centers(int(lbound(cell_centers(:,2))),2) + (1.d0-r)*cell_centers(int(ubound(cell_centers(:,2))),2)
        point_loc(2)= temp(1)
        point_loc_cc=gc2cc_vect(point_loc)
        
        write(*,*) "point_loc", point_loc
        write(*,*) "in cartesian random point", point_loc_cc
else
        point_loc_cc=gc2cc_vect(plam_phi/180*pi)
        write(*,*) "plam_phi, plam_phi/180*pi, point_loc",plam_phi, plam_phi/180*pi
endif 

CALL FIND_TRIANGLE(vertex_cc, vertex_of_cell, point_loc_cc, cell_cc, cell_of_triangle)
! calculate weights for barycentric interpolation
if (is_in) then 
        weights = bic_weights(point_loc_cc,         & 
          cell_cc(int(cell_of_triangle(ii,1)),1:3), &
          cell_cc(int(cell_of_triangle(ii,2)),1:3), &
          cell_cc(int(cell_of_triangle(ii,3)),1:3))
end if 

do i =1,3
if (dabs(weights(i)) .le. 1.d-8) then 
 weights(i) =0.d0
 end if
enddo

!check if sum of weights is equal to 1
write(6,*) "sum(weights)-1.0d0=", sum(weights)-1.d0,"weights=", weights

CALL INTERPOLATE_FIELD("U")
CALL INTERPOLATE_FIELD("V")
CALL INTERPOLATE_FIELD("HSURF")
#if !defined(NOMPI)
  CALL stop_mpi()
#endif

contains 
!----------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------
  FUNCTION gc2cc (p_pos, ncells)  RESULT(p)
         integer:: ncells, i
         REAL(kind=dp) :: p(ncells,3)                   ! Cart. coordinates
         REAL(kind=dp), INTENT(in) :: p_pos(ncells, 2)     ! geo. coordinates (lon,lat)
         REAL (kind=dp) :: z_clt(ncells), z_slt(ncells)
         REAL (kind=dp) :: z_sln(ncells), z_cln(ncells)
        ! z_clt = DCOS(p_pos(2))
        ! p_x(1:3) = (/ DCOS(p_pos(1))*z_clt, DSIN(p_pos(1))*z_clt, DSIN(p_pos(2)) /)
        
        do i =1, ncells 
                z_sln(i) = dsin(p_pos(i,1))
                z_cln(i) = dcos(p_pos(i,1))
                z_slt(i) = dsin(p_pos(i,2))
                z_clt(i) = dcos(p_pos(i,2))

                p(i,1) = z_cln(i)*z_clt(i)
                p(i,2) = z_sln(i)*z_clt(i)
                p(i,3) = z_slt(i)
        enddo 
  END FUNCTION gc2cc
!---------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------
  FUNCTION gc2cc_vect(p_pos)  RESULT(p)
         integer::  i
         REAL(kind=dp) :: p(3)                   ! Cart. coordinates
         REAL(kind=dp), INTENT(in) :: p_pos( 2)     ! geo. coordinates (lon,lat)
         REAL (kind=dp) :: z_clt, z_slt
         REAL (kind=dp) :: z_sln, z_cln
        ! z_clt = DCOS(p_pos(2))
        ! p_x(1:3) = (/ DCOS(p_pos(1))*z_clt, DSIN(p_pos(1))*z_clt, DSIN(p_pos(2)) /)

      
                z_sln = dsin(p_pos(1))
                z_cln = dcos(p_pos(1))
                z_slt = dsin(p_pos(2))
                z_clt = dcos(p_pos(2))

                p(1) = z_cln*z_clt
                p(2) = z_sln*z_clt
                p(3) = z_slt
  
  END FUNCTION gc2cc_vect
!---------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
FUNCTION inside_triangle(v, v1,v2,v3)
    LOGICAL :: inside_triangle
    REAL(kind=dp),       INTENT(IN)     :: v(1:3)          !< query point
    REAL(kind=dp),  INTENT(IN)     :: v1(1:3),v2(1:3),v3(1:3)      !< vertex longitudes/latitudes
! local variables
    LOGICAL       :: c1,c2,c3

    c1  = ccw_spherical(v1,v2,v)
    c2  = ccw_spherical(v2,v3,v)
    c3  = ccw_spherical(v3,v1,v)
   inside_triangle = ((      c1) .AND. (      c2) .AND. (      c3)) .OR. &
      &               ((.NOT. c1) .AND. (.NOT. c2) .AND. (.NOT. c3))
  END FUNCTION inside_triangle
!---------------------------------------------------------------------------------------------------------------
!       barycentric weights / coordinate transformation
!       All inputs in cartesian coordinates
FUNCTION bic_weights(p, v1, v2, v3) RESULT(w)
        REAL(kind=dp),       INTENT(IN)     :: p(1:3)          !< query point
        REAL(kind=dp),  INTENT(IN)     :: v1(1:3),v2(1:3),v3(1:3)      !< vertex longitudes/latitudes
        REAL(kind=dp):: w(3), a2
        a2 = v1(1) * (v2(2) - v3(2)) + v2(1) * (v3(2) - v1(2)) + v3(1) * (v1(2) - v2(2))

        w(1) = 1/a2 * (&
                &(v2(1)*v3(2) - v3(1)*v2(2)) +&
                &(v2(2) - v3(2)) * p(1) +&
                &(v3(1) - v2(1)) * p(2) )
        w(2) = 1/a2 * (&
                &(v3(1)*v1(2) - v1(1)*v3(2)) +&
                &(v3(2) - v1(2)) * p(1) +&
                &(v1(1) - v3(1)) * p(2) )
        w(3) = 1/a2 * (&
                &(v1(1)*v2(2) - v2(1)*v1(2)) +&
                &(v1(2) - v2(2)) * p(1) +&
                &(v2(1) - v1(1)) * p(2) )

        END FUNCTION

!---------------------------------------------------------------------------------------------------------------
!
! ! --------------------------------------------------------------------
!  !> Locates a point relative to a directed arc. Let v1, v2, and v3 be
!  !  distinct points on the sphere, and denote by v1 -> v2 the
!  !  geodesic connecting v1 and v2 and directed toward v2. This test
!  !  determines which of the two hemispheres defined by v1 -> v2
!  !  contains v3, or, in other words, if we "turn left" when going
!  !  from v1 to v2 to v3.
  FUNCTION ccw_spherical(v1,v2,v3)
    LOGICAL :: ccw_spherical
    REAL(kind=dp), INTENT(IN)  :: v1(1:3),v2(1:3),v3(1:3)
    REAL(kind=dp) :: ccw

    ! det(v1,v2,v3) = <v1 x v2, v3> = | v1 x v2 | cos(a) 
    !  
    ! where a is the angle between v3 and the normal to the plane
    ! defined by v1 and v2.
!    
    ccw =           v3(1)*(v1(2)*v2(3) - v2(2)*v1(3)) &
      &        -    v3(2)*(v1(1)*v2(3) - v2(1)*v1(3)) &
      &        +    v3(3)*(v1(1)*v2(2) - v2(1)*v1(2)) 
!
!    ! we apply a static error of
!    !   | e - e'| <= 3*2^-48
!    ! to decide if a floating-point evaluation e' of an expression e
!    ! has the correct sign, see Section 2.2 of
!    !
!    ! Burnikel, C.; Funke, S. & Seel, M. 
!    ! "Exact geometric computation using Cascading"
!    ! International Journal of Computational Geometry & Applications, 
!    ! World Scientific, 2001, 11, 245-266
    IF (DABS(ccw) <= 1.1d-14) THEN
      ccw_spherical = ccw_spherical_q128(v1,v2,v3)
    ELSE
      ccw_spherical = ccw <= 0.d0
    END IF
  END FUNCTION ccw_spherical
!! --------------------------------------------------------------------
!  !> Locates a point relative to a directed arc. Let v1, v2, and v3 be
!  !  distinct points on the sphere, and denote by v1 -> v2 the
!  !  geodesic connecting v1 and v2 and directed toward v2. This test
!  !  determines which of the two hemispheres defined by v1 -> v2
!  !  contains v3, or, in other words, if we "turn left" when going
!  !  from v1 to v2 to v3.
FUNCTION ccw_spherical_q128(v1,v2,v3)
    LOGICAL :: ccw_spherical_q128
    REAL(kind=dp), INTENT(IN)  :: v1(1:3),v2(1:3),v3(1:3)
    REAL(kind=dp) :: v1_x, v1_y, v1_z,v2_x, v2_y, v2_z,v3_x, v3_y, v3_z
!
!    ! det(v1,v2,v3) = <v1 x v2, v3> = | v1 x v2 | cos(a) 
!    !  
!    ! where a is the angle between v3 and the normal to the plane
!    ! defined by v1 and v2.
!    
    v1_x = v1(1)
    v1_y = v1(2)
    v1_z = v1(3)
    v2_x = v2(1)
    v2_y = v2(2)
    v2_z = v2(3)
    v3_x = v3(1)
    v3_y = v3(2)
    v3_z = v3(3)
!
    ccw_spherical_q128 = v3_x*(v1_y*v2_z - v2_y*v1_z) &
      &             -    v3_y*(v1_x*v2_z - v2_x*v1_z) &
      &             +    v3_z*(v1_x*v2_y - v2_x*v1_y)  <= 0.d0
  END FUNCTION ccw_spherical_q128
!------------------------------------------------------------------------------------------------------
  subroutine nf90_err(status)
    integer, intent (in) :: status
     if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop 'Stopped'
      end if
  end subroutine nf90_err
!------------------------------------------------------------------------------------------------------
  SUBROUTINE READ_GRID()
  USE precisions
  use GlobalParam
  IMPLICIT NONE
  
  ! open/read ICON grid definition 
  ! -----------------------------------------
  !read number of cells 
  status =nf90_inq_dimid(ncfileID, "cell", dimid)
  call nf90_err(status)
  status =nf90_inquire_dimension(ncfileID, dimid,dimname, ncells)
  call nf90_err(status)

  status =nf90_inq_dimid(ncfileID, "vertex", dimid)
  call nf90_err(status)
  status =nf90_inquire_dimension(ncfileID, dimid, dimname, nvertex)
  call nf90_err(status)

  status =nf90_inq_dimid(ncfileID, "cell_delaunay", dimid)
  call nf90_err(status)
  status =nf90_inquire_dimension(ncfileID, dimid,dimname, ncells_delaunay)
  call nf90_err(status)
  write(6,*) "ncells_delaunay", ncells_delaunay, "ncells", ncells, "nvertex", nvertex
  return
  END SUBROUTINE
!------------------------------------------------------------------------------------------------------
 SUBROUTINE GET_VAR_ID(lon, lat)
  USE precisions
  USE GlobalParam
  IMPLICIT NONE
  real(kind=dp):: lon(ncells), lat(ncells)
   ! Get the varid of the data variable, based on its name.
  status =nf90_inq_varid(ncfileID, "lon_cell_centre", varid)
  call nf90_err(status)
  status =nf90_get_var(ncfileID, varid, lon)
  call nf90_err(status)
  status =nf90_inq_varid(ncfileID, "lat_cell_centre", varid)
  call nf90_err(status)
  status =nf90_get_var(ncfileID, varid, lat)
  call nf90_err(status)

  return
  END SUBROUTINE
!------------------------------------------------------------------------------------------------------
SUBROUTINE GET_VLON_VLAT(LON,LAT)
  USE precisions
  USE GlobalParam
  IMPLICIT NONE
  real(kind=dp):: lon(nvertex), lat(nvertex)

   status =nf90_inq_varid(ncfileID, "vlon", varid)
  call nf90_err(status)
  status =nf90_get_var(ncfileID, varid, lon)
  call nf90_err(status)
  status =nf90_inq_varid(ncfileID, "vlat", varid)
  call nf90_err(status)
  status =nf90_get_var(ncfileID, varid, lat)
  call nf90_err(status)

return
END SUBROUTINE
!------------------------------------------------------------------------------------------------------
SUBROUTINE GET_VERTEX_OF_CELL()
  USE precisions
  USE GlobalParam
  IMPLICIT NONE
  real(kind=dp):: vertex_of_cell(ncells,3)

  status =nf90_inq_varid(ncfileID, "vertex_of_cell", varid)
  call nf90_err(status)
  status =nf90_get_var(ncfileID, varid, vertex_of_cell)
  call nf90_err(status)

return
END SUBROUTINE
!------------------------------------------------------------------------------------------------------
SUBROUTINE GET_CELL_OF_TRIANGLE(cell_of_triangle)
  USE precisions
  USE GlobalParam
  IMPLICIT NONE
  real(kind=dp):: cell_of_triangle(ncells_delaunay,3)

  status =nf90_inq_varid(ncfileID, "cc_delaunay", varid)
  call nf90_err(status)
  status =nf90_get_var(ncfileID, varid, cell_of_triangle)
  call nf90_err(status)

return
END SUBROUTINE
!------------------------------------------------------------------------------------------------------
SUBROUTINE FIND_TRIANGLE(vertex_cc, vertex_of_cell, point_loc_cc, cell_cc, cell_of_triangle)
  USE precisions
  USE GlobalParam
  IMPLICIT NONE
  real(kind=dp):: vertex_of_cell(ncells,3), vertex_cc(nvertex,3)
  real(kind=dp):: cell_of_triangle(ncells_delaunay,3), cell_cc(ncells,3)
  real(kind=dp):: point_loc_cc(3)
  REAL(kind=dp):: v1(3), v2(3), v3(3)
  v1=0.d0; v2=0.d0;v3=0.d0
 ncellsloop: do ii =1,ncells
  is_in = inside_triangle(point_loc_cc,&
               & vertex_cc(int(vertex_of_cell(ii,1)),1:3),&
               & vertex_cc(int(vertex_of_cell(ii,2)),1:3),&
               & vertex_cc(int(vertex_of_cell(ii,3)),1:3) )
  if (is_in) then
    write(6,*) "ncells loop:is_in", is_in, "ii=", ii
    exit ncellsloop
  endif

enddo ncellsloop

! find triangulation triangle containing point (simple loop)
ntriloop: do ii =1, ncells_delaunay

  v1=cell_cc(int(cell_of_triangle(ii,1)),1:3)
  v2=cell_cc(int(cell_of_triangle(ii,2)),1:3)
  v3=cell_cc(int(cell_of_triangle(ii,3)),1:3)

  is_in=  inside_triangle(point_loc_cc,v1,v2,v3)
  if (is_in) then
    write(*,*) "found containing triangle ii=", ii
    exit ntriloop
  endif
enddo ntriloop
return
END SUBROUTINE
!------------------------------------------------------------------------------------------------------
SUBROUTINE INTERPOLATE_FIELD(SHORT_NAME)
  USE precisions
  USE GlobalParam
  use eccodes
  use netcdf
  IMPLICIT NONE
  CHARACTER (LEN=*):: SHORT_NAME
  real(kind=dp), allocatable   :: FIELD(:)
  REAL(kind=dp) :: fp
  call grib_open_file(rfile,in_filename, open_mode,istat)
ifield=0         ! counter for number of fields found
  ! Loop over all grib fields of the grib file:
  gribfields: do
    ifield=ifield+1

    ! GET NEXT FIELDS

    call grib_new_from_file(rfile,igrib,istat)

    call grib_get(igrib, "shortName",sn, istat)
    if (sn ==SHORT_NAME) then
        exit gribfields
    end if

 end do gribfields
 call grib_get_size(igrib,'values',arsize,istat)

allocate(FIELD(arsize))

call grib_get(igrib,'values',FIELD,istat)


write(6,*) SHORT_NAME, maxval(FIELD), minval(FIELD)

fp =0
do j =1,3
 fp =fp+ FIELD(int(cell_of_triangle(ii,j)))*weights(j)
 write(*,*)  "int(cell_of_triangle(ii,j))", int(cell_of_triangle(ii,j)), SHORT_NAME, FIELD(int(cell_of_triangle(ii,j)))

enddo
! write(6,*) "hsurf:", hsurf(int(cell_of_triangle(ii,1))), hsurf(int(cell_of_triangle(ii,2))), hsurf(int(cell_of_triangle(ii,3)))
 write(*,*) "INTERPOLATED VALUE OF ",SHORT_NAME, fp
return
END SUBROUTINE
!------------------------------------------------------------------------------------------------------
END PROGRAM interpol_delaunay

