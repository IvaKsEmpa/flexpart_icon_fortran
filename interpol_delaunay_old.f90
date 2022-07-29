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

PROGRAM interpol_delaunay
  USE precisions 
  use eccodes
  use netcdf

  IMPLICIT NONE

  ! local constants
  CHARACTER (LEN=*), PARAMETER :: grid_filename   = "/input/ICON/ivme/grid/ICON-1E_DOM01_tri.nc"
  CHARACTER (LEN=*), PARAMETER :: in_filename     = "/project/ivme/MCH-1/icon-art-BRM/icon_scripts_output_code/icon_art_output/ICON-ART-OEM_1504Unsrt_DOM01_00140000.grb"
  REAL(kind=dp) :: point_loc(2), point_loc_cc(3), weights(3)
  REAL(kind=dp),          PARAMETER :: plam_phi(2)     = (/ 8.d0, 47.d0 /)
  CHARACTER (LEN=26):: dimname
  REAL (kind=dp), PARAMETER ::  pi        = 3.14159265358979323846264338327950288d0
  REAL (kind=dp), PARAMETER ::  pi_180    = pi/180.0d0

  ! local variables:
  integer :: ncfileID, status, nn_tri, n, iret, istat, ifield, arsize
  INTEGER :: i,j, index, min_node_idx(3), jc, jb, jc_n, jb_n, varid, dimid 
  INTEGER  ::  neighbor(2,2), interpol_cell(3,2) 
  INTEGER :: jc_v, jb_v, min_loc_vertex, cells_of_vertex(6), jc_min, jb_min
  REAL(kind=dp) :: distance(3), v, r, v1(3), v2(3), v3(3)
  REAL(kind=dp) :: A(2),B(2),C(2),AB(3),AC(3),AP(3),BC(3), PC(3), BP(3)
  REAL(kind=dp) :: vertex1(2), vertex2(2), vertex3(2)
  REAL(kind=dp) :: weights_u,weights_v,weights_w,fa,fb,fc,fP
  REAL(kind=dp) :: AREA_ABC,AREA_ABP,AREA_CBP,AREA_CPA
  real(kind=dp), allocatable:: lon(:), lat(:), cell_centers(:,:), vertex(:,:), vertex_of_cell(:,:)
  real(kind=dp), allocatable:: cell_cc(:,:), vertex_cc(:,:)
  LOGICAL       :: is_in
  INTEGER :: rfile, igrib
  CHARACTER(LEN=10), PARAMETER :: open_mode='r'
  CHARACTER(LEN=129) ::sn 

  integer    ::  errstat, ncells, ncells_delaunay, ii
  real(kind=dp)   ::  config_pole(2) = (/ -180.0d0,  90.0d0 /)
  integer    ::  config_nxpoints = 100
  integer    ::  config_nypoints = 100
  real(kind=dp)   ::  config_corner(2) = (/ 8.0d0, 47.0d0 /) 
  real(kind=dp)   ::  config_dx = 0.01d0
  real(kind=dp)   ::  config_dy = 0.01d0
  real(kind=dp), allocatable   :: hsurf(:), zsfc_in(:,:), zsfc_out(:,:), rfield1D(:), rfield2D(:,:), cell_of_triangle(:,:)
  distance =0.d0; v1=0.d0;v2=0.d0;v3=0.d0
  is_in=.False.
#if !defined(NOMPI)
  CALL start_mpi()
#endif
  ! open/read ICON grid definition 
  ! -----------------------------------------
  status = nf90_open(grid_filename, nf90_nowrite, ncfileID)
  !read number of cells 
  status =nf90_inq_dimid(ncfileID, "cell", dimid)
  status =nf90_inquire_dimension(ncfileID, dimid,dimname, ncells)

  status =nf90_inq_dimid(ncfileID, "cell_delaunay", dimid)
  status =nf90_inquire_dimension(ncfileID, dimid,dimname, ncells_delaunay)
  write(6,*) "ncells_delaunay", ncells_delaunay, "ncells", ncells, "dimid", dimid

  allocate(lon(ncells), lat(ncells), cell_centers(ncells, 2), vertex(ncells,2 ), vertex_of_cell(ncells,3),&
&  cell_of_triangle(ncells_delaunay,3))
  allocate(cell_cc(ncells,3), vertex_cc(ncells,3))
  ! Get the varid of the data variable, based on its name.
  status =nf90_inq_varid(ncfileID, "lon_cell_centre", varid) 
  status =nf90_get_var(ncfileID, varid, lon)

  status =nf90_inq_varid(ncfileID, "lat_cell_centre", varid)
  status =nf90_get_var(ncfileID, varid, lat)
  ! cell=(lon,lat)
  cell_centers(:,1)=lon(:)
  cell_centers(:,2)=lat(:)

  status =nf90_inq_varid(ncfileID, "vlon", varid)
  status =nf90_get_var(ncfileID, varid, lon)

  status =nf90_inq_varid(ncfileID, "vlat", varid)
  status =nf90_get_var(ncfileID, varid, lat)
  ! cell=(lon,lat)
  vertex(:,1)=lon(:)
  vertex(:,2)=lat(:)

  status =nf90_inq_varid(ncfileID, "vertex_of_cell", varid)
  status =nf90_get_var(ncfileID, varid, vertex_of_cell)
  
  status =nf90_inq_varid(ncfileID, "cc_delaunay", varid)
  status =nf90_get_var(ncfileID, varid, cell_of_triangle)
  nn_tri=size(cell_of_triangle, 1)
  write(6,*) "nn_tri", nn_tri
  status = nf90_close(ncfileID)
  ! convert to cartesian 
  vertex_cc = gc2cc(vertex, ncells)
  cell_cc = gc2cc(cell_centers, ncells)
  
  ! pick up random location from domain

  call random_seed
  call random_number(r)
write(*,*) "r", r
 point_loc(1)=cell_centers((int(r*size(cell_centers(:,1)))+1 ), 1)
 point_loc(2)=cell_centers((int(r*size(cell_centers(:,2)))+1), 2)
 point_loc_cc=gc2cc_vect(point_loc)
 !do i=1,2
 !plam_phi(i) = plam_phi(i)*pi_180
 !enddo
 point_loc_cc=gc2cc_vect(plam_phi)
 write(*,*) "point_loc ", plam_phi
 write(*,*) "in cartesian point_loc_cc point ", point_loc_cc

! find triangulation triangle containing point (simple loop)

! ncellsloop: do ii =1,ncells
 !       is_in = inside_triangle(point_loc_cc,&
  !             & vertex_cc(int(vertex_of_cell(ii,1)),1:3),&
   !            & vertex_cc(int(vertex_of_cell(ii,2)),1:3),&
    !           & vertex_cc(int(vertex_of_cell(ii,3)),1:3) )

     !   if (is_in) then 
      !          write(6,*) "ncells loop:is_in", is_in, "ii=", ii
 ! exit ncellsloop
 ! endif
! enddo ncellsloop

 ! find triangulation triangle containing point (simple loop)
ntriloop: do ii =1, nn_tri
v1=cell_cc(int(cell_of_triangle(ii,1)),1:3)
v2=cell_cc(int(cell_of_triangle(ii,2)),1:3)
v3=cell_cc(int(cell_of_triangle(ii,3)),1:3)

is_in=  inside_triangle(point_loc_cc,v1,v2,v3)

! call MY_BYC_WEIGHTS(point_loc_cc,V1,V2,V3, weights_u, weights_v, weights_w)
 !      if (dabs(weights_u+weights_v+weights_w-1.d0).le. 1.d-8) then
 !        write(*,*) "found containing triangle with my weights ii=", ii,"weights=", weights_u, weights_v,weights_w
  !     exit ntriloop  
   !    endif 
       if (is_in) then 
        write(*,*) "found containing triangle ii=", ii
  exit ntriloop
  endif 
 enddo ntriloop
write(6,*) "found containing triangle number ii=", ii
! calculate weights for barycentric interpolation
if (is_in) then 
        weights = bic_weights(point_loc_cc, cell_cc(int(cell_of_triangle(ii,1)),1:3),&
             & cell_cc(int(cell_of_triangle(ii,2)),1:3),&
              & cell_cc(int(cell_of_triangle(ii,3)),1:3))
 
    V1=cell_cc(int(cell_of_triangle(ii,1)),1:3)
     V2=cell_cc(int(cell_of_triangle(ii,2)),1:3)
     V3=cell_cc(int(cell_of_triangle(ii,3)),1:3)
     
     call MY_BYC_WEIGHTS(point_loc_cc,V1,V2,V3, weights_u, weights_v, weights_w)
write(*,*) "my weights:", weights_u, weights_v, weights_w, "sum:", weights_u+weights_v+weights_w
end if 

do i =1,3
if (dabs(weights(i)) .le. 1.d-6) then 
 weights(i) =0.d0
 end if
enddo

!check if sum of weights is equal to 1
write(6,*) "sum(weights)-1.0d0=", sum(weights)-1.d0,"weights=", weights
call grib_open_file(rfile,in_filename, open_mode,istat)

ifield=0         ! counter for number of fields found
  ! Loop over all grib fields of the grib file:
  gribfields: do

    ifield=ifield+1

    ! GET NEXT FIELDS

    call grib_new_from_file(rfile,igrib,istat)

    call grib_get(igrib, "shortName",sn, istat)

    if (sn =="HSURF") then 
        exit gribfields
    end if 
!    write(6,*) "sn", sn 

 end do gribfields 
 call grib_get_size(igrib,'values',arsize,istat)

allocate(hsurf(arsize)) 

call grib_get(igrib,'values',hsurf,istat)


write(6,*) "hsurf=", maxval(hsurf), minval(hsurf)

fp =0
do j =1,3 
 fp =fp+ hsurf(int(cell_of_triangle(ii,j)))*weights(j)
 write(*,*)  "int(cell_of_triangle(ii,j))", int(cell_of_triangle(ii,j)), "hsurf", hsurf(int(cell_of_triangle(ii,j)))

enddo 
! write(6,*) "hsurf:", hsurf(int(cell_of_triangle(ii,1))), hsurf(int(cell_of_triangle(ii,2))), hsurf(int(cell_of_triangle(ii,3)))
 write(*,*) "fp:", fp
#if !defined(NOMPI)
  CALL stop_mpi()
#endif

contains 

        !----------------------------------------------------------------------------------------------------------------
        SUBROUTINE TRIANGLE_AREA(A,B,C, AREA)
          USE precisions
          real(kind=dp)   :: A,B,C, AREA, S
          S =0.5D0*(A+B+C)
          AREA=DSQRT(S*(S-A)*(S-B)*(S-C))
        END SUBROUTINE
!-----------------------------------------------------------------------------------------------------------------
  FUNCTION gc2cc (p_pos, ncells)  RESULT(p)
          USE precisions
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
          USE precisions
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
 !> type declaration: single point in R^3
!  TYPE t_point
 !   USE precisions
  !  REAL(kind=dp)       :: x,y,z, ps
   ! INTEGER        :: gindex                      !< global index
 ! CONTAINS
  !  PROCEDURE :: norm2          => point_norm2
  !  PROCEDURE :: spherical_dist => point_spherical_dist
 ! END TYPE t_point
!-------------------------------------------------------------------------------------------------------------------
FUNCTION inside_triangle(v, v1,v2,v3)
    USE precisions
    LOGICAL :: inside_triangle
    REAL(kind=dp),       INTENT(IN)     :: v(3)          !< query point
    REAL(kind=dp),  INTENT(IN)     :: v1(3),v2(3),v3(3)     !< vertex longitudes/latitudes
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
        USE precisions
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

 SUBROUTINE MY_BYC_WEIGHTS(P,V1,V2,V3, weights_u, weights_v, weights_w)
 USE precisions
 REAL(kind=dp) :: AB(3),AC(3),AP(3),BC(3), PC(3), BP(3)
 REAL(kind=dp) :: P(3), V1(3),V2(3), V3(3)
 REAL(kind=dp) ::weights_u, weights_v, weights_w
        
 AB(1)=V2(1)-V1(1)
 AB(2)=V2(2)-V1(2)
 AB(3)=V2(3)-V1(3)

 AC(1)=V3(1)-V1(1)
 AC(2)=V3(2)-V1(2)
 AC(3)=V3(3)-V1(3)

 BC(1)=V3(1)-V2(1)
 BC(2)=V3(2)-V2(2)
 BC(3)=V3(3)-V2(3)

 AP(1)=P(1)-V1(1)
 AP(2)=P(2)-V1(2)
 AP(3)=P(3)-V1(3)

 BP(1)=P(1)-V2(1)
 BP(2)=P(2)-V2(2)
 BP(3)=P(3)-V2(3)

 PC(1)=P(1)-V3(1)
 PC(2)=P(2)-V3(2)
 PC(3)=P(3)-V3(3)
 !weights_u =NORM2(cross(AB,AP)/cross(AB,AC))
 CALL TRIANGLE_AREA(NORM2(AB),NORM2(AC),NORM2(BC),AREA_ABC)
 CALL TRIANGLE_AREA(NORM2(AB),NORM2(AP),NORM2(BP),AREA_ABP)
 CALL TRIANGLE_AREA(NORM2(BC),NORM2(BP),NORM2(PC),AREA_CBP)
 CALL TRIANGLE_AREA(NORM2(AC),NORM2(AP),NORM2(PC), AREA_CPA)
 !write(*,*) "AREA_ABC, AREA_ABP,AREA_CBP, AREA_CPA", AREA_ABC, AREA_ABP,AREA_CBP, AREA_CPA
 if (AREA_ABC .gt. 0.d0) then 
 weights_u=AREA_CBP/AREA_ABC
 weights_v=AREA_CPA/AREA_ABC
 weights_w=AREA_ABP/AREA_ABC
 else 
 weights_u=999999.d0
 weights_v=999999.d0
 weights_w=999999.d0
 
 endif
 END SUBROUTINE 
!---------------------------------------------------------------------------------------------------------------
!
! ! --------------------------------------------------------------------
!  !> Locates a point relative to a directed arc. Let v1, v2, and v3 be
!  !  distinct points on the sphere, and denote by v1 -> v2 the
!  !  geodesic connecting v1 and v2 and directed toward v2. This test
!  !  determines which of the two hemispheres defined by v1 -> v2
!  !  contains v3, or, in other words, if we "turn left" when going
!  !  from v1 to v2 to v3.
PURE  FUNCTION ccw_spherical(v1,v2,v3)
    USE precisions
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
PURE FUNCTION ccw_spherical_q128(v1,v2,v3)
    USE precisions
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

END PROGRAM interpol_delaunay

