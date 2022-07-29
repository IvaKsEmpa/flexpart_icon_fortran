! ----------------------------------------------------------------------
! EXAMPLE04
!
! Example program for the DWD ICON tools library
! ----------------------------------------------------------------------
!
! o3/2022 - K.Ivanova, EMPA
!
! This example program shows
! - how to load ICON grids
! - how to translate (lam/phi) -> grid index
! - how to get cell center and vertex coordinates
!
! This example program requires the DWD ICON tools libraries,
! which can be created by
!   cd ../../icontools/; make lce_intel; cd ../example/programs
!
! A short description of the ICON grid structure can be found in
! Chapter 2 and in Section 3.1 of the ICON database documentation,
! see
!   http://www.dwd.de/SharedDocs/downloads/DE/modelldokumentationen/nwv/icon/icon_dbbeschr_aktuell.pdf
!
! ----------------------------------------------------------------------
!
PROGRAM interpol_01

#if !defined(NOMPI)
  USE mo_mpi,               ONLY: start_mpi, stop_mpi
#endif
  USE mo_utilities,         ONLY: wp, t_geographical_coordinates, pi_180, blk_no, idx_no, & 
                                    nproma
  USE mo_remap_sync,        ONLY: t_gather, allocate_gather_c
  USE mo_metadata_types,    ONLY: t_file_metadata
  USE mo_remap_grid_types,  ONLY: t_grid, t_grid_metadata
  USE mo_remap_grid_query,  ONLY: init_query_algorithm
  USE mo_kdtree,            ONLY: kdtree_query
  !USE  kdtree_dist
  USE mo_icongpi,           ONLY: icongpi_setup, icongpi_finalize
  use mo_remap_config,      ONLY: GRID_TYPE_ICON, GRID_TYPE_REGULAR_LOC, MAX_NSTENCIL_RBF_SCALAR
  use mo_remap_io,          ONLY: load_grid, open_file
  use mo_cdi,               ONLY: GRID_LONLAT
  use mo_remap_input,       ONLY: read_field2D
  use mo_remap_intp_types,  ONLY: t_intp_data, allocate_intp_data
  use mo_remap_intp,        ONLY: interpolate_c_2D
  use mo_remap_weights_nnb, ONLY: prepare_intp_nnb
  use mo_remap_weights_rbf, ONLY: prepare_intp_rbf_scalar_stencil

  use netcdf, only: nf90_open, nf90_close, nf90_nowrite

  IMPLICIT NONE

  ! local constants
  CHARACTER (LEN=*), PARAMETER :: grid_filename   = "/input/ICON/ivme/grid/ICON-1E_DOM01.nc"
  CHARACTER (LEN=*), PARAMETER :: in_filename     = "/project/ivme/MCH-1/icon-art-BRM/icon_scripts_output_code/icon_art_output/ICON-ART-OEM_1504Unsrt_DOM01_00140000.grb"
  REAL(wp),          PARAMETER :: plam_phi(2)     = (/ 8._wp, 47._wp /) ! search point
  TYPE(t_geographical_coordinates) :: p, cell_center, vertex

  ! local variables:
  integer                          :: ncfileID, status
  TYPE (t_file_metadata)           :: file_metadata
  INTEGER                          :: i,j, index, min_node_idx(3), jc, jb, jc_n, jb_n
 INTEGER                           ::  neighbor(2,2), interpol_cell(3,2) 
  INTEGER                          :: jc_v, jb_v, min_loc_vertex, cells_of_vertex(6), jc_min, jb_min
  REAL(wp)                         :: min_dist, distance(3), v
  REAL(wp)                         :: A(2),B(2),C(2),AB(2),AC(2),AP(2),BC(2), PC(2), BP(2)
  REAL(wp)                         :: vertex1(2), vertex2(2), vertex3(2)
  REAL(wp)                         :: weights_u,weights_v,weights_w,fa,fb,fc,fP
  REAL(wp)                         :: AREA_ABC,AREA_ABP,AREA_CBP,AREA_CPA
  TYPE (t_grid)                    :: gridA, gridB
  TYPE (t_grid_metadata)           :: output_grid
  type (t_intp_data)               :: intp_data_nnb_B, intp_data_rbf_B

  integer    ::  errstat
  real(wp)   ::  config_pole(2) = (/ -180._wp,  90._wp /)
  integer    ::  config_nxpoints = 100
  integer    ::  config_nypoints = 100
  real(wp)   ::  config_corner(2) = (/ 8._wp, 46._wp /) 
  real(wp)   ::  config_dx = 0.01
  real(wp)   ::  config_dy = 0.01
  real(wp), allocatable   ::  zsfc_in(:,:), zsfc_out(:,:), rfield1D(:), rfield2D(:,:)
  type(t_gather) :: gather_c
  distance =0.d0

#if !defined(NOMPI)
  CALL start_mpi()
#endif

  ! open/read ICON grid definition 
  ! -----------------------------------------
  status = nf90_open(grid_filename, nf90_nowrite, ncfileID)
  CALL load_grid(gridA, GRID_TYPE_ICON, rank0=0, opt_fileID=ncfileID)
  status = nf90_close(ncfileID)
  CALL init_query_algorithm(gridA, kdtree=.TRUE.)
 ! icondelaunay -g grid_filename --vtk tri.vtk -v  
  ! define a new output grid
  ! -----------------------------------------
  output_grid%grid_type = GRID_LONLAT
  output_grid%cell_type = 4
  output_grid%nxpoints = config_nxpoints
  output_grid%nypoints = config_nypoints
  output_grid%north_pole(:) = config_pole(:)
  ALLOCATE(output_grid%xvals(output_grid%nxpoints), output_grid%yvals(output_grid%nypoints))
  output_grid%xvals = (/ ( config_corner(1) + config_dx*(i-1), i=1,output_grid%nxpoints ) /)
  output_grid%yvals = (/ ( config_corner(2) + config_dy*(i-1), i=1,output_grid%nypoints ) /)
  CALL load_grid(gridB, GRID_TYPE_REGULAR_LOC, rank0=0, opt_grid_metadata=output_grid)

  ! allocate/prepare interpolation data 
  ! -----------------------------------------
  CALL allocate_intp_data(intp_data_nnb_B, gridB%p_patch%nblks_c, gridB%p_patch%npromz_c, 1)
  CALL allocate_intp_data(intp_data_rbf_B, gridB%p_patch%nblks_c, gridB%p_patch%npromz_c, &
    MAX_NSTENCIL_RBF_SCALAR)
  ! first nnb, which is required for getting rbf weights
  CALL prepare_intp_nnb(gridA, gridB, intp_data_nnb_B)
  CALL prepare_intp_rbf_scalar_stencil(gridA, gridB, intp_data_nnb_B, intp_data_rbf_B)


  ! allocate 2D fields
  ! -----------------------------------------
  allocate(zsfc_in(nproma, gridA%p_patch%nblks_c), stat=errstat)
  allocate(zsfc_out(nproma, gridB%p_patch%nblks_c), stat=errstat)

  ! read from ICON output file
  ! -----------------------------------------
  ! prepare
  CALL allocate_gather_c("cell comm. pattern gridA_cov", gridA, 0, gather_c)
  allocate(rfield2D(nproma, gather_c%nblks), stat=errstat)
  allocate(rfield1D(gather_c%total_dim), stat=errstat)
  ! open data file
  CALL open_file(in_filename, GRID_TYPE_ICON, file_metadata, rank0=0)
  ! read a field
  CALL read_field2D(file_metadata, "HSURF", 6, gather_c, rfield1D, rfield2D, zsfc_in)

  ! interpolation to regular grid
  ! -----------------------------------------
  CALL interpolate_c_2D(zsfc_in, zsfc_out, intp_data_rbf_B, 0._wp)

  ! some field stats
  ! -----------------------------------------
  write(*,*) shape(zsfc_in)
  write(*,*) minval(zsfc_in), maxval(zsfc_in)

  write(*,*) shape(zsfc_out)
  write(*,*) minval(zsfc_out), maxval(zsfc_out)

  ! TODO: write output to file


  ! query point, give field value at nearest grid cell and neighbors
  ! -----------------------------------------------

  ! translation  lam/phi    ->  index of nearest cell
  p%lon = plam_phi(1) * pi_180
  p%lat = plam_phi(2) * pi_180
  CALL kdtree_query(gridA%kdtree, p, min_dist, min_node_idx) ! global index
  index = min_node_idx(3)
  jc = min_node_idx(1) ! translate global index into pair (index/block)
  jb = min_node_idx(2)

  cell_center = gridA%p_patch%cells%center(jc,jb)
  WRITE (0,*) "lam/pi      ->  nearest cell index = ", index
  WRITE (0,'(2(a,F6.3))') " grid index  ->  lam/phi = ", cell_center%lon, ", ", cell_center%lat
  write (*,*) "Distance:", min_dist, "min_node_idx:", min_node_idx
  write (*,*) "zsfc_in(jc, jb):", zsfc_in(jc, jb)
 
 ! print the cell neighbors' coordinates
  DO i=1,3
    jc_n = gridA%p_patch%cells%neighbor_idx(jc,jb,i)
    jb_n = gridA%p_patch%cells%neighbor_blk(jc,jb,i)
    cell_center = gridA%p_patch%cells%center(jc_n,jb_n)
    if (i==1) then 
            A=(/cell_center%lon, cell_center%lat/)
            fa= zsfc_in(jc_n, jb_n)
     
    elseif (i==2) then 
            B=(/cell_center%lon, cell_center%lat/)
            fb=zsfc_in(jc_n, jb_n)

    elseif (i==3) then 
            C=(/cell_center%lon, cell_center%lat/)
            fc=zsfc_in(jc_n, jb_n)
    endif
    WRITE (0,'(a,i0,2(a,F6.3),a,F10.2)') " neighbor_center ", i, " : ", &
        cell_center%lon, ", ", cell_center%lat, ", ", zsfc_in(jc_n, jb_n)
  END DO
  WRITE (6,*) "the wrong A,B,C,plam_phi=",A,B,C,plam_phi
  ! print the cell vertex coordinates
  DO i=1,3
    jc_v = gridA%p_patch%cells%vertex_idx(jc,jb,i)
    jb_v = gridA%p_patch%cells%vertex_blk(jc,jb,i)
    vertex = gridA%p_patch%verts%vertex(jc_v,jb_v)
    if (i==1) vertex1 = (/vertex%lon, vertex%lat/)
    if (i==2) vertex2 = (/vertex%lon, vertex%lat/)
    if (i==3) vertex3 = (/vertex%lon, vertex%lat/)
    WRITE (0,'(a,i0,2(a,F6.3))') " vertex ", i, " : ", vertex%lon, ", ", vertex%lat
    distance(i)=SQRT((plam_phi(1)-vertex%lon)**2.d0 +(plam_phi(2)-vertex%lat)**2.d0 )
    WRITE (0,'(a,i0,2(a,F6.3))') " distance to vertex", i, " : ", distance(i)
  END DO
  min_loc_vertex =MINLOC(distance, dim=1)
  jc_min = gridA%p_patch%cells%vertex_idx(jc,jb,min_loc_vertex)
  jb_min =gridA%p_patch%cells%vertex_blk(jc,jb,min_loc_vertex)
  WRITE (6,*) "the location of closest vertex of the cell to the particle", min_loc_vertex
  cells_of_vertex=gridA%p_patch%verts%cell_idx(jc_min, jb_min, :)
  do i=1, 3
    jc_n = gridA%p_patch%cells%neighbor_idx(jc,jb,i)
    jb_n = gridA%p_patch%cells%neighbor_blk(jc,jb,i)
    write(6,*) "jc_n, jb_n", jc_n, jb_n
    if (ANY(cells_of_vertex==jc_n)) then
           neighbor(i,:)=(/jc_n, jb_n/) 
    end if  
  end do  
  write(6,*) "neighbor", neighbor
  DO j=1,2
   DO i=1,3
    jc_n = gridA%p_patch%cells%neighbor_idx(neighbor(j,1),neighbor(j,2),i)
    jb_n = gridA%p_patch%cells%neighbor_blk(neighbor(j,1),neighbor(j,2),i)
    if (ANY(cells_of_vertex==jc_n) .and. (jc_n .ne. jc))  then
           interpol_cell(j,:)=(/jc_n, jb_n/)
    end if

    end do
  end do 
  interpol_cell(3,:)=(/jc, jb/)
  write(6,*) "interpolation cells", interpol_cell

   do i =1,3
    cell_center = gridA%p_patch%cells%center(interpol_cell(i,1),interpol_cell(i,2))
    if (i==1) then
            A=(/cell_center%lon, cell_center%lat/)
            fa= zsfc_in(jc_n, jb_n)

    elseif (i==2) then
            B=(/cell_center%lon, cell_center%lat/)
            fb=zsfc_in(jc_n, jb_n)

    elseif (i==3) then
            C=(/cell_center%lon, cell_center%lat/)
            fc=zsfc_in(jc_n, jb_n)
    endif
    end do 
 WRITE (6,*) "The right A,B,C,plam_phi=",A,B,C,plam_phi

 ! WRITE (6,*) "cells_of_vertex", gridA%p_patch%edges%cell_idx(min_loc_vertex,jc,jb)
 WRITE (6,*) "cells_of_vertex", gridA%p_patch%verts%cell_idx(jc_min, jb_min, :)
!calculate the weights
 AB=(/B(1)-A(1), B(2)-A(2)/)
 AC=(/C(1)-A(1), C(2)-A(2)/)
 BC=(/C(1)-B(1), C(2)-B(2)/)
 AP=(/PLAM_PHI(1)-A(1), PLAM_PHI(2)-A(2)/) 
 BP=(/PLAM_PHI(1)-B(1), PLAM_PHI(2)-B(2)/)
 PC=(/PLAM_PHI(1)-C(1), PLAM_PHI(2)-C(2)/)
 !weights_u =NORM2(cross(AB,AP)/cross(AB,AC))
 CALL TRIANGLE_AREA(NORM2(AB),NORM2(AC),NORM2(BC),AREA_ABC)
 CALL TRIANGLE_AREA(NORM2(AB),NORM2(AP),NORM2(BP),AREA_ABP)
 CALL TRIANGLE_AREA(NORM2(BC),NORM2(BP),NORM2(PC),AREA_CBP)
 CALL TRIANGLE_AREA(NORM2(AC),NORM2(AP),NORM2(PC), AREA_CPA)
 weights_u=AREA_CBP/AREA_ABC
 weights_v=AREA_CPA/AREA_ABC
 weights_w=AREA_ABP/AREA_ABC
 WRITE (6,*) "weights=", weights_u, weights_v, weights_w, weights_u+weights_v+weights_w
! f(P)=uf(A)+vf(B)+wf(C)
fP=weights_u*fa+weights_v*fb+weights_w*fc
WRITE (6,*) "interpolated value at point location", fP
#if !defined(NOMPI)
  CALL stop_mpi()
#endif
contains 
        SUBROUTINE TRIANGLE_AREA(A,B,C, AREA)
          real(wp)   :: A,B,C, AREA, S
          S =0.5D0*(A+B+C)
          AREA=SQRT(S*(S-A)*(S-B)*(S-C))
        END SUBROUTINE

    PURE FUNCTION gc2cc (p_pos)  RESULT(p_x)
         REAL(wp) :: p_x(3)                   ! Cart. coordinates
         REAL(wp), INTENT(in) :: p_pos(2)     ! geo. coordinates (lon,lat)
         REAL (wp) :: z_clt
         z_clt = COS(p_pos(2))
         p_x(1:3) = (/ COS(p_pos(1))*z_clt, SIN(p_pos(1))*z_clt, SIN(p_pos(2)) /)
  END FUNCTION gc2cc
END PROGRAM interpol_01

