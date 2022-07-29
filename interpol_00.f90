! ----------------------------------------------------------------------
! EXAMPLE04
!
! Example program for the DWD ICON tools library
! ----------------------------------------------------------------------
!
! 10/2015 - F. Prill, DWD
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
  CHARACTER (LEN=*), PARAMETER :: grid_filename   = "/input/ICON/ivme/ICON-1E_DOM01.nc"
  CHARACTER (LEN=*), PARAMETER :: in_filename     = "/project/ivme/MCH-1/icon-art-BRM/icon_output/ICON-ART-OEM_1504Unsrt_DOM01_00140000.grb"
  REAL(wp),          PARAMETER :: plam_phi(2)     = (/ 8._wp, 47._wp /) ! search point
  TYPE(t_geographical_coordinates) :: p, cell_center, vertex

  ! local variables:
  integer                          :: ncfileID, status
  TYPE (t_file_metadata)           :: file_metadata
  INTEGER                          :: i, index, min_node_idx(3), jc, jb, jc_n, jb_n
  REAL(wp)                         :: min_dist

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
  

#if !defined(NOMPI)
  CALL start_mpi()
#endif

  ! open/read ICON grid definition 
  ! -----------------------------------------
  status = nf90_open(grid_filename, nf90_nowrite, ncfileID)
  CALL load_grid(gridA, GRID_TYPE_ICON, rank0=0, opt_fileID=ncfileID)
  status = nf90_close(ncfileID)
  CALL init_query_algorithm(gridA, kdtree=.TRUE.)
  
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
  write (*,*) "Distance:", min_dist
  write (*,*) zsfc_in(jc, jb)

  ! print the cell neighbors' coordinates
  DO i=1,3
    jc_n = gridA%p_patch%cells%neighbor_idx(jc,jb,i)
    jb_n = gridA%p_patch%cells%neighbor_blk(jc,jb,i)
    cell_center = gridA%p_patch%cells%center(jc_n,jb_n)
    WRITE (0,'(a,i0,2(a,F6.3),a,F10.2)') " neighbor_center ", i, " : ", &
        cell_center%lon, ", ", cell_center%lat, ", ", zsfc_in(jc_n, jb_n)
  END DO

#if !defined(NOMPI)
  CALL stop_mpi()
#endif

END PROGRAM interpol_01
