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
PROGRAM example04

#if !defined(NOMPI)
  USE mo_mpi,               ONLY: start_mpi, stop_mpi
#endif
  USE mo_utilities,         ONLY: wp, t_geographical_coordinates, pi_180, blk_no, idx_no
  USE mo_metadata_types,    ONLY: t_file_metadata
  USE mo_remap_grid_types,  ONLY: t_grid
  USE mo_kdtree,            ONLY: kdtree_query
  USE mo_icongpi,           ONLY: icongpi_setup, icongpi_finalize

  IMPLICIT NONE

  ! local constants
  CHARACTER (LEN=*), PARAMETER :: grid_filename   = "/input/ICON/ivme/grid/ICON-1E_DOM01.nc"
  REAL(wp),          PARAMETER :: plam_phi(2)     = (/ 8._wp, 47._wp /) ! search point

  ! local variables:
  TYPE (t_grid)                    :: grid
  TYPE (t_file_metadata)           :: grid_file
  TYPE(t_geographical_coordinates) :: p, cell_center, vertex
  INTEGER                          :: index, jc, jb, i, jc_v, jb_v, min_node_idx(3)
  REAL(wp)                         :: min_dist

#if !defined(NOMPI)
  CALL start_mpi()
#endif

  WRITE (0,*) "point ", plam_phi
  CALL icongpi_setup(0, grid_filename, grid, grid_file)

  ! translation  lam/phi    ->  index of nearest cell
  p%lon = plam_phi(1) * pi_180
  p%lat = plam_phi(2) * pi_180
  CALL kdtree_query(grid%kdtree, p, min_dist, min_node_idx) ! global index
  index = min_node_idx(3)
  WRITE (0,*) "lam/pi      ->  nearest cell index = ", index

  ! print the cell center coordinates
  jc = idx_no(index) ! translate global index into pair (index/block)
  jb = blk_no(index)
  cell_center = grid%p_patch%cells%center(jc,jb)
  WRITE (0,'(2(a,F6.3))') " grid index  ->  lam/phi = ", cell_center%lon, ", ", cell_center%lat

  ! print the cell vertex coordinates
  DO i=1,3
    jc_v = grid%p_patch%cells%vertex_idx(jc,jb,i)
    jb_v = grid%p_patch%cells%vertex_blk(jc,jb,i)
    vertex = grid%p_patch%verts%vertex(jc_v,jb_v)
    WRITE (0,'(a,i0,2(a,F6.3))') " vertex ", i, " : ", vertex%lon, ", ", vertex%lat
  END DO

  CALL icongpi_finalize(grid, grid_file)

#if !defined(NOMPI)
  CALL stop_mpi()
#endif

END PROGRAM example04
