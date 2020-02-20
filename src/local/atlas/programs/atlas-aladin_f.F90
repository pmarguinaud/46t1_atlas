
use atlas_module, only : atlas_library, atlas_log, atlas_LambertRegionalGrid, &
                       & atlas_Grid

implicit none

integer, parameter :: ndglg = 32
integer :: nloeng (ndglg) = &
[ 20, 30, 40, 48, 54, 60, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
& 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 60, 54, 48, 40, 30, 20 ]

type (atlas_LambertRegionalGrid) :: grid

integer, parameter :: Nx = 64, Ny = 64;
integer, parameter :: Nux = 53, Ncx = 8, Nuy = 53, Ncy = 8;

real (8), parameter :: LaDInDegrees = 67.36, LoVInDegrees = 26.64;
real (8), parameter :: DxInMetres = 2500, DyInMetres = 2500;
real (8), parameter :: Latin1InDegrees = 67.36, Latin2InDegrees = 67.36;

integer :: ix, iy

call atlas_library%initialise()

grid = atlas_LambertRegionalGrid (Nx, Ny, -Nux / 2 * DxInMetres, -Nuy / 2 * DyInMetres, DxInMetres, DyInMetres, &
                                & LoVInDegrees, LaDInDegrees, Latin1InDegrees, Latin2InDegrees)

print *, " size = ", grid%size ()

do iy = 1, grid%ny ()
  do ix = 1, grid%nx (iy)
    write (*, '(2F20.10)') grid%lonlat (ix, iy)
  enddo
enddo


!call testGrid (grid)

call grid%final ()

call atlas_library%finalise()

contains

subroutine testGrid (grid)

type (atlas_Grid) :: grid

end

end 

