
use atlas_module, only : atlas_library, atlas_log, atlas_ReducedGaussianGrid, &
                       & atlas_Grid

implicit none

integer, parameter :: ndglg = 32
integer :: nloeng (ndglg) = &
[ 20, 30, 40, 48, 54, 60, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
& 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 60, 54, 48, 40, 30, 20 ]

type (atlas_ReducedGaussianGrid) :: grid

real (8), parameter :: latitudeOfStretchingPoleInDegrees = 46.4688;
real (8), parameter :: longitudeOfStretchingPoleInDegrees = 2.57831;
real (8), parameter :: stretchingFactor = 2.4;

integer :: ix, iy

call atlas_library%initialise()


grid = atlas_ReducedGaussianGrid (nloeng, [longitudeOfStretchingPoleInDegrees, &
                                & latitudeOfStretchingPoleInDegrees], stretchingFactor)

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

