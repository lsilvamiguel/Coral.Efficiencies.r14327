C commonblock containing information of hit cells
C ncell   number of hit cells in current event
C e       energy of a cell
C iaddr   index of cell in array of calibrations, ...
C x       x position of the cell (in cm)
C y       y position of the cell (in cm)
C err     error on energy measurement
C time    time of signal in cell

      common/cells/ ncell,
     $              e(numc),
     $              iaddr(numc),
     $              x(numc),
     $              y(numc),
     $              err(numc),
     $              time(numc)

      integer ncell,
     $        iaddr

      real e,
     $     x,
     $     y,
     $     err,
     $     time
