# Parameters of the calorimeters cell types (class CellType)
# ===============================================================
# cell  name       Xsize   Ysize   Zsize   RadLen  NuclLen  DispPhE 
# ---------------------------------------------------------------

  cell  gams_cell  38.40   38.40  450.00   10      10        1       0.5 0.5 0.5 0.5
  cell  pwo_cell   19.20   19.20  250.00   10      10        1       0.5 0.5 0.5 0.5

#         Parameters of the calorimeter GAMS
# ============================================================
# calo  name   Xcenter  Ycenter   Zcenter
# -----------------------------------------------------------
  calo  GAMS    0.00     0.00    1000.00

#        Parameters  GAMS  matrixes
# ============================================================
#        CalorimeterName MatrixName    CellType    Nhoriz    Nvert  Xcenter  Ycenter  Zcenter Nholes
# -----------------------------------------------------------
  calm     GAMS          gams_matrix  gams_cell       32      32     0.000    0.000    0.000    1     14  17  14  17
  calm     GAMS          pwo_matrix   pwo_cell         8       8     0.000    0.000 -100.000    1      3   4   3   4 
# 
#
