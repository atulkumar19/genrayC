
c      the data for the profiles of the absorbed power and current drive.
c      PARAMETER (NRA) is set in param.i
c      It has namelist data only

c-----from namelist /tokamak/
      integer
     &         NR, NRgrid, NZgrid

      common/onetwo_nml/
     &         NR, NRgrid, NZgrid

c---------------------------------------------------------------------------------
c    NR is the maximal number of points in the small radius direction 
c      for the
c      calculations of the power and current drive radial profiles 
c      in each radial bin.
c    NRA= maximal value of NR. It is set in param.i
c    YuP[Nov-2014] Added:
c    NRgrid and NZgrid == (R,Z) grid for 'local' power deposition profiles
c     (without averaging over flux surface, which may not exist)
c---------------------------------------------------------------------------------
