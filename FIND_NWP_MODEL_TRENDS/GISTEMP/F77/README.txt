For more information about the datasets in this directory, please see

    http://data.giss.nasa.gov/gistemp/

This directory contains the basic data files that are used to create
the various maps of temperature anomalies or trends:

1: SBBX.Tsurf1200  (44 MB) surf air temp (1880-present) 1200km smoothing
2: SBBX.Tsurf250   (19 MB) surf air temp (1880-present)  250km smoothing
3a: SBBX.SSTHadR2   (35 MB) sea surf temp. data (1880-present) - used until 11/2012
3b: SBBX.ERSST      (35 MB) sea surf temp. data (1880-present) - used until  6/2015
3c: SBBX.ERSSTv4    (35 MB) sea surf temp. data (1880-present) - used until  7/2017
3d: SBBX.ERSSTv5    (35 MB) sea surf temp. data (1880-present) - currently used

All data are anomalies, the base period depends on the file and may
also change from location to location. Each file starts with a header
record followed by 8000 data records, each containing the full time
series for one of the 8000 equal area boxes (SBBX: subgrid boxes.)

The files are binary, unformatted, sequential access, containing
4-byte reals and integers (IEEE standard, normal byte ordering) **)
They were created by a FORTRAN program on an AIX Unix system, hence
each record starts and ends with a 4-byte integer that represents
the length of the record in bytes (NOT counting those 8 bytes).

The remaining details of the format of these files may be obtained
by looking at or using the program SBBX_to_1x1.f or mkTsMap.f .

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ The new conversion-to-netcdf program                              + 
+ ------------------------------------                              +
+ The following script produced the available netcdf files:         +
+                                                                   +    
+ ln SBBX.ERSSTv5 ERSST_DATA                                        +
+ ln SBBX1880.Ts.GHCN.CL.PA.1200 TS1200_DATA                        +
+ ln SBBX1880.Ts.GHCN.CL.PA.250  TS250_DATA                         +
+ gfortran -fconvert=big-endian SBBX_to_ncint.f -o SBBX_to_nc.exe \ +
+          -L/usr/lib64 -lnetcdf -lnetcdff -I/usr/include           +
+ ./SBBX_to_nc.exe 1200 5 1880 2016 9 gistemp1200_ERSST.nc          +
+ ./SBBX_to_nc.exe  250 0 1880 2016 9 gistemp250.nc                 +
+                                                                   +  
+ The arguments are in order:                                       +
+ smoothing radius            250 or 1200                           +
+ ocean index                 1=Hadley/Reynolds 3=ERSST 4=ERSSTv4   +
+ first year of series        1880 or later             5=ERSSTv5   +
+ last year month of series                                         +
+ name of output file                                               +
+                                                                   + 
+ Any other parameter changes have to be made in SBBX_to_ncint.f    +
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

To use the other programs
=========================
1) rename (ln or mv) the appropriate input files:
   Tsurf file (1  or  2)  -->  TS_DATA
   SST file   (3a,..,3d)  -->  SST_DATA

2) fortran compile the programs (fortran90) with appropriate options
   (command and options are machine/systems dependent), e.g.
      f90 SBBX_to_1x1.f -o SBBX_to_1x1
      f90 mkTsMap.f     -o mkTsMap

3) run the programs typing  "SBBX_to_1x1"  or "mkTsMap < I"

Running SBBX_to_1x1: (obsolete - see next section)
====================
  The program asks for a month, year, ocean_flag
  that flag decides which files are actually used:
    ocean_flag = 0    only  TS_DATA is used
    ocean_flag = 1    only SST_DATA is used
    ocean_flag = 2    both files are used, in that case:
    if land data are present from a station that lies within
         rland=100km of the center of the box, use TS_DATA,
    else use SST_DATA for that box (if present).
  The program finds the appropriate anomaly (the base period
  1951-1980 may be changed by modifying one statement near the
  beginning of the code) and replicates it onto a 1x1 degree grid.
  Finally, an appropriate title and the data are written out
  creating the file dT1x1MAP .

Running mkTsMap < I (more flexible)
===================================
  Create the parameter text file I
  It contains 11 numbers determining:
   1 - task   0=time series of anomalies 1=anomaly 2=mean_change ("trend") ...
           -N<0, same as task=0 but creating only 1 record every N-th year
           (e.g. -1 may be used to create a Jan-series or a NH-Summer-series)
   2 - type   mmLL LL=length_of_period(01=mon,03=seas,12=ann), mm=beg.month
   3 - period:first_year; if LL > 1 use year in which the first period ends
   4 - period:last_year;  if LL > 1 use year in which the last period ends
         if type=1203 and the first period=Dec1950-Feb1951: first_year=1951
                    and if the last period=Dec1960-Feb1961:  last_year=1961
       that same rule applies to the base period below (its type is also mmLL)
   5 - base period:first year  or (if task=2) % of good data needed (e.g. 66)
   6 - base period:last  year  or (if task=2) max.starting year for data, i.e.
       if non-missing data start after that year, the trend is set to "missing"
   7 - rland : if =9999,<0,=100  then ocean_flag=0,1,2 (see:Running SBBX_to_1x1)
   8 - IM = number of output grid boxes in a latitude zone
   9 - JM = number of output grid boxes in the South-North direction
  10 - offI = (Western edge of box 1 minus IDL)/(West-East extent of box)
  11 - dlat = South-North extent of non-polar box in degrees
       if dlat<0, treat polar boxes around a pole as a single box
       Note: dlat=180./JM if all boxes have the same latitudinal width
             else: the width of a pole box is 1/2 of 180-dlat*(JM-2)

  Below are some examples of I_files, mostly creating time series:

  Example 1: Itest (contains 1 line:)
-1 1880 2

  Example 2: I_Ts_mon_anom_series (monthly mean anomalies, base 1951-80)
0  0101  1880 2013  1951 1980  9999.  180 90 0. 2.

  Example 3: I_LOTI_mon_anom_series (monthly mean anomalies, base 1951-80)
0  0101  1880 2007  1951 1980   100.  180 90 0. 2.

  Example 4: I_Ts_ann_anom_series (annual mean anomalies, base 1951-80)
0  0112  1880 2007  1951 1980  9999.  180 90 0. 2.

  Example 5: I_Ts_Jan_anom_series (January anomalies, base 1951-80)
-1  0101  1880 2008  1951 1980  9999.  180 90 0. 2.

  Example 6: I_Ts_Jan2007_anom (single map, Jan 2007 anomaly)
1 0101  2007 2007  1951 1980  9999.  180 90 0. 2.

  Example 7: I_Ts_Jan2001-2010_anom (single map, mean Jan 2001-10 anomaly)
1 0101  2001 2010  1951 1980  9999.  180 90 0. 2.

  Individual maps are easier to create interactively using the web utilities.

  Output file(s):
  ---------------
     out_data      is a binary file (1 or more records)
               If out_data is a sigle record, two additional files are created:
     grid.txt      out_data reformatted as text file
     zonal.txt     table of zonal means (text file)

  Each record of out_data is of the form  "title,array" with
        character*80 title ; real*4 array(1:im,1:jm)

  and may be read using:

        open(1,file='out_data',form='unformatted')

        do
          read (1,end=100) title,array
          write (*,*) titles
       !  bulk of program:
        end do

  100   stop
        end

  The output GRID is rectangular (i:West to East, j:South to North),
  where the pole boxes might be of a different latitudinal size than
  the other boxes (but North and South pole boxes having the same size).
     Examples:
  if offI=0, Western edges of i=1 boxes lie on the international date line;
  if offI=-.5, centers of i=1 boxes lie on the international date line;
  if dlat=180./JM, all boxes have the same latitudinal extent;
  if dlat=180./(JM-1), pole boxes are half boxes.

Problems
========
For questions and remarks, please contact     Reto A. Ruedy
by sending email to                           rruedy@giss.nasa.gov

If you have problems reading the basic binary files, use Itest (Example 1
above) and send me the output of  "SBBX_to_1x1 < Itest" .

**) If you are working with little-endian chips and your compiler
cannot handle the conversion, you may use the attached program
(swap.f) to convert the 3 basic files: rename each file in turn DAT,
apply the program, rename the output file DAT_REV to your liking.

swap.f :

C**** swap.f ; to compile:  f77 swap.f  (or: f90 swap.f)
C****
C**** Converts big_endian to little_endian byte arrangement
C**** for the 3 basic gistemp input files
C****
      parameter (monmx=12*(2200-1700))

      integer info(8)
      real  tin(monmx)
      character*4 infoc(8),tinc(monmx)
      equivalence (infoc(1),info(1)),(tin(1),tinc(1))
      character*80 title

      write(*,*) 'rename (link) your file to DAT'
      write(*,*) 'the converted file will be DAT_REV'

c**** Try assuming rec.lengths are measured in 4-byte words
      open(8,file='DAT',form='unformatted',access='direct',recl=1)
      read(8,rec=1,err=10) len
      go to 20
c**** if the section above leads to a crash, delete it

c**** looks like, rec.lengths are measured in bytes
   10 close (8)
      open(8,file='DAT',form='unformatted',access='direct',recl=4)
      read(8,rec=1) len

   20 if(len.eq.112) stop 'data already big_endian'

c**** just to make sure:
      read(8,rec=30) len1
      if(len.ne.len1) then  ! let's hope we did not yet try recl=4
        close (8)
        open(8,file='DAT',form='unformatted',access='direct',recl=4)
      end if

C****
C**** Read, convert, and save the header record
C****
!!    read(8) info,title   !! but do it word by word
      do n=2,9
         read(8,rec=n) info(n-1)
      end do
      call swap_bytes_4 (infoc,8)
      do n=10,29
        read(8,rec=n) title(1+4*(n-10):4+4*(n-10))
      end do
      nr=30+2

      open(9,file='DAT_REV',form='unformatted',access='sequential')
      write(9) info,title

      Mnow=info(1)       ! length of first data record (months)
C****
C**** Loop over subboxes - read, convert, save
C****
      do n=1,8000
C**** Read in time series TIN/TINO of monthly means: land/Ocean data
        call sread (8,tin(1),tinc(1),Mnow,next,nr)
        Mnow=next  ! Mnow/next: length of current/next time series
      end do
      stop
      end

      SUBROUTINE SREAD (NDISK,A,ac,LEN, LNEXT,nr)
      real A(LEN),dstn
      character*4 ac(len),dstnc,nsc(6)
      integer ns(6)
      equivalence (dstn,dstnc),(ns,nsc)

!!    READ(NDISK) LNEXT, NS, DSTN, A
      read(ndisk,rec=nr) lnext
      do n=1,6
        read(ndisk,rec=nr+n) ns(n)
      end do
      read(ndisk,rec=nr+7) dstn
      nr=nr+8
      do n=1,len
      read(ndisk,rec=nr) a(n)
      nr=nr+1
      end do
      nr=nr+2 ! skip over record end and beg. markers

      call swap_bytes_4 (lnext,1)
      call swap_bytes_4 (nsc,6)
      call swap_bytes_4 (dstnc,1)
      call swap_bytes_4 (ac,len)
      write(9) LNEXT, NS, DSTN, A
      return
      end

      subroutine swap_bytes_4( c, ndim )
      integer n,ndim
      character*1 c(4,ndim),temp
C**** Reverses bytes in array c (integer or real*4) with ndim elements
        do n=1,ndim
          temp   = c(1,n)
          c(1,n) = c(4,n)
          c(4,n) = temp
          temp   = c(2,n)
          c(2,n) = c(3,n)
          c(3,n) = temp
        end do
      return
      end
