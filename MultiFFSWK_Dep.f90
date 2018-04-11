! rewritie the ffswk.f90 for multi-mode
! use the thought of Deng Kai

! Shenjian Zhang
! July. 2017

! need to change: tstart(gcomr) and fresnel(wn)

! compile: $ ifort -132 -o MultiFFSWK_Dep 



module constant_para
integer, parameter :: Maxsrec = 200    ! max layer for kernels
integer, parameter :: Laydim = 500    ! dimension of model
integer, parameter :: Maxdat=2500    ! dimension of taper(max value)
integer, parameter :: NDISP=1.0E+5    ! dimension of disoersion(?)
integer, parameter :: Maxaniso=6    ! number of kernels(dimension of anisotropic)
integer, parameter :: ModeIndexMinTor = 1    ! min N-order for Toridal
integer, parameter :: ModeIndexMaxTor = 500    ! max N-order for T
integer, parameter :: ModeIndexMinSph = 1    ! min N-order for Spherical
integer, parameter :: ModeIndexMaxSph = 500    ! max N-order for S
integer, parameter :: MaxMode = ModeIndexMaxTor    ! max N-order
complex, parameter :: cc=(0,1)    ! imaginary unit
real, parameter :: pi = 3.14159265
real, parameter :: twopi = 3.14159265*2
real, parameter :: RR0 = 6371000.0
real, parameter :: D2K = 6371.0*twopi/360    ! Degree to Km, Dist(km)=R0*Delta(rad)
real, parameter :: R2D = 360.0/twopi    ! Radian to Degree
real, parameter :: Scale_k = (10**9)*(10**8)    ! convert kernel unit to 1E-8 1/km^3(?)
integer Maxrec    ! number of layers from the depath_max to the surface
real Depth_max    ! maximun depth to record kernels
!============NEW==================
integer target_index    ! index of the target interface in model file
integer tidxm     ! index of the interface, from 1 to maxrec
!=================================
end module constant_para

module Fourier_trans
use constant_para
integer, parameter :: NPTS=2**11    ! sampling points
real, parameter :: DT=1.953125    ! time interval
real, parameter :: DF=1.0/DT/float(NPTS)    ! freq interval
real, parameter :: DOmega=TWOPI*DF
!real, parameter :: LWindow=80.0    ! Window length?
!integer, parameter :: NPOINTS=NINT(LWindow/DT)+1
integer, parameter :: MaxF=NINT(70.0*0.001/DF)    ! 20mHZ, number of freq
end module Fourier_trans

module eigen_functions
use constant_para
use Fourier_trans
real rs, av(laydim), ah(laydim), bv(laydim), bh(laydim)    ! rs: source depth
real rho(laydim),ani(laydim)    
real radius(Laydim)
real(kind=4),dimension(MaxMode,0:MaxF) :: vphaser,vphase0,gcomr,fl,wn,qinv,omega    ! read from eigf.love.0
real(kind=4),dimension(MaxMode,0:MaxF,Maxsrec) :: UU,DU,VV,DV,WW,DW    ! IMPORTANT
real(kind=4),dimension(2) :: UU1,DU1,VV1,DV1,WW1,DW1    ! UU1(1):below, UU1(2):above
real(kind=4),dimension(3) :: UU2,DU2,VV2,DV2,WW2,DW2    ! UU2(1):below, UU2(2):above, UU2(3):maxrec(surface)
integer NLay, LS00, MaxLay
integer is_Love, is_UVW
character(len=70) Dir_eigen
end module eigen_functions

module synthetic_ray
use constant_para
use Fourier_trans
real scola,slong,rcola,rlong,fla,flo,dist_act,distrad
real xmom(6), scale
real freq_main    ! New, target frequency
integer imhz_main    ! New, index of main frequency
real tarray(MaxDat)
real tstart, wintt
real tcut1,tcut2    ! New, for time window
real Taper(0:NPTS-1)    ! New, (cosine) taper function,h_j(t)
integer is_cut    ! New
integer MaxLen,index_start    ! Max length of the window, for tempt taper functon
complex s_omega(0:NPTS-1)    ! New, synthetic response: s(\omega)
complex delta_s_omega(0:NPTS-1)
complex s_tp_omega_cj(0:NPTS-1)    ! s_j(w), used in the last step
real  wseisbot(0:NPTS-1)    ! denominator term of the kernel: s(w)*s_j(w)
real wkn(0:NPTS-1)
real wkn_main    ! result of target frequency
end module synthetic_ray



program ffswk_multimode

use constant_para
use Fourier_trans
use eigen_functions
use synthetic_ray

!implicit none
character(len=70) file_kn
real(kind=8) g2e(3,3), e2g(3,3)    ! rotation matrix
real(kind=8) xn,yn,zn,r,th,ph,ts,ps,tr,pr,rad,colat_r,alon_r,xo,yo,zo
real steplon,cla,clo,rlo,gcarc,az12,az21,dist
real istart_lon, max_lon_steps, isteplon, alon_d,gamm
real target_depth, target_radius
integer fresnel_number,frac
integer iunit_rev,iunit_meca,iunit_kern_dep
integer i
character(len=3) stemp

steplon = 1.0    ! grid of the path
is_cut = 1    ! as an input??????
! read input file
! ===================NEW==========================
write(*,*)'Target Depth =  km'
read(*,*)target_depth
target_depth = target_depth * 1000.0    ! in meter
!===================================================
write(*,*)'Main frequency =  mHz'
read(*,*)freq_main
imhz_main = NINT(freq_main*0.001/DF)

write(*,*)'Love(1) or Rayleigh(2)'
read(*,*)is_Love
if(is_Love .eq. 1) is_UVW = 3
if(is_Love .eq. 2) is_UVW = 1

write(*,*)'eigenfunctions dir (e.g. eigfs)'
read(*,'(A)')dir_eigen

write(*,*)'source_long source_lat depth'
read(*,*)slong,sla,dep
if(slong .lt. 0) slong=slong+360

write(*,*)'moment tensor xmom(6)'
read(*,*)xmom

write(*,*)'receiver_long receiver_lat'
read(*,*)rlong, rla
if(rlong .lt. 0) rlong=rlong+360

write(*,*)'length of measurement window in seconds'
read(*,*)wintt

write(*,*)'Fresnel number (3 or 2) frac(20 40)='
read(*,*)fresnel_number, frac
fresnel_number = fresnel_number*2

write(*,*)'maximum depth(km) to record kernels(800km)'
read(*,*)depth_max
depth_max = depth_max*1000

dir_eigen = adjustl(dir_eigen)    ! delete space at beginning

iunit_rev = 601
iunit_meca = 602

iunit_kern_dep = 701

! open output files
open(iunit_kern_dep,file="kern/multi.kernel.dep")

! set eigenfunctions set (3-D matrix)
call read_eigf()

! ===================NEW====================
target_index = -1
target_radius = RR0-target_depth    ! in meters
do i = 1, maxlay
    if(abs(radius(i) - target_radius) .lt. 1.0E-5) then
        target_index = i    ! the index of the target interface
        exit
    endif
end do
if (target_index .lt. 0) then
    write(*,*)'Bad target depth,abort'
    stop
end if
tidxm = maxrec-(maxlay-target_index)
!==========================================
! convert vphase to km/s
vphase0 = vphaser/1000.0

! scale moment tensor to avoid numerical overflow
scale = 1E-10
do i=1,6
    scale = max(abs(xmom(i)),scale)
end do
do i=1,6
    xmom(i) = xmom(i) / scale
    xmom(i) = xmom(i) * 1.0E20
end do

! convert latitude to co-latitude
scola = 90.0d0-sla
rcola = 90.0d0-rla

! source depth
rs = RR0-dep*1000

! recalculate distance: dist
cla = 90-scola
clo = slong
rla = 90-rcola
rlo = rlong
call distaz(cla,clo,rla,rlo,gcarc,az12,az21)
write(*,*)'recalculated epicentral distance = ',gcarc
dist_act = gcarc    ! gcarc: dist. in degree

open(iunit_rev,file='kern/example.sr')
write(iunit_rev,*)slong,(90-scola), '15 0  4  5 S'
write(iunit_rev,*)rlong,(90-rcola), '15 0  4  5 R'
close(iunit_rev)

open(iunit_meca,file='kern/meca.example')
write(iunit_meca,2)slong,(90-scola), dep, xmom, slong, (90-scola)
2 format(3(F12.2),6(E12.4)," 1 ",2(F12.2))
close(iunit_meca)

write(*,*)'GMT output files: kern/example.sr'

! convert to radian
ts = scola*pi/180
tr = rcola*pi/180
ps = slong*pi/180
pr = rlong*pi/180

! rotation matrix: g2e, e2g
call euler(ts,ps,tr,pr,g2e,e2g)

! taper and time window
dist = dist_act*pi/180
!tstart = dist_act*D2K*1000/gcomr(3,imhz_main)-0.5*wintt    ! gcomr(3,imhz_main): choose N=2 for Group Vel.
tstart = 605-0.5*wintt
tcut1 = tstart
tcut2 = tstart+wintt    ! time window
MaxLen = nint(wintt/dt)
write(*,*)'measurement window'
write(*,*)tcut1,'s  ',tcut2,'s '
write(*,*)
! calculate the taper-function h_j(t), consine-type
call calc_taper()    ! get: Taper(0:NPTS-1)
! ==========IMPORTANT===========
! calculate synthetic response
call calc_synthetic_ray()
! ==========IMPORTANT===========
! calculate scattered reponse
! 1. along-ray direction
istart_lon = 0    ! start lon in Euler sys
max_lon_steps = (gcarc/steplon)+0.5
! =========along ray loop========
do  isteplon = istart_lon, max_lon_steps
    alon_d = (isteplon-0.5)*steplon
    alon_r = alon_d*pi/180
! 2. cross ray loop
    gamm = abs(sin(gcarc*pi/180)/(sin(alon_r)*sin(gcarc*pi/180-alon_r)))
    fresnel = sqrt(2*pi/wn(2,imhz_main)/gamm)*180/pi    ! width of fresnel zone in degree
    ! frac: number of grid points per fresnel zone
    steplat = fresnel/frac    ! frac is from input file
    isteplat_max = frac*fresnel_number
    half = fresnel_number/2
    ! ================ cross-ray loop =================
    do  isteplat = 1, isteplat_max
        ! from -half of zone to +half of zone
        alat_d = (-1)*half*fresnel + steplat*float(isteplat-1)
        colat_r = pi/2-alat_d*pi/180
        rad = 1    ! radius for rotate

        call rtp2xyz(rad,colat_r,alon_r,xo,yo,zo)    ! r-theta-phi to x-y-z
        call rotate(xo,yo,zo,e2g,xn,yn,zn)    ! use rotate matrix e2g to rotate
        call xyz2rtp(xn,yn,zn,r,th,ph)    ! x-y-z to r-theta-phi

        t2 = pi/2-th
        p2 = ph
        ! flo,fla: geographic longitude and latitude of scatter point, for
        ! output
        fla = t2*180/pi
        flo = p2*180/pi
        if (flo .lt. 0)  flo = flo+360
        ! ========= IMPORTANT ===========
        call scatter(alon_d,alat_d)
        ! output kernels
        dep_k = target_depth/1000.0    ! depth in km
        write(iunit_kern_dep,1)flo,fla,wkn_main*scale_k,dep_k    ! alpha_v
1 format(F10.3,F10.3,E12.4,F10.2)
    end do     ! end of cross-ray loop
end do    ! end of along-ray loop

print *, 'kernel unit is 10^-8 km^3/s'

close(iunit_kern_dep)

stop


end program ffswk_multimode
