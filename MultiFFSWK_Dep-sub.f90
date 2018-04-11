! subroutines of MultiFFSWK_Dep.f90
! Use the thought of Deng Kai

! Shenjian Zhang
! July 2017



subroutine read_eigf()

use constant_para
use Fourier_trans
use eigen_functions

integer nlayr(NDISP), ls00r(NDISP)
character(len=70) file_layer, file_eigf, file_model    ! output files of pre-do
integer imHzIndx    ! imHz of last frquency

WW(:,:,:) = 0
DW(:,:,:) = 0
VV(:,:,:) = 0
DV(:,:,:) = 0
UU(:,:,:) = 0
DU(:,:,:) = 0
vphaser(:,:) = 0
gcomr(:,:) = 0
fl(:,:) = 0
wn(:,:) = 0
omega(:,:) = 0
qinv(:,:) = 0

! read model file to get maxrec and model para.
file_model = dir_eigen(1:len_trim(dir_eigen))//'/model'

! maxrec is the number of layers from depth_max to surface
open(10,file=file_model)
read(10,*)    ! first line is no use
read(10,*)ifanis,tref,ifdeck
read(10,*)maxlay,nic,noc,rx
maxrec = 0
do i = 1,maxlay
    read(10,*)ind,radius(i),rho(i),av(i),bv(i),Qk,Qmu,ah(i),bh(i),ani(i)
    if(radius(i) .ge. (6371000-depth_max)) maxrec = maxrec+1
end do
close(10)

if(is_Love .eq. 1) then
    file_layer = dir_eigen(1:len_trim(dir_eigen))//'/layer.love.0'
    file_eigf = dir_eigen(1:len_trim(dir_eigen))//'/eigf.love.0'
end if
if(is_Love .eq. 2) then
    file_layer = dir_eigen(1:len_trim(dir_eigen))//'/layer.rayleigh.0'
    file_eigf = dir_eigen(1:len_trim(dir_eigen))//'/eigf.rayleigh.0'
end if

! read layer file to get maxdisp(the number of all modes)
open(10,file=file_layer)
j = 1
read(10,*)
    40 read(10,91,end=41)w4,acr,gcom,nlayr(j),ls00r(j),maxl
    91 format(1X,E11.4,1X,F12.2,1X,F12.2,1X,I5,I5,I5)
    j = j+1
    if(j .eq. NDISP) then
        write(*,*)'dimension overflow, abort!'
        stop
    end if
goto 40
41 close(10)
maxdisp = j-1
! read eigf file to get eigenfunction set
open(10,file=file_eigf,form='unformatted')
! out cycle: mode, inner cycle: layer
imHzIndx = -1    ! record of imHz of last l-order
wmhzIndx = 100    ! the smallest distance from wmhz to imhz*df
do k = 1,maxdisp
    read(10)nord,l,vc,wmhz,tcom,gu,qmod,wdiff
    imHz = NINT(wmhz*0.001/df)
! ===================================
    nord = nord+1    ! nord=1 means N-order is zero!! Only consider Tor
! ===================================
!    if(imHz .eq. imHzIndx) then    ! ensure we record the nearest l-order of each frequency
!        if(abs(wmhz*0.001-imHz*df) .ge. abs(wmhzIndx*0.001-imHz*df)) then
! if the freq. of new mode is father than the record one, then abandon it.
!            do kk = 1,nlayr(k)
!                read(10)
!            end do
!            cycle    ! imHz of new l-order is the same as last one
!        else
!            wmhzIndx = wmhz    ! the nearest freq.
!        end if
!    else
!        imHzIndx = imHz
!        wmhzIndx = 100    ! update
!    end if
    vphaser(nord,imHz) = vc
    gcomr(nord,imHz) = gu
    fl(nord,imHz) = float(l)
    w4 = TWOPI*wmhz*0.001
    wn(nord,imHz) = w4*RR0/vphaser(nord,imHz)
    omega(nord,imHz) = w4
    qinv(nord,imHz) = qmod
    nlay = nlayr(k)
    ls00 = ls00r(k)
    do i = 1,nlay
        il = i-(maxlay-ls00+1)+maxrec
        if(il .gt. 0) then
            if(is_Love .eq. 1) then
                read(10)WW0,DW0
                WW(nord,imHz,il)=WW0*w4*RR0*sqrt( fl(nord,imHz)*(fl(nord,imHz)+1)/(vc*gu) )
                DW(nord,imHz,il)=DW0*w4*RR0*sqrt( fl(nord,imHz)*(fl(nord,imHz)+1)/(vc*gu) )
            end if
            if(is_Love .eq. 2) then
                read(10)UU0,DU0,VV0,DV0
                VV(nord,imHz,il)=VV0*w4*RR0*sqrt( fl(nord,imHz)*(fl(nord,imHz)+1)/(vc*gu) )
                DV(nord,imHz,il)=DV0*w4*RR0*sqrt( fl(nord,imHz)*(fl(nord,imHz)+1)/(vc*gu) )
                UU(nord,imHz,il)=UU0*w4*RR0*sqrt( 1.0/(vc*gu) )
                DU(nord,imHz,il)=DU0*w4*RR0*sqrt( 1.0/(vc*gu) )
            end if
        else
            read(10)    ! layer does not match
        end if
    end do
end do
close(10)
end



subroutine distaz(lat1,lon1,lat2,lon2,distdeg,az12deg,az21deg)
use constant_para

implicit double precision (a-h)
real lat1,lon1,lat2,lon2,distdeg,az12deg,az21deg

slat1 = sin(lat1*pi/180)
clat1 = cos(lat1*pi/180)
slon1 = sin(lon1*pi/180)
clon1 = cos(lon1*pi/180)

slat2 = sin(lat2*pi/180)
clat2 = cos(lat2*pi/180)
slon2 = sin(lon2*pi/180)
clon2 = cos(lon2*pi/180)

a1 = clat1 * clon1
b1 = clat1 * slon1
g1 = slat1 * clon1
h1 = slat1 * slon1

a2 = clat2 * clon2
b2 = clat2 * slon2
g2 = slat2 * clon2
h2 = slat2 * slon2

c1 = a1 * a2 + b1 * b2 + slat1 * slat2    ! c1=vector(s) point vector(r)
if(abs(c1) .lt. 0.94) then
    distrad = acos(c1)
else
    if(c1 .ge. 0) then
        distrad = asin(sqrt((a1 - a2) * (a1 - a2) + (b1 - b2) * (b1 - b2) + (slat1 - slat2) * (slat1 - slat2)) / 2.) * 2.
    else
        distrad = acos(sqrt((a1 + a2) * (a1 + a2) + (b1 + b2) * (b1 + b2) + (slat1 + slat2) * (slat1 + slat2)) / 2.) * 2.
    end if
end if

distdeg = distrad * R2D    ! distant in degree

c3 = (a2 - slon1) * (a2 - slon1) + (b2 + clon1) * (b2 + clon1) + slat2 * slat2 - 2.
c4 = (a2 - g1) * (a2 - g1) + (b2 - h1) * (b2 - h1) + (slat2 + clat1) * (slat2 + clat1) - 2.
az12rad = atan2(c3, c4)
az12deg = az12rad * R2D    ! t


c5 = (a1 - slon2) * (a1 - slon2) + (b1 + clon2) * (b1 + clon2) + slat1 * slat1 - 2.
c6 = (a1 - g2) * (a1 - g2) + (b1 - h2) * (b1 - h2) + (slat1 + clat2) * (slat1 + clat2) - 2.
az21rad = atan2(c5, c6)
az21deg = az21rad * R2D

if(az12deg.lt. 0)  az12deg = az12deg + 360
if(az21deg.lt. 0)  az21deg = az21deg + 360

! ccw from south -180 -- 180     
az12deg = 180 - az12deg
az21deg = 180 - az21deg

if(az12deg.lt. 0)  az12deg = az12deg + 360
if(az21deg.lt. 0)  az21deg = az21deg + 360

return
end



subroutine euler(ts,ps,tr,pr,g2e,e2g)

implicit double precision (a-h,o-z)
dimension g2e(3,3),e2g(3,3)    ! define matrix

pi2 = asin(1.0d0)    ! pi/2
pi = 2.0d0*pi2
twopi = 2.0d0*pi
rpd = pi/180.0d0
dpr = 1.0d0/rpd

call cart(ts,ps,xi,yi,zi)    ! spherical to Cartesian
call cart(tr,pr,xr,yr,zr)
! pole(xk,yk,zk)=(xi,yi,zi)x(xr,yr,zr), axis z'
xk = yi*zr-yr*zi
yk = xr*zi-xi*zr
zk = xi*yr-xr*yi
! call cross(xi,yi,zi,xr,yr,zr,xk,yk,zk)    ! use subroutine to calculate cross
! product
dk = dsqrt(xk*xk+yk*yk+zk*zk)    ! normalization
xk = xk/dk
yk = yk/dk
zk = zk/dk
! (xj,yj,zj)=(xk,yk,zk)*(xi,yi,zi)    ! axis y'
xj = yk*zi-yi*zk
yj = xi*zk-xk*zi
zj = xk*yi-xi*yk
dj = dsqrt(xj*xj+yj*yj+zj*zj)
xj = xj/dj
yj = yj/dj
zj = zj/dj

! rotation matrix cos(x',x)
g2e(1,1) = xi
g2e(2,1) = xj
g2e(3,1) = xk
g2e(1,2) = yi
g2e(2,2) = yj
g2e(3,2) = yk
g2e(1,3) = zi
g2e(2,3) = zj
g2e(3,3) = zk

e2g(1,1) = xi
e2g(2,1) = yi
e2g(3,1) = zi
e2g(1,2) = xj
e2g(2,2) = yj
e2g(3,2) = zj
e2g(1,3) = xk
e2g(2,3) = yk
e2g(3,3) = zk

return
end



subroutine cart(thet,phi,x,y,z)
implicit real*8(a-h,o-z)

s = dsin(thet)
x = s*cos(phi)
y = s*sin(phi)
z = dcos(thet)

return
end



subroutine cross(sx,sy,sz,rx,ry,rz,px,py,pz)
implicit real*8(a-h,o-z)

px = sy*rz-sz*ry
py = sz*rx-sx*rz
pz = sx*ry-sy*rx

return
end



subroutine calc_taper()
use synthetic_ray

implicit none
integer k
real TmpTaper(MaxLen)    ! MaxLen: the index length of window,global

Taper = 0
index_start = floor(tstart/dt)    ! both are defined in module synthetic_ray
do k=0,MaxLen
    TmpTaper(k) = 1-cos(TWOPI*k/MaxLen)
    TmpTaper(k) = TmpTaper(k)/2.0
    Taper(k+index_start) = TmpTaper(k)
end do

end



! =============IMPORTANT==================
subroutine calc_synthetic_ray()
use constant_para
use Fourier_trans
use eigen_functions
use synthetic_ray

implicit none
complex radiat0,exptmp,seis0,wseis0
real omega_tmp,wn_tmp
real cla,clo,rla,rlo,distdeg,az12,az21
real tarrival
integer i,imhz,nn,maslv
complex s_t(0:NPTS-1)
complex s_tapered_omega(0:NPTS-1)
real s_tapered_t(0:NPTS-1)

s_omega(:) = cmplx(0,0)    ! synthetic response s(\omege)
omega_tmp = 0
distrad = dist_act * pi / 180
cla = 90-scola
clo = slong
rla = 90-rcola
rlo = rlong
call distaz(cla,clo,rla,rlo,distdeg,az12,az21)  
! out cycle: frequency, inner cycle: N-order
do imhz = 0,MaxF
    omega_tmp = TWOPI*DF*imhz    ! \omega_n, discrete frequency
    radiat0 = cmplx(0.0,0.0)
    exptmp = cmplx(0.0,0.0)
    seis0 = cmplx(0.0,0.0)
    wseis0 = cmplx(0.0,0.0)
    wn_tmp = 0
   ! s_omega = cmplx(0.0,0.0)
    do nn = 1,MaxMode
        wn_tmp = wn(nn,imhz)
        if(wn_tmp .lt. 1.0E-5) cycle    ! for nn-order, no l-order with freq(imhz)
        if(gcomr(nn,imhz) .le. 0.0) cycle ! ?????
!        tarrival = RR0*distrad / gcomr(nn,imhz)    ! arrival time of this mode
        tarrival = 605.0
        if(is_cut .eq. 1 .and. (tarrival .lt. tcut1 .or. tarrival .gt. tcut2)) cycle
        ! ensure the arrival of this mode is in the time window
        ! ==== source term ====
        call radiation(nn,imhz,az12,radiat0)
        ! ==== phase, attenuation and spread term ====
        maslv = 0    ! pole passage
        exptmp = exp((-1)*cc*(wn_tmp*distrad+pi/4-pi/2*maslv))
        exptmp = exptmp*exp(-0.5*distrad*RR0*omega(nn,imhz)*qinv(nn,imhz)/gcomr(nn,imhz))
        seis0 = exptmp/sqrt(8*pi*wn(nn,imhz)*abs(sin(distrad)))
        ! ==== receiver term ====
        if(is_UVW .eq. 1) seis0 = seis0*UU(nn,imhz,maxrec)
        if(is_UVW .eq. 2) seis0 = -cc*seis0*VV(nn,imhz,maxrec)
        if(is_UVW .eq. 3) seis0 = cc*seis0*WW(nn,imhz,maxrec)
        wseis0 = radiat0*seis0
        s_omega(imhz) = s_omega(imhz)+wseis0    ! attention to imhz
    end do
end do
! ensure the s(t) is a real time series.
do i = 1,MaxF
    s_omega(NPTS-i) = conjg(s_omega(i))
end do
! implement taper function in the frequency domain
! s_tapered_omega = s_omega * h_j(w), * means convolution
call convlt(s_omega,s_tapered_omega)
! s_j(w) \time conjg(s_j(w))
do i = 0,NPTS-1
    s_tp_omega_cj(i) = conjg(s_tapered_omega(i))    ! s^*_j(w)
    wseisbot(i) = s_tapered_omega(i)*s_tp_omega_cj(i)
end do

return
end



subroutine radiation(nn,imhz,ang,radiat)
use eigen_functions
use synthetic_ray

implicit none
integer i,is1,is2,i_1,i_2,M,nn,imhz
complex e1,e2,pt1,pt2,pt3,pt4,pt5,radiat     ! pt = part
real source_factor,theta,ang
real Mrr,Mtt,Mpp,Mrp,Mrt,Mtp
real wn_tmp,omega_r
real WWS,DWS,UUS,DUS,VVS,DVS
e1 = cmplx(0,1)
e2 = cmplx(1,0)
Mrr = xmom(1)
Mtt = xmom(2)
Mpp = xmom(3)
Mrt = xmom(4)
Mrp = xmom(5)
Mtp = xmom(6)
wn_tmp = wn(nn,imhz)
omega_r = omega(nn,imhz)

! find source location
do i = 1,maxlay-1
    if(rs .gt. radius(i) .and. rs .le. radius(i+1)) then
        is1 = i
        is2 = i+1
        source_factor = (rs-radius(is1))/(radius(is2)-radius(is1))
    end if
end do

i_1 = is1-maxlay+maxrec
i_2 = is2-maxlay+maxrec
WWS = WW(nn,imhz,i_1) + source_factor*(WW(nn,imhz,i_2)-WW(nn,imhz,i_1))
DWS = DW(nn,imhz,i_1) + source_factor*(DW(nn,imhz,i_2)-DW(nn,imhz,i_1))
UUS = UU(nn,imhz,i_1) + source_factor*(UU(nn,imhz,i_2)-UU(nn,imhz,i_1))
DUS = DU(nn,imhz,i_1) + source_factor*(DU(nn,imhz,i_2)-DU(nn,imhz,i_1))
VVS = VV(nn,imhz,i_1) + source_factor*(VV(nn,imhz,i_2)-VV(nn,imhz,i_1))
DVS = DV(nn,imhz,i_1) + source_factor*(DV(nn,imhz,i_2)-DV(nn,imhz,i_1))

! radiation term
! compute radiation term
M = 1    ! M is the index n in fomula in Zhou, 2004
theta = ang * pi/180    ! az(from south,ccw) in rad, source take-off angle of the reference ray
pt1 = (Mrr*DUS+(Mtt+Mpp)*(1/rs)*(UUS-0.5*wn_tmp*VVS)) * e1
pt2 = (-1)**M*(DVS-(1/rs)*VVS+wn_tmp*(1/rs)*UUS)*(Mrp*sin(theta)+Mrt*cos(theta)) * e2
pt3 = (-1)*wn_tmp*(1/rs)*VVS*(Mtp*sin(2*theta)+0.5*(Mtt-Mpp)*cos(2*theta)) * e1
pt4 = (-1)**M*(DWS-(1/rs)*WWS)*(Mrt*sin(theta)-Mrp*cos(theta)) * e2
pt5 = (-1)*wn_tmp*(1/rs)*WWS*(0.5*(Mtt-Mpp)*sin(2*theta)-Mtp*cos(2*theta)) * e1
radiat = (-1)*(pt1+pt2+pt3+pt4+pt5)/omega_r
radiat = radiat

return
end


subroutine rtp2xyz(r,t,p,x,y,z)
implicit double precision (a-h,o-z)

pi2 = asin(1.0d0)
if (t .ne. pi2) then
    st = dsin(t)
    ct = dcos(t)
else
    st = 1.0d0    ! to avoid cos(theta) .ne. 0
    ct = 0.0d0
end if
sp = dsin(p)
cp = dcos(p)
x = r*st*cp
y = r*st*sp
z = r*ct

return
end



subroutine xyz2rtp(x,y,z,r,t,p)
implicit double precision(a-h,o-z)

pi2 = asin(1.0d0)
r = dsqrt(x*x+y*y+z*z)
t = dacos(z/r)
if (x .eq. 0.0d0) then
    p = pi2
    return
end if
p = datan2(y,x)

return
end



subroutine rotate(xo,yo,zo,rotm,xn,yn,zn)
implicit real*8 (a-h,o-z)
dimension rotm(3,3)

xn = xo*rotm(1,1)+yo*rotm(1,2)+zo*rotm(1,3)
yn = xo*rotm(2,1)+yo*rotm(2,2)+zo*rotm(2,3)
zn = xo*rotm(3,1)+yo*rotm(3,2)+zo*rotm(3,3)

return
end



subroutine scatter(blo,bla)
use eigen_functions
use synthetic_ray
complex radiatS,exptmp1,exptmp2,seis
real wn1,wn2,omega_tmp
real tarrival
real akern    ! NEW, for depth
complex K_omega(0:NPTS-1),K_tapered_omega(0:NPTS-1)
complex wkntop(0:NPTS-1)

! initialization
wkn(:) = 0
wkn_main = 0
delta_s_omega(:) = cmplx(0,0)

do imhz = 0,MaxF
    omega_tmp = TWOPI*imhz*DF
! mode from source to scatter    
    do nn1 = 1,MaxMode
        UU1(:) = 0
        DU1(:) = 0
        VV1(:) = 0
        DV1(:) = 0
        WW1(:) = 0
        DW1(:) = 0
        wn1 = wn(nn1,imhz)
        if(wn1 .lt. 1.0E-5) cycle
        if(gcomr(nn1,imhz) .le. 0.0) cycle
        ! source to scatter
        cla = 90-scola
        clo = slong
        rla = fla
        rlo = flo
        call distaz(cla,clo,rla,rlo,distdeg1,az12,az21)
        call dis1(0,blo,distdeg1,dist1)    ! ???
        ! source term
        call radiation(nn1,imhz,az12,radiatS)
        ! first phase term
        exptmp1 = exp( (-1)*cc*(wn1*dist1+pi/4) )    ! ignore polar passage
        exptmp1 = exptmp1 * exp( -0.5*dist1*RR0*omega(nn1,imhz)*qinv(nn1,imhz)/gcomr(nn1,imhz) )
        exptmp1 = exptmp1 / sqrt(8*pi*wn1*abs(sin(dist1)))
        WW1(1:2) = WW(nn1,imhz,tidxm:tidxm+1)    ! only consider the target interface
        DW1(1:2) = DW(nn1,imhz,tidxm:tidxm+1)     ! tidxm: the index from 1 to maxrec
        UU1(1:2) = UU(nn1,imhz,tidxm:tidxm+1)
        DU1(1:2) = DU(nn1,imhz,tidxm:tidxm+1)
        VV1(1:2) = VV(nn1,imhz,tidxm:tidxm+1)
        DV1(1:2) = DV(nn1,imhz,tidxm:tidxm+1)

! mode from scatter to receiver   
        do nn2 = 1,MaxMode
            UU2(:) = 0
            DU2(:) = 0
            VV2(:) = 0
            DV2(:) = 0
            WW2(:) = 0
            DW2(:) = 0
            wn2 = wn(nn2,imhz)
            if(wn2 .lt. 1.0E-5) cycle
            if(gcomr(nn2,imhz) .le. 0.0) cycle
            ! scatter to receiver
            cla = fla
            clo = flo
            rla = 90-rcola
            rlo = rlong
            call distaz(cla,clo,rla,rlo,distdeg2,bz12,bz21)
            call dis2(0,blo,distdeg2,dist2)    ! delta" in rad
            call angcal(blo,bla,sphi,ang2)    ! scattering angle and receiver arrival angle  
            sphi = sphi*pi/180
            ang2 = ang2*pi/180
            tarrival = RR0*dist1/gcomr(nn1,imhz)+RR0*dist2/gcomr(nn2,imhz)
            if (iscut .eq. 1 .and. (tarrival .lt. tcut1 .or. tarrival .gt. tcut2)) cycle
            ! second phase term
            exptmp2 = exp( (-1)*cc*(wn2*dist2+pi/4) )    ! ignore polar passage
            exptmp2 = exptmp2 * exp( -0.5*dist2*RR0*omega(nn2,imhz)*qinv(nn2,imhz)/gcomr(nn2,imhz) )
            exptmp2 = exptmp2 / sqrt(8*pi*wn2*abs(sin(dist2)))
            WW2(1:2) = WW(nn2,imhz,tidxm:tidxm+1)
            DW2(1:2) = DW(nn2,imhz,tidxm:tidxm+1)
            UU2(1:2) = UU(nn2,imhz,tidxm:tidxm+1)
            DU2(1:2) = DU(nn2,imhz,tidxm:tidxm+1)
            VV2(1:2) = VV(nn2,imhz,tidxm:tidxm+1)
            DV2(1:2) = DV(nn2,imhz,tidxm:tidxm+1)
            WW2(3) = WW(nn2,imhz,Maxrec)
            DW2(3) = DW(nn2,imhz,Maxrec)
            UU2(3) = UU(nn2,imhz,Maxrec)
            DU2(3) = DU(nn2,imhz,Maxrec)
            VV2(3) = VV(nn2,imhz,Maxrec)
            DV2(3) = DV(nn2,imhz,Maxrec)

            ! seis = source term \time phase term           
            if(is_UVW .eq. 1) seis = radiatS*exptmp1*exptmp2
            if(is_UVW .eq. 2) seis = -cc*radiatS*exptmp1*exptmp2
            if(is_UVW .eq. 3) seis = cc*radiatS*exptmp1*exptmp2
            call layer_coeff(omega_tmp,wn1,wn2,sphi,ang2,akern)    ! akern: the interaction term
            
            delta_s_omega(imhz) = delta_s_omega(imhz)+akern*seis
           
        end do    ! N-order from scatter to receiver
    end do    ! N-order from source to scatter
end do    ! frequency
! ========================New at 2017,1,17=============
! K_j(w) = K(w) convolute H_j(w)
K_omega = cmplx(0.0,0.0)    ! K(w) series
K_tapered_omega = cmplx(0.0,0.0)    !K_j(w) series
wkntop = cmplx(0.0,0.0)
do i = 0,MaxF
    K_omega(i) = delta_s_omega(i) 
end do
! ensure real time series
do i = 1,MaxF
    K_omega(NPTS-i) = conjg(K_omega(i))
end do
! K_tapered_omega = K_omega * H_j(w)
call convlt(K_omega,K_tapered_omega)
do i = 0,NPTS-1
    wkntop(i) = K_tapered_omega(i)*s_tp_omega_cj(i)    ! K_j(w) \time s^*_j(w)
    wkn(i) = (-1.0d0) * aimag(wkntop(i)/wseisbot(i))    ! phase kernel
!   wkn(i,ly,ianiso) = 1.0d0*real(wkntop(i)/wseisbot(i))   ! Amp kernel
end do
wkn_main = wkn(imhz_main)    ! target result

return
end



subroutine layer_coeff(omega_tmp,wn1,wn2,sphi,ang2,akern)
use eigen_functions
use synthetic_ray

real sigma1_a, sigma1_b, sigma2_a, sigma2_b, sigma_total
integer tidx
real RD

tidx = target_index   ! absolute index
RD = radius(tidx)
! Love-Love
if(is_love .eq. 1) then
    sigma1_b = -rho(tidx)*omega_tmp*omega_tmp*WW1(1)*WW2(1)*cos(sphi)
    sigma1_b = sigma1_b+rho(tidx)*bh(tidx)*bh(tidx)*wn1*wn2*WW1(1)*WW2(1)*cos(2*sphi)/RD/RD
    sigma1_b = sigma1_b+rho(tidx)*bv(tidx)*bv(tidx)*(DW1(1)-WW1(1)/RD)*(DW2(1)-WW2(1)/RD)*cos(sphi)

    sigma1_a = -rho(tidx+1)*omega_tmp*omega_tmp*WW1(2)*WW2(2)*cos(sphi)
    sigma1_a = sigma1_a+rho(tidx+1)*bh(tidx+1)*bh(tidx+1)*wn1*wn2*WW1(2)*WW2(2)*cos(2*sphi)/RD/RD
    sigma1_a = sigma1_a+rho(tidx+1)*bv(tidx+1)*bv(tidx+1)*(DW1(2)-WW1(2)/RD)*(DW2(2)-WW2(2)/RD)*cos(sphi)
    
    sigma2_b = DW2(1)*(DW1(1)-WW1(1)/RD)*cos(sphi)
    sigma2_b = sigma2_b+DW1(1)*(DW2(1)-WW2(1)/RD)*cos(sphi)
    sigma2_b = -sigma2_b*rho(tidx)*bv(tidx)*bv(tidx)

    sigma2_a = DW2(2)*(DW1(2)-WW1(2)/RD)*cos(sphi)
    sigma2_a = sigma2_a+DW1(2)*(DW2(2)-WW2(2)/RD)*cos(sphi)
    sigma2_a = -sigma2_a*rho(tidx+1)*bv(tidx+1)*bv(tidx+1)

    sigma_total = (sigma1_a+sigma2_a)-(sigma1_b+sigma2_b)
end if
! Rayleigh-Rayleigh 
if(is_love .eq. 2) then
    sigma1_b = -rho(tidx)*omega_tmp*omega_tmp*(UU1(1)*UU2(1)+VV1(1)*VV2(1)*cos(sphi))
    sigma1_b = sigma1_b+rho(tidx)*av(tidx)*av(tidx)*DU1(1)*DU2(1)
    sigma1_b = sigma1_b+rho(tidx)*ah(tidx)*ah(tidx)*(2*UU1(1)-wn1*VV1(1))*(2*UU2(1)-wn2*VV2(1))/RD/RD
    sigma1_b = sigma1_b+rho(tidx)*ani(tidx)*(ah(tidx)*ah(tidx)-2*bv(tidx)*bv(tidx))&
&   *(DU2(1)*(2*UU1(1)-wn1*VV1(1))+DU1(1)*(2*UU2(1)-wn2*VV2(1)))/RD
    sigma1_b = sigma1_b+rho(tidx)*bh(tidx)*bh(tidx)*&
&   (wn1*wn2*VV1(1)*VV2(1)*cos(2*sphi)-(2*UU1(1)-wn1*VV1(1))*(2*UU2(1)-wn2*VV2(1)))/RD/RD
    sigma1_b = sigma1_b+rho(tidx)*bv(tidx)*bv(tidx)*(DV1(1)-VV1(1)/RD+wn1*UU1(1)/RD)&
&   *(DV2(1)-VV2(1)/RD+wn2*UU2(1)/RD)*cos(sphi)

    sigma1_a = -rho(tidx+1)*omega_tmp*omega_tmp*(UU1(2)*UU2(2)+VV1(2)*VV2(2)*cos(sphi))
    sigma1_a = sigma1_a+rho(tidx+1)*av(tidx+1)*av(tidx+1)*DU1(2)*DU2(2)
    sigma1_a = sigma1_a+rho(tidx+1)*ah(tidx+1)*ah(tidx+1)*(2*UU1(2)-wn1*VV1(2))*(2*UU2(2)-wn2*VV2(2))/RD/RD
    sigma1_a = sigma1_a+rho(tidx+1)*ani(tidx+1)*(ah(tidx+1)*ah(tidx+1)-2*bv(tidx+1)*bv(tidx+1))&
&   *(DU2(2)*(2*UU1(2)-wn1*VV1(2))+DU1(2)*(2*UU2(2)-wn2*VV2(2)))/RD
    sigma1_a = sigma1_a+rho(tidx+1)*bh(tidx+1)*bh(tidx+1)*&
&   (wn1*wn2*VV1(2)*VV2(2)*cos(2*sphi)-(2*UU1(2)-wn1*VV1(2))*(2*UU2(2)-wn2*VV2(2)))/RD/RD
    sigma1_a = sigma1_a+rho(tidx+1)*bv(tidx+1)*bv(tidx+1)*(DV1(2)-VV1(2)/RD+wn1*UU1(2)/RD)&
&   *(DV2(2)-VV2(2)/RD+wn2*UU2(2)/RD)*cos(sphi)

    sigma2_b = -rho(tidx)*av(tidx)*av(tidx)*2*DU1(1)*DU2(1)
    sigma2_b = sigma2_b-ani(tidx)*rho(tidx)*(ah(tidx)*ah(tidx)-2*bv(tidx)*bv(tidx))*&
&   (DU2(1)*(2*UU1(1)-wn1*VV1(1))/RD+DU1(1)*(2*UU2(1)-wn2*VV2(1))/RD)
    sigma2_b = sigma2_b-rho(tidx)*bv(tidx)*bv(tidx)*DV2(1)*(DV1(1)-VV1(1)/RD+wn1*UU1(1)/RD)*cos(sphi)
    sigma2_b = sigma2_b-rho(tidx)*bv(tidx)*bv(tidx)*DV1(1)*(dv2(1)-VV2(1)/RD+wn2*UU2(1)/RD)*cos(sphi)

    sigma2_a = -rho(tidx+1)*av(tidx+1)*av(tidx+1)*2*DU1(2)*DU2(2)
    sigma2_a = sigma2_a-ani(tidx+1)*rho(tidx+1)*(ah(tidx+1)*ah(tidx+1)-2*bv(tidx+1)*bv(tidx+1))*&
&   (DU2(2)*(2*UU1(2)-wn1*VV1(2))/RD+DU1(2)*(2*UU2(2)-wn2*VV2(2))/RD)
    sigma2_a = sigma2_a-rho(tidx+1)*bv(tidx+1)*bv(tidx+1)*DV2(2)*(DV1(2)-VV1(2)/RD+wn1*UU1(2)/RD)*cos(sphi)
    sigma2_a = sigma2_a-rho(tidx+1)*bv(tidx+1)*bv(tidx+1)*DV1(2)*(DV2(2)-VV2(2)/RD+wn2*UU2(2)/RD)*cos(sphi)

    sigma_total = (sigma1_a+sigma2_a)-(sigma1_b+sigma2_b)

! Love-Rayleigh
!    sigma1_b = 0.0
!    sigma1_a = 0.0
!    sigma2_b = 0.0
!    sigma2_a = 0.0

!    sigma1_b = rho(tidx)*omega_tmp*omega_tmp*VV2(1)*WW1(1)*sin(sphi)
!    sigma1_b = sigma1_b-rho(tidx)*bh(tidx)*bh(tidx)/RD/RD*wn1*wn2*VV2(1)*WW1(1)*sin(2*sphi)
!    sigma1_b = sigma1_b-rho(tidx)*bv(tidx)*bv(tidx)*(DV2(1)-VV2(1)/RD+wn2*UU2(1)/RD)*&
!&              (DW1(1)-WW1(1)/RD)*sin(sphi)
!
!    sigma1_a = rho(tidx+1)*omega_tmp*omega_tmp*VV2(2)*WW1(2)*sin(sphi)
!    sigma1_a = sigma1_a-rho(tidx+1)*bh(tidx+1)*bh(tidx+1)/RD/RD*wn1*wn2*VV2(2)*WW1(2)*sin(2*sphi)
!    sigma1_a = sigma1_a-rho(tidx+1)*bv(tidx+1)*bv(tidx+1)*(DV2(2)-VV2(2)/RD+wn2*UU2(2)/RD)*&
!&              (DW1(2)-WW1(2)/RD)*sin(sphi)

!    sigma2_b = rho(tidx)*bv(tidx)*bv(tidx)*DV2(1)*(DW1(1)-WW1(1)/RD)*sin(sphi)
!    sigma2_b = sigma2_b+rho(tidx)*bv(tidx)*bv(tidx)*DW1(1)*(DV2(1)-VV2(1)/RD+wn2*UU2(1)/RD)*sin(sphi)
    
!    sigma2_a = rho(tidx+1)*bv(tidx+1)*bv(tidx+1)*DV2(2)*(DW1(2)-WW1(2)/RD)*sin(sphi)
!    sigma2_a = sigma2_a+rho(tidx+1)*bv(tidx+1)*bv(tidx+1)*DW1(2)*(DV2(2)-VV2(2)/RD+wn2*UU2(2)/RD)*sin(sphi)

!    sigma_total = sigma_total+(sigma1_a+sigma2_a)-(sigma1_b+sigma2_b)
end if

! receiver term, ang2-pi/2 = ksi'-ksi"       
if(is_uvw .eq. 1) receiver_term = UU2(3)
if(is_uvw .eq. 2) receiver_term = VV2(3)*sin(ang2) - WW2(3)*cos(ang2)
if(is_uvw .eq. 3) receiver_term = VV2(3)*cos(ang2) + WW2(3)*sin(ang2)

! kernels
akern = sigma_total * receiver_term


return
end



subroutine dis1(is_major,blo,distdeg1,dist1)
use constant_para
use synthetic_ray

dist1 = distdeg1*pi/180    ! delta' in rad
if(is_major .eq. 1 .or. is_major .eq. 3) tlong = (-1)*blo
if(is_major.eq.0.or.is_major.eq.2) tlong = blo

if(tlong.gt.180.and.tlong.lt.360) dist1 =  twopi - dist1
if(tlong.gt.360.and.tlong.lt.540)  dist1 = twopi + dist1
if(tlong.gt.540.and.tlong.lt.720)  dist1 = 2*twopi - dist1

return
end



subroutine dis2(is_major,blo,distdeg2,dist2)
use constant_para
use synthetic_ray

dist2 = distdeg2*pi/180    ! delta" in rad
if(is_major.eq.1.or.is_major.eq.3) tlong = dist_act - (-1)*blo
if(is_major.eq.0.or.is_major.eq.2) tlong = dist_act - blo

if(tlong.gt.180.and.tlong.lt.360) dist2 =  twopi - dist2
if(tlong.gt.360.and.tlong.lt.540)  dist2 = twopi + dist2
if(tlong.gt.540.and.tlong.lt.720)  dist2 = 2*twopi - dist2

return
end


subroutine angcal(blo,bla,sphi,ang2)
use constant_para
use synthetic_ray

! All calculation is under Euler system
cla = 0
clo = 0
rla = bla
rlo = blo
call distaz(cla,clo,rla,rlo,dist1,az12,az21)    ! from source to scatter
cla = bla
clo = blo
rla = 0
rlo = dist_act    ! defined in module synthetic_ray
call distaz(cla,clo,rla,rlo,distdeg2,bz12,bz21)    ! from scatter to receiver

ak1 = az21+180
ak2 = bz12

sphi = ak2-ak1    ! ccw sphi>0
if(sphi .ge. 180) sphi = sphi-360
if(sphi .le. -180) sphi = sphi+360

ang2 = bz21-180    ! receiver arrival angle

return
end

subroutine convlt(f_w,g_w)
! g_w = f_w convolute Taper_w
use synthetic_ray
complex g_w(0:NPTS-1),f_w(0:NPTS-1)
complex f_t(0:NPTS-1),g_t(0:NPTS-1)
 
! f(w) -> f(t)
call fft(f_w,NPTS,1,f_t)
! g(t) = f(t) \time h(t)
g_t = (0,0)
do i=0,NPTS-1
    g_t(i) = f_t(i) * Taper(i)
end do
! g(t) -> g(w)
call fft(g_t,NPTS,-1,g_w)

return
end



