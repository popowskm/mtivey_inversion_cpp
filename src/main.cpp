#include <valarray>
#include <cmath>
#include <complex>
#include <numeric>
#include <algorithm>
#include <iostream>

std::valarray<double> magfd(int date, int itype, double alt, double colat, double elong);

void prints(std::valarray<std::valarray<double>> f3d)
{
    for (int i = 0; i < f3d.size(); i++)
    {
        for (int j = 0; j < f3d[i].size(); j++)
        {
            std::cout << f3d[i][j] << ' ';
        }
        std::cout << "\n";
    }
}

void prints1(std::valarray<double> f3d)
{
    for (int i = 0; i < f3d.size(); i++)
    {
        std::cout << f3d[i] << " ";
    }
    std::cout << "\n";
}

std::valarray<double> nskew(double yr, double rlat, double rlon, double zobs, double slin, double sdec, double sdip, bool opts)
{
    // NSKEW - Compute skewness parameter and amplitude factor
    //  following Schouten (1971)
    //  Computes GEOCENTRIC DIPOLE unless given
    //  declination and dip of magnetization
    // Usage:
    //    [theta,ampfac]=nskew(yr,rlat,rlon,zobs,slin)
    // or
    //    [theta,ampfac]=nskew(yr,rlat,rlon,zobs,slin,sdec,sdip)
    //
    //  Input variables:
    //   yr : decimal year of survey
    //   rlat,rlon : regional latitude, longitude in decimal degrees
    //   zobs : level of observation in km above sealevel
    //   slin : strike of lineations normal to profile (+cw degrees from north)
    //   sdec,sdip : magnetization declination,inclination
    // Output variables
    //   theta : phase angle
    //   ampfac : amplitude factor
    // Calls <magfd>
    //
    // Maurice A. Tivey February 3, 1993
    //    checked April 1996
    // MATLAB 5 Nov 1998
    // See skew.m for more general calculation
    //---------------------------------------------------------
    double rad = atan(1.0) * 4 / 180;
    // get unit valarrays
    double colat = 90.0 - rlat;
    std::valarray<double> y = magfd(abs(yr), 1, zobs, colat, rlon);
    // compute skewness parameter
    double bx = y[0];
    double by = y[1];
    double bz = y[2];
    double bh = sqrt(pow(bx, 2) + pow(by, 2));
    double decl1 = atan2(by, bx) / rad;
    double incl1 = atan2(bz, bh) / rad;
    if (yr > 0){
        printf(" EARTH' 'S MAGNETIC FIELD DIRECTION:\n");
        printf(" //10.3f = MAGNETIC DECLINATION ( STRIKE, CW FROM N )\n", decl1);
        printf(" //10.4f = MAGNETIC INCLINATION ( DIP, POS DOWN )\n", incl1);
    }
    if (opts){
        //  NOTE FOR GEOCENTRIC DIPOLE TAN(INC)=2*TAN(LAT)
        //if abs(sdec) > 0. | abs(sdip) > 0.
        if (yr > 0) {
            printf(" NON-GEOCENTRIC MAGNETIZATION valarray SPECIFIED:\n");
            printf(" //10.4f = DESIRED MAGNETIZATION DECLINATION (+CW FROM N)\n", sdec);
            printf(" //10.4f = DESIRED MAGNETIZATION INCLINATION (+DN)\n", sdip);
        } 
    }else {
        sdip = atan2(2. * sin(rlat * rad), cos(rlat * rad)) / rad;
        sdec = 0;
        if (yr > 0) {
            printf(" GEOCENTRIC MAGNETIZATION valarray SPECIFIED:\n");
            printf(" //10.4f = GEOCENTRIC DIPOLE INCLINATION \n", sdip);
            printf(" //10.3f = GEOCENTRIC DECLINATION ASSUMED\n", sdec);
        }
    }
    // compute phase and amplitude factors
    double ra1 = incl1 * rad;
    double rb1 = (decl1 - slin) * rad;
    double ra2 = sdip * rad;
    double rb2 = (sdec - slin) * rad;
    // compute phase and amplitude factors
    double inclm = atan2(tan(ra2), sin(rb2));
    double inclf = atan2(tan(ra1), sin(rb1));
    double ampfac = ((sin(ra2)) * (sin(ra1))) / ((sin(inclm)) * (sin(inclf)));
    double theta = (inclm / rad) + (inclf / rad) - 180.;
    if (theta <= -360){
        theta = theta + 360;
    }else{
        if (theta >= 360){
        theta = theta - 360;
        }
    }
    // compute unit valarrays for a check
    std::valarray<double> hatm;
    std::valarray<double>hatb;
    hatm[0] = (cos(sdip * rad) * sin((sdec - slin) * rad));
    hatm[1] = (cos(sdip * rad) * cos((sdec - slin) * rad));
    hatm[2] = (-sin(sdip * rad));
    hatb[0] = (cos(incl1 * rad) * sin((decl1 - slin) * rad));
    hatb[1] = (cos(incl1 * rad) * cos((decl1 - slin) * rad));
    hatb[2] = (-sin(incl1 * rad));
    //
    if (yr > 0) {
        printf("  //10.6f //10.6f //10.6f = MAGNETIZATION UNIT valarray\n", hatm[0], hatm[1], hatm[2]);
        printf("  //10.6f //10.6f //10.6f = AMBIENT FIELD UNIT valarray\n", hatb[0], hatb[1], hatb[2]);
        printf("  COMPONENTS ARE (X,Y,Z=ALONG, ACROSS PROFILE, AND UP\n\n");
    }

    std::valarray<double> a({theta, ampfac});
    return a;
}

//  MAGFD
//  Function to compute Earths magnetic field
//  and components: x,y,z,t for a given latitude
//  and longitude, date and altitude.
//  Uses MATLAB MAT files sh1900.mat to sh2015.mat in 5 yr
//  intervals.
//
//  DATE = date of survey (decimal years)
//  ITYPE=1 for geodetic coordinates (usual case)
//  ITYPE=0 for geocentric coordinates
//  alt = (for ITYPE=1) altitude of survey relative to sealevel (km +ve up)
//  alt = (for ITYPE=0) radial distance from center of earth in km
//  colat=90-latitude (decimal degrees)
//  elong=longitude of survey (decimal degrees)
//
//  Output array out contains components x,y,z,t in nanoteslas
//   x north component
//   y east component
//   z vertical component +ve down
//   t total field magnitude
//
//  Usage: out=magfd(DATE,ITYPE,alt,colat,elong);
//
//  ref: IAGA, Division V, Working Group VMOD,
//   International Geomagnetic Reference Field: the 12th generation,
//   Earth Planets and Space, 67 (79), 2014. doi:10.1186/s40623-015-0228-9
//
// March 1997
// Mod Dec 1999 (add igrf2000 and y2k compliance
// Mod Nov 2000 (use up to degree 10 sh coefficients)
// Mod Apr 2005 added 2005 coeffs
// Mod Sep 2006 some clean up and info added
// Mod 2010 coeffs Ref: IGRF11 Finlay et al., 2010
// Mod Jan 2015 coeffs Ref: IGRF12 IAGA V-MOD Working Group
// Mod Jun 2017 uses coefficients thru degree 13
//
// http://deeptow.whoi.edu/matlab.html
// Copyright: Maurice A. Tivey, 2017
// Woods Hole Oceanographic Institution
std::valarray<double> magfd(int date, int itype, double alt, double colat, double elong)
{
    // Initialize IGRFYEAR as 2015
    int igrfyear = 2015;
    std::valarray<int> dgrf;
    int j = 0;
    for (int i = 1000; i < 2015; i += 5)
    {
        dgrf[j] = i;
        j++;
    }
    std::string igrffile = "sh" + std::to_string(igrfyear);

    // simple switch if printout needed
    // negative DATE means don't print out
    int pl = 0;
    if (date < 0)
        pl = 1;
    date = abs(date);

    // Matlab scripts load several .mat files depending on the date with hardcoded variables
    // TODO: Determine the best way to supply the data to the program -- CSV?

    // WARNING: TEMPORARILY HARDCODED VARIABLES
    std::valarray<int> agh({-31543, -2298, 5922, -677, 2905, -1061, 924, 1121, 1022, -1469, -330, 1256, 3, 572, 523, 876, 628, 195, 660, -69, -361, -210, 134, -75, -184, 328, -210, 264, 53, 5, -33, -86, -124, -16, 3, 63, 61, -9, -11, 83, -217, 2, -58, -35, 59, 36, -90, -69, 70, -55, -45, 0, -13, 34, -10, -41, -1, -21, 28, 18, -12, 6, -22, 11, 8, 8, -4, -14, -9, 7, 1, -13, 2, 5, -9, 16, 5, -5, 8, -18, 8, 10, -20, 1, 14, -11, 5, 12, -3, 1, -2, -2, 8, 2, 10, -1, -2, -1, 2, -3, -4, 2, 2, 1, -5, 2, -2, 6, 6, -4, 4, 0, 0, -2, 2, 4, 2, 0, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
    std::valarray<double> dgh({16.6000000000000, 12.8000000000000, -20, -13.8000000000000, 2.20000000000000, -17.4000000000000, -1, -8, 4.20000000000000, -5.60000000000000, 4.40000000000000, 0.200000000000000, 1.80000000000000, -8.60000000000000, -15, 0.200000000000000, 0, 3, -7, 0.800000000000000, 1, 2.60000000000000, -3.80000000000000, -1.40000000000000, 0, -0.200000000000000, 0, -2, 2.20000000000000, -1.80000000000000, 2, -0.200000000000000, 2.80000000000000, 3.80000000000000, 2, 1.40000000000000, 0.400000000000000, -0.200000000000000, 1.80000000000000, -2, 1.60000000000000, -0.400000000000000, -0.800000000000000, -1.20000000000000, 0.200000000000000, 0, 0.600000000000000, 2.40000000000000, 0, -1.60000000000000, 2.20000000000000, -0.200000000000000, 0.200000000000000, 0.400000000000000, 0.800000000000000, 1.20000000000000, 0.600000000000000, -0.200000000000000, 0, -0.200000000000000, -0.200000000000000, -0.400000000000000, -0.400000000000000, 0.400000000000000, 0.200000000000000, 0.200000000000000, -1, -0.400000000000000, 0.200000000000000, 0.400000000000000, -0.400000000000000, -0.200000000000000, 1.20000000000000, 0.600000000000000, 0.400000000000000, -0.200000000000000, -1.40000000000000, 0, -0.200000000000000, 1.20000000000000, 0, 0, 0, 0.400000000000000, 0, 0.400000000000000, 0.200000000000000, -0.200000000000000, 0.200000000000000, -0.800000000000000, -0.200000000000000, 0.200000000000000, -0.200000000000000, 0.600000000000000, -0.600000000000000, -0.600000000000000, -0.200000000000000, -0.400000000000000, 0.200000000000000, 0, -0.400000000000000, -0.200000000000000, 0, -0.200000000000000, 0.200000000000000, 0.200000000000000, 0.200000000000000, -0.200000000000000, 0, -0.200000000000000, -0.200000000000000, -0.200000000000000, 0.200000000000000, 0, 0.400000000000000, -0.400000000000000, -0.400000000000000, -0.200000000000000, 0, -0.200000000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
    int base = 1990;
    int i = 199;
    double t = 0;
    // WARNING ENDS

    double d2r = std::atan(1.0) * 4 / 180; // pi/180
    double r = alt;
    double slat = cos(colat * d2r);
    double clat = sin(colat * d2r);
    std::valarray<double> cl({cos(elong * d2r)});
    std::valarray<double> sl({sin(elong * d2r)});
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double cd = 1.0;
    double sd = 0.0;
    int l = 1;
    int m = 1;
    int n = 0;
    double re = 6371.2; // Earth's mean radius

    if (itype == 1)
    { // CONVERSION FROM GEODETIC TO GEOCENTRIC COORDINATES
        //a2    = 40680925.;  // squared semi major axis
        //b2    = 40408588.;  // squared semi minor axis
        // WGS84
        double a2 = 40680631.6; //6378.137^2;  // squared semi major axis
        double b2 = 40408296.0; //6356.7523^2;  // squared semi minor axis
        double one = a2 * clat * clat;
        double two = b2 * slat * slat;
        double three = one + two;
        double four = sqrt(three);
        r = sqrt(alt * (alt + 2.0 * four) + (a2 * one + b2 * two) / three);
        cd = (alt + four) / r;
        sd = (a2 - b2) / four * slat * clat / r;
        one = slat;
        slat = slat * cd - clat * sd;
        clat = clat * cd + one * sd;
    }
    double ratio = re / r;

    std::valarray<double> p({2.0 * slat, 2.0 * clat, 4.5 * slat * slat - 1.5, sqrt(27) * clat * slat});
    std::valarray<double> q({-clat, slat, -3.0 * clat * slat, sqrt(3) * (slat * slat - clat * clat)});

    double nmax = 13; // Max number of harmonic degrees , 13

    double npq = (nmax * (nmax + 3)) / 2;

    double fn = 0; // GUESS - in matlab fn allowed to be potentially 0/undefined
    double rr;
    for (int k = 0; k < npq; k++)
    {
        if (n < m)
        {
            m = 0;
            n = n + 1;
            rr = pow(ratio, (n + 2));
            fn = n;
        }

        double fm = m;
        if (k >= 4)
        {
            if ((m - n) == 0)
            {
                double one = sqrt(1.0 - 0.5 / fm);
                double j = k - n - 1;
                p[k] = ((1.0 + 1.0 / fm) * one * clat * p[j]);
                q[k] = (one * (clat * q[j] + slat / fm * p[j]));
                sl[m] = (sl[m - 2] * cl[0] + cl[m - 2] * sl[0]);
                cl[m] = (cl[m - 2] * cl[0] - sl[m - 2] * sl[0]);
            }
            else
            {
                double one = sqrt(fn * fn - fm * fm);
                double two = sqrt(pow((fn - 1.0), 2) - fm * fm) / one;
                double three = (2.0 * fn - 1.0) / one;
                i = k - n;
                double j = k - 2 * n + 1;
                p[k] = ((fn + 1.0) * (three * slat / fn * p[i] - two / (fn - 1.0) * p[j]));
                q[k] = (three * (slat * q[i] - clat / fn * p[i]) - two * q[j]);
            }
        }
        //
        //     SYNTHESIS OF x, y AND z IN GEOCENTRIC COORDINATES
        //
        double one = (agh[l - 1] + dgh[l - 1] * t) * rr;

        if (m == 0)
        {
            x = x + one * q[k];
            z = z - one * p[k];
            l = l + 1;
        }
        else
        {
            double two = (agh[l] + dgh[l] * t) * rr;
            double three = one * cl[m - 1] + two * sl[m - 1];
            x = x + three * q[k];
            z = z - three * p[k];
            if (clat > 0)
            {
                y = y + (one * sl[m - 1] - two * cl[m - 1]) * fm * p[k] / ((fn + 1.0) * clat);
            }
            else
            {
                y = y + (one * sl[m - 1] - two * cl[m - 1]) * q[k] * slat;
            }
            l = l + 2;
        }
        m = m + 1;
    }

    //  CONVERSION TO COORDINATE SYSTEM SPECIFIED BY ITYPE

    double one = x;
    x = x * cd + z * sd;
    z = z * cd - one * sd;
    t = sqrt(x * x + y * y + z * z);

    std::valarray<double> b({x, y, z, t});
    return b;
}

// INV3D Calculate magnetization from magnetic field and bathymetry for a map.
// Assumes constant thickness source layer whose upper bound is bathymetry
// Use the Parker & Huestis [1974] Fourier inversion approach.
//
// Usage: m3d=inv3d(f3d,h,wl,ws,rlat,
//                      rlon,yr,zobs,thick,slin,dx,dy,sdec,sdip);
//   or for geocentric dipole
//     m3d=inv3d(f3d,h,wl,ws,rlat,rlon,yr,zobs,thick,slin,dx,dy);
//
// Input arrays:
//    f3d 	magnetic field (nT)
//    h 		bathymetry (km +ve up)
//    wl		filter long wavelength cutoff (km)
//    ws		filter short wavelength cutoff (km)
//    rlat 	latitude of survey area dec. deg.
//    rlon 	longitude of survey area dec. deg.
//    yr 	year of survey (dec. year)
//    slin 	azimuth of lineations (deg)
//    zobs 	observation level (+km up)
//    thick 	thickness of source layer (km)
//    slin	azimuth of grid (degrees) hard wired to 0
//    dx 	x grid spacing  (km)
//    dy 	y grid spacing  (km)
//    sdec	declination of magnetization (optional)
//    sdip	inclination of magnetization (optional)
// Output array:
//    m3d	magnetization (A/m)
//
//
//
// Maurice A. Tivey  MATLAB August 27 1992
// MAT May  5 1995
// MAT Mar 1996 (new igrf)
// calls <syn3d,magfd,nskew,bpass3d>
std::valarray<std::valarray<double>> inv3d(std::valarray<std::valarray<double>> f3d, std::valarray<std::valarray<double>> h, float wl, float ws, float rlat, float rlon, int yr, float zobs, float thick, float slin, float dx, float dy, float sdec, float sdip)
{
    //error
    std::valarray<std::valarray<double>> a;
    // parameters defined
    const std::complex<double> i_math(0.0, 1.0);
    const double pi = std::atan(1.0) * 4;
    const double rad = pi / 180; // conversion radians to degrees
    const int mu = 100;          // conversion factor to nT
    // changeable parameters
    int nterms = 20;
    int nitrs = 20;
    float tol = 0.0001;
    float tolmag = 0.0001;
    int flag = 0;
    int xmin = 0;

    printf("       3D MAGNETIC INVERSE MODEL\n");
    printf("                  INV3D\n");
    printf("        Constant thickness layer\n");
    printf("Version : 2/24/2015\n");

    printf(" Zobs= //12.5f\n Rlat= //12.5f Rlon= //12.5f\n", zobs, rlat, rlon);
    printf(" Yr= //12.5f\n", yr);
    printf(" Thick= //12.5f\n", thick);
    printf(" slin = //12.6f\n", slin);
    printf(" Nterms,Tol //6.0f //10.5f \n", nterms, tol);

    if (h[0][0] == NULL || f3d[0][0] == NULL)
    {
        printf("Bathy and field arrays must have values\n");
        return a;
    }
    int ny = f3d.size();
    int nx = f3d[0].size();
    // print out the input files header
    if ((h.size() != ny) || (h[0].size() != nx))
    {
        printf(" bathy and field arrays must be of the same length\n");
        return a;
    }
    printf(" READ //6.0f x //6.0f matrix by columns \n", nx, ny);
    printf(" DX,DY= //10.3f //10.3f XMIN,YMIN= //10.3f  //10.3f\n", dx, dy, xmin, xmin);

    // remove mean from input field
    // double total = 0;
    // std::for_each(f3d.begin(), f3d.end(), [&total](std::valarray<double> v) { std::for_each(v.begin(), v.end(), [&total](double &d) { total += d; }); });

    // const double mnf3d = total / (f3d.size() * f3d[0].size());
    // std::for_each(f3d.begin(), f3d.end(), [mnf3d](std::valarray<double> &v) { std::for_each(v.begin(), v.end(), [mnf3d](double &d) { d -= mnf3d; }); });
    double total = 0;
    printf("%5.5f", f3d.sum());
    std::for_each(std::begin(f3d), std::end(f3d), [&total](std::valarray<double> &v) { total += v.sum();});
    const double mnf3d = total / (f3d.size() * f3d[0].size());
    printf("Remove mean of //10.3f from field \n", mnf3d);

    double colat = 90. - rlat;
    std::valarray<double> y = magfd(yr, 1, zobs, colat, rlon);

    double bx = y[0];
    double by = y[1];
    double bz = y[2];
    double bh = sqrt(pow(bx, 2) + pow(by, 2));
    double decl1 = atan2(by, bx) / rad;
    double incl1 = atan2(bz, bh) / rad;
    double theta;
    double ampfac;
    std::valarray<double> skew_array;
    if (rlat == 90)
    { // rtp anomaly
        incl1 = 90;
        decl1 = 0;
        sdip = 90;
        sdec = 0;
        printf("Assume an RTP anomaly\n");
    }
    else
    {
        if (abs(sdec) > 0. || abs(sdip) > 0.)
        {
            skew_array = nskew(yr, rlat, rlon, zobs, slin, sdec, sdip, true); //skew_array = [theta, ampfac]
        }
        else
        {
            skew_array = nskew(yr, rlat, rlon, zobs, slin, 0, 0, false);
            sdip = atan2(2. * sin(rlat * rad), cos(rlat * rad)) / rad;
            sdec = 0;
        }
        theta = skew_array[0];
        ampfac= skew_array[1];
    }

    slin=0; // slin is forced to zero
    double ra1=incl1*rad;
    double rb1=(decl1-slin)*rad;
    // rb1=(slin-decl1)*rad;
    double ra2=sdip*rad;
    double rb2=(sdec-slin)*rad;
    // rb2=(slin-sdec)*rad;

    // make wave number array
    // ni=1/nx;
    double nx2=nx/2;
    double nx2plus=nx2+1;
    // x=-.5:ni:.5-ni;
    // ni=1/ny;
    double ny2=ny/2;
    double ny2plus=ny2+1;
    // y=-.5:ni:.5-ni;
    // X=ones(size(y))'*x;
    // Y=y'*ones(size(x));
    // k=2*pi*sqrt(X.^2+Y.^2);  // wavenumber array
    // k= fftshift(k);
    // compute another way
    double dkx=pi/(nx*dx);
    double dky=pi/(ny*dy);

    //TODO
    // std::valarray<double> kx=(-nx2:nx2-1).*dkx;
    // ky=(-ny2:ny2-1).*dky;
    // X=ones(size(ky))'*kx;
    // Y=ky'*ones(size(kx));
    // k= 2*sqrt(X.^2+Y.^2);  // wavenumber array
    // k=fftshift(k);

    return a;
}

int main()
{
    std::valarray<std::valarray<double>> f3d;
    for (int i = 0; i < 10; i++)
    {
        std::valarray<double> temp;
        for (int j = 0; j < 10; j++)
        {
            temp[j] = (i + j);
        }
        f3d[i] = (temp);
    }
    std::valarray<std::valarray<double>> h;
    for (int i = 0; i < 10; i++)
    {
        std::valarray<double> temp;
        for (int j = 0; j < 10; j++)
        {
            temp[j] = (j);
        }
        h[i] = (temp);
    }
    float wl = 0;
    float ws = 0;
    float rlat = 10;
    float rlon = 0;
    int yr = 1990;
    float zobs = 0;
    float thick = 0;
    float slin = 0;
    float dx = 0;
    float dy = 0;

    // Optional values, default assumes geocentric dipole hypothesis
    float sdec = 0;
    float sdip = 0;
    prints(f3d);
    prints(h);
    inv3d(f3d, h, wl, ws, rlat, rlon, yr, zobs, thick, slin, dx, dy, sdec, sdip);
    return 0;
}
