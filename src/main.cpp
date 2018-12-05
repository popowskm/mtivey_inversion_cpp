#include <vector>
#include <cmath>
#include <complex>
#include <numeric>
#include <algorithm>
#include <iostream>

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
std::vector<double> magfd(int date, int itype, double alt, double colat, double elong)
{
    // Initialize IGRFYEAR as 2015
    int igrfyear = 2015;
    std::vector<int> dgrf;
    for (int i = 1000; i < 2015; i += 5)
    {
        dgrf.push_back(i);
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
    std::vector<int> agh({-31543, -2298, 5922, -677, 2905, -1061, 924, 1121, 1022, -1469, -330, 1256, 3, 572, 523, 876, 628, 195, 660, -69, -361, -210, 134, -75, -184, 328, -210, 264, 53, 5, -33, -86, -124, -16, 3, 63, 61, -9, -11, 83, -217, 2, -58, -35, 59, 36, -90, -69, 70, -55, -45, 0, -13, 34, -10, -41, -1, -21, 28, 18, -12, 6, -22, 11, 8, 8, -4, -14, -9, 7, 1, -13, 2, 5, -9, 16, 5, -5, 8, -18, 8, 10, -20, 1, 14, -11, 5, 12, -3, 1, -2, -2, 8, 2, 10, -1, -2, -1, 2, -3, -4, 2, 2, 1, -5, 2, -2, 6, 6, -4, 4, 0, 0, -2, 2, 4, 2, 0, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
    std::vector<double> dgh({16.6000000000000, 12.8000000000000, -20, -13.8000000000000, 2.20000000000000, -17.4000000000000, -1, -8, 4.20000000000000, -5.60000000000000, 4.40000000000000, 0.200000000000000, 1.80000000000000, -8.60000000000000, -15, 0.200000000000000, 0, 3, -7, 0.800000000000000, 1, 2.60000000000000, -3.80000000000000, -1.40000000000000, 0, -0.200000000000000, 0, -2, 2.20000000000000, -1.80000000000000, 2, -0.200000000000000, 2.80000000000000, 3.80000000000000, 2, 1.40000000000000, 0.400000000000000, -0.200000000000000, 1.80000000000000, -2, 1.60000000000000, -0.400000000000000, -0.800000000000000, -1.20000000000000, 0.200000000000000, 0, 0.600000000000000, 2.40000000000000, 0, -1.60000000000000, 2.20000000000000, -0.200000000000000, 0.200000000000000, 0.400000000000000, 0.800000000000000, 1.20000000000000, 0.600000000000000, -0.200000000000000, 0, -0.200000000000000, -0.200000000000000, -0.400000000000000, -0.400000000000000, 0.400000000000000, 0.200000000000000, 0.200000000000000, -1, -0.400000000000000, 0.200000000000000, 0.400000000000000, -0.400000000000000, -0.200000000000000, 1.20000000000000, 0.600000000000000, 0.400000000000000, -0.200000000000000, -1.40000000000000, 0, -0.200000000000000, 1.20000000000000, 0, 0, 0, 0.400000000000000, 0, 0.400000000000000, 0.200000000000000, -0.200000000000000, 0.200000000000000, -0.800000000000000, -0.200000000000000, 0.200000000000000, -0.200000000000000, 0.600000000000000, -0.600000000000000, -0.600000000000000, -0.200000000000000, -0.400000000000000, 0.200000000000000, 0, -0.400000000000000, -0.200000000000000, 0, -0.200000000000000, 0.200000000000000, 0.200000000000000, 0.200000000000000, -0.200000000000000, 0, -0.200000000000000, -0.200000000000000, -0.200000000000000, 0.200000000000000, 0, 0.400000000000000, -0.400000000000000, -0.400000000000000, -0.200000000000000, 0, -0.200000000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
    int base = 1990;
    int i = 199;
    int t = 0;
    // WARNING ENDS

    double d2r = std::atan(1.0) * 4 / 180; // pi/180
    double R = alt;
    double slat = cos(colat * d2r);
    double clat = sin(colat * d2r);
    std::vector<double> cl(cos(elong * d2r));
    std::vector<double> sl(sin(elong * d2r));
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double cd = 1.0;
    double sd = 0.0;
    int l = 1;
    int m = 1;
    int n = 0;
    double re = 6371.2; // Earth's mean radius

    // Warning: Guess. Matlab code allowed r to be undefined/0 if it wasn't defined in the if.
    // Can divide by 0 in matlab - produces inf.
    double r = 1;
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
        double cd = (alt + four) / R;
        double sd = (a2 - b2) / four * slat * clat / R;
        one = slat;
        double slat = slat * cd - clat * sd;
        double clat = clat * cd + one * sd;
    }
    double ratio = re / r;

    std::vector<double> p({2.0 * slat, 2.0 * clat, 4.5 * slat * slat - 1.5, sqrt(27) * clat * slat});
    std::vector<double> q({-clat, slat, -3.0 * clat * slat, sqrt(3) * (slat * slat - clat * clat)});

    double nmax = 13; // Max number of harmonic degrees , 13

    double npq = (nmax * (nmax + 3)) / 2;

    double fn = 0;// GUESS - in matlab fn allowed to be potentially 0/undefined 
    double rr;
    for (int k = 0; k < npq; k++)
    {
        printf("Test %i\n", k);
        if (n < m)
        {
            m = 0;
            n = n + 1;
            rr = pow(ratio, (n + 2));
            fn = n;
        }

        double fm = m;
        printf("Nested ifs \n");
        if (k >= 5)
        {
            if ((m - n) == 0)
            {
                double one = sqrt(1.0 - 0.5 / fm);
                double j = k - n - 1;
                p[k] = (1.0 + 1.0 / fm) * one * clat * p[j];
                q[k] = one * (clat * q[j] + slat / fm * p[j]);
                sl[m] = sl[m - 2] * cl[0] + cl[m - 2] * sl[0];
                cl[m] = cl[m - 2] * cl[0] - sl[m - 2] * sl[0];
            }
            else
            {
                double one = sqrt(fn * fn - fm * fm);
                double two = sqrt(pow((fn - 1.0), 2) - fm * fm) / one;
                double three = (2.0 * fn - 1.0) / one;
                i = k - n;
                double j = k - 2 * n + 1;
                p[k] = (fn + 1.0) * (three * slat / fn * p[i] - two / (fn - 1.0) * p[j]);
                q[k] = three * (slat * q[i] - clat / fn * p[i]) - two * q[j];
            }
        }
        //
        //     SYNTHESIS OF x, y AND z IN GEOCENTRIC COORDINATES
        //
        double one = (agh[l-1] + dgh[l-1] * t) * rr;

        printf("Synthesis\n");
        if(m == 0){ 
            x = x + one * q[k];
            z = z - one * p[k];
            l = l + 1;
        }
        else{
            double two = (agh[l] + dgh[l] * t) * rr;
            double three = one * cl[m-1] + two * sl[m-1];
            x = x + three * q[k];
            z = z - three * p[k];
            if(clat > 0){
                y = y + (one * sl[m-1] - two * cl[m-1]) * fm * p[k] / ((fn + 1.0) * clat);
            }
            else{
                y = y + (one * sl[m-1] - two * cl[m-1]) * q[k] * slat;
            }
            l = l + 2;
        }
        m = m + 1;
    }

    std::vector<double> b({1,2,3,4});
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
std::vector<std::vector<double>> inv3d(std::vector<std::vector<double>> f3d, std::vector<std::vector<double>> h, float wl, float ws, float rlat, float rlon, int yr, float zobs, float thick, float slin, float dx, float dy, float sdec, float sdip)
{
    //error
    std::vector<std::vector<double>> a;
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

    printf(" Zobs= %12.5f\n Rlat= %12.5f Rlon= %12.5f\n", zobs, rlat, rlon);
    printf(" Yr= %12.5f\n", yr);
    printf(" Thick= %12.5f\n", thick);
    printf(" slin = %12.6f\n", slin);
    printf(" Nterms,Tol %6.0f %10.5f \n", nterms, tol);

    if (h.empty() || f3d.empty())
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
    printf(" READ %6.0f x %6.0f matrix by columns \n", nx, ny);
    printf(" DX,DY= %10.3f %10.3f XMIN,YMIN= %10.3f  %10.3f\n", dx, dy, xmin, xmin);

    // remove mean from input field
    // double mnf3d=std::accumulate(f3d.begin(), f3d.end(), 0.0)/f3d.size();
    double total = 0;
    std::for_each(f3d.begin(), f3d.end(), [&total](std::vector<double> v) { std::for_each(v.begin(), v.end(), [&total](double &d) { total += d; }); });

    const double mnf3d = total / (f3d.size() * f3d[0].size());
    std::for_each(f3d.begin(), f3d.end(), [mnf3d](std::vector<double> &v) { std::for_each(v.begin(), v.end(), [mnf3d](double &d) { d -= mnf3d; }); });
    printf("Remove mean of %10.3f from field \n", mnf3d);

    double colat = 90. - rlat;
    std::vector<double> y = magfd(yr, 1, zobs, colat, rlon);

    for (int i = 0; i < f3d.size(); i++)
    {
        for (int j = 0; j < f3d[i].size(); j++)
        {
            std::cout << f3d[i][j] << ' ';
        }
        std::cout << "\n";
    }
    return a;
}

int main(int argc, char const *argv[])
{
    std::vector<std::vector<double>> f3d;
    for (int i = 0; i < 10; i++)
    {
        std::vector<double> temp;
        for (int j = 0; j < 10; j++)
        {
            temp.push_back(i + j);
        }
        f3d.push_back(temp);
    }
    std::vector<std::vector<double>> h;
    for (int i = 0; i < 10; i++)
    {
        std::vector<double> temp;
        for (int j = 0; j < 10; j++)
        {
            temp.push_back(j);
        }
        h.push_back(temp);
    }
    float wl = 0;
    float ws = 0;
    float rlat = 0;
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
    inv3d(f3d, h, wl, ws, rlat, rlon, yr, zobs, thick, slin, dx, dy, sdec, sdip);
    return 0;
}
