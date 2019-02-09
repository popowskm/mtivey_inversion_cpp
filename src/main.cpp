#include <iomanip>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <sstream>
#include <complex>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fftw3.h>

#define pi 3.141592653589793
std::vector<double> magfd(int date, int itype, double alt, double colat, double elong);

std::complex<double> twodmax(std::vector<std::vector<std::complex<double>>> base){
    std::vector<std::complex<double>> temp;
    for(auto vec: base){
        temp.push_back(*std::max_element(vec.begin(), vec.end(), [](std::complex<double> first, std::complex<double> second) {return (first.real() < second.real());}));
    }
    std::complex<double> max = *std::max_element(temp.begin(), temp.end(), [](std::complex<double> first, std::complex<double> second) {return (first.real() < second.real());});
    return max;
}

std::complex<double> twodmin(std::vector<std::vector<std::complex<double>>> base){
    std::vector<std::complex<double>> temp;
    for(auto vec: base){
        temp.push_back(*std::min_element(vec.begin(), vec.end(), [](std::complex<double> first, std::complex<double> second) {return (first.real() < second.real());}));
    }
    std::complex<double> min = *std::min_element(temp.begin(), temp.end(), [](std::complex<double> first, std::complex<double> second) {return (first.real() < second.real());});
    return min;
}

double realtwodmax(std::vector<std::vector<double>> base){
    std::vector<double> temp;
    for(auto vec: base){
        temp.push_back(*std::max_element(vec.begin(), vec.end()));
    }
    double max = *std::max_element(temp.begin(), temp.end());
    return max;
}

double realtwodmin(std::vector<std::vector<double>> base){
    std::vector<double> temp;
    for(auto vec: base){
        temp.push_back(*std::min_element(vec.begin(), vec.end()));
    }
    double min = *std::min_element(temp.begin(), temp.end());
    return min;
}

void writes(std::vector<std::vector<std::complex<double>>> a, std::string b){
    std::ofstream file(b);
    for(int i = 0; i < a.size(); i++){
        auto temps = a[i].size();
        for(int j = 0; j < a[i].size()-1; j++){
            file << std::setprecision(16) << a[i][j].real();
            if(a[i][j].imag()>=0) file << std::setprecision(16) << '+';
            file << std::setprecision(16) << a[i][j].imag() << "i,";
        }
        file << std::setprecision(16) << a[i][a.size()-1].real();
        if(a[i][a.size()-1].imag()>=0) file << std::setprecision(16) << '+';
        file << std::setprecision(16) << a[i][a.size()-1].imag() << "i" << std::endl;
    }
    file.close();
}

void reads(std::string filename, std::vector<std::vector<double>> &v){
    std::ifstream file(filename);
    int l = 0;

    while (file) {
        l++;
        std::string s;
        if (!getline(file, s)) break;
        if (s[0] != '#') {
            std::istringstream ss(s);
            std::vector<double> record;
 
            while (ss) {
                std::string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    record.push_back(stod(line));
                }
                catch (const std::invalid_argument e) {
                    std::cout << "NaN found in file " << filename << " line " << l
                         << std::endl;
                    e.what();
                }
            }
 
            v.push_back(record);
        }
    }
 
    if (!file.eof()) {
        std::cerr << "Could not read file " << filename << "\n";
    }
 
    file.close();
}

void prints(std::vector<std::vector<double>> f3d)
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

//MATLAB FUNCTION REPLICATION
//Rotates halfway in x and y
std::vector<std::vector<double>> fftshift(std::vector<std::vector<double>> a)
{
    std::vector<std::vector<double>>::const_iterator first = a.begin() + ceil(a.size() / 2);
    std::vector<std::vector<double>>::const_iterator last = a.end();
    std::vector<std::vector<double>> temp(first, last);
    for (auto ptr = a.begin(); ptr < first; ptr++)
    {
        temp.push_back(*ptr);
    }

    for (std::vector<double> &j : temp)
    {
        std::vector<double>::const_iterator first2 = j.begin() + ceil(j.size() / 2);
        std::vector<double>::const_iterator last2 = j.end();
        std::vector<double> temp2(first2, last2);
        for (auto ptr = j.begin(); ptr < first2; ptr++)
        {
            temp2.push_back(*ptr);
        }
        j = temp2;
    }

    return temp;
}

//Phase angle
std::complex<double> angle(std::complex<double> a){
    return atan2(a.imag(), a.real());
}

std::vector<std::vector<std::complex<double>>> fftshift_complex(std::vector<std::vector<std::complex<double>>> a)
{
    std::vector<std::vector<std::complex<double>>>::const_iterator first = a.begin() + ceil(a.size() / 2);
    std::vector<std::vector<std::complex<double>>>::const_iterator last = a.end();
    std::vector<std::vector<std::complex<double>>> temp(first, last);
    for (auto ptr = a.begin(); ptr < first; ptr++)
    {
        temp.push_back(*ptr);
    }

    for (std::vector<std::complex<double>> &j : temp)
    {
        std::vector<std::complex<double>>::const_iterator first2 = j.begin() + ceil(j.size() / 2);
        std::vector<std::complex<double>>::const_iterator last2 = j.end();
        std::vector<std::complex<double>> temp2(first2, last2);
        for (auto ptr = j.begin(); ptr < first2; ptr++)
        {
            temp2.push_back(*ptr);
        }
        j = temp2;
    }

    return temp;
}

void fft2(std::vector<std::vector<std::complex<double>>> &input, std::vector<std::vector<std::complex<double>>> &output){
    fftw_complex *in, *out;
    fftw_plan p;

    int size = input.size(), size2 = input[0].size();
    
    in  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size * size2);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size * size2);
    p   = fftw_plan_dft_2d(size, size2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    int i = 0;
    for (int x = 0; x < size; x++)
    {
        for (int z = 0; z < size2; z++)
        {
            in[i][0] = std::real(input[x][z]);
            in[i][1] = std::imag(input[x][z]);
            i++;
        }
    }

    fftw_execute(p);
    i = 0;
    for (int x = 0; x < size; x++)
    {
        for (int z = 0; z < size2; z++)
        {
            output[x][z] = std::complex<double>(out[i][0], out[i][1]);
            i++;
        }
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

void ifft2(std::vector<std::vector<std::complex<double>>> &input, std::vector<std::vector<std::complex<double>>> &output){
    fftw_complex *in, *out;
    fftw_plan p;

    int size = input.size(), size2 = input[0].size();
    
    in  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size * size2);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size * size2);
    p   = fftw_plan_dft_2d(size, size2, in, out, FFTW_BACKWARD, FFTW_PATIENT);

    int i = 0;
    for (int x = 0; x < size; x++)
        for (int z = 0; z < size2; z++)
        {
            in[i][0] = std::real(input[x][z]);
            in[i][1] = std::imag(input[x][z]);
            i++;
        }
    fftw_execute(p);

    i = 0;
    for (int x = 0; x < size; x++){
        for (int z = 0; z < size2; z++)
        {
            output[x][z] = std::complex<double>(out[i][0], out[i][1]);
            i++;
        }
    }
}


double nfac(double N){
    // NFAC - calculate N-factorial
    //
    //  Maurice A. Tivey 30-Oct-90

    double nsum=1;
    for (int i=1; i<N+1;i++){
        nsum=nsum*i;
    }
    return nsum;
}

std::vector<std::vector<double>> bpass3d(double nnx, double nny, double dx, double dy, double wlong, double wshort){
// BPASS3D set up bandpass filter weights in 2 dimensions
// using a cosine tapered filter
// Usage:  wts3d=bpass3d(nnx,nny,dx,dy,wlong,wshort);
//
// Maurice A. Tivey MATLAB March 1996
// MAT Jun 2006
// Calls <>
//---------------------------------------------------
    double twopi=pi*2;
    double dk1=2*pi/((nnx-1)*dx);
    double dk2=2*pi/((nny-1)*dy);
    // calculate wavenumber array
    double nx2=nnx/2;
    double nx2plus=nx2+1;
    double ny2=nny/2;
    double ny2plus=ny2+1;
    double dkx=2*pi/(nnx*dx);
    double dky=2*pi/(nny*dy);
    std::vector<double> kx;
    for(int i = -nx2; i < nx2; i++){
        kx.push_back(i*dkx);
    }
    std::vector<double> ky;
    for(int i = -ny2; i < ny2; i++){
        ky.push_back(i*dky);
    }
    std::vector<std::vector<double>> X(ky.size(), ky);
    std::vector<std::vector<double>> Y;
    for (auto num : kx)
    {
        std::vector<double> temp(kx.size(), num);
        Y.push_back(temp);
    }
    std::vector<std::vector<double>> k;
    for (int i = 0; i < X.size(); i++)
    {
        std::vector<double> temp;
        for (int j = 0; j < X[0].size(); j++)
        {
            temp.push_back(sqrt(pow(X[i][j], 2) + pow(Y[i][j], 2)));
        }
        k.push_back(temp);
    } // wavenumber array
    k=fftshift(k);

//
    if(wshort==0)
        wshort = (dx*2 > dy*2) ? dx*2: dy*2;
    if(wlong==0)
        wlong = (nnx*dx < nny*dy) ? nnx*dx: nny*dy; 

    double klo=twopi/wlong;
    double khi=twopi/wshort;
    double khif=0.5*khi;
    double klof=2*klo;
    double dkl=klof-klo;
    double dkh=khi-khif;
    printf(" BPASS3D\n SET UP BANDPASS WEIGHTS ARRAY :\n");
    printf(" HIPASS COSINE TAPER FROM K= %10.6f TO K= %10.6f\n",klo,klof);
    printf(" LOPASS COSINE TAPER FROM K= %10.6f TO K= %10.6f\n",khif,khi);
    printf(" DK1,DK2= %10.4f  %10.4f\n",dk1,dk2);

    double wl1=1000;
    double  wl2=1000;
    if(klo>0) wl1=twopi/klo; 
    if(klof >0) wl2=twopi/klof; 
    double wl3=twopi/khif;
    double wl4=twopi/khi;
    double wnx=twopi/(dk1*(nnx-1)/2);
    double wny=twopi/(dk2*(nny-1)/2);

    printf("IE BANDPASS OVER WAVELENGTHS\n");
    printf("   INF CUT-- %8.3f --TAPER-- %8.3f (PASS) %8.3f --TAPER--%8.3f\n",wl1,wl2,wl3,wl4);
    printf("   --  CUT TO NYQUIST X,Y= %8.3f  %8.3f\n",wnx,wny);
    double nnx2=nnx/2+1;
    double nny2=nny/2+1;
    std::vector<std::vector<double>> wts(k.size(), std::vector<double>(k[0].size(), 0));  // initialise to zero
    for (int i=0; i<nny;i++){
        for (int j=0; j<nnx;j++){
            if (k[i][j]>klo){ 
                if (k[i][j]<khi) wts[i][j]=1; 
            }
        }
    }
    for (int i=0;i<nny;i++){
        for (int j=0;j<nnx;j++){
            if (k[i][j]>klo){ 
                if (k[i][j]<klof) wts[i][j]=wts[i][j]*(1-cos(pi*(k[i][j]-klo)/dkl))/2;
            }
            if (k[i][j]>khif){ 
                if (k[i][j]<khi) wts[i][j]=wts[i][j]*(1-cos(pi*(khi-k[i][j])/dkh))/2;
            }
        }
    }
    return wts;
}


std::vector<double> nskew(double yr, double rlat, double rlon, double zobs, double azim, double sdec, double sdip, bool opts)
{
    // NSKEW - Compute skewness parameter and amplitude factor
    //  following Schouten (1971)
    //  Computes GEOCENTRIC DIPOLE unless given
    //  declination and dip of magnetization
    // Usage:
    //    [theta,ampfac]=nskew(yr,rlat,rlon,zobs,azim)
    // or
    //    [theta,ampfac]=nskew(yr,rlat,rlon,zobs,azim,sdec,sdip)
    //
    //  Input variables:
    //   yr : decimal year of survey
    //   rlat,rlon : regional latitude, longitude in decimal degrees
    //   zobs : level of observation in km above sealevel
    //   azim : strike of lineations normal to profile (+cw degrees from north)
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
    // get unit vectors
    double colat = 90.0 - rlat;
    std::vector<double> y = magfd(abs(yr), 1, zobs, colat, rlon);
    // compute skewness parameter
    double bx = y[0];
    double by = y[1];
    double bz = y[2];
    double bh = sqrt(pow(bx, 2) + pow(by, 2));
    double decl1 = atan2(by, bx) / rad;
    double incl1 = atan2(bz, bh) / rad;
    if (yr > 0)
    {
        printf(" EARTH' 'S MAGNETIC FIELD DIRECTION:\n");
        printf(" %10.3f = MAGNETIC DECLINATION ( STRIKE, CW FROM N )\n", decl1);
        printf(" %10.4f = MAGNETIC INCLINATION ( DIP, POS DOWN )\n", incl1);
    }
    if (opts)
    {
        //  NOTE FOR GEOCENTRIC DIPOLE TAN(INC)=2*TAN(LAT)
        //if abs(sdec) > 0. | abs(sdip) > 0.
        if (yr > 0)
        {
            printf(" NON-GEOCENTRIC MAGNETIZATION VECTOR SPECIFIED:\n");
            printf(" %10.4f = DESIRED MAGNETIZATION DECLINATION (+CW FROM N)\n", sdec);
            printf(" %10.4f = DESIRED MAGNETIZATION INCLINATION (+DN)\n", sdip);
        }
    }
    else
    {
        sdip = atan2(2. * sin(rlat * rad), cos(rlat * rad)) / rad;
        sdec = 0;
        if (yr > 0)
        {
            printf(" GEOCENTRIC MAGNETIZATION VECTOR SPECIFIED:\n");
            printf(" %10.4f = GEOCENTRIC DIPOLE INCLINATION \n", sdip);
            printf(" %10.3f = GEOCENTRIC DECLINATION ASSUMED\n", sdec);
        }
    }
    // compute phase and amplitude factors
    double ra1 = incl1 * rad;
    double rb1 = (decl1 - azim) * rad;
    double ra2 = sdip * rad;
    double rb2 = (sdec - azim) * rad;
    // compute phase and amplitude factors
    double inclm = atan2(tan(ra2), sin(rb2));
    double inclf = atan2(tan(ra1), sin(rb1));
    double ampfac = ((sin(ra2)) * (sin(ra1))) / ((sin(inclm)) * (sin(inclf)));
    double theta = (inclm / rad) + (inclf / rad) - 180.;
    if (theta <= -360)
    {
        theta = theta + 360;
    }
    else
    {
        if (theta >= 360)
        {
            theta = theta - 360;
        }
    }
    // compute unit vectors for a check
    std::vector<double> hatm;
    std::vector<double> hatb;
    hatm.push_back(cos(sdip * rad) * sin((sdec - azim) * rad));
    hatm.push_back(cos(sdip * rad) * cos((sdec - azim) * rad));
    hatm.push_back(-sin(sdip * rad));
    hatb.push_back(cos(incl1 * rad) * sin((decl1 - azim) * rad));
    hatb.push_back(cos(incl1 * rad) * cos((decl1 - azim) * rad));
    hatb.push_back(-sin(incl1 * rad));
    //
    if (yr > 0)
    {
        printf("  %10.6f %10.6f %10.6f = MAGNETIZATION UNIT VECTOR\n", hatm[0], hatm[1], hatm[2]);
        printf("  %10.6f %10.6f %10.6f = AMBIENT FIELD UNIT VECTOR\n", hatb[0], hatb[1], hatb[2]);
        printf("  COMPONENTS ARE (X,Y,Z=ALONG, ACROSS PROFILE, AND UP\n\n");
    }

    std::vector<double> a({theta, ampfac});
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
    std::vector<std::vector<double>> agh_prime;
    reads("agh",agh_prime);
    std::vector<std::vector<double>> dgh_prime;
    reads("dgh",dgh_prime);
    std::vector<double> agh = agh_prime[0];
    std::vector<double> dgh = dgh_prime[0];
    int base = 1990;
    int i = 199;
    double t = 0;
    // WARNING ENDS

    double d2r = std::atan(1.0) * 4 / 180; // pi/180
    double r = alt;
    double slat = cos(colat * d2r);
    double clat = sin(colat * d2r);
    std::vector<double> cl({cos(elong * d2r)});
    std::vector<double> sl({sin(elong * d2r)});
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

    std::vector<double> p({2.0 * slat, 2.0 * clat, 4.5 * slat * slat - 1.5, sqrt(27) * clat * slat});
    std::vector<double> q({-clat, slat, -3.0 * clat * slat, sqrt(3) * (slat * slat - clat * clat)});

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
                p.push_back((1.0 + 1.0 / fm) * one * clat * p[j]);
                q.push_back(one * (clat * q[j] + slat / fm * p[j]));
                sl.push_back(sl[m - 2] * cl[0] + cl[m - 2] * sl[0]);
                cl.push_back(cl[m - 2] * cl[0] - sl[m - 2] * sl[0]);
            }
            else
            {
                double one = sqrt(fn * fn - fm * fm);
                double two = sqrt(pow((fn - 1.0), 2) - fm * fm) / one;
                double three = (2.0 * fn - 1.0) / one;
                i = k - n;
                double j = k - 2 * n + 1;
                p.push_back((fn + 1.0) * (three * slat / fn * p[i] - two / (fn - 1.0) * p[j]));
                q.push_back(three * (slat * q[i] - clat / fn * p[i]) - two * q[j]);
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

    std::vector<double> b({x, y, z, t});
    return b;
}

// INV3D Calculate magnetization from magnetic field and bathymetry for a map.
// Assumes constant thickness source layer whose upper bound is bathymetry
// Use the Parker & Huestis [1974] Fourier inversion approach.
//
// Usage: m3d=inv3d(f3d,h,wl,ws,rlat,
//                      rlon,yr,zobs,thick,azim,dx,dy,sdec,sdip);
//   or for geocentric dipole
//     m3d=inv3d(f3d,h,wl,ws,rlat,rlon,yr,zobs,thick,azim,dx,dy);
//
// Input arrays:
//    f3d 	magnetic field (nT)
//    h 		bathymetry (km +ve up)
//    wl		filter long wavelength cutoff (km)
//    ws		filter short wavelength cutoff (km)
//    rlat 	latitude of survey area dec. deg.
//    rlon 	longitude of survey area dec. deg.
//    yr 	year of survey (dec. year)
//    azim 	azimuth of lineations (deg)
//    zobs 	observation level (+km up)
//    thick 	thickness of source layer (km)
//    azim	azimuth of grid (degrees) hard wired to 0
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
std::vector<std::vector<double>> inv3da(std::vector<std::vector<double>> f3d, std::vector<std::vector<double>> h, double wl, double ws, double rlat, double rlon, double yr, double zobs, std::vector<std::vector<double>> thick, double azim, double dx, double dy, double sdec, double sdip)
{
    //error
    std::vector<std::vector<double>> a;
    // parameters defined
    const std::complex<double> i_math(0.0, 1.0);
    const double rad = pi / 180; // conversion radians to degrees
    const int mu = 100;          // conversion factor to nT
    // changeable parameters
    int nterms = 40;
    int nitrs = 40;
    double tol = 0.0001;
    double tolmag = 0.01;
    int flag = 0;
    int xmin = 0;

    printf("       3D MAGNETIC INVERSE MODEL\n");
    printf("                  INV3DA\n");
    printf("        Variable thickness layer\n");

    printf(" Zobs= %12.5f\n Rlat= %12.5f Rlon= %12.5f\n", zobs, rlat, rlon);
    printf(" Yr= %12.5f\n", yr);
    printf(" Thick= %12.5f\n", thick);
    printf(" azim = %12.6f\n", azim);
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
    
    const double mnf3d = -1.7763568394002505e-15;
    std::for_each(f3d.begin(), f3d.end(), [mnf3d](std::vector<double> &v) { std::for_each(v.begin(), v.end(), [mnf3d](double &d) { d -= mnf3d; }); });
    printf("Remove mean of %10.3f from field \n", mnf3d);

    double colat = 90. - rlat;
    std::vector<double> y = magfd(yr, 1, zobs, colat, rlon);

    double bx = y[0];
    double by = y[1];
    double bz = y[2];
    double bh = sqrt(pow(bx, 2) + pow(by, 2));
    double decl1 = atan2(by, bx) / rad;
    double incl1 = atan2(bz, bh) / rad;
    double theta;
    double ampfac;
    std::vector<double> skew_array;

    if (abs(sdec) > 0. || abs(sdip) > 0.)
    {
        skew_array = nskew(yr, rlat, rlon, zobs, azim, sdec, sdip, true); //skew_array = [theta, ampfac]
    }
    else
    {
        skew_array = nskew(yr, rlat, rlon, zobs, azim, 0, 0, false);
        sdip = atan2(2.0 * sin(rlat * rad), cos(rlat * rad)) / rad;
        sdec = 0;
    }
    theta = skew_array[0];
    ampfac = skew_array[1];

    azim = 0; // azim is forced to zero
    double ra1 = incl1 * rad;
    double rb1 = (decl1 - azim) * rad;
    // rb1=(azim-decl1)*rad;
    double ra2 = sdip * rad;
    double rb2 = (sdec - azim) * rad;
    // rb2=(azim-sdec)*rad;

    // make wave number array
    // ni=1/nx;
    double nx2 = nx / 2;
    double nx2plus = nx2 + 1;
    // x=-.5:ni:.5-ni;
    // ni=1/ny;
    double ny2 = ny / 2;
    double ny2plus = ny2 + 1;
    double dkx = pi / (nx * dx);
    double dky = pi / (ny * dy);

    std::vector<double> kx;
    for (int i = nx2 * -1; i < nx2; i++)
    {
        kx.push_back(i * dkx);
    }
    std::vector<double> ky;
    for (int i = ny2 * -1; i < ny2; i++)
    {
        ky.push_back(i * dky);
    }
    std::vector<std::vector<double>> X(ky.size(), ky);
    std::vector<std::vector<double>> Y;
    for (auto num : kx)
    {
        std::vector<double> temp(kx.size(), num);
        Y.push_back(temp);
    }
    std::vector<std::vector<double>> k;
    for (int i = 0; i < X.size(); i++)
    {
        std::vector<double> temp;
        for (int j = 0; j < X[0].size(); j++)
        {
            temp.push_back(2 * (sqrt(pow(X[i][j], 2) + pow(Y[i][j], 2))));
        }
        k.push_back(temp);
    } // wavenumber array

    k = fftshift(k);

    std::vector<std::vector<std::complex<double>>> ob;
    std::vector<std::vector<std::complex<double>>> om;
    for(int i = 0; i<X.size();i++){
        std::vector<std::complex<double>> temp1;
        std::vector<std::complex<double>> temp2;
        for(int j = 0; j<X[i].size();j++){
            temp1.push_back(sin(ra1) + i_math * cos(ra1) * sin(atan2(Y[i][j], X[i][j]) + rb1));
            temp2.push_back(sin(ra2) + i_math * cos(ra2) * sin(atan2(Y[i][j], X[i][j]) + rb2));
        }
        ob.push_back(temp1);
        om.push_back(temp2);
    }
    std::vector<std::vector<std::complex<double>>> o;
    for(int i = 0; i<X.size();i++){
        std::vector<std::complex<double>> temp;
        for(int j = 0; j<X[i].size();j++){
            temp.push_back(ob[i][j]*om[i][j]);
        }
        o.push_back(temp);
    }
    
    o = fftshift_complex(o);
    std::vector<std::vector<std::complex<double>>> amp;
    for(auto a: o){
        std::vector<std::complex<double>> temp;
        for(auto b: a){
            temp.push_back(abs(b));
        }
        amp.push_back(temp);
    }
    // amplitude factor
    std::vector<std::vector<std::complex<double>>> phase;
    for(int i = 0; i<ob.size(); i++){
        std::vector<std::complex<double>> temp;
        for(int j = 0; j<ob[i].size(); j++){
            temp.push_back(exp(i_math * (angle(ob[i][j]) + angle(om[i][j]))));
        }
        phase.push_back(temp);
    }
    phase = fftshift_complex(phase);
    // phase angle
    double math_constant = 2 * pi * mu;
    // calculate base layer 
    std::vector<std::vector<std::complex<double>>> g;
    for(int i = 0; i < h.size(); i++){
        std::vector<std::complex<double>> temp;
        for(int j = 0; j < h[0].size(); j++){
            temp.push_back(-(abs(h[i][j])+thick[i][j]));
        }
        g.push_back(temp);
    }
    //shift zero level of bathy
    double hmax = realtwodmax(h);
    double hmin= realtwodmin(h);
    std::complex<double> gmax = twodmax(g);
    std::complex<double> gmin = twodmin(g);
    double conv=1;
    printf(" %10.3f %10.3f = MIN, MAX OBSERVED BATHY\n",hmin,hmax);
    printf("CONVERT BATHY (M OR KM), +DOWN or +UP)\n");
    printf("TO BATHY (KM, +UP)\n");
    double shift=hmax;
    double hwiggl=abs(hmax-gmin)/2;
    double zup=zobs-shift;
    printf(" SHIFT ZERO OF BATHY WILL BE %8.3f\n",shift);
    printf("THIS IS OPTIMUM FOR INVERSION.\n");
    printf("NOTE OBSERVATIONS ARE %8.3f KM ABOVE BATHY\n",zup);
    printf("ZOBS=%8.3f ZUP=%8.3f\n",zobs,zup);

    printf("%8.3f = HWIGGL, DISTANCE TO MID-LINE OF BATHY\n",hwiggl);
    printf("THIS IS OPTIMUM ZERO LEVEL FOR FORWARD PROBLEM\n");

    // bathy zero placed halfway between extremes
    // this is optimum for summation but not for iteration
    // which needs zero at highest point of bathy
    for(auto &i: h){
        for(auto &j: i){
            j = j - shift + hwiggl;
        }
    }
    for(auto &vec: g){
        for(auto &val: vec){
            val=val-shift+hwiggl;
        }
    }
    // set up bandpass filter
    std::vector<std::vector<double>> wts=bpass3d(nx,ny,dx,dy,wl,ws);
    // do eterm
    std::vector<std::vector<double>> dexpz;
    for(auto a: k){
        std::vector<double> temp;
        for(auto b: a){
            temp.push_back(exp(b*zup));
        }
        dexpz.push_back(temp);
    }
    std::vector<std::vector<double>> dexpw;
    for(auto a: k){
        std::vector<double> temp;
        for(auto b: a){
            temp.push_back(exp(-b*hwiggl));
        }
        dexpw.push_back(temp);
    }

    // take fft of observed magnetic field and initial m3d
    std::vector<std::vector<std::complex<double>>> m3d(ny, std::vector<std::complex<double>>(nx,0)); // make an initial guess of 0 for m3d
    std::vector<std::vector<std::complex<double>>> f3d_2;
    for(auto a: f3d){
        std::vector<std::complex<double>> temp;
        for(auto b: a){
                temp.push_back(b);
        }
        f3d_2.push_back(temp);
    }
    std::vector<std::vector<std::complex<double>>> F(ny, std::vector<std::complex<double>>(nx,0));
    fft2(f3d_2, F);
    std::vector<std::vector<std::complex<double>>> HG(h.size(), std::vector<std::complex<double>>(h[0].size(), 0));
    for(int i = 0; i < h.size(); i++){
        for(int j = 0; j < h[0].size(); j++){
            HG[i][j] = h[i][j] - g[i][j];
        }
    }

    int intsum=0;
    std::vector<std::vector<std::complex<double>>> mlast(ny, std::vector<std::complex<double>>(nx,0));
    std::vector<std::vector<std::complex<double>>> lastm3d(ny, std::vector<std::complex<double>>(nx,0));
    std::vector<std::vector<std::complex<double>>> B;
    
    for(int i =0; i< F.size(); i++){
        std::vector<std::complex<double>> temp;
        for(int j = 0; j< F[i].size(); j++){
            temp.push_back((F[i][j]*dexpz[i][j])/(math_constant*amp[i][j]*phase[i][j]));

        }
        B.push_back(temp);
    }
 
    B[0][0]=0;
    //
    printf(" CONVERGENCE :\n");
    printf(" ITER  MAX_PERTURB  #_TERMS  AVG ERR  \n");
    int nkount;
    double first1;
    double erpast;
    double errmax;

    for (int iter = 0; iter<nitrs; iter++){
        // summation loop, start with n = 2
        std::vector<std::vector<std::complex<double>>> sum(ny, std::vector<std::complex<double>>(nx,0));
        for (nkount = 2;nkount < nterms + 1; nkount++){
            int n=nkount;
            std::vector<std::vector<std::complex<double>>> MH(m3d.size(), std::vector<std::complex<double>>(m3d[0].size(), 0));
            std::vector<std::vector<std::complex<double>>> m3d2;
            for(int i = 0; i<m3d.size(); i++){
                std::vector<std::complex<double>> temp;
                for(int j = 0; j < m3d[i].size(); j++){
                    temp.push_back(m3d[i][j]*(pow(h[i][j],n)-pow(g[i][j], n)));
                }
                m3d2.push_back(temp);
            }
            fft2(m3d2, MH);
            for(int i = 0; i<sum.size(); i++){
                for(int j=0; j<sum[i].size(); j++){
                    sum[i][j] += dexpw[i][j]*(pow(k[i][j],n-1)/nfac(n))*MH[i][j];
                }
            }
            errmax = abs(twodmax(sum));
        }
        // transform to get new solution
        std::vector<std::vector<std::complex<double>>> M;
        for(int i = 0; i < B.size(); i++){
            std::vector<std::complex<double>> temp;
            for(int j = 0; j< B[i].size(); j++){
                temp.push_back(B[i][j]-(sum[i][j]));
            }
            M.push_back(temp);
        }
        // filter before transforming to ensure no blow ups
        M[0][0]=0;
        for(int i = 0; i < B.size(); i++){
            for(int j = 0; j< B[i].size(); j++){
                mlast[i][j]=(M[i][j]/thick[i][j])*wts[i][j];
            }
        }
        // writes(M, "M");
        ifft2(mlast, m3d);
        // writes(m3d, "m3d");
        // do convergence test
        for(auto &vec: m3d){
            for(auto &val: vec){
                val = val/static_cast<double>(m3d.size()*m3d[0].size());
            }
        }
        errmax=0;
        std::vector<std::vector<std::complex<double>>> s1(ny, std::vector<std::complex<double>>(nx,0));
        std::vector<std::vector<std::complex<double>>> s2(ny, std::vector<std::complex<double>>(nx,0));
        std::vector<std::vector<std::complex<double>>> dif;
        for(int i = 0; i < B.size(); i++){
            std::vector<std::complex<double>> temp;
            for(int j = 0; j< B[i].size(); j++){
                temp.push_back(abs(lastm3d[i][j]-m3d[i][j]));
            }
            dif.push_back(temp);
        }
        for(int i = 0; i < B.size(); i++){
            std::vector<std::complex<double>> temp;
            for(int j = 0; j< B[i].size(); j++){
                s2[i][j]=s2[i][j]+dif[i][j]*dif[i][j];
            }
        }
        std::complex<double> difmax = dif[0][0].imag() + dif[0][0].real();
        std::complex<double> difsum = 0;
        for(auto a: dif){
            for(auto b: a){
                difsum+=b;
                if(b.real() > difmax.real()) difmax = b;
            }
        }
        if (errmax-difmax.real() < 0) errmax=difmax.real();
        
        lastm3d=m3d;
        
        std::complex<double> avg=difsum/(dif.size()*dif[0].size()*1.0);
        //  rms=sqrt(s2/(nx*ny) - avg^2);
        if (iter==0){ 
            first1=errmax+1e-10; 
            erpast=errmax;
        }
        if (errmax > erpast){
            flag=1;  // set the flag to show diverging solution
            break;
        }
        erpast=errmax;
        // test for errmax less than tolerance
        if (errmax < tolmag){
            flag=0;
            break;
        }
        printf("%3.0i, %10.4e, %6.0i ",iter,errmax,nkount-1);
        printf(" %10.4e",avg);
        std::cout<<std::endl;
    }  // end of iteration loop

    if (flag == 1){ 
        printf(" I would be quitting now error < tolerance ");
    }
    else{
        printf(" RESTORE ORIGINAL ZERO LEVEL\n");
        printf(" SHIFT ZERO LEVEL OF BATHY BY %8.3f\n",shift);
        for(auto &a: h){
            for(auto &b: a){
                b = b+shift-hwiggl;
            }
        }
        for(auto &vec: g){
            for(auto &val: vec){
                val = val + shift - hwiggl;
            }
        }
        std::vector<std::vector<double>> greal;
        for(auto vec: g){
            std::vector<double> temp;
            for(auto val: vec){
                temp.push_back(val.real());
            }
            greal.push_back(temp);
        }
    }
    //
    std::vector<std::vector<double>> returns;
    writes(m3d, "final");
    for(auto a: m3d){
        std::vector<double> temp;
        for(auto b: a){
            temp.push_back(b.real());
        }
        returns.push_back(temp);
    }
    return returns;
}

int main()
{
    std::vector<std::vector<double>> f3d;
    reads("f3d", f3d);
    std::vector<std::vector<double>> h;
    reads("h", h);
    std::vector<std::vector<double>> other;
    reads("other", other);
    std::vector<std::vector<double>> thick;
    reads("thick", thick);


    double wl = other[0][0];
    double ws = other[0][1];
    double rlat = other[0][2];
    double rlon = other[0][3];
    double yr = other[0][4];
    double zobs = other[0][5];
    //TODO: thick is array in inv3da
    double azim = other[0][6];
    double dx = other[0][7];
    double dy = other[0][8];


    // Optional values, default assumes geocentric dipole hypothesis
    double sdec = 0;
    double sdip = 0;
    inv3da(f3d, h, wl, ws, rlat, rlon, yr, zobs, thick, azim, dx, dy, sdec, sdip);
    return 0;
}
