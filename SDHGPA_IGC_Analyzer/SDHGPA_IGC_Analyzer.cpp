//
//  main.cpp
//  SDHGPA_IGC_Analyzer
//
//  Created by Jeffrey H. Brown on 09/26/2017.
//  Copyright (c) 2017 Jeffrey H. Brown. All rights reserved.
//

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <unistd.h>
#include <vector>

// Geodesic library for WGS84
#include <GeographicLib/Geodesic.hpp>

// Bump this whenever a substantial change is made
#define PROGRAM_NAME   "SDHGPA_IGC_Analyzer"
#define VERSION_STRING "1.0.4"


inline double max(double a, double b) {return a > b ? a : b;}
inline double min(double a, double b) {return a < b ? a : b;}

constexpr double cylrad = 400.0;    // radius of cylinder around start point in meters

constexpr double meterspermile = 1609.344L;
constexpr double rptmin = 3.141592653589793238462643383279502884195L / 10800000.0L; // radians per thousandth of a minute

// The FAI Sphere distance formula is based on an earth radius of exactly 6371 km.
constexpr double metersperradian2 = 6371000.0 * 2.0;

typedef int_fast32_t  pt_index_t;   // type for indexing track points
typedef int_least32_t lat_t;        // we store latitude and longitude in thousandths of minutes
typedef int_least32_t lon_t;
typedef std::pair<lat_t, lon_t> latlon;   // latitude and longitude

static bool FAIsphere = false;
static bool WGS84 = true;
static bool KMLfile = false;
static bool verbose = false;
static bool useWGS84 = true;
std::ofstream kmlstream;

using namespace GeographicLib;
const Geodesic& geod = Geodesic::WGS84();

// This will return the distance between two points in meters given lat/lon in thousandths of minutes
double Distance(const latlon& p1, const latlon& p2)
{
    if (useWGS84)
    {
        // Use the Geodesic library to get the WGS84 distance
        // Input is in thousandths of minutes, so divide by 60000 to get degrees
        double s12;
        geod.Inverse(p1.first/60000.0, p1.second/60000.0, p2.first/60000.0, p2.second/60000.0, s12);
        return s12;
    }
    else
    {
        // Use the small angle formula for great circle distance on a sphere
        double lat1 = p1.first*rptmin;
        double lon1 = p1.second*rptmin;
        double lat2 = p2.first*rptmin;
        double lon2 = p2.second*rptmin;
        double coslat1 = std::cos(lat1);
        double coslat2 = std::cos(lat2);
        double sinlatd = std::sin((lat1 - lat2)/2.0);
        double sinlond = std::sin((lon1 - lon2)/2.0);
        double d = metersperradian2 * std::asin(std::sqrt(coslat1*coslat2*sinlond*sinlond + sinlatd*sinlatd));
        return d;
    }
}

// This is just a nice function for writing out lat/lon in a nice degrees and decimal minutes format
static std::ostream& operator<<(std::ostream& op, const latlon& ll)
{
    int latdeg = ll.first/60000;        // keep this signed
    int londeg = ll.second/60000;       // keep this signed
    lat_t lat = std::abs(ll.first);
    lon_t lon = std::abs(ll.second);
    int latthous = lat % 1000;
    int lonthous = lon % 1000;
    int latmin = (lat/1000) % 60;
    int lonmin = (lon/1000) % 60;
    
	char	buf[60];
	sprintf(buf, "%2d,%2d.%03d   %3d,%2d.%03d", latdeg, latmin, latthous, londeg, lonmin, lonthous);
	op << buf;
	return op;
}

// This is the flight geometry analyzer
static void AnalyzePath(const std::vector<latlon>& pts)
{
    // Build all point-point distances first, but while we're at it, do max straight line distance
    size_t npts = pts.size();
    double* dist = new double [npts*npts];
#define DIST(i,j) dist[(i)*npts+(j)]
    // Straight Line
    double  maxDSL = 0.0;
	pt_index_t	imaxSL = 0;
	pt_index_t	jmaxSL = 0;
    for (pt_index_t i=0; i<npts; i++)
    {
        DIST(i,i) = 0.0;
        for (pt_index_t j=i+1; j<npts; j++)
        {
            double d1 = DIST(i,j) = Distance(pts[i], pts[j]);
            if (d1 > maxDSL)
            {
                maxDSL = d1;
                imaxSL = i;
                jmaxSL = j;
            }
        }
    }

    // Find the max distance between any two consecutive points, then double it
    double  max2con2 = 0.0;
    for (pt_index_t i = 1; i<npts; i++)
    {
        double d = DIST(i-1,i);
        if (max2con2 < d)
            max2con2 = d;
    }
    max2con2 *= 2.0;
    
    // Build list of best first leg for any given first turnpoint, and
    // list of best last leg for any given last turnpoint
    pt_index_t* bestF = new pt_index_t [npts];
    pt_index_t* bestL = new pt_index_t [npts];
    for (pt_index_t i=0; i<npts; i++)
    {
        double max = 0.0;
        pt_index_t jmax = i;
        for (pt_index_t j=0; j<i; j++)
        {
            double d1 = DIST(j,i);
            if (d1 > max)
            {
                max = d1;
                jmax = j;
            }
        }
        bestF[i] = jmax;

        max = 0.0;
        jmax = i;
        for (pt_index_t j=i+1; j<npts; j++)
        {
            double d1 = DIST(i,j);
            if (d1 > max)
            {
                max = d1;
                jmax = j;
            }
        }
        bestL[i] = jmax;
    }

    // Out & Return
    double  maxDOR = 0.0;
	size_t	imaxOR = 0;
	size_t	jmaxOR = 0;
	size_t	lmaxOR = 0;
    double  penaltyOR = 0.0;
    // FAI Triangle
    double	maxDT = 0.0;
    double	maxDT28 = 0.0;
	size_t	imaxT = 0;
	size_t	jmaxT = 0;
	size_t	kmaxT = 0;
	size_t	lmaxT = 0;
    double  penaltyT = 0.0;
    // 4 segment (like Leonardo)
    double	maxD4 = maxDSL; // maxD4 will be >= maxDSL, so might as well start out that way
    maxD4 = max(maxD4,(DIST(0,npts/3) + DIST(npts/3,2*npts/3) + DIST(2*npts/3,npts-1)));
    maxD4 = max(maxD4,(DIST(0,npts/4) + DIST(npts/4,npts/2) + DIST(npts/2,3*npts/4) + DIST(3*npts/4,npts-1)));
	size_t	imax4 = 0;
	size_t	jmax4 = 0;
	size_t	kmax4 = 0;
	size_t	lmax4 = 0;
	size_t	mmax4 = 0;
	size_t	i, j, k, l;
    long nTimesInnerLoopOR = 0;
    long nTimesInnerLoopT = 0;
    long nTimesInnerLoop4 = 0;
    long nskip = 0;
    long nskipboth = 0;
    
    for (i=0; i<npts-2; i++)
    {
        // 4 segment
        // Divide into sum of two segments before i, two segments after i
        double bestf = DIST(bestF[i],i);
        double bestl = DIST(i,bestL[i]);
        double b = bestf + min(2.0*bestf,maxDSL);
        double e = bestl + min(2.0*bestl,maxDSL);
        if (b + e < maxD4)
            nskipboth++;
        else
        {
            // We do the after part first, because it has better cache performance
            double max4a = 0.0;
            size_t jmaxa = 0;
            size_t kmaxa = 0;
            for (j=i+1; j<npts-1; j++)
            {
                nTimesInnerLoop4++;
                double d3 = DIST(i,j);
                k = bestL[j];
                double d4 = DIST(j,k);
                if ((d3 + d4) > max4a)
                {
                    max4a = d3 + d4;
                    jmaxa = j;
                    kmaxa = k;
                }
            }
            // If the after 2-segment max is so short that there is no way to exceed maxD4, then we can skip the before part
            if ((b+max4a) < maxD4)
            {
                nskip++;
            }
            else
            {
                double max4b = 0.0;
                size_t jmaxb = 0;
                size_t kmaxb = 0;
                for (j=1; j<i; j++)
                {
                    nTimesInnerLoop4++;
                    double d1 = DIST(j,i);
                    k = bestF[j];
                    double d2 = DIST(k,j);
                    if ((d1 + d2) > max4b)
                    {
                        max4b = d1 + d2;
                        jmaxb = j;
                        kmaxb = k;
                    }
                }
                if ((max4b + max4a) > maxD4)
                {
                    maxD4 = max4b + max4a;
                    imax4 = kmaxb;
                    jmax4 = jmaxb;
                    kmax4 = i;
                    lmax4 = jmaxa;
                    mmax4 = kmaxa;
                }
            }
        }
        
        // Out & Return
        double maxpen = 0.0;        // the max penalized distance for this i
        double maxunpen = 0.0;      // the max unpenalized distance for this i
        for (l=npts-1; l>i+1; l--)  // Going in reverse is an optimization
        {
            // See if we returned close enough
            double gap = DIST(i,l);
            if (gap <= cylrad)   // 400 meters
            {
                for (j=i+1; j<l-1; j++)
                {
                    nTimesInnerLoopOR++;
                    double d1 = DIST(i,j);
                    if (maxunpen < 2.0*d1)
                        maxunpen = 2.0*d1;
                    double d2 = DIST(j,l);
                    double dOR = d1 + min(d1,d2);
                    if (maxpen < dOR)
                        maxpen = dOR;
                    if (maxDOR < dOR)
                    {
                        maxDOR = dOR;
                        imaxOR = i;
                        jmaxOR = j;
                        lmaxOR = l;
                        penaltyOR = d1 - min(d1,d2);
                    }
                }
                if ((maxpen == maxunpen) || (maxunpen<=maxDOR))
                    break;
            }
        }
        
        // FAI Triangle
        maxpen = 0.0;
        maxunpen = 0.0;
        for (l=npts-1; l>i+2; l--)  // Going in reverse is an optimization
        {
            // See if we returned close enough
            double gap = DIST(i,l);
            if (gap <= cylrad)   // 400 meters
            {
//                std::cout << "i, l: " << i << ", " << l << std::endl;
                for (j=i+1; j<l-1; j++)
                {
                    double d1 = DIST(i,j);
                    
                    if (d1 < maxDT28) continue;             //it won't be >= 28% of a larger triangle
                    for (k=j+1; k<l; k++)
                    {
                        nTimesInnerLoopT++;
                        
                        double d2 = DIST(j,k);
                        double d3 = DIST(i,k);              // This is the "simple" third leg
                        double dunpen = d1 + d2 + d3;
                        double d28 = 0.28L*dunpen;
                        if (d1 < d28) continue;
                        if (d2 < d28) continue;
                        if (d3 < d28) continue;
                        if (maxunpen < dunpen)
                            maxunpen = dunpen;
                        double dp = min(d3,DIST(k,l));      // This is the "penalized" third leg
                        double dt = d1 + d2 + dp;
                        if (maxpen < dt)
                            maxpen = dt;
                        if (dt > maxDT)
                        {
                            maxDT = dt;
                            maxDT28 = d28;
                            imaxT = i;
                            jmaxT = j;
                            kmaxT = k;
                            lmaxT = l;
                            penaltyT = d3 - dp;
                        }
                        double kinc = (maxDT-dt)/max2con2;
                        if (kinc > 2.0)
                            k += (size_t)kinc - 1;
                    }
                }
                if ((maxpen == maxunpen) || (maxunpen<=maxDT))
                    break;
            }
        }
    }
    if (penaltyOR == 0.0)
    {
        // Find the earliest final point that also has zero penalty
        double d1 = DIST(imaxOR,jmaxOR);
        for (l = jmaxOR+1; l<lmaxOR; l++)
        {
            if ((DIST(jmaxOR,l) >= d1) && (DIST(imaxOR,l) <= cylrad))
            {
                lmaxOR = l;
                break;
            }
        }
    }
    if (penaltyT == 0.0)
    {
        // Find the earliest final point that also has zero penalty
        double d3 = DIST(imaxT,kmaxT);
        for (l = kmaxT+1; l<lmaxT; l++)
        {
            if ((DIST(kmaxT,l) >= d3) && (DIST(imaxT,l) <= cylrad))
            {
                lmaxT = l;
                break;
            }
        }
    }
    std::cout << std::fixed;
    pt_index_t jfl = bestL[0];
    double maxfl = DIST(0,jfl);
    std::cout << std::setprecision(4);
    std::cout << "Max from launch:  " << maxfl/1000.0 << " km, " << maxfl/meterspermile << " mi" << std::endl;
    std::cout << "  A       " << pts[0] << std::endl;
    std::cout << "  F       " << pts[jfl] << std::endl;
    if (KMLfile)
    {
        kmlstream <<
        "      <Placemark>\n"
        "        <name>Max From Launch</name>\n"
        "        <description>" << std::setprecision(4) << maxfl/1000.0 << " km, " << maxfl/meterspermile << " mi ";
        if (useWGS84)
            kmlstream << "WGS84";
        else
            kmlstream << "FAI Sphere";
        kmlstream << "</description>\n"
        "        <styleUrl>#MaxFromLaunch</styleUrl>\n"
        "        <LineString>\n"
        "          <tessellate>1</tessellate>\n"
        "          <altitudeMode>clampToGround</altitudeMode>\n"
        "          <coordinates>\n" << std::setprecision(6) <<
        pts[0].second/60000.0 << "," << pts[0].first/60000.0 << "\n" <<
        pts[jfl].second/60000.0 << "," << pts[jfl].first/60000.0 << "\n" <<
        "          </coordinates>\n"
        "        </LineString>\n"
        "      </Placemark>\n";
    }
    std::cout << std::setprecision(4);
    std::cout << "Straight Line:  " << maxDSL/1000.0 << " km, " << maxDSL/meterspermile << " mi" << std::endl;
    std::cout << "  A       " << pts[imaxSL] << std::endl;
    std::cout << "  F       " << pts[jmaxSL] << std::endl;
    if (KMLfile)
    {
        kmlstream <<
        "      <Folder>\n"
        "        <name>Straight Line</name>\n"
        "        <description>" << std::setprecision(4) << maxDSL/1000.0 << " km, " << maxDSL/meterspermile << " mi</description>\n"
        "        <open>1</open>\n"
        "        <Placemark>\n"
        "          <name>Start</name>\n"
        "          <visibility>0</visibility>\n"
        "          <Point>\n"
        "            <coordinates>\n" << std::setprecision(6) <<
        pts[imaxSL].second/60000.0 << "," << pts[imaxSL].first/60000.0 << "\n" <<
        "            </coordinates>\n"
        "          </Point>\n"
        "        </Placemark>\n"
        "        <Placemark>\n"
        "          <styleUrl>#StraightLine</styleUrl>\n"
        "          <LineString>\n"
        "            <tessellate>1</tessellate>\n"
        "            <altitudeMode>clampToGround</altitudeMode>\n"
        "            <coordinates>\n" << std::setprecision(6) <<
        pts[imaxSL].second/60000.0 << "," << pts[imaxSL].first/60000.0 << "\n" <<
        pts[jmaxSL].second/60000.0 << "," << pts[jmaxSL].first/60000.0 << "\n" <<
        "            </coordinates>\n"
        "          </LineString>\n"
        "        </Placemark>\n"
        "      </Folder>\n";
    }
    std::cout << std::setprecision(4);
    std::cout << "Out & Return:   " << maxDOR/1000.0 << " km, " << maxDOR/meterspermile << " mi" << ", penalty " << std::setprecision(1) << penaltyOR << " m" << std::endl;
    std::cout << "  A       " << pts[imaxOR] << std::endl;
    std::cout << "  B       " << pts[jmaxOR] << std::endl;
    std::cout << "  F       " << pts[lmaxOR] << std::endl;
    if (KMLfile)
    {
        kmlstream <<
        "      <Folder>\n"
        "        <name>Out And Return</name>\n"
        "        <description>" << std::setprecision(4) << maxDOR/1000.0 << " km, " << maxDOR/meterspermile << " mi<br/>penalty " << std::setprecision(1) << penaltyOR << " m</description>\n"
        "        <open>1</open>\n"
        "        <Placemark>\n"
        "          <name>Start</name>\n"
        "          <visibility>0</visibility>\n"
        "          <Point>\n"
        "            <coordinates>\n" << std::setprecision(6) <<
        pts[imaxOR].second/60000.0 << "," << pts[imaxOR].first/60000.0 << "\n" <<
        "            </coordinates>\n"
        "          </Point>\n"
        "        </Placemark>\n"
        "        <Placemark>\n"
        "          <styleUrl>#OutAndReturn</styleUrl>\n"
        "          <LineString>\n"
        "            <tessellate>1</tessellate>\n"
        "            <altitudeMode>clampToGround</altitudeMode>\n"
        "            <coordinates>\n" <<
        pts[imaxOR].second/60000.0 << "," << pts[imaxOR].first/60000.0 << "\n" <<
        pts[jmaxOR].second/60000.0 << "," << pts[jmaxOR].first/60000.0 << "\n" <<
        pts[lmaxOR].second/60000.0 << "," << pts[lmaxOR].first/60000.0 << "\n" <<
        "            </coordinates>\n"
        "          </LineString>\n"
        "        </Placemark>\n"
        "      </Folder>\n";
    }
    std::cout << std::setprecision(4);
    std::cout << "FAI Triangle:   " << maxDT/1000.0 << " km, " << maxDT/meterspermile << " mi" << ", penalty " << std::setprecision(1) << penaltyT << " m" << std::endl;
    std::cout << "  A       " << pts[imaxT] << std::endl;
    std::cout << "  B       " << pts[jmaxT] << std::endl;
    std::cout << "  C       " << pts[kmaxT] << std::endl;
    std::cout << "  F       " << pts[lmaxT] << std::endl;
    double p1 = DIST(imaxT,jmaxT)*100.0/maxDT;
    double p2 = DIST(jmaxT,kmaxT)*100.0/maxDT;
    double p3 = DIST(imaxT,kmaxT)*100.0/maxDT;
    std::cout << std::setprecision(4);
    std::cout << "  %s      " << p1 << "%, " << p2 << "%, " << p3 << "%" << std::endl;
    if (KMLfile)
    {
        kmlstream <<
        "      <Folder>\n"
        "        <name>FAI Triangle</name>\n"
        "        <description>" << std::setprecision(4) << maxDT/1000.0 << " km, " << maxDT/meterspermile << " mi<br/>penalty " << std::setprecision(1) << penaltyT << " m<br/>" << std::setprecision(4) << p1 << "%, " << p2 << "%, " << p3 << "%</description>\n"
        "        <open>1</open>\n"
        "        <Placemark>\n"
        "          <name>Start</name>\n"
        "          <visibility>0</visibility>\n"
        "          <Point>\n"
        "            <coordinates>\n" << std::setprecision(6) <<
        pts[imaxT].second/60000.0 << "," << pts[imaxT].first/60000.0 << "\n" <<
        "            </coordinates>\n"
        "          </Point>\n"
        "        </Placemark>\n"
        "        <Placemark>\n"
        "          <styleUrl>#FAITriangle</styleUrl>\n"
        "          <LineString>\n"
        "            <tessellate>1</tessellate>\n"
        "            <altitudeMode>clampToGround</altitudeMode>\n"
        "            <coordinates>\n" << std::setprecision(6) <<
        pts[kmaxT].second/60000.0 << "," << pts[kmaxT].first/60000.0 << "\n" <<
        pts[imaxT].second/60000.0 << "," << pts[imaxT].first/60000.0 << "\n" <<
        pts[jmaxT].second/60000.0 << "," << pts[jmaxT].first/60000.0 << "\n" <<
        pts[kmaxT].second/60000.0 << "," << pts[kmaxT].first/60000.0 << "\n" <<
        pts[lmaxT].second/60000.0 << "," << pts[lmaxT].first/60000.0 << "\n" <<
        "            </coordinates>\n"
        "          </LineString>\n"
        "        </Placemark>\n"
        "      </Folder>\n";
    }
    std::cout << std::setprecision(4);
    std::cout << "4 segment:     " << maxD4/1000.0 << " km, " << maxD4/meterspermile << " mi" << std::endl;
    std::cout << "  A       " << pts[imax4] << std::endl;
    std::cout << "  B       " << pts[jmax4] << std::endl;
    std::cout << "  C       " << pts[kmax4] << std::endl;
    std::cout << "  D       " << pts[lmax4] << std::endl;
    std::cout << "  F       " << pts[mmax4] << std::endl;
    if (KMLfile)
    {
        kmlstream <<
        "      <Folder>\n"
        "        <name>Four Segment</name>\n"
        "        <description>" << std::setprecision(4) << maxD4/1000.0 << " km, " << maxD4/meterspermile << " mi</description>\n"
        "        <open>1</open>\n"
        "        <Placemark>\n"
        "          <name>Start</name>\n"
        "          <visibility>0</visibility>\n"
        "          <Point>\n"
        "            <coordinates>\n" << std::setprecision(6) <<
        pts[imax4].second/60000.0 << "," << pts[imax4].first/60000.0 << "\n" <<
        "            </coordinates>\n"
        "          </Point>\n"
        "        </Placemark>\n"
        "        <Placemark>\n"
        "          <styleUrl>#FourSegment</styleUrl>\n"
        "          <LineString>\n"
        "            <tessellate>1</tessellate>\n"
        "            <altitudeMode>clampToGround</altitudeMode>\n"
        "            <coordinates>\n" << std::setprecision(6) <<
        pts[imax4].second/60000.0 << "," << pts[imax4].first/60000.0 << "\n" <<
        pts[jmax4].second/60000.0 << "," << pts[jmax4].first/60000.0 << "\n" <<
        pts[kmax4].second/60000.0 << "," << pts[kmax4].first/60000.0 << "\n" <<
        pts[lmax4].second/60000.0 << "," << pts[lmax4].first/60000.0 << "\n" <<
        pts[mmax4].second/60000.0 << "," << pts[mmax4].first/60000.0 << "\n" <<
        "            </coordinates>\n"
        "          </LineString>\n"
        "        </Placemark>\n"
        "      </Folder>\n";
    }
    if (verbose)
    {
        std::cout << "n Times Through Inner Loop OR: " << nTimesInnerLoopOR << std::endl;
        std::cout << "n Times Through Inner Loop T : " << nTimesInnerLoopT << std::endl;
        std::cout << "n Times Through Inner Loop 4 : " << nTimesInnerLoop4 << std::endl;
        std::cout << "n skip                       : " << nskip << std::endl;
        std::cout << "n skip both                  : " << nskipboth << std::endl;
    }

    delete[] dist;
    delete[] bestF;
    delete[] bestL;
}

void OpenKML(const std::string& filename, const std::string& date, const std::string& pilot)
{
    std::string kmlfilename(filename + ".kml");
    kmlstream.open(kmlfilename);
    kmlstream <<
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n"
    "  <Document>\n"
    "    <name>" << filename << "</name>\n"
    "    <open>1</open>\n"
    "    <description>" << date << "<br/>" << pilot << "</description>\n"
    "    <Style id=\"FlightPath\">\n"
    "      <LineStyle>\n"
    "        <color>7f0000ff</color>\n"
    "        <width>2</width>\n"
    "      </LineStyle>\n"
    "      <PolyStyle>\n"
    "        <color>3f0000ff</color>\n"
    "      </PolyStyle>\n"
    "    </Style>\n"
    "    <Style id=\"MaxFromLaunch\">\n"
    "      <LineStyle>\n"
    "        <color>ffffffff</color>\n"
    "      </LineStyle>\n"
    "    </Style>\n"
    "    <Style id=\"StraightLine\">\n"
    "      <LineStyle>\n"
    "        <color>ffffff00</color>\n"
    "      </LineStyle>\n"
    "    </Style>\n"
    "    <Style id=\"OutAndReturn\">\n"
    "      <LineStyle>\n"
    "        <color>ffff00ff</color>\n"
    "      </LineStyle>\n"
    "    </Style>\n"
    "    <Style id=\"FAITriangle\">\n"
    "      <LineStyle>\n"
    "        <color>ff00ffff</color>\n"
    "      </LineStyle>\n"
    "    </Style>\n"
    "    <Style id=\"FourSegment\">\n"
    "      <LineStyle>\n"
    "        <color>ff00ff00</color>\n"
    "      </LineStyle>\n"
    "    </Style>\n"
    "    <Folder>\n"
    "      <name>SDHGPA HG XC Contest</name>\n"
    "      <open>1</open>\n";
    kmlstream << std::fixed;
}

void StartKMLFlightPath()
{
    kmlstream <<
    "      <Placemark>\n"
    "        <name>Flight Path</name>\n"
    "        <styleUrl>#FlightPath</styleUrl>\n"
    "        <LineString>\n"
    "          <extrude>1</extrude>\n"
    "          <tessellate>1</tessellate>\n"
    "          <altitudeMode>absolute</altitudeMode>\n"
    "          <coordinates>\n";
}

void FinishKMLFlightPath()
{
    kmlstream <<
    "          </coordinates>\n"
    "        </LineString>\n"
    "      </Placemark>\n";
}

void CloseKML()
{
    kmlstream <<
    "    </Folder>\n"
    "  </Document>\n"
    "</kml>" << std::endl;
    kmlstream.close();
}

void AnalyzeIGCFile(const std::string& filename)
{
    std::cout << filename << std::endl;
    std::vector<latlon> latlons;
    
    std::ifstream source;
    source.open(filename);
    std::string line;
    latlon lastpt;
    int lastalt(INT_MIN);
    bool KMLB = false;
    std::string pilot;
    std::stringstream date;
    while (std::getline(source, line))
    {
//        std::cout << '\n' << line << std::endl;

        // Handle the date if found
        if (line.find("HFDTE")==0)
        {
            int month = 0;
            int day = 0;
            int year = 0;
            if (std::sscanf(line.c_str(),"HFDT%*6[DATE:]%2d%2d%2d%*3s", &day, &month, &year) == 3)
            {
                date << "UTC Date: "
                << std::setw(2) << std::setfill('0') << month << '/'
                << std::setw(2) << std::setfill('0') << day << '/'
                << std::setw(2) << std::setfill('0') << year;
                std::cout << date.str() << '\n';
            }
            else
                std::cerr << "Date record found but could not be parsed: " << line << std::endl;
        }
        
        // Handle the pilot if found
        if (line.find("HFPLT")==0)
        {
            pilot = line.c_str()+5;
            std::cout << pilot << std::endl;
        }
        
        // Handle the GPS Datum if found
        if (line.find("HFDTM")==0)
        {
            std::cout << "GPS Datum" << line.c_str()+line.find(":") << std::endl;
        }
        
        // Handle the B records
        // TBD Need to trim off pre-launch and post-landing points
        if (line.find("B")==0)
        {
            int timestamp = 0;
            int latdeg = 0;
            int latmin = 0;
            char latNS = '\0';
            int londeg = 0;
            int lonmin = 0;
            char lonEW = '\0';
            int alt = 0;
            if (std::sscanf(line.c_str(),"B%6d%2d%5d%c%3d%5d%c%*c%5d", &timestamp, &latdeg, &latmin, &latNS, &londeg, &lonmin, &lonEW, &alt) == 8)
            {
                if (alt == 0)
                {
                    // Some IGCs record GPS altitude but not baro altitude, so see if we can use GPS altitude
                    int gpsalt = 0;
                    if (std::sscanf(line.c_str(),"B%*6d%*2d%*5d%*c%*3d%*5d%*c%*c%*5d%5d", &gpsalt) == 1)
                    {
                        alt = gpsalt;
                    }
                }
                lat_t lattmin = latdeg*60000 + latmin;     // latitude in thousandths of minutes
                lon_t lontmin = londeg*60000 + lonmin;     // longitude in thousandths of minutes
                if (latNS == 'S') lattmin = -lattmin;
                if (lonEW == 'W') lontmin = -lontmin;
                // Add this point if it is different from the last (or if the list is empty)
                latlon thispt(lattmin, lontmin);
                if ((thispt != lastpt) || (alt != lastalt) || (latlons.size() == 0))
                {
                    if (KMLfile)
                    {
                        if (!KMLB)
                        {
                            OpenKML(filename, date.str(), pilot);
                            StartKMLFlightPath();
                            KMLB = true;
                        }
                        kmlstream << std::setprecision(6) << lontmin/60000.0 << "," << lattmin/60000.0 << "," << alt << '\n';
                    }
                    latlons.push_back(thispt);
                    lastpt = thispt;
                    lastalt = alt;
                }
            }
            else
                std::cerr << "B record found but could not be parsed: " << line << std::endl;
        }
    }
    if (KMLB)
        FinishKMLFlightPath();
    std::cout << latlons.size() << " B records found" << std::endl;
    std::cout << "Launch  " << latlons[0] << std::endl;
    std::cout << "Landing " << latlons[latlons.size()-1] << std::endl;
    if (FAIsphere)
    {
        useWGS84 = false;
        std::cout << "\nAnalysis using FAI Sphere" << std::endl;
        AnalyzePath(latlons);
    }
    if (WGS84)
    {
        useWGS84 = true;
        std::cout << "\nAnalysis using WGS84 ellipsoid" << std::endl;
        AnalyzePath(latlons);
    }
    if (KMLfile)
        CloseKML();
}

int main(int argc, char * argv[])
{
    int optcnt = 0;
    int opt;
    while ((opt = getopt(argc,argv,"fwkdvh")) != -1)
    {
        optcnt++;
        switch(opt)
        {
            case 'f': FAIsphere = true; break;
            case 'w': WGS84 = true; break;
            case 'k': KMLfile = true; break;
            case 'd': verbose = true; break;
            case 'v': std::cout << PROGRAM_NAME << " " << VERSION_STRING << std::endl; break;
            case 'h':
            case '?': std::cout << "usage is\n"
                "-f : Use FAI Sphere for distance calculations\n"
                "-w : Use WGS84 for distance calculations (default if FAI Sphere is not specified)\n"
                "-k : Produce KML file for each input file\n"
                "-d : print debug info (verbose)\n"
                "-v : print version\n"
                "-h : help (usage)" << std::endl; return 0;
            default: return 1;
        }
    }
    if (!FAIsphere) WGS84 = true;

//    AnalyzeIGCFile("/Users/jeff/Documents/From T43/Flights copy/2011-09-11-BRA-008-01.igc");
//    AnalyzeIGCFile("/Users/jeff/Documents/From T43/Flights/20090927GPSVAR.igc");
//    AnalyzeIGCFile("/Users/jeff/Documents/From T43/Flights copy/38V_20030831GPSVAR.igc");
//    AnalyzeIGCFile("/Users/jeff/Downloads/2017-07-27-XCS-AAA-01.igc");
//    AnalyzeIGCFile("/Users/jeff/Downloads/77RA0OV1.igc");
//    AnalyzeIGCFile("/Users/jeff/Downloads/79HA0OV1.igc");
//    AnalyzeIGCFile("/Users/jeff/Downloads/79HA0OV1-JB.igc");
//    AnalyzeIGCFile("/Users/jeff/Downloads/2017-06-08-XFH-000-01.igc");
//    AnalyzeIGCFile("/Users/jeff/devel/JHB/Ball2IGC/IGC Files From Ball2IGC/2003-07-27-XBV-000-01.IGC");
//    AnalyzeIGCFile("/Users/jeff/devel/JHB/Ball2IGC/IGC Files From Ball2IGC/2003-08-10-XBV-000-01.IGC");
//    AnalyzeIGCFile("/Users/jeff/devel/JHB/Ball2IGC/IGC Files From Ball2IGC/2003-08-31-XBV-000-01.IGC");
//    AnalyzeIGCFile("/Users/jeff/Downloads/2003-07-09-XBV-000-01.IGC");
//    AnalyzeIGCFile("/Users/jeff/Library/Developer/Xcode/DerivedData/Ball2IGC-dbgimcvfckizjhejwkjkbzdxdivt/Build/Products/Debug/1997-05-18-XBV-000-01.IGC");
    if (verbose)
    {
        for (int i=1; i< argc; i++) {
            printf("arg%d=%s\n", i, argv[i]);
        }
    }
    for (int i=optcnt+1; i< argc; i++) {
        AnalyzeIGCFile(argv[i]);
    }
    return 0;
}

