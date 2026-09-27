
//----------------------------------------------------------------
// Using Gauss-Legendre quadrature in theta direction
// and   trapezoidal rule in phi direction (from Second Euler-Maclaurin summation formula, we can see that
// this method gives expolential convergence for periodic function)
//----------------------------------------------------------------
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include <cstdlib>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <math.h>
#include <map.h>
#endif
#include <mpi.h>

#include "misc.h"
#include "cgh.h"
#include "Parallel.h"
#include "surface_integral.h"
#include "fadmquantites_bssn.h"
#include "getnpem2.h"
#include "getnp4.h"
#include "parameters.h"
#include "fmisc.h"

#define PI M_PI

namespace {
struct BHFastPoint {
  Block *b;
  int det, th, ph;
  double x[3];
  int inds[3];
  double coef[6 * ghost_width];
};

struct BHFastCache {
  bool configured, ready, valid;
  int lev, maxl, modes, ndet, ord, nth, nph, nm;
  double rmax, dr;
  cgh *gh;
  std::vector<double> radii, cphi, sphi, thetaA, thetaB;
  std::vector<int> ml, mm;
  std::vector<BHFastPoint> pts;
  std::vector<double> fourier, local, global;
  BHFastCache():configured(false),ready(false),valid(false),lev(-1),maxl(0),
    modes(0),ndet(0),ord(0),nth(0),nph(0),nm(0),rmax(0),dr(0),gh(NULL){}
};
static std::map<const surface_integral*,BHFastCache> bhfc;

static int bh_mode_count(int L) {
  int n=0; for(int l=2;l<=L;++l) n+=2*l+1; return n;
}

static bool bh_read_cfg(int &L,int &nd,double &rmax,double &dr) {
  std::map<string,string>::iterator it=parameters::str_par.find("inputpar");
  if(it==parameters::str_par.end()) return false;
  ifstream in(it->second.c_str(),ifstream::in); if(!in.good()) return false;
  bool a=false,b=false,c=false,d=false; const int LEN=256; char line[LEN];
  string str,g,k,v; int si;
  while(in.good()) {
    in.getline(line,LEN); str=line;
    int st=misc::parse_parts(str,g,k,v,si); if(st<=0||g!="ABE") continue;
    if(k=="Max mode l"){L=atoi(v.c_str());a=true;}
    else if(k=="detector number"){nd=atoi(v.c_str());b=true;}
    else if(k=="farest detector position"){rmax=atof(v.c_str());c=true;}
    else if(k=="detector distance"){dr=atof(v.c_str());d=true;}
  }
  return a&&b&&c&&d;
}

struct BHBounds { Block *b; double lo[3],hi[3]; };
static bool bh_owns(const BHBounds &bb,const double p[3],const double h[3]) {
  for(int d=0;d<3;++d)
    if(p[d]-bb.lo[d] < -0.5*h[d] || p[d]-bb.hi[d] > 0.5*h[d]) return false;
  return true;
}

static bool bh_prepare_interp(BHFastPoint &p,Block *b,const double x[3],
                              int ord,int sym) {
  if(ord<2 || ord>2*ghost_width) return false;
  double xa[2*ghost_width], ya[2*ghost_width];
  for(int q=0;q<ord;++q) xa[q]=(double)q;

  for(int d=0;d<3;++d) {
    if(!b->X[d] || b->shape[d]<ord) return false;
    const double h=b->X[d][1]-b->X[d][0];
    int ci=(int)((x[d]-b->X[d][0])/h+0.4)+1; // Fortran 1-based
    int cb=ci-ord/2+1;
    int ct=cb+ord-1;
    int cmin=1, cmax=b->shape[d];
#ifdef Cell
    if(sym==2 && d<2 && fabs(b->X[d][0])<h) cmin=-ord/2+1;
    if(sym!=0 && d==2 && fabs(b->X[d][0])<h) cmin=-ord/2+1;
#else
    if(sym==2 && d<2 && fabs(b->X[d][0])<h) cmin=-ord/2+2;
    if(sym!=0 && d==2 && fabs(b->X[d][0])<h) cmin=-ord/2+2;
#endif
    if(cb<cmin){cb=cmin;ct=cb+ord-1;}
    if(ct>cmax){ct=cmax;cb=ct+1-ord;}

    double cx;
    if(cb>0) {
      cx=(x[d]-b->X[d][cb-1])/h;
    } else {
#ifdef Cell
      cx=(x[d]+b->X[d][-cb])/h;
#else
      cx=(x[d]+b->X[d][1-cb])/h;
#endif
    }

    p.inds[d]=cb-1; // f_global_interpind expects C-style start index
    for(int basis=0;basis<ord;++basis) {
      for(int q=0;q<ord;++q) ya[q]=(q==basis)?1.0:0.0;
      double xq=cx, yq=0.0, dy=0.0;
      int oq=ord;
      f_polint(xa,ya,xq,yq,dy,oq);
      p.coef[d*ord+basis]=yq;
    }
  }
  return true;
}

static size_t bh_fidx(const BHFastCache &c,int d,int t,int m,int z) {
  return ((((size_t)d*c.nth+t)*c.nm+m)*4+z);
}

static bool bh_build(BHFastCache &c,int lev,cgh *gh,var *rp,var *ip,
                     int sym,int nth,int nph,int ntot,double *nx,double *ny,double *nz,
                     double *ct,double *qw,double dphi,int rank) {
  if(!c.configured) {
    if(!bh_read_cfg(c.maxl,c.ndet,c.rmax,c.dr)) return false;
    c.modes=bh_mode_count(c.maxl); c.nm=c.maxl+1;
    std::map<string,int>::iterator oi=parameters::int_par.find("surface interpolation order");
    c.ord=(oi==parameters::int_par.end())?2*ghost_width:oi->second;
    c.radii.resize(c.ndet);
    for(int d=0;d<c.ndet;++d) c.radii[d]=c.rmax-d*c.dr;
    c.configured=true;
  }
  if(sym!=1||!gh||lev<0||lev>=gh->levels||gh->grids[lev]!=1||
     c.ord<2||c.ord>2*ghost_width||(c.ord&1)||c.maxl<2||c.ndet<1) return false;
  Patch *pa=gh->PatL[lev]?gh->PatL[lev]->data:NULL; if(!pa) return false;
  c.lev=lev;c.gh=gh;c.nth=nth;c.nph=nph;
  c.ml.clear();c.mm.clear();
  for(int l=2;l<=c.maxl;++l) for(int m=-l;m<=l;++m){c.ml.push_back(l);c.mm.push_back(m);}
  c.cphi.resize((size_t)nph*c.nm); c.sphi.resize((size_t)nph*c.nm);
  for(int j=0;j<nph;++j){double ph=(j+0.5)*dphi;for(int m=0;m<c.nm;++m){
    c.cphi[(size_t)j*c.nm+m]=cos(m*ph); c.sphi[(size_t)j*c.nm+m]=sin(m*ph);}}
  c.thetaA.resize((size_t)c.modes*nth); c.thetaB.resize((size_t)c.modes*nth);
  for(int q=0;q<c.modes;++q){int l=c.ml[q],m=c.mm[q];double norm=sqrt((2.0*l+1)/(4.0*PI));
    for(int i=0;i<nth;++i){double tp=norm*misc::Wigner_d_function(l,m,2,ct[i]);
      double tm=norm*misc::Wigner_d_function(l,m,2,-ct[i]);
#ifdef GaussInt
      double w=qw[i]*dphi;
#else
      double w=dphi;
#endif
      c.thetaA[(size_t)q*nth+i]=(tp+tm)*w;
      c.thetaB[(size_t)q*nth+i]=(tp-tm)*w;
    }}
  double h[3]={pa->getdX(0),pa->getdX(1),pa->getdX(2)};
  std::vector<BHBounds> bs;
  for(MyList<Block>*bp=pa->blb;bp;bp=bp->next){Block *b=bp->data;BHBounds bb;bb.b=b;
    for(int d=0;d<3;++d){
#ifdef Vertex
      bb.lo[d]=feq(b->bbox[d],pa->bbox[d],h[d]/2)?b->bbox[d]+pa->lli[d]*h[d]:b->bbox[d]+(ghost_width-0.5)*h[d];
      bb.hi[d]=feq(b->bbox[3+d],pa->bbox[3+d],h[d]/2)?b->bbox[3+d]-pa->uui[d]*h[d]:b->bbox[3+d]-(ghost_width-0.5)*h[d];
#else
      bb.lo[d]=feq(b->bbox[d],pa->bbox[d],h[d]/2)?b->bbox[d]+pa->lli[d]*h[d]:b->bbox[d]+ghost_width*h[d];
      bb.hi[d]=feq(b->bbox[3+d],pa->bbox[3+d],h[d]/2)?b->bbox[3+d]-pa->uui[d]*h[d]:b->bbox[3+d]-ghost_width*h[d];
#endif
    }bs.push_back(bb);if(bp==pa->ble)break;}
  int world=1;MPI_Comm_size(MPI_COMM_WORLD,&world);
  c.pts.clear();c.pts.reserve(((size_t)ntot*c.ndet+world-1)/world*2+64);
  for(int d=0;d<c.ndet;++d){double r=c.radii[d];
    for(int n=0;n<ntot;++n){double x[3]={r*nx[n],r*ny[n],r*nz[n]};Block *own=NULL;
      for(size_t k=0;k<bs.size();++k)if(bh_owns(bs[k],x,h)){own=bs[k].b;break;}
      if(!own)return false;if(own->rank!=rank)continue;
      BHFastPoint p;p.b=own;p.det=d;p.th=n/nph;p.ph=n-p.th*nph;
      p.x[0]=x[0];p.x[1]=x[1];p.x[2]=x[2];
      if(!bh_prepare_interp(p,own,x,c.ord,sym)) return false;
      c.pts.push_back(p);
    }}
  c.fourier.assign((size_t)c.ndet*nth*c.nm*4,0.0);
  c.local.assign((size_t)c.ndet*c.modes*2,0.0);
  c.global.assign((size_t)c.ndet*c.modes*2,0.0);
  c.valid=true;return true;
}

static void bh_eval(BHFastCache &c,var *rp,var *ip) {
  std::fill(c.fourier.begin(),c.fourier.end(),0.0);
  std::fill(c.local.begin(),c.local.end(),0.0);
  for(size_t z=0;z<c.pts.size();++z){BHFastPoint &p=c.pts[z];
    double rr=0.0,ii=0.0;
    int sym=1, sst=-1;
    double xx=p.x[0],yy=p.x[1],zz=p.x[2];
    f_global_interpind(p.b->shape,p.b->X[0],p.b->X[1],p.b->X[2],
                       p.b->fgfs[rp->sgfn],rr,xx,yy,zz,c.ord,rp->SoA,sym,
                       p.inds,p.coef,sst);
    xx=p.x[0];yy=p.x[1];zz=p.x[2];
    f_global_interpind(p.b->shape,p.b->X[0],p.b->X[1],p.b->X[2],
                       p.b->fgfs[ip->sgfn],ii,xx,yy,zz,c.ord,ip->SoA,sym,
                       p.inds,p.coef,sst);
    size_t tr=(size_t)p.ph*c.nm;
    for(int m=0;m<c.nm;++m){double cs=c.cphi[tr+m],sn=c.sphi[tr+m];
      double *f=&c.fourier[bh_fidx(c,p.det,p.th,m,0)];
      f[0]+=rr*cs;f[1]+=rr*sn;f[2]+=ii*cs;f[3]+=ii*sn;}
  }
  for(int d=0;d<c.ndet;++d){double *out=&c.local[(size_t)d*c.modes*2];double r=c.radii[d];
    for(int q=0;q<c.modes;++q){int m=c.mm[q],ma=abs(m);double ss=m<0?-1.0:1.0,ar=0,ai=0;
      const double *A=&c.thetaA[(size_t)q*c.nth],*B=&c.thetaB[(size_t)q*c.nth];
      for(int t=0;t<c.nth;++t){const double *f=&c.fourier[bh_fidx(c,d,t,ma,0)];
        double rc=f[0],rs=ss*f[1],ic=f[2],is=ss*f[3];
        ar+=A[t]*rc+B[t]*is;ai+=B[t]*ic-A[t]*rs;}
      out[q]=ar*r;out[c.modes+q]=ai*r;
    }}
  MPI_Reduce(c.local.data(),c.global.data(),(int)c.local.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  c.ready=true;
}
}
//|============================================================================
//| Constructor
//|============================================================================

surface_integral::surface_integral(int iSymmetry) : Symmetry(iSymmetry)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
  int N = 40;
  // read parameter from file
  {
    const int LEN = 256;
    char pline[LEN];
    string str, sgrp, skey, sval;
    int sind;
    char pname[50];
    {
      map<string, string>::iterator iter = parameters::str_par.find("inputpar");
      if (iter != parameters::str_par.end())
      {
        strcpy(pname, (iter->second).c_str());
      }
      else
      {
        cout << "Error inputpar" << endl;
        exit(0);
      }
    }
    ifstream inf(pname, ifstream::in);
    if (!inf.good() && myrank == 0)
    {
      cout << "Can not open parameter file " << pname << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 1; inf.good(); i++)
    {
      inf.getline(pline, LEN);
      str = pline;

      int status = misc::parse_parts(str, sgrp, skey, sval, sind);
      if (status == -1)
      {
        cout << "error reading parameter file " << pname << " in line " << i << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      else if (status == 0)
        continue;

      if (sgrp == "SurfaceIntegral")
      {
        if (skey == "number of points for quarter sphere")
          N = atoi(sval.c_str());
      }
    }
    inf.close();
  }
  //|-----number of points for whole [0,pi] x [0,2pi]
  N_phi   = 4 * N;   // for simplicity, we require this number must be 4*N
  N_theta = 2 * N;   //                                                2*N

  if (myrank == 0)
  {
    cout << "-----------------------------------------------------------------------" << endl;
#ifdef GaussInt
    cout << " spherical integration for wave form extraction with Gauss method      " << endl;
#else
    cout << " spherical integration for wave form extraction with mid point method  " << endl;
#endif
    cout << " N_phi   = " << N_phi   << endl;
    cout << " N_theta = " << N_theta << endl;
    cout << "-----------------------------------------------------------------------" << endl;
  }

#ifdef GaussInt
  //  weight function cover all of [0,pi]
  arcostheta = new double[N_theta];
  wtcostheta = new double[N_theta];

  // note: theta in [0,pi/2], upper half sphere, corresponds to 1 < costheta < 0
  misc::gaulegf(-1.0, 1.0, arcostheta, wtcostheta, N_theta);
  // due to symmetry, I need first half array corresponds to upper sphere, note these two arrays must match each other
  misc::inversearray(arcostheta, N_theta);
  misc::inversearray(wtcostheta, N_theta);
#endif

  if (Symmetry == 2)
  {
    N_phi = N_phi / 4;
    N_theta = N_theta / 2;
    dphi = PI / (2.0 * N_phi);
    dcostheta = 1.0 / N_theta;
    factor = 8;
  }
  else if (Symmetry == 1)
  {
    N_theta = N_theta / 2;
    dphi = 2.0 * PI / N_phi;
    dcostheta = 1.0 / N_theta;
    factor = 2;
  }
  else if (Symmetry == 0)
  {
    dphi = 2.0 * PI / N_phi;
    dcostheta = 2.0 / N_theta;
    factor = 1;
  }
  else if (myrank == 0)
  {
    cout << "surface_integral::surface_integral: not supported Symmetry setting!" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

#ifndef GaussInt
  //  weight function cover all of [0,pi]
  arcostheta = new double[N_theta];
#endif
  n_tot = N_theta * N_phi;
  nx_g = new double[n_tot];
  ny_g = new double[n_tot];
  nz_g = new double[n_tot];

  int n = 0;
  double costheta, sintheta, ph;

  for (int i = 0; i < N_theta; ++i)
  {
#ifndef GaussInt
    arcostheta[i] = 1.0 - (i + 0.5) * dcostheta;
#endif
    costheta = arcostheta[i];
    sintheta = sqrt(1.0 - costheta * costheta);

    for (int j = 0; j < N_phi; ++j)
    {
      ph = (j + 0.5) * dphi;
      // normal vector respect to the constant R sphere
      nx_g[n] = sintheta * cos(ph);
      ny_g[n] = sintheta * sin(ph);
      nz_g[n] = costheta;
      n++;
    }
  }
}

//|============================================================================
//| Destructor
//|============================================================================
surface_integral::~surface_integral()
{
  delete[] nx_g;
  delete[] ny_g;
  delete[] nz_g;
  delete[] arcostheta;
#ifdef GaussInt
  delete[] wtcostheta;
#endif
}
void surface_integral::begin_wave_benchmark(int lev, cgh *GH, var *Rpsi4, var *Ipsi4,
                                            int field_generation)
{
  (void)field_generation;
  BHFastCache &c=bhfc[this];
  c.ready=false;
  if(!c.valid || c.gh!=GH || c.lev!=lev)
    c.valid=bh_build(c,lev,GH,Rpsi4,Ipsi4,Symmetry,N_theta,N_phi,n_tot,
                     nx_g,ny_g,nz_g,arcostheta,
#ifdef GaussInt
                     wtcostheta,
#else
                     NULL,
#endif
                     dphi,myrank);
  if(c.valid) {
    bh_eval(c,Rpsi4,Ipsi4);
    if(field_generation==1) {
      std::vector<double> fast=c.global;
      std::vector<double> rr(c.modes),ii(c.modes);
      double max_abs=0.0,max_rel=0.0; int bad_det=-1,bad_mode=-1,bad_part=-1;
      c.ready=false;
      for(int d=0; d<c.ndet; ++d) {
        surf_Wave(c.radii[d],lev,GH,Rpsi4,Ipsi4,2,c.maxl,c.modes,
                  rr.data(),ii.data(),NULL);
        if(myrank==0) {
          const double *q=&fast[(size_t)d*c.modes*2];
          for(int m=0;m<c.modes;++m) {
            double vals[2]={rr[m],ii[m]};
            double refs[2]={q[m],q[c.modes+m]};
            for(int part=0;part<2;++part) {
              double da=fabs(vals[part]-refs[part]);
              double dr=da/std::max(1.0e-300,fabs(vals[part]));
              if(da>max_abs){max_abs=da;bad_det=d;bad_mode=m;bad_part=part;}
              if(dr>max_rel)max_rel=dr;
            }
          }
        }
      }
      if(myrank==0)
        cerr<<"BH_DEBUG fast-vs-legacy max_abs="<<setprecision(17)<<max_abs
            <<" max_rel="<<max_rel<<" det="<<bad_det<<" mode="<<bad_mode
            <<" part="<<(bad_part==0?"R":"I")<<endl;
      c.ready=true;
    }
  }
}
void surface_integral::end_wave_benchmark()
{
  bhfc[this].ready=false;
}
//|----------------------------------------------------------------
//  spin weighted spinw component of psi4, general routine
//  l takes from spinw to maxl; m takes from -l to l
//|----------------------------------------------------------------
void surface_integral::surf_Wave(double rex, int lev, cgh *GH, var *Rpsi4, var *Ipsi4,
                                 int spinw, int maxl, int NN, double *RP, double *IP,
                                 monitor *Monitor) // NN is the length of RP and IP
{
  BHFastCache &fc=bhfc[this];
  if(fc.valid && fc.ready && spinw==2 && maxl==fc.maxl && NN==fc.modes &&
     lev==fc.lev && GH==fc.gh) {
    int d=-1; double tol=1.0e-12*std::max(1.0,fabs(rex));
    for(int k=0;k<fc.ndet;++k) if(fabs(rex-fc.radii[k])<=tol){d=k;break;}
    if(d>=0) {
      if(myrank==0) {
        const double *q=&fc.global[(size_t)d*fc.modes*2];
        memcpy(RP,q,sizeof(double)*NN);
        memcpy(IP,q+fc.modes,sizeof(double)*NN);
      } else {
        memset(RP,0,sizeof(double)*NN);
        memset(IP,0,sizeof(double)*NN);
      }
      return;
    }
  }
  if (myrank == 0 && GH->grids[lev] != 1)
    if (Monitor->outfile)
      Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
    else
      cout << "WARNING: surface integral on multipatches" << endl;

  const int InList = 2;

  MyList<var> *DG_List = new MyList<var>(Rpsi4);
  DG_List->insert(Ipsi4);

  int n;
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  GH->PatL[lev]->data->Interp_Points(DG_List, n_tot, pox, shellf, Symmetry);

  int mp, Lp, Nmin, Nmax;

  mp = n_tot / cpusize;
  Lp = n_tot - cpusize * mp;

  if (Lp > myrank)
  {
    Nmin = myrank * mp + myrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = myrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  //|~~~~~> Integrate the dot product of Dphi with the surface normal.

  double *RP_out, *IP_out;
  RP_out = new double[NN];
  IP_out = new double[NN];

  for (int ii = 0; ii < NN; ii++)
  {
    RP_out[ii] = 0;
    IP_out[ii] = 0;
  }
  // theta part
  double costheta, thetap;
  double cosmphi, sinmphi;

  int i, j;
  int lpsy = 0;
  if (Symmetry == 0)
    lpsy = 1;
  else if (Symmetry == 1)
    lpsy = 2;
  else if (Symmetry == 2)
    lpsy = 8;

  double psi4RR, psi4II;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1
    j = n - i * N_phi;

    int countlm = 0;
    for (int pl = spinw; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        for (int lp = 0; lp < lpsy; lp++)
        {
          switch (lp)
          {
          case 0: //+++ (theta, phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            psi4RR = shellf[InList * n];
            psi4II = shellf[InList * n + 1];
            break;
          case 1: //++- (pi-theta, phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            psi4RR = Rpsi4->SoA[2] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * shellf[InList * n + 1];
            break;
          case 2: //+-+ (theta, 2*pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            psi4RR = Rpsi4->SoA[1] * shellf[InList * n];
            psi4II = Ipsi4->SoA[1] * shellf[InList * n + 1];
            break;
          case 3: //+-- (pi-theta, 2*pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            psi4RR = Rpsi4->SoA[2] * Rpsi4->SoA[1] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * Ipsi4->SoA[1] * shellf[InList * n + 1];
            break;
          case 4: //-++ (theta, pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[0] * shellf[InList * n + 1];
            break;
          case 5: //-+- (pi-theta, pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[2] * Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * Ipsi4->SoA[0] * shellf[InList * n + 1];
            break;
          case 6: //--+ (theta, pi+phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[1] * Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[1] * Ipsi4->SoA[0] * shellf[InList * n + 1];
            break;
          case 7: //--- (pi-theta, pi+phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[2] * Rpsi4->SoA[1] * Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * Ipsi4->SoA[1] * Ipsi4->SoA[0] * shellf[InList * n + 1];
          }

          thetap = sqrt((2 * pl + 1.0) / 4.0 / PI) * misc::Wigner_d_function(pl, pm, spinw, costheta); // note the variation from -2 to 2
#ifdef GaussInt
          // wtcostheta is even function respect costheta
          RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi) * wtcostheta[i];
          IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi) * wtcostheta[i];
#else
          RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi);
          IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
        }
        countlm++; // no sanity check for countlm and NN which should be noted in the input parameters
      }
  }

  for (int ii = 0; ii < NN; ii++)
  {
#ifdef GaussInt
    RP_out[ii] = RP_out[ii] * rex * dphi;
    IP_out[ii] = IP_out[ii] * rex * dphi;
#else
    RP_out[ii] = RP_out[ii] * rex * dphi * dcostheta;
    IP_out[ii] = IP_out[ii] * rex * dphi * dcostheta;
#endif
  }
  //|------+  Communicate and sum the results from each processor.

  MPI_Allreduce(RP_out, RP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(IP_out, IP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //|------= Free memory.

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  delete[] RP_out;
  delete[] IP_out;
  DG_List->clearList();
}
void surface_integral::surf_Wave(double rex, int lev, cgh *GH, var *Rpsi4, var *Ipsi4,
                                 int spinw, int maxl, int NN, double *RP, double *IP,
                                 monitor *Monitor, MPI_Comm Comm_here) // NN is the length of RP and IP
{
  //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"start surface_integral::surf_Wave");

  int lmyrank;
  MPI_Comm_rank(Comm_here, &lmyrank);
  if (lmyrank == 0 && GH->grids[lev] != 1)
    if (Monitor->outfile)
      Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
    else
      cout << "WARNING: surface integral on multipatches" << endl;

  const int InList = 2;

  MyList<var> *DG_List = new MyList<var>(Rpsi4);
  DG_List->insert(Ipsi4);

  int n;
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  //    misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"before Interp_Points");

  GH->PatL[lev]->data->Interp_Points(DG_List, n_tot, pox, shellf, Symmetry, Comm_here);

  //    misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"after Interp_Points");

  int mp, Lp, Nmin, Nmax;

  int cpusize_here;
  MPI_Comm_size(Comm_here, &cpusize_here);

  mp = n_tot / cpusize_here;
  Lp = n_tot - cpusize_here * mp;

  if (Lp > lmyrank)
  {
    Nmin = lmyrank * mp + lmyrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = lmyrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  //|~~~~~> Integrate the dot product of Dphi with the surface normal.

  double *RP_out, *IP_out;
  RP_out = new double[NN];
  IP_out = new double[NN];

  for (int ii = 0; ii < NN; ii++)
  {
    RP_out[ii] = 0;
    IP_out[ii] = 0;
  }
  // theta part
  double costheta, thetap;
  double cosmphi, sinmphi;

  int i, j;
  int lpsy = 0;
  if (Symmetry == 0)
    lpsy = 1;
  else if (Symmetry == 1)
    lpsy = 2;
  else if (Symmetry == 2)
    lpsy = 8;

  double psi4RR, psi4II;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1
    j = n - i * N_phi;

    int countlm = 0;
    for (int pl = spinw; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        for (int lp = 0; lp < lpsy; lp++)
        {
          switch (lp)
          {
          case 0: //+++ (theta, phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            psi4RR = shellf[InList * n];
            psi4II = shellf[InList * n + 1];
            break;
          case 1: //++- (pi-theta, phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            psi4RR = Rpsi4->SoA[2] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * shellf[InList * n + 1];
            break;
          case 2: //+-+ (theta, 2*pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            psi4RR = Rpsi4->SoA[1] * shellf[InList * n];
            psi4II = Ipsi4->SoA[1] * shellf[InList * n + 1];
            break;
          case 3: //+-- (pi-theta, 2*pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            psi4RR = Rpsi4->SoA[2] * Rpsi4->SoA[1] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * Ipsi4->SoA[1] * shellf[InList * n + 1];
            break;
          case 4: //-++ (theta, pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[0] * shellf[InList * n + 1];
            break;
          case 5: //-+- (pi-theta, pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[2] * Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * Ipsi4->SoA[0] * shellf[InList * n + 1];
            break;
          case 6: //--+ (theta, pi+phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[1] * Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[1] * Ipsi4->SoA[0] * shellf[InList * n + 1];
            break;
          case 7: //--- (pi-theta, pi+phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[2] * Rpsi4->SoA[1] * Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * Ipsi4->SoA[1] * Ipsi4->SoA[0] * shellf[InList * n + 1];
          }

          thetap = sqrt((2 * pl + 1.0) / 4.0 / PI) * misc::Wigner_d_function(pl, pm, spinw, costheta); // note the variation from -2 to 2
#ifdef GaussInt
          // wtcostheta is even function respect costheta
          RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi) * wtcostheta[i];
          IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi) * wtcostheta[i];
#else
          RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi);
          IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
        }
        countlm++; // no sanity check for countlm and NN which should be noted in the input parameters
      }
  }

  for (int ii = 0; ii < NN; ii++)
  {
#ifdef GaussInt
    RP_out[ii] = RP_out[ii] * rex * dphi;
    IP_out[ii] = IP_out[ii] * rex * dphi;
#else
    RP_out[ii] = RP_out[ii] * rex * dphi * dcostheta;
    IP_out[ii] = IP_out[ii] * rex * dphi * dcostheta;
#endif
  }
  //|------+  Communicate and sum the results from each processor.

  MPI_Allreduce(RP_out, RP, NN, MPI_DOUBLE, MPI_SUM, Comm_here);
  MPI_Allreduce(IP_out, IP, NN, MPI_DOUBLE, MPI_SUM, Comm_here);

  //|------= Free memory.

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  delete[] RP_out;
  delete[] IP_out;
  DG_List->clearList();
}
//|----------------------------------------------------------------
//  for shell patch
//|----------------------------------------------------------------
void surface_integral::surf_Wave(double rex, int lev, ShellPatch *GH, var *Rpsi4, var *Ipsi4,
                                 int spinw, int maxl, int NN, double *RP, double *IP,
                                 monitor *Monitor) // NN is the length of RP and IP
{
  const int InList = 2;

  MyList<var> *DG_List = new MyList<var>(Rpsi4);
  DG_List->insert(Ipsi4);

  int n;
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  GH->Interp_Points(DG_List, n_tot, pox, shellf, Symmetry);

  int mp, Lp, Nmin, Nmax;

  mp = n_tot / cpusize;
  Lp = n_tot - cpusize * mp;

  if (Lp > myrank)
  {
    Nmin = myrank * mp + myrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = myrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  //|~~~~~> Integrate the dot product of Dphi with the surface normal.

  double *RP_out, *IP_out;
  RP_out = new double[NN];
  IP_out = new double[NN];

  for (int ii = 0; ii < NN; ii++)
  {
    RP_out[ii] = 0;
    IP_out[ii] = 0;
  }
  // theta part
  double costheta, thetap;
  double cosmphi, sinmphi;

  int i, j;
  int lpsy = 0;
  if (Symmetry == 0)
    lpsy = 1;
  else if (Symmetry == 1)
    lpsy = 2;
  else if (Symmetry == 2)
    lpsy = 8;

  double psi4RR, psi4II;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1
    j = n - i * N_phi;

    int countlm = 0;
    for (int pl = spinw; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        for (int lp = 0; lp < lpsy; lp++)
        {
          switch (lp)
          {
          case 0: //+++ (theta, phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            psi4RR = shellf[InList * n];
            psi4II = shellf[InList * n + 1];
            break;
          case 1: //++- (pi-theta, phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            psi4RR = Rpsi4->SoA[2] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * shellf[InList * n + 1];
            break;
          case 2: //+-+ (theta, 2*pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            psi4RR = Rpsi4->SoA[1] * shellf[InList * n];
            psi4II = Ipsi4->SoA[1] * shellf[InList * n + 1];
            break;
          case 3: //+-- (pi-theta, 2*pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            psi4RR = Rpsi4->SoA[2] * Rpsi4->SoA[1] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * Ipsi4->SoA[1] * shellf[InList * n + 1];
            break;
          case 4: //-++ (theta, pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[0] * shellf[InList * n + 1];
            break;
          case 5: //-+- (pi-theta, pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[2] * Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * Ipsi4->SoA[0] * shellf[InList * n + 1];
            break;
          case 6: //--+ (theta, pi+phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[1] * Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[1] * Ipsi4->SoA[0] * shellf[InList * n + 1];
            break;
          case 7: //--- (pi-theta, pi+phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            psi4RR = Rpsi4->SoA[2] * Rpsi4->SoA[1] * Rpsi4->SoA[0] * shellf[InList * n];
            psi4II = Ipsi4->SoA[2] * Ipsi4->SoA[1] * Ipsi4->SoA[0] * shellf[InList * n + 1];
          }

          thetap = sqrt((2 * pl + 1.0) / 4.0 / PI) * misc::Wigner_d_function(pl, pm, spinw, costheta); // note the variation from -2 to 2
#ifdef GaussInt
          // wtcostheta is even function respect costheta
          RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi) * wtcostheta[i];
          IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi) * wtcostheta[i];
#else
          RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi);
          IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
        }
        countlm++; // no sanity check for countlm and NN which should be noted in the input parameters
      }
  }

  for (int ii = 0; ii < NN; ii++)
  {
#ifdef GaussInt
    RP_out[ii] = RP_out[ii] * rex * dphi;
    IP_out[ii] = IP_out[ii] * rex * dphi;
#else
    RP_out[ii] = RP_out[ii] * rex * dphi * dcostheta;
    IP_out[ii] = IP_out[ii] * rex * dphi * dcostheta;
#endif
  }
  //|------+  Communicate and sum the results from each processor.

  MPI_Allreduce(RP_out, RP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(IP_out, IP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //|------= Free memory.

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  delete[] RP_out;
  delete[] IP_out;
  DG_List->clearList();
}
//|----------------------------------------------------------------
//  for shell patch
//  for EM wave specially symmetric case
//|----------------------------------------------------------------
void surface_integral::surf_Wave(double rex, int lev, ShellPatch *GH,
                                 var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
                                 var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
                                 int spinw, int maxl, int NN, double *RP, double *IP,
                                 monitor *Monitor) // NN is the length of RP and IP
{
  const int InList = 13;

  MyList<var> *DG_List = new MyList<var>(Ex);
  DG_List->insert(Ey);
  DG_List->insert(Ez);
  DG_List->insert(Bx);
  DG_List->insert(By);
  DG_List->insert(Bz);
  DG_List->insert(chi);
  DG_List->insert(gxx);
  DG_List->insert(gxy);
  DG_List->insert(gxz);
  DG_List->insert(gyy);
  DG_List->insert(gyz);
  DG_List->insert(gzz);

  int n;
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  GH->Interp_Points(DG_List, n_tot, pox, shellf, Symmetry);

  int mp, Lp, Nmin, Nmax;

  mp = n_tot / cpusize;
  Lp = n_tot - cpusize * mp;

  if (Lp > myrank)
  {
    Nmin = myrank * mp + myrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = myrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  //|~~~~~> Integrate the dot product of Dphi with the surface normal.

  double *RP_out, *IP_out;
  RP_out = new double[NN];
  IP_out = new double[NN];

  for (int ii = 0; ii < NN; ii++)
  {
    RP_out[ii] = 0;
    IP_out[ii] = 0;
  }
  // theta part
  double costheta, thetap;
  double cosmphi, sinmphi;

  int i, j;
  int lpsy = 0;
  if (Symmetry == 0)
    lpsy = 1;
  else if (Symmetry == 1)
    lpsy = 2;
  else if (Symmetry == 2)
    lpsy = 8;

  double psi4RR, psi4II;
  double px, py, pz;
  double pEx, pEy, pEz, pBx, pBy, pBz;
  double pchi, pgxx, pgxy, pgxz, pgyy, pgyz, pgzz;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1
    j = n - i * N_phi;

    int countlm = 0;
    for (int pl = spinw; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        for (int lp = 0; lp < lpsy; lp++)
        {
          px = pox[0][n];
          py = pox[1][n];
          pz = pox[2][n];
          pEx = shellf[InList * n];
          pEy = shellf[InList * n + 1];
          pEz = shellf[InList * n + 2];
          pBx = shellf[InList * n + 3];
          pBy = shellf[InList * n + 4];
          pBz = shellf[InList * n + 5];
          pchi = shellf[InList * n + 6];
          pgxx = shellf[InList * n + 7];
          pgxy = shellf[InList * n + 8];
          pgxz = shellf[InList * n + 9];
          pgyy = shellf[InList * n + 10];
          pgyz = shellf[InList * n + 11];
          pgzz = shellf[InList * n + 12];
          switch (lp)
          {
          case 0: //+++ (theta, phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            break;
          case 1: //++- (pi-theta, phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            pz = -pz;
            pEz = -pEz;
            pBx = -pBx;
            pBy = -pBy;
            pgxz = -pgxz;
            pgyz = -pgyz;
            break;
          case 2: //+-+ (theta, 2*pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            py = -py;
            pEy = -pEy;
            pBx = -pBx;
            pBz = -pBz;
            pgxy = -pgxy;
            pgyz = -pgyz;
            break;
          case 3: //+-- (pi-theta, 2*pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            py = -py;
            pz = -pz;
            pEz = -pEz;
            pBz = -pBz;
            pgxz = -pgxz;
            pEy = -pEy;
            pBy = -pBy;
            pgxy = -pgxy;
            break;
          case 4: //-++ (theta, pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            px = -px;
            pEx = -pEx;
            pBy = -pBy;
            pBz = -pBz;
            pgxy = -pgxy;
            pgxz = -pgxz;
            break;
          case 5: //-+- (pi-theta, pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            pz = -pz;
            px = -px;
            pEz = -pEz;
            pBz = -pBz;
            pgyz = -pgyz;
            pEx = -pEx;
            pBx = -pBx;
            pgxy = -pgxy;
            break;
          case 6: //--+ (theta, pi+phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            px = -px;
            py = -py;
            pEx = -pEx;
            pBx = -pBx;
            pgxz = -pgxz;
            pEy = -pEy;
            pBy = -pBy;
            pgyz = -pgyz;
            break;
          case 7: //--- (pi-theta, pi+phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            px = -px;
            py = -py;
            pz = -pz;
            pEx = -pEx;
            pEy = -pEy;
            pEz = -pEz;
          }

          f_getnpem2_point(px, py, pz, pchi, pgxx, pgxy, pgxz, pgyy, pgyz, pgzz, pEx, pEy, pEz, pBx, pBy, pBz,
                           psi4RR, psi4II);
          thetap = sqrt((2 * pl + 1.0) / 4.0 / PI) * misc::Wigner_d_function(pl, pm, spinw, costheta); // note the variation from -2 to 2

          //	 find back the one
          pchi = pchi + 1;
#ifdef GaussInt
          // wtcostheta is even function respect costheta
          RP_out[countlm] = RP_out[countlm] + thetap / pchi / pchi * (psi4RR * cosmphi + psi4II * sinmphi) * wtcostheta[i];
          IP_out[countlm] = IP_out[countlm] + thetap / pchi / pchi * (psi4II * cosmphi - psi4RR * sinmphi) * wtcostheta[i];
#else
          RP_out[countlm] = RP_out[countlm] + thetap / pchi / pchi * (psi4RR * cosmphi + psi4II * sinmphi);
          IP_out[countlm] = IP_out[countlm] + thetap / pchi / pchi * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
        }
        countlm++; // no sanity check for countlm and NN which should be noted in the input parameters
      }
  }

  for (int ii = 0; ii < NN; ii++)
  {
#ifdef GaussInt
    RP_out[ii] = RP_out[ii] * rex * dphi;
    IP_out[ii] = IP_out[ii] * rex * dphi;
#else
    RP_out[ii] = RP_out[ii] * rex * dphi * dcostheta;
    IP_out[ii] = IP_out[ii] * rex * dphi * dcostheta;
#endif
  }
  //|------+  Communicate and sum the results from each processor.

  MPI_Allreduce(RP_out, RP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(IP_out, IP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //|------= Free memory.

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  delete[] RP_out;
  delete[] IP_out;
  DG_List->clearList();
}
//|----------------------------------------------------------------
//  for shell patch
//  for EM wave specially symmetric case
//  unify for phi1 and phi2
//|----------------------------------------------------------------
void surface_integral::surf_Wave(double rex, int lev, ShellPatch *GH,
                                 var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
                                 var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
                                 int spinw, int maxl, int NN, double *RP, double *IP,
                                 monitor *Monitor,
                                 void (*funcs)(double &, double &, double &,
                                               double &, double &, double &, double &, double &, double &, double &,
                                               double &, double &, double &, double &, double &, double &,
                                               double &, double &)) // NN is the length of RP and IP
{
  const int InList = 13;

  MyList<var> *DG_List = new MyList<var>(Ex);
  DG_List->insert(Ey);
  DG_List->insert(Ez);
  DG_List->insert(Bx);
  DG_List->insert(By);
  DG_List->insert(Bz);
  DG_List->insert(chi);
  DG_List->insert(gxx);
  DG_List->insert(gxy);
  DG_List->insert(gxz);
  DG_List->insert(gyy);
  DG_List->insert(gyz);
  DG_List->insert(gzz);

  int n;
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  GH->Interp_Points(DG_List, n_tot, pox, shellf, Symmetry);

  double *RP_out, *IP_out;
  RP_out = new double[NN];
  IP_out = new double[NN];

  for (int ii = 0; ii < NN; ii++)
  {
    RP_out[ii] = 0;
    IP_out[ii] = 0;
  }

#if 0
// for debug    
  if(myrank==0)
  {
    double costheta, thetap;
    double cosmphi,sinmphi;

    int i,j;
    int lpsy=0;
         if( Symmetry == 0 )     lpsy=1;
    else if( Symmetry == 1 )     lpsy=2;
    else if( Symmetry == 2 )     lpsy=8;

    double psi4RR,psi4II;
    double px,py,pz;
    double pEx,pEy,pEz,pBx,pBy,pBz;
    double pchi,pgxx,pgxy,pgxz,pgyy,pgyz,pgzz;
    for( n = 0; n <= n_tot-1; n++) 
     {
//       need round off always	     
        i = int(n/N_phi); // int(1.723) = 1, int(-1.732) = -1
        j = n - i * N_phi;
        
	for(int lp=0;lp<lpsy;lp++)
	{
         px = pox[0][n];		 
         py = pox[1][n];		 
         pz = pox[2][n];		 
	 pEx = shellf[InList*n  ];
	 pEy = shellf[InList*n+1];
	 pEz = shellf[InList*n+2];
	 pBx = shellf[InList*n+3];
	 pBy = shellf[InList*n+4];
	 pBz = shellf[InList*n+5];
	 pchi = shellf[InList*n+6];
	 pgxx = shellf[InList*n+7];
	 pgxy = shellf[InList*n+8];
	 pgxz = shellf[InList*n+9];
	 pgyy = shellf[InList*n+10];
	 pgyz = shellf[InList*n+11];
	 pgzz = shellf[InList*n+12];
 	 switch(lp)
	 {
	  case 1:  //++- (pi-theta, phi)
	  pz = -pz;
	  pEz = -pEz;
	  pBx = -pBx;
	  pBy = -pBy;
	  pgxz = -pgxz;
	  pgyz = -pgyz;
	  break;
	  case 2:  //+-+ (theta, 2*pi-phi)
	  py = -py;
	  pEy = -pEy;
	  pBx = -pBx;
	  pBz = -pBz;
	  pgxy = -pgxy;
	  pgyz = -pgyz;
	  break;
	  case 3:  //+-- (pi-theta, 2*pi-phi)
	  py = -py;
	  pz = -pz;
	  pEz = -pEz;
	  pBz = -pBz;;
	  pgxz = -pgxz;
	  pEy = -pEy;
	  pBy = -pBy;
	  pgxy = -pgxy;
	  break;
	  case 4:  //-++ (theta, pi-phi)
	  px = -px;
	  pEx = -pEx;
	  pBy = -pBy;
	  pBz = -pBz;
	  pgxy = -pgxy;
	  pgxz = -pgxz;
	  break;
	  case 5:  //-+- (pi-theta, pi-phi)
	  pz = -pz;
	  px = -px;
	  pEz = -pEz;
	  pBz = -pBz;
	  pgyz = -pgyz;
	  pEx = -pEx;
	  pBx = -pBx;
	  pgxy = -pgxy;
	  break;
	  case 6:  //--+ (theta, pi+phi)
	  px = -px;
	  py = -py;
	  pEx = -pEx;
	  pBx = -pBx;
	  pgxz = -pgxz;
	  pEy = -pEy;
	  pBy = -pBy;
	  pgyz = -pgyz;
	  break;
	  case 7:  //--- (pi-theta, pi+phi)
	  px = -px;
	  py = -py;
	  pz = -pz;
	  pEx = -pEx;
	  pEy = -pEy;
	  pEz = -pEz;
	 }
	  
	 funcs(px,py,pz,pchi,pgxx,pgxy,pgxz,pgyy,pgyz,pgzz,pEx,pEy,pEz,pBx,pBy,pBz,
			  psi4RR,psi4II);
//	 if(n==0 || n==N_phi/2-1 || n==N_phi/2 || n==N_phi-1 ||
//	    n==N_phi*(N_theta-1)+0 || n==N_phi*(N_theta-1)+N_phi/2-1 || n==N_phi*(N_theta-1)+N_phi/2 || n==N_phi*(N_theta-1)+N_phi-1)
//	 cout<<px<<","<<py<<","<<pz<<","<<pchi<<","<<pgxx<<","<<pgxy<<","<<pgxz<<","<<pgyy<<","<<pgyz<<","<<pgzz<<","<<pEx<<","
//	     <<pEy<<","<<pEz<<","<<pBx<<","<<pBy<<","<<pBz<<","<<psi4RR<<","<<psi4II<<endl<<endl;

// find back the one
        pchi = pchi+1;

	int countlm=0;
	for(int pl=spinw;pl<maxl+1;pl++)
          for(int pm=-pl;pm<pl+1;pm++)
	  {
 	 switch(lp)
	 {
	  case 0:  //+++ (theta, phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi = sin(pm * (j+0.5) * dphi);
	  break;
	  case 1:  //++- (pi-theta, phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi = sin(pm * (j+0.5) * dphi);
	  break;
	  case 2:  //+-+ (theta, 2*pi-phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi =-sin(pm * (j+0.5) * dphi);
	  break;
	  case 3:  //+-- (pi-theta, 2*pi-phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi =-sin(pm * (j+0.5) * dphi);
	  break;
	  case 4:  //-++ (theta, pi-phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (PI - (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI - (j+0.5) * dphi));
	  break;
	  case 5:  //-+- (pi-theta, pi-phi)
          costheta = -arcostheta[i];
 	  cosmphi = cos(pm * (PI - (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI - (j+0.5) * dphi));
	  break;
	  case 6:  //--+ (theta, pi+phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (PI + (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI + (j+0.5) * dphi));
	  break;
	  case 7:  //--- (pi-theta, pi+phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (PI + (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI + (j+0.5) * dphi));
	 }
   	    thetap = sqrt((2*pl+1.0)/4.0/PI)*misc::Wigner_d_function(pl,pm,spinw,costheta); //note the variation from -2 to 2

#ifdef GaussInt
// wtcostheta is even function respect costheta
            RP_out[countlm] = RP_out[countlm] + thetap/pchi/pchi * (psi4RR * cosmphi + psi4II * sinmphi)*wtcostheta[i];
  	    IP_out[countlm] = IP_out[countlm] + thetap/pchi/pchi * (psi4II * cosmphi - psi4RR * sinmphi)*wtcostheta[i];
	    if(pl==2 && pm==0) cout<<countlm+1<<","<<RP_out[countlm] * rex * dphi<<endl;
#else	 
            RP_out[countlm] = RP_out[countlm] + thetap/pchi/pchi * (psi4RR * cosmphi + psi4II * sinmphi);
  	    IP_out[countlm] = IP_out[countlm] + thetap/pchi/pchi * (psi4II * cosmphi - psi4RR * sinmphi);
#endif	 
  	    countlm++;  //no sanity check for countlm and NN which should be noted in the input parameters
	  }
	}
//        if(Symmetry == 2) MPI_Abort(MPI_COMM_WORLD,1);
     }
     MPI_Abort(MPI_COMM_WORLD,1);
  }
#else
  int mp, Lp, Nmin, Nmax;

  mp = n_tot / cpusize;
  Lp = n_tot - cpusize * mp;

  if (Lp > myrank)
  {
    Nmin = myrank * mp + myrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = myrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  // theta part
  double costheta, thetap;
  double cosmphi, sinmphi;

  int i, j;
  int lpsy = 0;
  if (Symmetry == 0)
    lpsy = 1;
  else if (Symmetry == 1)
    lpsy = 2;
  else if (Symmetry == 2)
    lpsy = 8;

  double psi4RR, psi4II;
  double px, py, pz;
  double pEx, pEy, pEz, pBx, pBy, pBz;
  double pchi, pgxx, pgxy, pgxz, pgyy, pgyz, pgzz;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1
    j = n - i * N_phi;

    int countlm = 0;
    for (int pl = spinw; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        for (int lp = 0; lp < lpsy; lp++)
        {
          px = pox[0][n];
          py = pox[1][n];
          pz = pox[2][n];
          pEx = shellf[InList * n];
          pEy = shellf[InList * n + 1];
          pEz = shellf[InList * n + 2];
          pBx = shellf[InList * n + 3];
          pBy = shellf[InList * n + 4];
          pBz = shellf[InList * n + 5];
          pchi = shellf[InList * n + 6];
          pgxx = shellf[InList * n + 7];
          pgxy = shellf[InList * n + 8];
          pgxz = shellf[InList * n + 9];
          pgyy = shellf[InList * n + 10];
          pgyz = shellf[InList * n + 11];
          pgzz = shellf[InList * n + 12];
          switch (lp)
          {
          case 0: //+++ (theta, phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            break;
          case 1: //++- (pi-theta, phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            pz = -pz;
            pEz = -pEz;
            pBx = -pBx;
            pBy = -pBy;
            pgxz = -pgxz;
            pgyz = -pgyz;
            break;
          case 2: //+-+ (theta, 2*pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            py = -py;
            pEy = -pEy;
            pBx = -pBx;
            pBz = -pBz;
            pgxy = -pgxy;
            pgyz = -pgyz;
            break;
          case 3: //+-- (pi-theta, 2*pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            py = -py;
            pz = -pz;
            pEz = -pEz;
            pBz = -pBz;
            pgxz = -pgxz;
            pEy = -pEy;
            pBy = -pBy;
            pgxy = -pgxy;
            break;
          case 4: //-++ (theta, pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            px = -px;
            pEx = -pEx;
            pBy = -pBy;
            pBz = -pBz;
            pgxy = -pgxy;
            pgxz = -pgxz;
            break;
          case 5: //-+- (pi-theta, pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            pz = -pz;
            px = -px;
            pEz = -pEz;
            pBz = -pBz;
            pgyz = -pgyz;
            pEx = -pEx;
            pBx = -pBx;
            pgxy = -pgxy;
            break;
          case 6: //--+ (theta, pi+phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            px = -px;
            py = -py;
            pEx = -pEx;
            pBx = -pBx;
            pgxz = -pgxz;
            pEy = -pEy;
            pBy = -pBy;
            pgyz = -pgyz;
            break;
          case 7: //--- (pi-theta, pi+phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            px = -px;
            py = -py;
            pz = -pz;
            pEx = -pEx;
            pEy = -pEy;
            pEz = -pEz;
          }

          funcs(px, py, pz, pchi, pgxx, pgxy, pgxz, pgyy, pgyz, pgzz, pEx, pEy, pEz, pBx, pBy, pBz,
                psi4RR, psi4II);
          thetap = sqrt((2 * pl + 1.0) / 4.0 / PI) * misc::Wigner_d_function(pl, pm, spinw, costheta); // note the variation from -2 to 2

          //	 find back the one
          pchi = pchi + 1;
#ifdef GaussInt
          // wtcostheta is even function respect costheta
          RP_out[countlm] = RP_out[countlm] + thetap / pchi / pchi * (psi4RR * cosmphi + psi4II * sinmphi) * wtcostheta[i];
          IP_out[countlm] = IP_out[countlm] + thetap / pchi / pchi * (psi4II * cosmphi - psi4RR * sinmphi) * wtcostheta[i];
#else
          RP_out[countlm] = RP_out[countlm] + thetap / pchi / pchi * (psi4RR * cosmphi + psi4II * sinmphi);
          IP_out[countlm] = IP_out[countlm] + thetap / pchi / pchi * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
        }
        countlm++; // no sanity check for countlm and NN which should be noted in the input parameters
      }
  }
#endif

  for (int ii = 0; ii < NN; ii++)
  {
#ifdef GaussInt
    RP_out[ii] = RP_out[ii] * rex * dphi;
    IP_out[ii] = IP_out[ii] * rex * dphi;
#else
    RP_out[ii] = RP_out[ii] * rex * dphi * dcostheta;
    IP_out[ii] = IP_out[ii] * rex * dphi * dcostheta;
#endif
  }
  //|------+  Communicate and sum the results from each processor.

  MPI_Allreduce(RP_out, RP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(IP_out, IP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //|------= Free memory.

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  delete[] RP_out;
  delete[] IP_out;
  DG_List->clearList();
}
//|----------------------------------------------------------------
//  for box
//  for EM wave specially symmetric case
//  unify for phi1 and phi2
//|----------------------------------------------------------------
void surface_integral::surf_Wave(double rex, int lev, cgh *GH,
                                 var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
                                 var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
                                 int spinw, int maxl, int NN, double *RP, double *IP,
                                 monitor *Monitor,
                                 void (*funcs)(double &, double &, double &,
                                               double &, double &, double &, double &, double &, double &, double &,
                                               double &, double &, double &, double &, double &, double &,
                                               double &, double &)) // NN is the length of RP and IP
{
  const int InList = 13;

  MyList<var> *DG_List = new MyList<var>(Ex);
  DG_List->insert(Ey);
  DG_List->insert(Ez);
  DG_List->insert(Bx);
  DG_List->insert(By);
  DG_List->insert(Bz);
  DG_List->insert(chi);
  DG_List->insert(gxx);
  DG_List->insert(gxy);
  DG_List->insert(gxz);
  DG_List->insert(gyy);
  DG_List->insert(gyz);
  DG_List->insert(gzz);

  int n;
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  GH->PatL[lev]->data->Interp_Points(DG_List, n_tot, pox, shellf, Symmetry);

  double *RP_out, *IP_out;
  RP_out = new double[NN];
  IP_out = new double[NN];

  for (int ii = 0; ii < NN; ii++)
  {
    RP_out[ii] = 0;
    IP_out[ii] = 0;
  }

#if 0
// for debug    
  if(myrank==0)
  {
    double costheta, thetap;
    double cosmphi,sinmphi;

    int i,j;
    int lpsy=0;
         if( Symmetry == 0 )     lpsy=1;
    else if( Symmetry == 1 )     lpsy=2;
    else if( Symmetry == 2 )     lpsy=8;

    double psi4RR,psi4II;
    double px,py,pz;
    double pEx,pEy,pEz,pBx,pBy,pBz;
    double pchi,pgxx,pgxy,pgxz,pgyy,pgyz,pgzz;
    for( n = 0; n <= n_tot-1; n++) 
     {
//       need round off always	     
        i = int(n/N_phi); // int(1.723) = 1, int(-1.732) = -1
        j = n - i * N_phi;
        
	for(int lp=0;lp<lpsy;lp++)
	{
         px = pox[0][n];		 
         py = pox[1][n];		 
         pz = pox[2][n];		 
	 pEx = shellf[InList*n  ];
	 pEy = shellf[InList*n+1];
	 pEz = shellf[InList*n+2];
	 pBx = shellf[InList*n+3];
	 pBy = shellf[InList*n+4];
	 pBz = shellf[InList*n+5];
	 pchi = shellf[InList*n+6];
	 pgxx = shellf[InList*n+7];
	 pgxy = shellf[InList*n+8];
	 pgxz = shellf[InList*n+9];
	 pgyy = shellf[InList*n+10];
	 pgyz = shellf[InList*n+11];
	 pgzz = shellf[InList*n+12];
 	 switch(lp)
	 {
	  case 1:  //++- (pi-theta, phi)
	  pz = -pz;
	  pEz = -pEz;
	  pBx = -pBx;
	  pBy = -pBy;
	  pgxz = -pgxz;
	  pgyz = -pgyz;
	  break;
	  case 2:  //+-+ (theta, 2*pi-phi)
	  py = -py;
	  pEy = -pEy;
	  pBx = -pBx;
	  pBz = -pBz;
	  pgxy = -pgxy;
	  pgyz = -pgyz;
	  break;
	  case 3:  //+-- (pi-theta, 2*pi-phi)
	  py = -py;
	  pz = -pz;
	  pEz = -pEz;
	  pBz = -pBz;;
	  pgxz = -pgxz;
	  pEy = -pEy;
	  pBy = -pBy;
	  pgxy = -pgxy;
	  break;
	  case 4:  //-++ (theta, pi-phi)
	  px = -px;
	  pEx = -pEx;
	  pBy = -pBy;
	  pBz = -pBz;
	  pgxy = -pgxy;
	  pgxz = -pgxz;
	  break;
	  case 5:  //-+- (pi-theta, pi-phi)
	  pz = -pz;
	  px = -px;
	  pEz = -pEz;
	  pBz = -pBz;
	  pgyz = -pgyz;
	  pEx = -pEx;
	  pBx = -pBx;
	  pgxy = -pgxy;
	  break;
	  case 6:  //--+ (theta, pi+phi)
	  px = -px;
	  py = -py;
	  pEx = -pEx;
	  pBx = -pBx;
	  pgxz = -pgxz;
	  pEy = -pEy;
	  pBy = -pBy;
	  pgyz = -pgyz;
	  break;
	  case 7:  //--- (pi-theta, pi+phi)
	  px = -px;
	  py = -py;
	  pz = -pz;
	  pEx = -pEx;
	  pEy = -pEy;
	  pEz = -pEz;
	 }
	  
	 funcs(px,py,pz,pchi,pgxx,pgxy,pgxz,pgyy,pgyz,pgzz,pEx,pEy,pEz,pBx,pBy,pBz,
			  psi4RR,psi4II);
//	 if(n==0 || n==N_phi/2-1 || n==N_phi/2 || n==N_phi-1 ||
//	    n==N_phi*(N_theta-1)+0 || n==N_phi*(N_theta-1)+N_phi/2-1 || n==N_phi*(N_theta-1)+N_phi/2 || n==N_phi*(N_theta-1)+N_phi-1)
//	 cout<<px<<","<<py<<","<<pz<<","<<pchi<<","<<pgxx<<","<<pgxy<<","<<pgxz<<","<<pgyy<<","<<pgyz<<","<<pgzz<<","<<pEx<<","
//	     <<pEy<<","<<pEz<<","<<pBx<<","<<pBy<<","<<pBz<<","<<psi4RR<<","<<psi4II<<endl<<endl;

// find back the one
        pchi = pchi+1;

	int countlm=0;
	for(int pl=spinw;pl<maxl+1;pl++)
          for(int pm=-pl;pm<pl+1;pm++)
	  {
 	 switch(lp)
	 {
	  case 0:  //+++ (theta, phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi = sin(pm * (j+0.5) * dphi);
	  break;
	  case 1:  //++- (pi-theta, phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi = sin(pm * (j+0.5) * dphi);
	  break;
	  case 2:  //+-+ (theta, 2*pi-phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi =-sin(pm * (j+0.5) * dphi);
	  break;
	  case 3:  //+-- (pi-theta, 2*pi-phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (j+0.5) * dphi);
	  sinmphi =-sin(pm * (j+0.5) * dphi);
	  break;
	  case 4:  //-++ (theta, pi-phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (PI - (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI - (j+0.5) * dphi));
	  break;
	  case 5:  //-+- (pi-theta, pi-phi)
          costheta = -arcostheta[i];
 	  cosmphi = cos(pm * (PI - (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI - (j+0.5) * dphi));
	  break;
	  case 6:  //--+ (theta, pi+phi)
          costheta = arcostheta[i];
	  cosmphi = cos(pm * (PI + (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI + (j+0.5) * dphi));
	  break;
	  case 7:  //--- (pi-theta, pi+phi)
          costheta = -arcostheta[i];
	  cosmphi = cos(pm * (PI + (j+0.5) * dphi));
	  sinmphi = sin(pm * (PI + (j+0.5) * dphi));
	 }
   	    thetap = sqrt((2*pl+1.0)/4.0/PI)*misc::Wigner_d_function(pl,pm,spinw,costheta); //note the variation from -2 to 2

#ifdef GaussInt
// wtcostheta is even function respect costheta
            RP_out[countlm] = RP_out[countlm] + thetap/pchi/pchi * (psi4RR * cosmphi + psi4II * sinmphi)*wtcostheta[i];
  	    IP_out[countlm] = IP_out[countlm] + thetap/pchi/pchi * (psi4II * cosmphi - psi4RR * sinmphi)*wtcostheta[i];
	    if(pl==2 && pm==0) cout<<countlm+1<<","<<RP_out[countlm] * rex * dphi<<endl;
#else	 
            RP_out[countlm] = RP_out[countlm] + thetap/pchi/pchi * (psi4RR * cosmphi + psi4II * sinmphi);
  	    IP_out[countlm] = IP_out[countlm] + thetap/pchi/pchi * (psi4II * cosmphi - psi4RR * sinmphi);
#endif	 
  	    countlm++;  //no sanity check for countlm and NN which should be noted in the input parameters
	  }
	}
//        if(Symmetry == 2) MPI_Abort(MPI_COMM_WORLD,1);
     }
     MPI_Abort(MPI_COMM_WORLD,1);
  }
#else
  int mp, Lp, Nmin, Nmax;

  mp = n_tot / cpusize;
  Lp = n_tot - cpusize * mp;

  if (Lp > myrank)
  {
    Nmin = myrank * mp + myrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = myrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  // theta part
  double costheta, thetap;
  double cosmphi, sinmphi;

  int i, j;
  int lpsy = 0;
  if (Symmetry == 0)
    lpsy = 1;
  else if (Symmetry == 1)
    lpsy = 2;
  else if (Symmetry == 2)
    lpsy = 8;

  double psi4RR, psi4II;
  double px, py, pz;
  double pEx, pEy, pEz, pBx, pBy, pBz;
  double pchi, pgxx, pgxy, pgxz, pgyy, pgyz, pgzz;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1
    j = n - i * N_phi;

    int countlm = 0;
    for (int pl = spinw; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        for (int lp = 0; lp < lpsy; lp++)
        {
          px = pox[0][n];
          py = pox[1][n];
          pz = pox[2][n];
          pEx = shellf[InList * n];
          pEy = shellf[InList * n + 1];
          pEz = shellf[InList * n + 2];
          pBx = shellf[InList * n + 3];
          pBy = shellf[InList * n + 4];
          pBz = shellf[InList * n + 5];
          pchi = shellf[InList * n + 6];
          pgxx = shellf[InList * n + 7];
          pgxy = shellf[InList * n + 8];
          pgxz = shellf[InList * n + 9];
          pgyy = shellf[InList * n + 10];
          pgyz = shellf[InList * n + 11];
          pgzz = shellf[InList * n + 12];
          switch (lp)
          {
          case 0: //+++ (theta, phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            break;
          case 1: //++- (pi-theta, phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            pz = -pz;
            pEz = -pEz;
            pBx = -pBx;
            pBy = -pBy;
            pgxz = -pgxz;
            pgyz = -pgyz;
            break;
          case 2: //+-+ (theta, 2*pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            py = -py;
            pEy = -pEy;
            pBx = -pBx;
            pBz = -pBz;
            pgxy = -pgxy;
            pgyz = -pgyz;
            break;
          case 3: //+-- (pi-theta, 2*pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            py = -py;
            pz = -pz;
            pEz = -pEz;
            pBz = -pBz;
            pgxz = -pgxz;
            pEy = -pEy;
            pBy = -pBy;
            pgxy = -pgxy;
            break;
          case 4: //-++ (theta, pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            px = -px;
            pEx = -pEx;
            pBy = -pBy;
            pBz = -pBz;
            pgxy = -pgxy;
            pgxz = -pgxz;
            break;
          case 5: //-+- (pi-theta, pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            pz = -pz;
            px = -px;
            pEz = -pEz;
            pBz = -pBz;
            pgyz = -pgyz;
            pEx = -pEx;
            pBx = -pBx;
            pgxy = -pgxy;
            break;
          case 6: //--+ (theta, pi+phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            px = -px;
            py = -py;
            pEx = -pEx;
            pBx = -pBx;
            pgxz = -pgxz;
            pEy = -pEy;
            pBy = -pBy;
            pgyz = -pgyz;
            break;
          case 7: //--- (pi-theta, pi+phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            px = -px;
            py = -py;
            pz = -pz;
            pEx = -pEx;
            pEy = -pEy;
            pEz = -pEz;
          }

          funcs(px, py, pz, pchi, pgxx, pgxy, pgxz, pgyy, pgyz, pgzz, pEx, pEy, pEz, pBx, pBy, pBz,
                psi4RR, psi4II);
          thetap = sqrt((2 * pl + 1.0) / 4.0 / PI) * misc::Wigner_d_function(pl, pm, spinw, costheta); // note the variation from -2 to 2

          //	 find back the one
          pchi = pchi + 1;
#ifdef GaussInt
          // wtcostheta is even function respect costheta
          RP_out[countlm] = RP_out[countlm] + thetap / pchi / pchi * (psi4RR * cosmphi + psi4II * sinmphi) * wtcostheta[i];
          IP_out[countlm] = IP_out[countlm] + thetap / pchi / pchi * (psi4II * cosmphi - psi4RR * sinmphi) * wtcostheta[i];
#else
          RP_out[countlm] = RP_out[countlm] + thetap / pchi / pchi * (psi4RR * cosmphi + psi4II * sinmphi);
          IP_out[countlm] = IP_out[countlm] + thetap / pchi / pchi * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
        }
        countlm++; // no sanity check for countlm and NN which should be noted in the input parameters
      }
  }
#endif

  for (int ii = 0; ii < NN; ii++)
  {
#ifdef GaussInt
    RP_out[ii] = RP_out[ii] * rex * dphi;
    IP_out[ii] = IP_out[ii] * rex * dphi;
#else
    RP_out[ii] = RP_out[ii] * rex * dphi * dcostheta;
    IP_out[ii] = IP_out[ii] * rex * dphi * dcostheta;
#endif
  }
  //|------+  Communicate and sum the results from each processor.

  MPI_Allreduce(RP_out, RP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(IP_out, IP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //|------= Free memory.

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  delete[] RP_out;
  delete[] IP_out;
  DG_List->clearList();
}
//|----------------------------------------------------------------
//  for null shell patch2
//|----------------------------------------------------------------
// rex is x instead of r
void surface_integral::surf_Wave(double rex, int lev, NullShellPatch2 *GH, var *Rpsi4, var *Ipsi4,
                                 int spinw, int maxl, int NN, double *RP, double *IP,
                                 monitor *Monitor) // NN is the length of RP and IP
// spinw 0 for scalar; 1 for electricmagnetic wave; 2 for gravitaitonal wave
// we always assume spinw >= 0
{
  const int InList = 2;

  MyList<var> *DG_List = new MyList<var>(Rpsi4);
  DG_List->insert(Ipsi4);

  int n;
  // since we used x instead of r, these global coordinates are fake
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  GH->Interp_Points_2D(DG_List, n_tot, pox, shellf, Symmetry);

  int mp, Lp, Nmin, Nmax;

  mp = n_tot / cpusize;
  Lp = n_tot - cpusize * mp;

  if (Lp > myrank)
  {
    Nmin = myrank * mp + myrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = myrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  //|~~~~~> Integrate the dot product of Dphi with the surface normal.

  double *RP_out, *IP_out;
  RP_out = new double[NN];
  IP_out = new double[NN];

  for (int ii = 0; ii < NN; ii++)
  {
    RP_out[ii] = 0;
    IP_out[ii] = 0;
  }
  // theta part
  double costheta, thetap;
  double cosmphi, sinmphi;

  int i, j;
  int lpsy = 0;
  if (Symmetry == 0)
    lpsy = 1;
  else if (Symmetry == 1)
    lpsy = 2;
  else if (Symmetry == 2)
    lpsy = 8;

  double psi4RR, psi4II;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1
    j = n - i * N_phi;

    int countlm = 0;
    for (int pl = spinw; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        for (int lp = 0; lp < lpsy; lp++)
        {
          switch (lp)
          {
          case 0: //+++ (theta, phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            psi4RR = shellf[InList * n];
            psi4II = shellf[InList * n + 1];
            break;
          case 1: //++- (pi-theta, phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            psi4RR = shellf[InList * n];
            psi4II = -shellf[InList * n + 1];
            break;
          case 2: //+-+ (theta, 2*pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            psi4RR = shellf[InList * n];
            psi4II = -shellf[InList * n + 1];
            break;
          case 3: //+-- (pi-theta, 2*pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            psi4RR = shellf[InList * n];
            psi4II = shellf[InList * n + 1];
            break;
          case 4: //-++ (theta, pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            psi4RR = shellf[InList * n];
            psi4II = -shellf[InList * n + 1];
            break;
          case 5: //-+- (pi-theta, pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            psi4RR = shellf[InList * n];
            psi4II = shellf[InList * n + 1];
            break;
          case 6: //--+ (theta, pi+phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            psi4RR = shellf[InList * n];
            psi4II = shellf[InList * n + 1];
            break;
          case 7: //--- (pi-theta, pi+phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            psi4RR = shellf[InList * n];
            psi4II = -shellf[InList * n + 1];
          }

          thetap = sqrt((2 * pl + 1.0) / 4.0 / PI) * misc::Wigner_d_function(pl, pm, spinw, costheta); // note the variation from -2 to 2
                                                                                                       // based on Eq.(41) of PRD 77, 024027 (2008)
#ifdef GaussInt
          // wtcostheta is even function respect costheta
          RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi) * wtcostheta[i];
          IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi) * wtcostheta[i];
#else
          RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi); // + is because \bar of \bar{Y^s_lm} in Eq.(40)
                                                                                              // of PRD 77, 024027 (2008)
          IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
        }
        countlm++; // no sanity check for countlm and NN which should be noted in the input parameters
      }
  }

  for (int ii = 0; ii < NN; ii++)
  {
// do not need multiply with rex for null shell
#ifdef GaussInt
    RP_out[ii] = RP_out[ii] * dphi;
    IP_out[ii] = IP_out[ii] * dphi;
#else
    RP_out[ii] = RP_out[ii] * dphi * dcostheta;
    IP_out[ii] = IP_out[ii] * dphi * dcostheta;
#endif
  }
  //|------+  Communicate and sum the results from each processor.

  MPI_Allreduce(RP_out, RP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(IP_out, IP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //|------= Free memory.

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  delete[] RP_out;
  delete[] IP_out;
  DG_List->clearList();
}
//|----------------------------------------------------------------
//  for null shell patch
//|----------------------------------------------------------------
// rex is x instead of r
void surface_integral::surf_Wave(double rex, int lev, NullShellPatch *GH, var *Rpsi4, var *Ipsi4,
                                 int spinw, int maxl, int NN, double *RP, double *IP,
                                 monitor *Monitor) // NN is the length of RP and IP
// spinw 0 for scalar; 1 for electricmagnetic wave; 2 for gravitaitonal wave
// we always assume spinw >= 0
{
  const int InList = 2;

  MyList<var> *DG_List = new MyList<var>(Rpsi4);
  DG_List->insert(Ipsi4);

  int n;
  // since we used x instead of r, these global coordinates are fake
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  GH->Interp_Points_2D(DG_List, n_tot, pox, shellf, Symmetry);

  int mp, Lp, Nmin, Nmax;

  mp = n_tot / cpusize;
  Lp = n_tot - cpusize * mp;

  if (Lp > myrank)
  {
    Nmin = myrank * mp + myrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = myrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  //|~~~~~> Integrate the dot product of Dphi with the surface normal.

  double *RP_out, *IP_out;
  RP_out = new double[NN];
  IP_out = new double[NN];

  for (int ii = 0; ii < NN; ii++)
  {
    RP_out[ii] = 0;
    IP_out[ii] = 0;
  }
  // theta part
  double costheta, thetap;
  double cosmphi, sinmphi;

  int i, j;
  int lpsy = 0;
  if (Symmetry == 0)
    lpsy = 1;
  else if (Symmetry == 1)
    lpsy = 2;
  else if (Symmetry == 2)
    lpsy = 8;

  double psi4RR, psi4II;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1
    j = n - i * N_phi;

    int countlm = 0;
    for (int pl = spinw; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        for (int lp = 0; lp < lpsy; lp++)
        {
          switch (lp)
          {
          case 0: //+++ (theta, phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            psi4RR = shellf[InList * n];
            psi4II = shellf[InList * n + 1];
            break;
          case 1: //++- (pi-theta, phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            psi4RR = shellf[InList * n];
            psi4II = -shellf[InList * n + 1];
            break;
          case 2: //+-+ (theta, 2*pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            psi4RR = shellf[InList * n];
            psi4II = -shellf[InList * n + 1];
            break;
          case 3: //+-- (pi-theta, 2*pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            psi4RR = shellf[InList * n];
            psi4II = shellf[InList * n + 1];
            break;
          case 4: //-++ (theta, pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            psi4RR = shellf[InList * n];
            psi4II = -shellf[InList * n + 1];
            break;
          case 5: //-+- (pi-theta, pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            psi4RR = shellf[InList * n];
            psi4II = shellf[InList * n + 1];
            break;
          case 6: //--+ (theta, pi+phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            psi4RR = shellf[InList * n];
            psi4II = shellf[InList * n + 1];
            break;
          case 7: //--- (pi-theta, pi+phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            psi4RR = shellf[InList * n];
            psi4II = -shellf[InList * n + 1];
          }

          thetap = sqrt((2 * pl + 1.0) / 4.0 / PI) * misc::Wigner_d_function(pl, pm, spinw, costheta); // note the variation from -2 to 2
                                                                                                       // based on Eq.(41) of PRD 77, 024027 (2008)
#ifdef GaussInt
          // wtcostheta is even function respect costheta
          RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi) * wtcostheta[i];
          IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi) * wtcostheta[i];
#else
          RP_out[countlm] = RP_out[countlm] + thetap * (psi4RR * cosmphi + psi4II * sinmphi); // + is because \bar of \bar{Y^s_lm} in Eq.(40)
                                                                                              // of PRD 77, 024027 (2008)
          IP_out[countlm] = IP_out[countlm] + thetap * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
        }
        countlm++; // no sanity check for countlm and NN which should be noted in the input parameters
      }
  }

  for (int ii = 0; ii < NN; ii++)
  {
// do not need multiply with rex for null shell
#ifdef GaussInt
    RP_out[ii] = RP_out[ii] * dphi;
    IP_out[ii] = IP_out[ii] * dphi;
#else
    RP_out[ii] = RP_out[ii] * dphi * dcostheta;
    IP_out[ii] = IP_out[ii] * dphi * dcostheta;
#endif
  }
  //|------+  Communicate and sum the results from each processor.

  MPI_Allreduce(RP_out, RP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(IP_out, IP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //|------= Free memory.

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  delete[] RP_out;
  delete[] IP_out;
  DG_List->clearList();
}
//|----------------------------------------------------
//|
//| ADM mass, linear momentum and angular momentum
//|
//|----------------------------------------------------
void surface_integral::surf_MassPAng(double rex, int lev, cgh *GH, var *chi, var *trK,
                                     var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
                                     var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
                                     var *Gmx, var *Gmy, var *Gmz,
                                     var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs, // temparay memory for mass^i
                                     double *Rout, monitor *Monitor)
{
  if (myrank == 0 && GH->grids[lev] != 1)
    if (Monitor && Monitor->outfile)
      Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
    else
      cout << "WARNING: surface integral on multipatches" << endl;

  double mass, px, py, pz, sx, sy, sz;

  MyList<Patch> *Pp = GH->PatL[lev];
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_admmass_bssn(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                       cg->fgfs[chi->sgfn], cg->fgfs[trK->sgfn],
                       cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                       cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                       cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                       cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                       Symmetry);
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }

  const int InList = 17;

  MyList<var> *DG_List = new MyList<var>(Sfx_rhs);
  DG_List->insert(Sfy_rhs);
  DG_List->insert(Sfz_rhs);
  DG_List->insert(chi);
  DG_List->insert(trK);
  DG_List->insert(gxx);
  DG_List->insert(gxy);
  DG_List->insert(gxz);
  DG_List->insert(gyy);
  DG_List->insert(gyz);
  DG_List->insert(gzz);
  DG_List->insert(Axx);
  DG_List->insert(Axy);
  DG_List->insert(Axz);
  DG_List->insert(Ayy);
  DG_List->insert(Ayz);
  DG_List->insert(Azz);

  int n;
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  // we have assumed there is only one box on this level,
  // so we do not need loop boxes
  GH->PatL[lev]->data->Interp_Points(DG_List, n_tot, pox, shellf, Symmetry);

  double Mass_out = 0;
  double ang_outx, ang_outy, ang_outz;
  double p_outx, p_outy, p_outz;
  ang_outx = ang_outy = ang_outz = 0.0;
  p_outx = p_outy = p_outz = 0.0;
  const double f1o8 = 0.125;

  int mp, Lp, Nmin, Nmax;

  mp = n_tot / cpusize;
  Lp = n_tot - cpusize * mp;

  if (Lp > myrank)
  {
    Nmin = myrank * mp + myrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = myrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  double Chi, Psi;
  double Gxx, Gxy, Gxz, Gyy, Gyz, Gzz;
  double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
  double TRK, axx, axy, axz, ayy, ayz, azz;
  double aupxx, aupxy, aupxz, aupyx, aupyy, aupyz, aupzx, aupzy, aupzz;
  int i;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1

    Chi = shellf[InList * n + 3]; // chi in fact
    TRK = shellf[InList * n + 4];
    Gxx = shellf[InList * n + 5] + 1.0;
    Gxy = shellf[InList * n + 6];
    Gxz = shellf[InList * n + 7];
    Gyy = shellf[InList * n + 8] + 1.0;
    Gyz = shellf[InList * n + 9];
    Gzz = shellf[InList * n + 10] + 1.0;
    axx = shellf[InList * n + 11];
    axy = shellf[InList * n + 12];
    axz = shellf[InList * n + 13];
    ayy = shellf[InList * n + 14];
    ayz = shellf[InList * n + 15];
    azz = shellf[InList * n + 16];

    Chi = 1.0 / (1.0 + Chi); // exp(4*phi)
    Psi = Chi * sqrt(Chi);   // Psi^6

// Chi^2 corresponds to metric determinant
// but this factor has been considered in f_admmass_bssn
#ifdef GaussInt
    // wtcostheta is even function respect costheta
    Mass_out = Mass_out + (shellf[InList * n] * nx_g[n] + shellf[InList * n + 1] * ny_g[n] + shellf[InList * n + 2] * nz_g[n]) * wtcostheta[i];
#else
    Mass_out = Mass_out + (shellf[InList * n] * nx_g[n] + shellf[InList * n + 1] * ny_g[n] + shellf[InList * n + 2] * nz_g[n]);
#endif

    gupzz = Gxx * Gyy * Gzz + Gxy * Gyz * Gxz + Gxz * Gxy * Gyz -
            Gxz * Gyy * Gxz - Gxy * Gxy * Gzz - Gxx * Gyz * Gyz;
    gupxx = (Gyy * Gzz - Gyz * Gyz) / gupzz;
    gupxy = -(Gxy * Gzz - Gyz * Gxz) / gupzz;
    gupxz = (Gxy * Gyz - Gyy * Gxz) / gupzz;
    gupyy = (Gxx * Gzz - Gxz * Gxz) / gupzz;
    gupyz = -(Gxx * Gyz - Gxy * Gxz) / gupzz;
    gupzz = (Gxx * Gyy - Gxy * Gxy) / gupzz;

    aupxx = gupxx * axx + gupxy * axy + gupxz * axz;
    aupxy = gupxx * axy + gupxy * ayy + gupxz * ayz;
    aupxz = gupxx * axz + gupxy * ayz + gupxz * azz;
    aupyx = gupxy * axx + gupyy * axy + gupyz * axz;
    aupyy = gupxy * axy + gupyy * ayy + gupyz * ayz;
    aupyz = gupxy * axz + gupyy * ayz + gupyz * azz;
    aupzx = gupxz * axx + gupyz * axy + gupzz * axz;
    aupzy = gupxz * axy + gupyz * ayy + gupzz * ayz;
    aupzz = gupxz * axz + gupyz * ayz + gupzz * azz;
    if (Symmetry == 0)
    {
#ifdef GaussInt
      // wtcostheta is even function respect costheta
      //  1/8\pi \int \psi^6 (y A^m_z - zA^m_y) dS_m
      ang_outx = ang_outx + f1o8 * Psi * (nx_g[n] * (pox[1][n] * aupxz - pox[2][n] * aupxy) + ny_g[n] * (pox[1][n] * aupyz - pox[2][n] * aupyy) + nz_g[n] * (pox[1][n] * aupzz - pox[2][n] * aupzy)) * wtcostheta[i];
      //  1/8\pi \int \psi^6 (z A^m_x - xA^m_z) dS_m
      ang_outy = ang_outy + f1o8 * Psi * (nx_g[n] * (pox[2][n] * aupxx - pox[0][n] * aupxz) + ny_g[n] * (pox[2][n] * aupyx - pox[0][n] * aupyz) + nz_g[n] * (pox[2][n] * aupzx - pox[0][n] * aupzz)) * wtcostheta[i];
      // 1/8\pi \int \psi^6 (x A^m_y - yA^m_x) dS_m
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx)) * wtcostheta[i];
#else
      //  1/8\pi \int \psi^6 (y A^m_z - zA^m_y) dS_m
      ang_outx = ang_outx + f1o8 * Psi * (nx_g[n] * (pox[1][n] * aupxz - pox[2][n] * aupxy) + ny_g[n] * (pox[1][n] * aupyz - pox[2][n] * aupyy) + nz_g[n] * (pox[1][n] * aupzz - pox[2][n] * aupzy));
      //  1/8\pi \int \psi^6 (z A^m_x - xA^m_z) dS_m
      ang_outy = ang_outy + f1o8 * Psi * (nx_g[n] * (pox[2][n] * aupxx - pox[0][n] * aupxz) + ny_g[n] * (pox[2][n] * aupyx - pox[0][n] * aupyz) + nz_g[n] * (pox[2][n] * aupzx - pox[0][n] * aupzz));
      // 1/8\pi \int \psi^6 (x A^m_y - yA^m_x) dS_m
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx));
#endif
    }
    else if (Symmetry == 1)
    {
#ifdef GaussInt
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx)) * wtcostheta[i];
#else
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx));
#endif
    }

    axx = Chi * (axx + Gxx * TRK / 3.0);
    axy = Chi * (axy + Gxy * TRK / 3.0);
    axz = Chi * (axz + Gxz * TRK / 3.0);
    ayy = Chi * (ayy + Gyy * TRK / 3.0);
    ayz = Chi * (ayz + Gyz * TRK / 3.0);
    azz = Chi * (azz + Gzz * TRK / 3.0);

    axx = axx - TRK;
    ayy = ayy - TRK;
    azz = azz - TRK;

    // 1/8\pi \int \psi^6 (K_mi - \delta_mi trK) dS^m: lower index linear momentum
    if (Symmetry == 0)
    {
#ifdef GaussInt
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz) * wtcostheta[i];
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz) * wtcostheta[i];
      p_outz = p_outz + f1o8 * Psi * (nx_g[n] * axz + ny_g[n] * ayz + nz_g[n] * azz) * wtcostheta[i];
#else
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz);
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz);
      p_outz = p_outz + f1o8 * Psi * (nx_g[n] * axz + ny_g[n] * ayz + nz_g[n] * azz);
#endif
    }
    else if (Symmetry == 1)
    {
#ifdef GaussInt
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz) * wtcostheta[i];
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz) * wtcostheta[i];
#else
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz);
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz);
#endif
    }
  }

  MPI_Allreduce(&Mass_out, &mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&ang_outx, &sx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ang_outy, &sy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ang_outz, &sz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&p_outx, &px, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&p_outy, &py, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&p_outz, &pz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#ifdef GaussInt
  mass = mass * rex * rex * dphi * factor;

  sx = sx * rex * rex * dphi * (1.0 / PI) * factor;
  sy = sy * rex * rex * dphi * (1.0 / PI) * factor;
  sz = sz * rex * rex * dphi * (1.0 / PI) * factor;

  px = px * rex * rex * dphi * (1.0 / PI) * factor;
  py = py * rex * rex * dphi * (1.0 / PI) * factor;
  pz = pz * rex * rex * dphi * (1.0 / PI) * factor;
#else
  mass = mass * rex * rex * dphi * dcostheta * factor;

  sx = sx * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  sy = sy * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  sz = sz * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;

  px = px * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  py = py * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  pz = pz * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
#endif

  Rout[0] = mass;
  Rout[1] = px;
  Rout[2] = py;
  Rout[3] = pz;
  Rout[4] = sx;
  Rout[5] = sy;
  Rout[6] = sz;

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  DG_List->clearList();
}
void surface_integral::surf_MassPAng(double rex, int lev, cgh *GH, var *chi, var *trK,
                                     var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
                                     var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
                                     var *Gmx, var *Gmy, var *Gmz,
                                     var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs, // temparay memory for mass^i
                                     double *Rout, monitor *Monitor, MPI_Comm Comm_here)
{
  int lmyrank;
  MPI_Comm_rank(Comm_here, &lmyrank);
  if (lmyrank == 0 && GH->grids[lev] != 1)
    if (Monitor && Monitor->outfile)
      Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
    else
      cout << "WARNING: surface integral on multipatches" << endl;

  double mass, px, py, pz, sx, sy, sz;

  MyList<Patch> *Pp = GH->PatL[lev];
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_admmass_bssn(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                       cg->fgfs[chi->sgfn], cg->fgfs[trK->sgfn],
                       cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                       cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                       cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                       cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                       Symmetry);
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }

  const int InList = 17;

  MyList<var> *DG_List = new MyList<var>(Sfx_rhs);
  DG_List->insert(Sfy_rhs);
  DG_List->insert(Sfz_rhs);
  DG_List->insert(chi);
  DG_List->insert(trK);
  DG_List->insert(gxx);
  DG_List->insert(gxy);
  DG_List->insert(gxz);
  DG_List->insert(gyy);
  DG_List->insert(gyz);
  DG_List->insert(gzz);
  DG_List->insert(Axx);
  DG_List->insert(Axy);
  DG_List->insert(Axz);
  DG_List->insert(Ayy);
  DG_List->insert(Ayz);
  DG_List->insert(Azz);

  int n;
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  // we have assumed there is only one box on this level,
  // so we do not need loop boxes
  GH->PatL[lev]->data->Interp_Points(DG_List, n_tot, pox, shellf, Symmetry, Comm_here);

  double Mass_out = 0;
  double ang_outx, ang_outy, ang_outz;
  double p_outx, p_outy, p_outz;
  ang_outx = ang_outy = ang_outz = 0.0;
  p_outx = p_outy = p_outz = 0.0;
  const double f1o8 = 0.125;

  int mp, Lp, Nmin, Nmax;

  int cpusize_here;
  MPI_Comm_size(Comm_here, &cpusize_here);

  mp = n_tot / cpusize_here;
  Lp = n_tot - cpusize_here * mp;

  if (Lp > lmyrank)
  {
    Nmin = lmyrank * mp + lmyrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = lmyrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  double Chi, Psi;
  double Gxx, Gxy, Gxz, Gyy, Gyz, Gzz;
  double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
  double TRK, axx, axy, axz, ayy, ayz, azz;
  double aupxx, aupxy, aupxz, aupyx, aupyy, aupyz, aupzx, aupzy, aupzz;
  int i;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1

    Chi = shellf[InList * n + 3]; // chi in fact
    TRK = shellf[InList * n + 4];
    Gxx = shellf[InList * n + 5] + 1.0;
    Gxy = shellf[InList * n + 6];
    Gxz = shellf[InList * n + 7];
    Gyy = shellf[InList * n + 8] + 1.0;
    Gyz = shellf[InList * n + 9];
    Gzz = shellf[InList * n + 10] + 1.0;
    axx = shellf[InList * n + 11];
    axy = shellf[InList * n + 12];
    axz = shellf[InList * n + 13];
    ayy = shellf[InList * n + 14];
    ayz = shellf[InList * n + 15];
    azz = shellf[InList * n + 16];

    Chi = 1.0 / (1.0 + Chi); // exp(4*phi)
    Psi = Chi * sqrt(Chi);   // Psi^6

// Chi^2 corresponds to metric determinant
// but this factor has been considered in f_admmass_bssn
#ifdef GaussInt
    // wtcostheta is even function respect costheta
    Mass_out = Mass_out + (shellf[InList * n] * nx_g[n] + shellf[InList * n + 1] * ny_g[n] + shellf[InList * n + 2] * nz_g[n]) * wtcostheta[i];
#else
    Mass_out = Mass_out + (shellf[InList * n] * nx_g[n] + shellf[InList * n + 1] * ny_g[n] + shellf[InList * n + 2] * nz_g[n]);
#endif

    gupzz = Gxx * Gyy * Gzz + Gxy * Gyz * Gxz + Gxz * Gxy * Gyz -
            Gxz * Gyy * Gxz - Gxy * Gxy * Gzz - Gxx * Gyz * Gyz;
    gupxx = (Gyy * Gzz - Gyz * Gyz) / gupzz;
    gupxy = -(Gxy * Gzz - Gyz * Gxz) / gupzz;
    gupxz = (Gxy * Gyz - Gyy * Gxz) / gupzz;
    gupyy = (Gxx * Gzz - Gxz * Gxz) / gupzz;
    gupyz = -(Gxx * Gyz - Gxy * Gxz) / gupzz;
    gupzz = (Gxx * Gyy - Gxy * Gxy) / gupzz;

    aupxx = gupxx * axx + gupxy * axy + gupxz * axz;
    aupxy = gupxx * axy + gupxy * ayy + gupxz * ayz;
    aupxz = gupxx * axz + gupxy * ayz + gupxz * azz;
    aupyx = gupxy * axx + gupyy * axy + gupyz * axz;
    aupyy = gupxy * axy + gupyy * ayy + gupyz * ayz;
    aupyz = gupxy * axz + gupyy * ayz + gupyz * azz;
    aupzx = gupxz * axx + gupyz * axy + gupzz * axz;
    aupzy = gupxz * axy + gupyz * ayy + gupzz * ayz;
    aupzz = gupxz * axz + gupyz * ayz + gupzz * azz;
    if (Symmetry == 0)
    {
#ifdef GaussInt
      // wtcostheta is even function respect costheta
      //  1/8\pi \int \psi^6 (y A^m_z - zA^m_y) dS_m
      ang_outx = ang_outx + f1o8 * Psi * (nx_g[n] * (pox[1][n] * aupxz - pox[2][n] * aupxy) + ny_g[n] * (pox[1][n] * aupyz - pox[2][n] * aupyy) + nz_g[n] * (pox[1][n] * aupzz - pox[2][n] * aupzy)) * wtcostheta[i];
      //  1/8\pi \int \psi^6 (z A^m_x - xA^m_z) dS_m
      ang_outy = ang_outy + f1o8 * Psi * (nx_g[n] * (pox[2][n] * aupxx - pox[0][n] * aupxz) + ny_g[n] * (pox[2][n] * aupyx - pox[0][n] * aupyz) + nz_g[n] * (pox[2][n] * aupzx - pox[0][n] * aupzz)) * wtcostheta[i];
      // 1/8\pi \int \psi^6 (x A^m_y - yA^m_x) dS_m
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx)) * wtcostheta[i];
#else
      //  1/8\pi \int \psi^6 (y A^m_z - zA^m_y) dS_m
      ang_outx = ang_outx + f1o8 * Psi * (nx_g[n] * (pox[1][n] * aupxz - pox[2][n] * aupxy) + ny_g[n] * (pox[1][n] * aupyz - pox[2][n] * aupyy) + nz_g[n] * (pox[1][n] * aupzz - pox[2][n] * aupzy));
      //  1/8\pi \int \psi^6 (z A^m_x - xA^m_z) dS_m
      ang_outy = ang_outy + f1o8 * Psi * (nx_g[n] * (pox[2][n] * aupxx - pox[0][n] * aupxz) + ny_g[n] * (pox[2][n] * aupyx - pox[0][n] * aupyz) + nz_g[n] * (pox[2][n] * aupzx - pox[0][n] * aupzz));
      // 1/8\pi \int \psi^6 (x A^m_y - yA^m_x) dS_m
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx));
#endif
    }
    else if (Symmetry == 1)
    {
#ifdef GaussInt
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx)) * wtcostheta[i];
#else
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx));
#endif
    }

    axx = Chi * (axx + Gxx * TRK / 3.0);
    axy = Chi * (axy + Gxy * TRK / 3.0);
    axz = Chi * (axz + Gxz * TRK / 3.0);
    ayy = Chi * (ayy + Gyy * TRK / 3.0);
    ayz = Chi * (ayz + Gyz * TRK / 3.0);
    azz = Chi * (azz + Gzz * TRK / 3.0);

    axx = axx - TRK;
    ayy = ayy - TRK;
    azz = azz - TRK;

    // 1/8\pi \int \psi^6 (K_mi - \delta_mi trK) dS^m: lower index linear momentum
    if (Symmetry == 0)
    {
#ifdef GaussInt
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz) * wtcostheta[i];
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz) * wtcostheta[i];
      p_outz = p_outz + f1o8 * Psi * (nx_g[n] * axz + ny_g[n] * ayz + nz_g[n] * azz) * wtcostheta[i];
#else
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz);
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz);
      p_outz = p_outz + f1o8 * Psi * (nx_g[n] * axz + ny_g[n] * ayz + nz_g[n] * azz);
#endif
    }
    else if (Symmetry == 1)
    {
#ifdef GaussInt
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz) * wtcostheta[i];
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz) * wtcostheta[i];
#else
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz);
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz);
#endif
    }
  }

  MPI_Allreduce(&Mass_out, &mass, 1, MPI_DOUBLE, MPI_SUM, Comm_here);

  MPI_Allreduce(&ang_outx, &sx, 1, MPI_DOUBLE, MPI_SUM, Comm_here);
  MPI_Allreduce(&ang_outy, &sy, 1, MPI_DOUBLE, MPI_SUM, Comm_here);
  MPI_Allreduce(&ang_outz, &sz, 1, MPI_DOUBLE, MPI_SUM, Comm_here);

  MPI_Allreduce(&p_outx, &px, 1, MPI_DOUBLE, MPI_SUM, Comm_here);
  MPI_Allreduce(&p_outy, &py, 1, MPI_DOUBLE, MPI_SUM, Comm_here);
  MPI_Allreduce(&p_outz, &pz, 1, MPI_DOUBLE, MPI_SUM, Comm_here);

#ifdef GaussInt
  mass = mass * rex * rex * dphi * factor;

  sx = sx * rex * rex * dphi * (1.0 / PI) * factor;
  sy = sy * rex * rex * dphi * (1.0 / PI) * factor;
  sz = sz * rex * rex * dphi * (1.0 / PI) * factor;

  px = px * rex * rex * dphi * (1.0 / PI) * factor;
  py = py * rex * rex * dphi * (1.0 / PI) * factor;
  pz = pz * rex * rex * dphi * (1.0 / PI) * factor;
#else
  mass = mass * rex * rex * dphi * dcostheta * factor;

  sx = sx * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  sy = sy * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  sz = sz * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;

  px = px * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  py = py * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  pz = pz * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
#endif

  Rout[0] = mass;
  Rout[1] = px;
  Rout[2] = py;
  Rout[3] = pz;
  Rout[4] = sx;
  Rout[5] = sy;
  Rout[6] = sz;

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  DG_List->clearList();
}
//|----------------------------------------------------------------
//  for shell patch
//|----------------------------------------------------------------
void surface_integral::surf_MassPAng(double rex, int lev, ShellPatch *GH, var *chi, var *trK,
                                     var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
                                     var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
                                     var *Gmx, var *Gmy, var *Gmz,
                                     var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs, // temparay memory for mass^i
                                     double *Rout, monitor *Monitor)
{
  if (lev != 0)
  {
    if (myrank == 0)
    {
      if (Monitor && Monitor->outfile)
        Monitor->outfile << "WARNING: shell surface integral not on level 0" << endl;
      else
        cout << "WARNING: shell surface integral not on level 0" << endl;
    }
    return;
  }

  double mass, px, py, pz, sx, sy, sz;

  MyList<ss_patch> *Pp = GH->PatL;
  while (Pp)
  {
    MyList<Block> *BL = Pp->data->blb;
    int fngfs = Pp->data->fngfs;
    while (BL)
    {
      Block *cg = BL->data;
      if (myrank == cg->rank)
      {
        f_admmass_bssn_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                          cg->fgfs[fngfs + ShellPatch::gx], cg->fgfs[fngfs + ShellPatch::gy], cg->fgfs[fngfs + ShellPatch::gz],
                          cg->fgfs[fngfs + ShellPatch::drhodx], cg->fgfs[fngfs + ShellPatch::drhody], cg->fgfs[fngfs + ShellPatch::drhodz],
                          cg->fgfs[fngfs + ShellPatch::dsigmadx], cg->fgfs[fngfs + ShellPatch::dsigmady], cg->fgfs[fngfs + ShellPatch::dsigmadz],
                          cg->fgfs[fngfs + ShellPatch::dRdx], cg->fgfs[fngfs + ShellPatch::dRdy], cg->fgfs[fngfs + ShellPatch::dRdz],
                          cg->fgfs[fngfs + ShellPatch::drhodxx], cg->fgfs[fngfs + ShellPatch::drhodxy], cg->fgfs[fngfs + ShellPatch::drhodxz],
                          cg->fgfs[fngfs + ShellPatch::drhodyy], cg->fgfs[fngfs + ShellPatch::drhodyz], cg->fgfs[fngfs + ShellPatch::drhodzz],
                          cg->fgfs[fngfs + ShellPatch::dsigmadxx], cg->fgfs[fngfs + ShellPatch::dsigmadxy], cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                          cg->fgfs[fngfs + ShellPatch::dsigmadyy], cg->fgfs[fngfs + ShellPatch::dsigmadyz], cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                          cg->fgfs[fngfs + ShellPatch::dRdxx], cg->fgfs[fngfs + ShellPatch::dRdxy], cg->fgfs[fngfs + ShellPatch::dRdxz],
                          cg->fgfs[fngfs + ShellPatch::dRdyy], cg->fgfs[fngfs + ShellPatch::dRdyz], cg->fgfs[fngfs + ShellPatch::dRdzz],
                          cg->fgfs[chi->sgfn], cg->fgfs[trK->sgfn],
                          cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                          cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                          cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                          cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                          Symmetry, Pp->data->sst);
      }
      if (BL == Pp->data->ble)
        break;
      BL = BL->next;
    }
    Pp = Pp->next;
  }

  const int InList = 17;

  MyList<var> *DG_List = new MyList<var>(Sfx_rhs);
  DG_List->insert(Sfy_rhs);
  DG_List->insert(Sfz_rhs);
  DG_List->insert(chi);
  DG_List->insert(trK);
  DG_List->insert(gxx);
  DG_List->insert(gxy);
  DG_List->insert(gxz);
  DG_List->insert(gyy);
  DG_List->insert(gyz);
  DG_List->insert(gzz);
  DG_List->insert(Axx);
  DG_List->insert(Axy);
  DG_List->insert(Axz);
  DG_List->insert(Ayy);
  DG_List->insert(Ayz);
  DG_List->insert(Azz);

  int n;
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  // we have assumed there is only one box on this level,
  // so we do not need loop boxes
  GH->Interp_Points(DG_List, n_tot, pox, shellf, Symmetry);

  double Mass_out = 0;
  double ang_outx, ang_outy, ang_outz;
  double p_outx, p_outy, p_outz;
  ang_outx = ang_outy = ang_outz = 0.0;
  p_outx = p_outy = p_outz = 0.0;
  const double f1o8 = 0.125;

  int mp, Lp, Nmin, Nmax;

  mp = n_tot / cpusize;
  Lp = n_tot - cpusize * mp;

  if (Lp > myrank)
  {
    Nmin = myrank * mp + myrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = myrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  double Chi, Psi;
  double Gxx, Gxy, Gxz, Gyy, Gyz, Gzz;
  double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
  double TRK, axx, axy, axz, ayy, ayz, azz;
  double aupxx, aupxy, aupxz, aupyx, aupyy, aupyz, aupzx, aupzy, aupzz;
  int i;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1

    Chi = shellf[InList * n + 3]; // chi in fact
    TRK = shellf[InList * n + 4];
    Gxx = shellf[InList * n + 5] + 1.0;
    Gxy = shellf[InList * n + 6];
    Gxz = shellf[InList * n + 7];
    Gyy = shellf[InList * n + 8] + 1.0;
    Gyz = shellf[InList * n + 9];
    Gzz = shellf[InList * n + 10] + 1.0;
    axx = shellf[InList * n + 11];
    axy = shellf[InList * n + 12];
    axz = shellf[InList * n + 13];
    ayy = shellf[InList * n + 14];
    ayz = shellf[InList * n + 15];
    azz = shellf[InList * n + 16];

    Chi = 1.0 / (1.0 + Chi); // exp(4*phi)
    Psi = Chi * sqrt(Chi);   // Psi^6
// Chi^2 corresponds to metric determinant
// but this factor has been considered in f_admmass_bssn
#ifdef GaussInt
    // wtcostheta is even function respect costheta
    Mass_out = Mass_out + (shellf[InList * n] * nx_g[n] + shellf[InList * n + 1] * ny_g[n] + shellf[InList * n + 2] * nz_g[n]) * wtcostheta[i];
#else
    Mass_out = Mass_out + (shellf[InList * n] * nx_g[n] + shellf[InList * n + 1] * ny_g[n] + shellf[InList * n + 2] * nz_g[n]);
#endif

    gupzz = Gxx * Gyy * Gzz + Gxy * Gyz * Gxz + Gxz * Gxy * Gyz -
            Gxz * Gyy * Gxz - Gxy * Gxy * Gzz - Gxx * Gyz * Gyz;
    gupxx = (Gyy * Gzz - Gyz * Gyz) / gupzz;
    gupxy = -(Gxy * Gzz - Gyz * Gxz) / gupzz;
    gupxz = (Gxy * Gyz - Gyy * Gxz) / gupzz;
    gupyy = (Gxx * Gzz - Gxz * Gxz) / gupzz;
    gupyz = -(Gxx * Gyz - Gxy * Gxz) / gupzz;
    gupzz = (Gxx * Gyy - Gxy * Gxy) / gupzz;

    aupxx = gupxx * axx + gupxy * axy + gupxz * axz;
    aupxy = gupxx * axy + gupxy * ayy + gupxz * ayz;
    aupxz = gupxx * axz + gupxy * ayz + gupxz * azz;
    aupyx = gupxy * axx + gupyy * axy + gupyz * axz;
    aupyy = gupxy * axy + gupyy * ayy + gupyz * ayz;
    aupyz = gupxy * axz + gupyy * ayz + gupyz * azz;
    aupzx = gupxz * axx + gupyz * axy + gupzz * axz;
    aupzy = gupxz * axy + gupyz * ayy + gupzz * ayz;
    aupzz = gupxz * axz + gupyz * ayz + gupzz * azz;
    if (Symmetry == 0)
    {
#ifdef GaussInt
      // wtcostheta is even function respect costheta
      //  1/8\pi \int \psi^6 (y A^m_z - zA^m_y) dS_m
      ang_outx = ang_outx + f1o8 * Psi * (nx_g[n] * (pox[1][n] * aupxz - pox[2][n] * aupxy) + ny_g[n] * (pox[1][n] * aupyz - pox[2][n] * aupyy) + nz_g[n] * (pox[1][n] * aupzz - pox[2][n] * aupzy)) * wtcostheta[i];
      //  1/8\pi \int \psi^6 (z A^m_x - xA^m_z) dS_m
      ang_outy = ang_outy + f1o8 * Psi * (nx_g[n] * (pox[2][n] * aupxx - pox[0][n] * aupxz) + ny_g[n] * (pox[2][n] * aupyx - pox[0][n] * aupyz) + nz_g[n] * (pox[2][n] * aupzx - pox[0][n] * aupzz)) * wtcostheta[i];
      // 1/8\pi \int \psi^6 (x A^m_y - yA^m_x) dS_m
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx)) * wtcostheta[i];
#else
      //  1/8\pi \int \psi^6 (y A^m_z - zA^m_y) dS_m
      ang_outx = ang_outx + f1o8 * Psi * (nx_g[n] * (pox[1][n] * aupxz - pox[2][n] * aupxy) + ny_g[n] * (pox[1][n] * aupyz - pox[2][n] * aupyy) + nz_g[n] * (pox[1][n] * aupzz - pox[2][n] * aupzy));
      //  1/8\pi \int \psi^6 (z A^m_x - xA^m_z) dS_m
      ang_outy = ang_outy + f1o8 * Psi * (nx_g[n] * (pox[2][n] * aupxx - pox[0][n] * aupxz) + ny_g[n] * (pox[2][n] * aupyx - pox[0][n] * aupyz) + nz_g[n] * (pox[2][n] * aupzx - pox[0][n] * aupzz));
      // 1/8\pi \int \psi^6 (x A^m_y - yA^m_x) dS_m
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx));
#endif
    }
    else if (Symmetry == 1)
    {
#ifdef GaussInt
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx)) * wtcostheta[i];
#else
      ang_outz = ang_outz + f1o8 * Psi * (nx_g[n] * (pox[0][n] * aupxy - pox[1][n] * aupxx) + ny_g[n] * (pox[0][n] * aupyy - pox[1][n] * aupyx) + nz_g[n] * (pox[0][n] * aupzy - pox[1][n] * aupzx));
#endif
    }

    axx = Chi * (axx + Gxx * TRK / 3.0);
    axy = Chi * (axy + Gxy * TRK / 3.0);
    axz = Chi * (axz + Gxz * TRK / 3.0);
    ayy = Chi * (ayy + Gyy * TRK / 3.0);
    ayz = Chi * (ayz + Gyz * TRK / 3.0);
    azz = Chi * (azz + Gzz * TRK / 3.0);

    axx = axx - TRK;
    ayy = ayy - TRK;
    azz = azz - TRK;

    // 1/8\pi \int \psi^6 (K_mi - \delta_mi trK) dS^m: lower index linear momentum
    if (Symmetry == 0)
    {
#ifdef GaussInt
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz) * wtcostheta[i];
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz) * wtcostheta[i];
      p_outz = p_outz + f1o8 * Psi * (nx_g[n] * axz + ny_g[n] * ayz + nz_g[n] * azz) * wtcostheta[i];
#else
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz);
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz);
      p_outz = p_outz + f1o8 * Psi * (nx_g[n] * axz + ny_g[n] * ayz + nz_g[n] * azz);
#endif
    }
    else if (Symmetry == 1)
    {
#ifdef GaussInt
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz) * wtcostheta[i];
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz) * wtcostheta[i];
#else
      p_outx = p_outx + f1o8 * Psi * (nx_g[n] * axx + ny_g[n] * axy + nz_g[n] * axz);
      p_outy = p_outy + f1o8 * Psi * (nx_g[n] * axy + ny_g[n] * ayy + nz_g[n] * ayz);
#endif
    }
  }

  MPI_Allreduce(&Mass_out, &mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&ang_outx, &sx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ang_outy, &sy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&ang_outz, &sz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&p_outx, &px, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&p_outy, &py, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&p_outz, &pz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#ifdef GaussInt
  mass = mass * rex * rex * dphi * factor;

  sx = sx * rex * rex * dphi * (1.0 / PI) * factor;
  sy = sy * rex * rex * dphi * (1.0 / PI) * factor;
  sz = sz * rex * rex * dphi * (1.0 / PI) * factor;

  px = px * rex * rex * dphi * (1.0 / PI) * factor;
  py = py * rex * rex * dphi * (1.0 / PI) * factor;
  pz = pz * rex * rex * dphi * (1.0 / PI) * factor;
#else
  mass = mass * rex * rex * dphi * dcostheta * factor;

  sx = sx * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  sy = sy * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  sz = sz * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;

  px = px * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  py = py * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
  pz = pz * rex * rex * dphi * dcostheta * (1.0 / PI) * factor;
#endif

  Rout[0] = mass;
  Rout[1] = px;
  Rout[2] = py;
  Rout[3] = pz;
  Rout[4] = sx;
  Rout[5] = sy;
  Rout[6] = sz;

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  DG_List->clearList();
}
//|----------------------------------------------------------------
//  do not discriminate box and shell
//  for Gravitational wave specially symmetric case
//|----------------------------------------------------------------
void surface_integral::surf_Wave(double rex, cgh *GH, ShellPatch *SH,
                                 var *chi, var *trK,
                                 var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
                                 var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
                                 var *chix, var *chiy, var *chiz,
                                 var *trKx, var *trKy, var *trKz,
                                 var *Axxx, var *Axxy, var *Axxz,
                                 var *Axyx, var *Axyy, var *Axyz,
                                 var *Axzx, var *Axzy, var *Axzz,
                                 var *Ayyx, var *Ayyy, var *Ayyz,
                                 var *Ayzx, var *Ayzy, var *Ayzz,
                                 var *Azzx, var *Azzy, var *Azzz,
                                 var *Gamxxx, var *Gamxxy, var *Gamxxz, var *Gamxyy, var *Gamxyz, var *Gamxzz,
                                 var *Gamyxx, var *Gamyxy, var *Gamyxz, var *Gamyyy, var *Gamyyz, var *Gamyzz,
                                 var *Gamzxx, var *Gamzxy, var *Gamzxz, var *Gamzyy, var *Gamzyz, var *Gamzzz,
                                 var *Rxx, var *Rxy, var *Rxz, var *Ryy, var *Ryz, var *Rzz,
                                 int spinw, int maxl, int NN, double *RP, double *IP,
                                 monitor *Monitor) // NN is the length of RP and IP
{
  const int InList = 62;

  MyList<var> *DG_List = new MyList<var>(chi);
  DG_List->insert(trK);
  DG_List->insert(gxx);
  DG_List->insert(gxy);
  DG_List->insert(gxz);
  DG_List->insert(gyy);
  DG_List->insert(gyz);
  DG_List->insert(gzz);
  DG_List->insert(Axx);
  DG_List->insert(Axy);
  DG_List->insert(Axz);
  DG_List->insert(Ayy);
  DG_List->insert(Ayz);
  DG_List->insert(Azz);
  DG_List->insert(chix);
  DG_List->insert(chiy);
  DG_List->insert(chiz);
  DG_List->insert(trKx);
  DG_List->insert(trKy);
  DG_List->insert(trKz);
  DG_List->insert(Axxx);
  DG_List->insert(Axxy);
  DG_List->insert(Axxz);
  DG_List->insert(Axyx);
  DG_List->insert(Axyy);
  DG_List->insert(Axyz);
  DG_List->insert(Axzx);
  DG_List->insert(Axzy);
  DG_List->insert(Axzz);
  DG_List->insert(Ayyx);
  DG_List->insert(Ayyy);
  DG_List->insert(Ayyz);
  DG_List->insert(Ayzx);
  DG_List->insert(Ayzy);
  DG_List->insert(Ayzz);
  DG_List->insert(Azzx);
  DG_List->insert(Azzy);
  DG_List->insert(Azzz);
  DG_List->insert(Gamxxx);
  DG_List->insert(Gamxxy);
  DG_List->insert(Gamxxz);
  DG_List->insert(Gamxyy);
  DG_List->insert(Gamxyz);
  DG_List->insert(Gamxzz);
  DG_List->insert(Gamyxx);
  DG_List->insert(Gamyxy);
  DG_List->insert(Gamyxz);
  DG_List->insert(Gamyyy);
  DG_List->insert(Gamyyz);
  DG_List->insert(Gamyzz);
  DG_List->insert(Gamzxx);
  DG_List->insert(Gamzxy);
  DG_List->insert(Gamzxz);
  DG_List->insert(Gamzyy);
  DG_List->insert(Gamzyz);
  DG_List->insert(Gamzzz);
  DG_List->insert(Rxx);
  DG_List->insert(Rxy);
  DG_List->insert(Rxz);
  DG_List->insert(Ryy);
  DG_List->insert(Ryz);
  DG_List->insert(Rzz);

  int n;
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[n_tot];
  for (n = 0; n < n_tot; n++)
  {
    pox[0][n] = rex * nx_g[n];
    pox[1][n] = rex * ny_g[n];
    pox[2][n] = rex * nz_g[n];
  }

  double *shellf;
  shellf = new double[n_tot * InList];

  SR_Interp_Points(DG_List, GH, SH, n_tot, pox, shellf);

  double *RP_out, *IP_out;
  RP_out = new double[NN];
  IP_out = new double[NN];

  for (int ii = 0; ii < NN; ii++)
  {
    RP_out[ii] = 0;
    IP_out[ii] = 0;
  }

  int mp, Lp, Nmin, Nmax;

  mp = n_tot / cpusize;
  Lp = n_tot - cpusize * mp;

  if (Lp > myrank)
  {
    Nmin = myrank * mp + myrank;
    Nmax = Nmin + mp;
  }
  else
  {
    Nmin = myrank * mp + Lp;
    Nmax = Nmin + mp - 1;
  }

  // theta part
  double costheta, thetap;
  double cosmphi, sinmphi;

  int i, j;
  int lpsy = 0;
  if (Symmetry == 0)
    lpsy = 1;
  else if (Symmetry == 1)
    lpsy = 2;
  else if (Symmetry == 2)
    lpsy = 8;

  double psi4RR, psi4II;
  double px, py, pz;
  double pchi, ptrK, pgxx, pgxy, pgxz, pgyy, pgyz, pgzz;
  double pAxx, pAxy, pAxz, pAyy, pAyz, pAzz;
  double pchix, pchiy, pchiz;
  double ptrKx, ptrKy, ptrKz;
  double pAxxx, pAxxy, pAxxz;
  double pAxyx, pAxyy, pAxyz;
  double pAxzx, pAxzy, pAxzz;
  double pAyyx, pAyyy, pAyyz;
  double pAyzx, pAyzy, pAyzz;
  double pAzzx, pAzzy, pAzzz;
  double pGamxxx, pGamxxy, pGamxxz, pGamxyy, pGamxyz, pGamxzz;
  double pGamyxx, pGamyxy, pGamyxz, pGamyyy, pGamyyz, pGamyzz;
  double pGamzxx, pGamzxy, pGamzxz, pGamzyy, pGamzyz, pGamzzz;
  double pRxx, pRxy, pRxz, pRyy, pRyz, pRzz;
  for (n = Nmin; n <= Nmax; n++)
  {
    //       need round off always
    i = int(n / N_phi); // int(1.723) = 1, int(-1.732) = -1
    j = n - i * N_phi;

    int countlm = 0;
    for (int pl = spinw; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        for (int lp = 0; lp < lpsy; lp++)
        {
          px = pox[0][n];
          py = pox[1][n];
          pz = pox[2][n];
          pchi = shellf[InList * n];
          ptrK = shellf[InList * n + 1];
          pgxx = shellf[InList * n + 2];
          pgxy = shellf[InList * n + 3];
          pgxz = shellf[InList * n + 4];
          pgyy = shellf[InList * n + 5];
          pgyz = shellf[InList * n + 6];
          pgzz = shellf[InList * n + 7];
          pAxx = shellf[InList * n + 8];
          pAxy = shellf[InList * n + 9];
          pAxz = shellf[InList * n + 10];
          pAyy = shellf[InList * n + 11];
          pAyz = shellf[InList * n + 12];
          pAzz = shellf[InList * n + 13];
          pchix = shellf[InList * n + 14];
          pchiy = shellf[InList * n + 15];
          pchiz = shellf[InList * n + 16];
          ptrKx = shellf[InList * n + 17];
          ptrKy = shellf[InList * n + 18];
          ptrKz = shellf[InList * n + 19];
          pAxxx = shellf[InList * n + 20];
          pAxxy = shellf[InList * n + 21];
          pAxxz = shellf[InList * n + 22];
          pAxyx = shellf[InList * n + 23];
          pAxyy = shellf[InList * n + 24];
          pAxyz = shellf[InList * n + 25];
          pAxzx = shellf[InList * n + 26];
          pAxzy = shellf[InList * n + 27];
          pAxzz = shellf[InList * n + 28];
          pAyyx = shellf[InList * n + 29];
          pAyyy = shellf[InList * n + 30];
          pAyyz = shellf[InList * n + 31];
          pAyzx = shellf[InList * n + 32];
          pAyzy = shellf[InList * n + 33];
          pAyzz = shellf[InList * n + 34];
          pAzzx = shellf[InList * n + 35];
          pAzzy = shellf[InList * n + 36];
          pAzzz = shellf[InList * n + 37];
          pGamxxx = shellf[InList * n + 38];
          pGamxxy = shellf[InList * n + 39];
          pGamxxz = shellf[InList * n + 40];
          pGamxyy = shellf[InList * n + 41];
          pGamxyz = shellf[InList * n + 42];
          pGamxzz = shellf[InList * n + 43];
          pGamyxx = shellf[InList * n + 44];
          pGamyxy = shellf[InList * n + 45];
          pGamyxz = shellf[InList * n + 46];
          pGamyyy = shellf[InList * n + 47];
          pGamyyz = shellf[InList * n + 48];
          pGamyzz = shellf[InList * n + 49];
          pGamzxx = shellf[InList * n + 50];
          pGamzxy = shellf[InList * n + 51];
          pGamzxz = shellf[InList * n + 52];
          pGamzyy = shellf[InList * n + 53];
          pGamzyz = shellf[InList * n + 54];
          pGamzzz = shellf[InList * n + 55];
          pRxx = shellf[InList * n + 56];
          pRxy = shellf[InList * n + 57];
          pRxz = shellf[InList * n + 58];
          pRyy = shellf[InList * n + 59];
          pRyz = shellf[InList * n + 60];
          pRzz = shellf[InList * n + 61];
          switch (lp)
          {
          case 0: //+++ (theta, phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            break;
          case 1: //++- (pi-theta, phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = sin(pm * (j + 0.5) * dphi);
            pz = -pz;
            pgxz = -pgxz;
            pgyz = -pgyz;
            pAxz = -pAxz;
            pAyz = -pAyz;
            pchiz = -pchiz;
            ptrKz = -ptrKz;
            pAxxz = -pAxxz;
            pAxyz = -pAxyz;
            pAxzx = -pAxzx;
            pAxzy = -pAxzy;
            pAyyz = -pAyyz;
            pAyzx = -pAyzx;
            pAyzy = -pAyzy;
            pAzzz = -pAzzz;
            pGamxxz = -pGamxxz;
            pGamxyz = -pGamxyz;
            pGamyxz = -pGamyxz;
            pGamyyz = -pGamyyz;
            pGamzxx = -pGamzxx;
            pGamzxy = -pGamzxy;
            pGamzyy = -pGamzyy;
            pGamzzz = -pGamzzz;
            pRxz = -pRxz;
            pRyz = -pRyz;
            break;
          case 2: //+-+ (theta, 2*pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            py = -py;
            pgxy = -pgxy;
            pgyz = -pgyz;
            pAxy = -pAxy;
            pAyz = -pAyz;
            pchiy = -pchiy;
            ptrKy = -ptrKy;
            pAxxy = -pAxxy;
            pAxyx = -pAxyx;
            pAxyz = -pAxyz;
            pAxzy = -pAxzy;
            pAyyy = -pAyyy;
            pAyzx = -pAyzx;
            pAyzz = -pAyzz;
            pAzzy = -pAzzy;
            pGamxxy = -pGamxxy;
            pGamxyz = -pGamxyz;
            pGamyxx = -pGamyxx;
            pGamyxz = -pGamyxz;
            pGamyyy = -pGamyyy;
            pGamyzz = -pGamyzz;
            pGamzxy = -pGamzxy;
            pGamzyz = -pGamzyz;
            pRxy = -pRxy;
            pRyz = -pRyz;
            break;
          case 3: //+-- (pi-theta, 2*pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (j + 0.5) * dphi);
            sinmphi = -sin(pm * (j + 0.5) * dphi);
            py = -py;
            pz = -pz;
            pgxy = -pgxy;
            pgxz = -pgxz;
            pAxy = -pAxy;
            pAxz = -pAxz;
            pchiy = -pchiy;
            pchiz = -pchiz;
            ptrKy = -ptrKy;
            ptrKz = -ptrKz;
            pAxxy = -pAxxy;
            pAxxz = -pAxxz;
            pAxyx = -pAxyx;
            pAxzx = -pAxzx;
            pAyyy = -pAyyy;
            pAyyz = -pAyyz;
            pAyzy = -pAyzy;
            pAyzz = -pAyzz;
            pAzzy = -pAzzy;
            pAzzz = -pAzzz;
            pGamxxy = -pGamxxy;
            pGamxxz = -pGamxxz;
            pGamyxx = -pGamyxx;
            pGamyyy = -pGamyyy;
            pGamyyz = -pGamyyz;
            pGamyzz = -pGamyzz;
            pGamzxx = -pGamzxx;
            pGamzyy = -pGamzyy;
            pGamzyz = -pGamzyz;
            pGamzzz = -pGamzzz;
            pRxy = -pRxy;
            pRxz = -pRxz;
            break;
          case 4: //-++ (theta, pi-phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            px = -px;
            pgxy = -pgxy;
            pgxz = -pgxz;
            pAxy = -pAxy;
            pAxz = -pAxz;
            pchix = -pchix;
            ptrKx = -ptrKx;
            pAxxx = -pAxxx;
            pAxyy = -pAxyy;
            pAxyz = -pAxyz;
            pAxzy = -pAxzy;
            pAxzz = -pAxzz;
            pAyyx = -pAyyx;
            pAyzx = -pAyzx;
            pAzzx = -pAzzx;
            pGamxxx = -pGamxxx;
            pGamxyy = -pGamxyy;
            pGamxyz = -pGamxyz;
            pGamxzz = -pGamxzz;
            pGamyxy = -pGamyxy;
            pGamyxz = -pGamyxz;
            pGamzxy = -pGamzxy;
            pGamzxz = -pGamzxz;
            pRxy = -pRxy;
            pRxz = -pRxz;
            break;
          case 5: //-+- (pi-theta, pi-phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI - (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI - (j + 0.5) * dphi));
            px = -px;
            pz = -pz;
            pgxy = -pgxy;
            pgyz = -pgyz;
            pAxy = -pAxy;
            pAyz = -pAyz;
            pchix = -pchix;
            pchiz = -pchiz;
            ptrKx = -ptrKx;
            ptrKz = -ptrKz;
            pAxxx = -pAxxx;
            pAxxz = -pAxxz;
            pAxyy = -pAxyy;
            pAxzx = -pAxzx;
            pAxzz = -pAxzz;
            pAyyx = -pAyyx;
            pAyyz = -pAyyz;
            pAyzy = -pAyzy;
            pAzzx = -pAzzx;
            pAzzz = -pAzzz;
            pGamxxx = -pGamxxx;
            pGamxxz = -pGamxxz;
            pGamxyy = -pGamxyy;
            pGamxzz = -pGamxzz;
            pGamyxy = -pGamyxy;
            pGamyyz = -pGamyyz;
            pGamzxx = -pGamzxx;
            pGamzxz = -pGamzxz;
            pGamzyy = -pGamzyy;
            pGamzzz = -pGamzzz;
            pRxy = -pRxy;
            pRyz = -pRyz;
            break;
          case 6: //--+ (theta, pi+phi)
            costheta = arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            px = -px;
            py = -py;
            pgxz = -pgxz;
            pgyz = -pgyz;
            pAxz = -pAxz;
            pAyz = -pAyz;
            pchix = -pchix;
            pchiy = -pchiy;
            ptrKx = -ptrKx;
            ptrKy = -ptrKy;
            pAxxx = -pAxxx;
            pAxxy = -pAxxy;
            pAxyx = -pAxyx;
            pAxyy = -pAxyy;
            pAxzz = -pAxzz;
            pAyyx = -pAyyx;
            pAyyy = -pAyyy;
            pAyzz = -pAyzz;
            pAzzx = -pAzzx;
            pAzzy = -pAzzy;
            pGamxxx = -pGamxxx;
            pGamxxy = -pGamxxy;
            pGamxyy = -pGamxyy;
            pGamxzz = -pGamxzz;
            pGamyxx = -pGamyxx;
            pGamyxy = -pGamyxy;
            pGamyyy = -pGamyyy;
            pGamyzz = -pGamyzz;
            pGamzxz = -pGamzxz;
            pGamzyz = -pGamzyz;
            pRxz = -pRxz;
            pRyz = -pRyz;
            break;
          case 7: //--- (pi-theta, pi+phi)
            costheta = -arcostheta[i];
            cosmphi = cos(pm * (PI + (j + 0.5) * dphi));
            sinmphi = sin(pm * (PI + (j + 0.5) * dphi));
            px = -px;
            py = -py;
            pz = -pz;
            pchix = -pchix;
            pchiy = -pchiy;
            pchiz = -pchiz;
            ptrKx = -ptrKx;
            ptrKy = -ptrKy;
            ptrKz = -ptrKz;
            pAxxx = -pAxxx;
            pAxxy = -pAxxy;
            pAxxz = -pAxxz;
            pAxyx = -pAxyx;
            pAxyy = -pAxyy;
            pAxyz = -pAxyz;
            pAxzx = -pAxzx;
            pAxzy = -pAxzy;
            pAxzz = -pAxzz;
            pAyyx = -pAyyx;
            pAyyy = -pAyyy;
            pAyyz = -pAyyz;
            pAyzx = -pAyzx;
            pAyzy = -pAyzy;
            pAyzz = -pAyzz;
            pAzzx = -pAzzx;
            pAzzy = -pAzzy;
            pAzzz = -pAzzz;
            pGamxxx = -pGamxxx;
            pGamxxy = -pGamxxy;
            pGamxxz = -pGamxxz;
            pGamxyy = -pGamxyy;
            pGamxyz = -pGamxyz;
            pGamxzz = -pGamxzz;
            pGamyxx = -pGamyxx;
            pGamyxy = -pGamyxy;
            pGamyxz = -pGamyxz;
            pGamyyy = -pGamyyy;
            pGamyyz = -pGamyyz;
            pGamyzz = -pGamyzz;
            pGamzxx = -pGamzxx;
            pGamzxy = -pGamzxy;
            pGamzxz = -pGamzxz;
            pGamzyy = -pGamzyy;
            pGamzyz = -pGamzyz;
            pGamzzz = -pGamzzz;
          }

          f_getnp4_point(px, py, pz, pchi, ptrK,
                         pgxx, pgxy, pgxz, pgyy, pgyz, pgzz,
                         pAxx, pAxy, pAxz, pAyy, pAyz, pAzz,
                         pchix, pchiy, pchiz,
                         ptrKx, ptrKy, ptrKz,
                         pAxxx, pAxxy, pAxxz,
                         pAxyx, pAxyy, pAxyz,
                         pAxzx, pAxzy, pAxzz,
                         pAyyx, pAyyy, pAyyz,
                         pAyzx, pAyzy, pAyzz,
                         pAzzx, pAzzy, pAzzz,
                         pGamxxx, pGamxxy, pGamxxz, pGamxyy, pGamxyz, pGamxzz,
                         pGamyxx, pGamyxy, pGamyxz, pGamyyy, pGamyyz, pGamyzz,
                         pGamzxx, pGamzxy, pGamzxz, pGamzyy, pGamzyz, pGamzzz,
                         pRxx, pRxy, pRxz, pRyy, pRyz, pRzz,
                         psi4RR, psi4II);

          thetap = sqrt((2 * pl + 1.0) / 4.0 / PI) * misc::Wigner_d_function(pl, pm, spinw, costheta); // note the variation from -2 to 2

          //	 find back the one
          pchi = pchi + 1;
#ifdef GaussInt
          // wtcostheta is even function respect costheta
          RP_out[countlm] = RP_out[countlm] + thetap / pchi / pchi * (psi4RR * cosmphi + psi4II * sinmphi) * wtcostheta[i];
          IP_out[countlm] = IP_out[countlm] + thetap / pchi / pchi * (psi4II * cosmphi - psi4RR * sinmphi) * wtcostheta[i];
#else
          RP_out[countlm] = RP_out[countlm] + thetap / pchi / pchi * (psi4RR * cosmphi + psi4II * sinmphi);
          IP_out[countlm] = IP_out[countlm] + thetap / pchi / pchi * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
        }
        countlm++; // no sanity check for countlm and NN which should be noted in the input parameters
      }
  }

  for (int ii = 0; ii < NN; ii++)
  {
#ifdef GaussInt
    RP_out[ii] = RP_out[ii] * rex * dphi;
    IP_out[ii] = IP_out[ii] * rex * dphi;
#else
    RP_out[ii] = RP_out[ii] * rex * dphi * dcostheta;
    IP_out[ii] = IP_out[ii] * rex * dphi * dcostheta;
#endif
  }
  //|------+  Communicate and sum the results from each processor.

  MPI_Allreduce(RP_out, RP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(IP_out, IP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //|------= Free memory.

  delete[] pox[0];
  delete[] pox[1];
  delete[] pox[2];
  delete[] shellf;
  delete[] RP_out;
  delete[] IP_out;
  DG_List->clearList();
}
//|----------------------------------------------------------------
//  do not discriminate box and shell
//|----------------------------------------------------------------
bool surface_integral::SR_Interp_Points(MyList<var> *VarList, cgh *GH, ShellPatch *SH,
                                        int NN, double **XX, double *Shellf)
{
  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double pox[3];
  for (int i = 0; i < NN; i++)
  {
    for (int j = 0; j < 3; j++)
      pox[j] = XX[j][i];
    int lev = GH->levels - 1;
    bool notfound = true;

    while (notfound)
    {
      if (lev < 0)
      {
        if (SH)
        {
          if (SH->Interp_One_Point(VarList, pox, Shellf + i * num_var, Symmetry))
          {
            return true;
          }
          if (myrank == 0)
            cout << "surface_integral::SR_Interp_Points point (" << pox[0] << "," << pox[1] << "," << pox[2] << ") is out of cgh and shell domain!" << endl;
        }
        else
        {
          if (myrank == 0)
            cout << "surface_integral::SR_Interp_Points: point (" << pox[0] << "," << pox[1] << "," << pox[2] << ") is out of cgh domain!" << endl;
        }
        return false;
      }
      MyList<Patch> *Pp = GH->PatL[lev];
      while (Pp)
      {
        if (Pp->data->Interp_ONE_Point(VarList, pox, Shellf + i * num_var, Symmetry))
        {
          notfound = false;
          break;
        }
        Pp = Pp->next;
      }
      lev--;
    }
  }
  return true;
}
