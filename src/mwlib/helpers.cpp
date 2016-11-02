/**
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2016  Michael Welter and Thierry Fredrich

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define INTERNAL
#if defined _WIN32 && !defined WIN32
#  define WIN32 _WIN32
#endif
#ifdef WIN32
    #define  WIN32_LEAN_AND_MEAN
    #ifndef _WIN32_WINNT
      #define _WIN32_WINNT 0x0400
    #endif
    #include <windows.h>
    #include <mmsystem.h>
    #include <Psapi.h>
    #include <shellapi.h>
//#ifdef DEBUG    
    #include <Dbghelp.h>
//#endif
    #undef min
    #undef max
#else
  #include<sys/time.h>
  #include<dirent.h>
  #include<sys/stat.h>
  #include<sys/resource.h> 
#endif

#include "helpers.h"
#include "drawing.h"
#include <stdio.h>
#include <memory.h>
#include <sstream>
#include <fstream>
#include "math_ext.h"
#include "random.h"

using std::ostream;


//void* operator new(size_t s) { return GeAlloc(s); }
//void* operator new[](size_t s) { return GeAlloc(s); }
//void operator delete(void *p) { GeFree(p); }
//void operator delete[](void *p) { GeFree(p); }

void _Fatal( const char *msg, const char *file, int line )
{
  fprintf(stderr,"Fatal Error: \"%s\" in %s:%i\n",msg,file,line);
  exit(0);
}


float GetRand()
{
  return ((float)(rand()))/RAND_MAX;
}


void _SetMem( void *ptr, int c, size_t n )
{
  memset( ptr, c, n );
}

void _CopyMem( const void *src, void *dst, size_t n )
{
  memcpy( dst, src, n );
}


namespace calc_strides
{

void first_dim_varies_fastest(int rank, const int *dims, int *strides)
{
  strides[0] = 1;
  for (int i=1; i<rank; ++i)
    strides[i] = dims[i-1] * strides[i-1];
}

void last_dim_varies_fastest(int rank, const int *dims, int *strides)
{
  strides[rank-1] = 1;
  for (int i=rank-2; i>=0; --i)
    strides[i] = dims[i+1] * strides[i+1];
}
  
}



void HPBToMatrix(const Float3 &w, Float3x4 &frame)
{
	//float cn,sn,ck,sk,cs,ss,cosksinn,sinksinn;

	//SinCos(w[0],ss,cs);
	//SinCos(w[1],sn,cn);
	//SinCos(w[2],sk,ck);

	//cosksinn=ck*sn;
	//sinksinn=sk*sn;
 //	
	//frame.X = Float3(ck*cs-sinksinn*ss,  ck*ss+sinksinn*cs, -sk*cn  );
	//frame.Y = Float3(sk*cs+cosksinn*ss,  sk*ss-cosksinn*cs  , ck*cn   );
	//frame.Z = Float3(-cn*ss           ,cn*cs            , sn      );

  float a,b,c,d,e,f,cosksinn,sinksinn;
	my::sincos<float>(w[0],d,c); //y = h
	my::sincos<float>(w[1],b,a); //x = p
	my::sincos<float>(w[2],f,e); //z = b
 // sinksinn = d * b;
	//cosksinn = c * b;
	//		frame.X=Float3(c * e + sinksinn * f, -c * f + sinksinn * e, d * a);
	//									frame.Y=Float3(a * f, a * e, -b);
	//									frame.Z=Float3(-d * e + cosksinn * f, d * f + cosksinn * e, c * a);
	cosksinn = e * b;
	sinksinn = f * b;
	frame.X = Float3(e * c - sinksinn * d, -f * a, e * d + sinksinn * c);
	frame.Y = Float3(f * c + cosksinn * d, e * a, f * d - cosksinn * c);
	frame.Z = Float3(-a * d, b, a * c);
}

Float3 operator^( const Float3x4 &M, const Float3 &x )
{
	Float3 y;
	y[0] = M.X[0]*x[0] + M.Y[0]*x[1] + M.Z[0]*x[2];
	y[1] = M.X[1]*x[0] + M.Y[1]*x[1] + M.Z[1]*x[2];
	y[2] = M.X[2]*x[0] + M.Y[2]*x[1] + M.Z[2]*x[2];
	return y;
}


const Float3x4 operator*(const Float3x4 &m1, const Float3x4 &m2)
{
  Float3x4 erg;
	erg.T[0] = m1.T[0] + m1.X[0] * m2.T[0] + m1.Y[0] * m2.T[1] + m1.Z[0] * m2.T[2];
	erg.T[1] = m1.T[1] + m1.X[1] * m2.T[0] + m1.Y[1] * m2.T[1] + m1.Z[1] * m2.T[2];
	erg.T[2] = m1.T[2] + m1.X[2] * m2.T[0] + m1.Y[2] * m2.T[1] + m1.Z[2] * m2.T[2];

	erg.X[0] = m1.X[0] * m2.X[0] + m1.Y[0] * m2.X[1] + m1.Z[0] * m2.X[2];
	erg.X[1] = m1.X[1] * m2.X[0] + m1.Y[1] * m2.X[1] + m1.Z[1] * m2.X[2];
	erg.X[2] = m1.X[2] * m2.X[0] + m1.Y[2] * m2.X[1] + m1.Z[2] * m2.X[2];

	erg.Y[0] = m1.X[0] * m2.Y[0] + m1.Y[0] * m2.Y[1] + m1.Z[0] * m2.Y[2];
	erg.Y[1] = m1.X[1] * m2.Y[0] + m1.Y[1] * m2.Y[1] + m1.Z[1] * m2.Y[2];
	erg.Y[2] = m1.X[2] * m2.Y[0] + m1.Y[2] * m2.Y[1] + m1.Z[2] * m2.Y[2];

	erg.Z[0] = m1.X[0] * m2.Z[0] + m1.Y[0] * m2.Z[1] + m1.Z[0] * m2.Z[2];
	erg.Z[1] = m1.X[1] * m2.Z[0] + m1.Y[1] * m2.Z[1] + m1.Z[1] * m2.Z[2];
	erg.Z[2] = m1.X[2] * m2.Z[0] + m1.Y[2] * m2.Z[1] + m1.Z[2] * m2.Z[2];
  return erg;
}

float Float3x4::Det() const
{
	return X[0]*Y[1]*Z[2] -
				 X[0]*Y[2]*Z[1] +
				 X[2]*Y[0]*Z[1] -
				 X[2]*Y[1]*Z[0] +
				 X[1]*Y[2]*Z[0] -
				 X[1]*Y[0]*Z[2];
}

void Inverse( const Float3x4 &M, Float3x4 &R )
{
	float det = M.Det();
	if( std::abs(det)<1.0e-10 ) return;
	det = 1.0f/det;

	R.X[0] = (M.Y[1] * M.Z[2] - M.Z[1] * M.Y[2]) * det;
	R.X[1] = (M.Z[1] * M.X[2] - M.X[1] * M.Z[2]) * det;
	R.X[2] = (M.X[1] * M.Y[2] - M.Y[1] * M.X[2]) * det;

	R.Y[0] = (M.Y[2] * M.Z[0] - M.Z[2] * M.Y[0]) * det;
	R.Y[1] = (M.Z[2] * M.X[0] - M.X[2] * M.Z[0]) * det;
	R.Y[2] = (M.X[2] * M.Y[0] - M.Y[2] * M.X[0]) * det;

	R.Z[0] = (M.Y[0] * M.Z[1] - M.Z[0] * M.Y[1]) * det;
	R.Z[1] = (M.Z[0] * M.X[1] - M.X[0] * M.Z[1]) * det;
	R.Z[2] = (M.X[0] * M.Y[1] - M.Y[0] * M.X[1]) * det;

	R.T[0] = -R.X[0]*M.T[0] - R.Y[0]*M.T[1] - R.Z[0]*M.T[2];
	R.T[1] = -R.X[1]*M.T[0] - R.Y[1]*M.T[1] - R.Z[1]*M.T[2];
	R.T[2] = -R.X[2]*M.T[0] - R.Y[2]*M.T[1] - R.Z[2]*M.T[2];
}

void Transpose( const Float3x4 &M, Float3x4 &R )
{
	R.T = Float3(0);
	R.X[1] = M.Y[0];
	R.X[2] = M.Z[0];
	R.Y[0] = M.X[1];
	R.Y[2] = M.Z[1];
	R.Z[0] = M.X[2];
	R.Z[1] = M.Y[2];
  R.X[0] = M.X[0];
  R.Y[1] = M.Y[1];
  R.Z[2] = M.Z[2];
}

void InverseTranspose(  const Float3x4 &M, Float3x4 &R )
{
	float det = M.Det();
	if( std::abs(det)<1.0e-10 ) return;
	det = 1.0f/det;

	R.X[0] = (M.Y[1] * M.Z[2] - M.Z[1] * M.Y[2]) * det;
	R.Y[0] = (M.Z[1] * M.X[2] - M.X[1] * M.Z[2]) * det;
	R.Z[0] = (M.X[1] * M.Y[2] - M.Y[1] * M.X[2]) * det;

	R.X[1] = (M.Y[2] * M.Z[0] - M.Z[2] * M.Y[0]) * det;
	R.Y[1] = (M.Z[2] * M.X[0] - M.X[2] * M.Z[0]) * det;
	R.Z[1] = (M.X[2] * M.Y[0] - M.Y[2] * M.X[0]) * det;

	R.X[2] = (M.Y[0] * M.Z[1] - M.Z[0] * M.Y[1]) * det;
	R.Y[2] = (M.Z[0] * M.X[1] - M.X[0] * M.Z[1]) * det;
	R.Z[2] = (M.X[0] * M.Y[1] - M.Y[0] * M.X[1]) * det;

	R.T = Float3(0);
}

const Float3x4 Inverse(const Float3x4 &M)
{
	Float3x4 R;
	Inverse(M,R);
	return R;
}

const Float3x4 InverseTranspose( const Float3x4 &M )
{
	Float3x4 MIT;
	InverseTranspose(M,MIT);
	return MIT;
}

void Diagonal( float s, Float3x4 &D )
{
	D.X = Float3(s,0,0);
	D.Y = Float3(0,s,0);
	D.Z = Float3(0,0,s);
	D.T = Float3(0);
}

Float3x4 OrthogonalSystem(const Float3 &Z)
{
  Float3x4 m;
  m.Z = Z.normalized();
  if(m.Z.x()*m.Z.x() + m.Z.y()*m.Z.y() > 1.0e-6) {
    m.Y = cross(m.Z,Float3(0,0,1)).normalized();
    m.X = cross(m.Y,m.Z).normalized();
  } else {
    m.Y = Float3(0,1,0);
    m.X = Float3(1,0,0);
  }
  return m;
}



/*
void Eigenvalues( const Float2x2 &m, Float2 &ev, Float2x2 &o )
{
  myAssert(m.v1[1]==m.v2[0]);
  const float tr  = m.v1[0]+m.v2[1];
  const float det = m.Det();
  const float cc = tr*tr*0.25f-det;
  if(cc>1.0e-19)
  {
    ev[1] = tr*0.5f+sqrtf(tr*tr*0.25f-det);
    ev[0] = tr*0.5f-sqrtf(tr*tr*0.25f-det);
    float dTheta;
    if ( m.v1[1] == 0.0f && m.v1[0] == m.v2[1] )
        dTheta = 0.0f;
    else
        dTheta = 0.5f*atan2f(-2.0f*m.v1[1], m.v2[1]-m.v1[0]);
    float cs = cosf(dTheta);
    float sn = sinf(dTheta);
    o.v1[0] = cs;
    o.v1[1] = sn;
    o.v2[0] = -sn;
    o.v2[1] = cs;
  }
  else
  {
    ev[0]=ev[1] = tr*0.5f;
    o.v1=Float2(1.0f,0.0f);
    o.v2=Float2(0.0f,1.0f);
  }
}*/


#define JACOBI_ROTATE(a,i,j,k,l) g=a[i][j]; h=a[k][l]; a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);
/*
void Eigenvalues( const Float3x4 &m, Float3 &d, Float3x4 &v )
{
  Float3x4 a(m);

  myAssert(Abs(a[0][1]-a[1][0])<1.0e-6);
  myAssert(Abs(a[0][2]-a[2][0])<1.0e-6);
  myAssert(Abs(a[1][2]-a[2][1])<1.0e-6);
  Float3 b,z(0);
  
  Diagonal(1.0f,v);
  b = d = Float3(a[0][0],a[1][1],a[2][2]);
  
  for(int i=0; i<10; ++i)
  {
    float sm = fabs(a[0][1])+fabs(a[0][2])+fabs(a[1][2]);
    if(sm==0.0f) return;
    float thres;
    if(i<4) 
      thres=0.2f*sm/(3*3);
    else
      thres=0.0f;
    for(int ip=0; ip<2; ++ip) for(int iq=ip+1; iq<3; ++iq)
    {
      float g = 100.0f*fabs(a[ip][iq]);
      //if (i>4 && fabs(d[ip])+g == fabs(d[ip]) &&
      //           fabs(d[iq])+g == fabs(d[iq]))
      //           a[ip][iq]=0.0f;
      //else if(fabs(a[ip][iq])>thres)     
      //{
        float h = d[iq]-d[ip];
        float t;
        if(fabs(h)+g==fabs(h))
          t = a[ip][iq]/h;
        else
        {
          float theta=0.5f*h/a[ip][iq];
          t = 1.0f/(fabs(theta)+sqrtf(1.0f+theta*theta));
          if(theta<0.0f) t = -t;
        }
        float c = 1.0f/sqrtf(1+t*t);
        float s = t*c;
        float tau = s/(1.0f+c);
        h = t*a[ip][iq];
        z[ip] -= h;
        z[iq] += h;
        d[ip] -= h;
        d[iq] += h;
        a[ip][iq]=a[iq][ip] = 0.0f;
        int j;
        for(j=0   ; j<ip; ++j) { JACOBI_ROTATE(a,j ,ip,j ,iq) }
        for(j=ip+1; j<iq; ++j) { JACOBI_ROTATE(a,ip,j ,j ,iq) }
        for(j=iq+1; j<3 ; ++j) { JACOBI_ROTATE(a,ip,j ,iq,j ) }
        //for(j=0   ; j<3 ; ++j) { JACOBI_ROTATE(v,j ,ip,j ,iq) }
        for(j=0   ; j<3 ; ++j) { JACOBI_ROTATE(v,ip,j,iq,j) }
      //}
    }
    for(int j=0; j<3; ++j)
    {
      b[j] += z[j];
      d[j]  = b[j];
      z[j]  = 0.0f;
    }

    //cout << "a = " << endl;
    //for(int y=0; y<3; ++y) 
    //{
    //  for(int x=0; x<3; ++x)
    //  {
    //    cout << a[y][x] << " ";
    //  }
    //  cout << endl;
    //}
  }
}

void TestEV()
{
  Float3x4 rot(0),scale(0),rot_transpose(0);
  scale[0][0]=3.0f; scale[1][1]=2.0f; scale[2][2]=1.0f;
  HPBToMatrix(Float3(Rad(45.f),0,Rad(138.0f)),rot);
  Transpose(rot,rot_transpose);
  Float3x4 m = rot*scale*rot_transpose;
  
  cout << "m: " << endl;
  cout << m.X <<endl<< m.Y <<endl<< m.Z << endl;
  cout << endl;

  cout << "rot: " << endl;
  cout << rot.X <<endl<< rot.Y <<endl<< rot.Z << endl;
  cout << endl;

  Float3 d;
  Float3x4 v;
  Eigenvalues(m,d,v);
  
  cout << "d: " << d << endl;
  cout << "v: " << endl;
  cout << v.X <<endl<< v.Y <<endl<< v.Z << endl;
  
  getchar();
}


Float2 LeastSquaresSolveNx2( const float* A, const float* b, const int N )
{
  Float2x2 ATA(0.f);
  Float2 ATb(0.f);
  
  for( int i=0; i<2; ++i )
  {
    for( int j=0; j<2; ++j )
    {
      for( int c=0; c<N; ++c )
      {
        ATA(i,j) += A[i+2*c]*A[j+2*c];
      }
    }   
    for( int c=0; c<N; ++c )
    {
      ATb[i] += A[i+2*c]*b[c];
    }
  }
  ATA.Invert();
  return ATA*ATb;  
}


Float3 LeastSquaresSolveNx3( const float* A, const float* b, const int N )
{
  Float3x4 ATA(0.f);
  Float3 ATb(0.f);
  
  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      for( int c=0; c<N; ++c )
      {
        ATA(i,j) += A[i+3*c]*A[j+3*c];
      }
    }   
    for( int c=0; c<N; ++c )
    {
      ATb[i] += A[i+3*c]*b[c];
    }
  }
  ATA = Inverse(ATA);
  return ATA^ATb;  
}*/



void TestRandom( Random &rand, const char* imgfilename )
{
  const uint s = 10;
  uint mcnt=0;
  std::vector<uint> field( s*s, 0 );  
  for( uint n=0; n<s*s*1; ++n )
  {
    //const uint x = rand.Get( s );
    //const uint y = rand.Get( s );
    const double fx = rand.Get01();
    const double fy = rand.Get01();
    int x = int(fx*s);
    int y = int(fy*s);
    ++(field[ y*s+x ]);
    mcnt = std::max( field[y*s+x], mcnt );
  }
  Image img( s, s );
  for( int y=0; y<s; ++y )
  {
    for( int x=0; x<s; ++x )
    {
      uint c =field[y*s+x];
      uchar uc = uchar((((float)c)/(mcnt))*255);
      img.SetColor( uc,uc,uc );
      img.DrawPixel( x,y );
    }
  }
  img.Write( imgfilename );
}

void TestRandom( const char* imgfilename )
{
  Random rand;
  TestRandom( rand, imgfilename );
}





#define M 397
#define MATRIX_A 0x9908b0dfUL
#define UMASK 0x80000000UL
#define LMASK 0x7fffffffUL
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))

void Random::Init(uint s)
{
    int j;
    state[0]= s & 0xffffffffUL;
    for (j=1; j<N; j++) {
        state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j);
        state[j] &= 0xffffffffUL;
    }
    left = 1;
}

void Random::NextState()
{
  uint *p=state;
  int j;

  left = N;
  //next = state;
  inext= 0;

  for (j=N-M+1; --j; p++)
      *p = p[M] ^ TWIST(p[0], p[1]);

  for (j=M; --j; p++)
      *p = p[M-N] ^ TWIST(p[0], p[1]);

  *p = p[M-N] ^ TWIST(p[0], state[0]);
}




/*

const float TMSNet::DENOM = 1.0f / (float)((long long)1 << 32);*/


//template< class T >
//void HistogramInfo<T>::Init( const MinMax<T> &_limits, uint _num_buckets, bool _bLogScale, T _eps )
//{
//  bLogScale = _bLogScale;
//  limits = _limits;
//  eps = _eps;
//  num_buckets = _num_buckets;
//  
//  if( bLogScale ) {
//    limits.min = LogTrafo( limits.min, eps );
//    limits.max = LogTrafo( limits.max, eps );
//  }  
//  size_bucket = (limits.max-limits.min)/num_buckets;  
//}
//
//template struct HistogramInfo<float>;
//template struct HistogramInfo<double>;


int RandomPick( Random& random, int cnt, double *prob )
{
  double r = random.Get01();
  double pp=0.0;
  for( int k=0; k<cnt; ++k )
  {
    pp += prob[k];
    if( r<pp ) return k;
  }
  myAssert( std::abs(pp-1.0)<1.0e-6 );
  return cnt-1;
}


double NormalizeSumNorm( double* x, int cnt )
{
	double n = 0.0;
	for( int i=0; i<cnt; ++i ) n+=std::abs(x[i]);
	if( n>0.0 ) {
		const double nn = 1.0/n;
		for( int i=0; i<cnt; ++i ) x[i]*=nn;
	}
  return n;
}



uint CalcPowOf2GE( uint i )
{
  if( i>(1u<<31) ) return 0;
  uint r = 1;
  while( r<i ) {
    r = r<<1;
  }
  return r;
}

int  CalcPowOf2GE( int i )
{
  if( i<0 ) return -CalcPowOf2GE(-i);
  if( i>(1u<<30) ) return 0;
  int r = 1;
  while( r<i ) {
    r = r<<1;
  }
  return r;
}





double MostSignificantImpl( double x, const uint b, double &_bn, uchar nk )
{
  myAssert( x >= 0.f );
  if( x == 0.0f ) { _bn = 0.f; return 0.f; }
  uint bn = 1;
  uint bbn = b;
  uint bbbn = 1;
  for( uchar i = 0; i<nk-1; ++i ) bbbn *= b;
  if( x < bbbn )
  {    
    while( x*bbn < bbbn && (bbn*b>bbn) ) // carefull to prevent overflow
    {
      bn *= b;
      bbn *= b;
    }
    _bn = 1.0f/(bbn);
    return double(uint(x*bbn))/bbn;
  }
  else
  {
    while( x >= bbn*bbbn && (bbn*b>bbn) )
    {
      bn *= b;
      bbn *= b;
    }
    _bn = bn;
    return double(uint(x/bn))*bn;
  }
}



double MostSignificantFloor( double x, const uint b, uchar nk )
{
  double bn;  
  if( x >= 0.f )
    return MostSignificantImpl( x, b, bn, nk );
  else 
  {
    double ret = -MostSignificantImpl( -x, b, bn, nk );
    return ret - bn;
  }
  return x;
}
double MostSignificantCeil( double x, const uint b, uchar nk )
{
  double bn;
  if( x >= 0.f )
  {
    double ret = MostSignificantImpl( x, b, bn, nk );
    return ret + bn;
  }
  else return -MostSignificantImpl( -x, b, bn, nk );
  return x;
}

void TestMostSignificant()
{
  Random rand;
  for( int i=0; i<100; ++i )
  {
    int n = 0;
    double x = pow( 10.0, my::lerp<double>( rand.Get01(), -5, 5 ) );
    //if( rand.Get01()<0.5f ) x = -x;
    double f = MostSignificantFloor( x, 10, 1 );
    double c = MostSignificantCeil( x, 10, 1 );
    myAssert( f <= x && c >= x );
    printf( "%010lf <= %010lf <= %010lf \n",f,x,c );
  }
  int c = getchar();
}



#if 1
#include<stdarg.h>


const string strprintf( const char* format, ... )
{
	va_list argp;
	char buffer[2048];
	va_start( argp, format );
#ifdef WIN32
	_vsnprintf( buffer, 2048, format, argp );
#else
	vsnprintf( buffer, 2048, format, argp );
#endif
	va_end( argp );
	return string( buffer );
}

const string RemoveExtension( const string &s )
{
  string::size_type pos = s.find_last_of('.');
  string::size_type dir_pos = s.find_last_of("\\/");
  if( pos != string::npos && (pos>dir_pos || dir_pos==string::npos)) {
    return s.substr( 0, pos );
  }
  else return s;
}

const string GetExtension( const string &s )
{
  string::size_type pos = s.find_last_of('.');
  if( pos != string::npos ) {
    return s.substr( pos+1, s.length()-pos-1 );
  }
  else return string();
}

int64 GetPathSep( const string &s )
{
  string::size_type pos0 = s.rfind('\\');
  string::size_type pos1 = s.rfind('/');
  if( pos0!=string::npos && pos1==string::npos ) {
    return pos0;
  }
  else if( pos1!=string::npos && pos0==string::npos ) {
    return pos1;
  }
  else if( pos0!=string::npos && pos1!=string::npos ) {
    return std::max( pos0, pos1 );
  }
  else return -1;
}



const string RemoveAllExtensions( const string &s )
{
  int64 psep = std::max<int64>(0,GetPathSep(s));
  string::size_type pos = s.find('.',psep);
  if( pos!=string::npos ) 
    return s.substr( 0, pos );
  else
    return s;
}


const string GetPath( const string &s )
{
  int64 p = GetPathSep( s );
  if( p>=0 ) return s.substr( 0, p );
  else return string("");
}

const string RemovePath( const string &s )
{
  int64 p = GetPathSep( s );
  if( p>=0 ) return s.substr( p+1, s.length()-p-1 );
  else return s;
}

const string JoinPaths( const string &a, const string &b )
{
#ifdef WIN32
  const char sep = '\\';
#else
  const char sep = '/';
#endif
  const int64 l = a.length();
  const int64 k = b.length();
  if( l > 0 )
  {
    int64 x = l-1;
    int64 y = 0;
    while( (a[x]=='\\' || a[x]=='/') && x>0 ) --x;
    while( (b[y]=='\\' || b[y]=='/') && y<k ) ++y;
    return a.substr( 0, x+1 ) + sep + b.substr( y, k );
  }
  else return b;
}

bool SplitPathBack( const string &s, string &a, string &b )
{
  int64 p = GetPathSep(s);
  if(p>=0) 
  {
    string tmp = s;
    a = tmp.substr( 0, p );
    b = tmp.substr( p+1, s.length()-p-1 );
    return true;
  }
  else
  {
    a = string();
    b = s;
    return false;
  }
}

bool   FileExists(const string &s )
{
  if(s.empty()) return false;
  FILE *fp = fopen( s.c_str(), "rb" );
  if( !fp ) return false;
  fclose(fp);
  return true;
}

bool   PathExists( const string &s, bool bCheckForDir )
{
#ifdef WIN32
  DWORD a = GetFileAttributes(s.c_str());
  return (a!=INVALID_FILE_ATTRIBUTES) && (!bCheckForDir || (a&FILE_ATTRIBUTE_DIRECTORY));
#else
  DIR* dir = opendir(s.c_str());
  if(!dir) return false;
  closedir(dir);
  return true;
#endif
}

bool _CreatePathRec( const string &s )
{
  if(s.empty() || PathExists(s)) return true;
  
  string sub = GetPath(s);
  if(!_CreatePathRec(sub)) return false;
  {
#ifdef WIN32
    BOOL b = CreateDirectory(s.c_str(),NULL);  
    return b==TRUE;
#else
    int res = mkdir(s.c_str(),S_IRWXU|S_IRWXG);
    return res==0;
#endif
  }
}

bool   CreatePath( const string &s )
{  
  if(PathExists(s,false)) return true;  
  return _CreatePathRec(s);
}

bool IsFileGZip( const string &s )
{
  if(s.empty()) return false;
  FILE *fp = fopen( s.c_str(), "rb" );
  if( !fp ) return false;
  uchar dd[10];
  size_t n = fread( dd, 1, 10, fp );
  fclose(fp);
  return n==10 && dd[0]==0x1f && dd[1]==0x8b;
}

const string FormatMemSize(size_t size)
{
  size_t kilos = size/1024;
  size_t megs  = kilos/1024;
  size  -=kilos*1024;
  kilos -=megs*1024;
  if( megs>0 ) {
    return strprintf("%uMB %ukB",megs,kilos);
  }
  else {
    return strprintf("%ukB %uBytes",kilos,size);
  }
}

inline uint Div( uint a, uint b, uint &rest )
{
  uint c = a/b;
  rest = a-c*b;
  return c;
}

inline uint Div( double a, double b, double &rest )
{
  double c=a/b;
  uint i = uint(c+0.5);
  rest = a-double(i)*b;
  return i;
}

const string FormatTimeSec( double s )
{
  uint d = Div( s, 24*60*60, s );
  uint h = Div( s, 60*60, s );
  uint m = Div( s, 60, s );
  if(d>0)
    return strprintf("%ud %uh %um %0.2lfs",d,h,m,s);
  else if(h>0)
    return strprintf("%uh %um %0.2lfs",h,m,s);
  else if(m>0)
    return strprintf("%um %0.2lfs",m,s);
  else
    return strprintf("%0.2lfs",s);
}




const string FormatTimeMs( uint ms )
{
  uint h = Div(ms,1000*60*60,ms);
  uint m = Div(ms,1000*60,ms);
  uint s = Div(ms,1000,ms);
  return strprintf("%uh %um %u.%us",h,m,s,ms);
}

bool strToInt( const char* s, int &x )
{
    std::istringstream is( s );
    is >> x;
    return !is.fail() && is.eof();
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//
//RefCounted::RefCounted() : cnt(0) 
//{
////#ifdef WIN32
////  NT_TIB* tib = (PNT_TIB)NtCurrentTeb();
////  isOnStack = this<=tib->StackBase && this>=tib->StackBase;
////#else
////
////#endif
//}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//void GeDebugOut( const char* format, ... )
//{
//  va_list argp;
//  va_start( argp, format );
//  vfprintf(stderr,format,argp);
//  fputc( '\n',stderr );
//  va_end( argp );
//}

// bool _CopyFile( const string &src, const string &dst )
// {
// #ifdef WIN32
//   return CopyFile( src.c_str(), dst.c_str(), FALSE );
// #else
// #error Implement Me
// #endif
// }

// bool _MoveFile( const std::string &src, const std::string &dst )
// {
// #ifdef WIN32
//   return MoveFile( src.c_str(), dst.c_str() );
// #else
// #error Implement Me
// #endif
// }

// bool _RemoveFile( const string &fn )
// {
// #ifdef WIN32
//   return DeleteFile( fn.c_str() );
// #else
// #error Implement Me
// #endif
// }

uint _GetTimer()
{
#ifdef WIN32
  return timeGetTime();
#else
  timeval tv;
  gettimeofday( &tv, NULL );
  return tv.tv_sec*1000+tv.tv_usec/1000;
#endif
}


void DebugOutString( const char* format, ... )
{
  va_list argp;  
  va_start( argp, format );
#ifdef WIN32
  char buffer[2048];
	_vsnprintf( buffer, 2048, format, argp );
	va_end( argp );
	fputs(buffer,stderr);
	OutputDebugString(buffer);  
#else
	vfprintf( stderr, format, argp );
	va_end( argp );
#endif	
}

#include <assert.h>

// had some trouble breaking into the debugger with default assert under windows
void _Assert( const char *msg, const char* file, int line )
{
#ifdef WIN32
  {
    const int SIZE = 1024;
    const int STACKSIZE = 8;
    TCHAR buffer[SIZE];
    
    _snprintf(buffer,SIZE,"\nAssertion \"%s\" failed in file %s (%i)\n",msg,file,line);
    DebugOutString(buffer);    

    bool bok = true;

    void* stack[STACKSIZE];
    GetCallstack(stack,STACKSIZE);

    for(int i=3; i<STACKSIZE && stack[i]!=NULL; ++i)
    {
      int line;
      char* file;
      if(GetLineFile(stack[i],line,file))
      {          
          _snprintf(buffer,SIZE-1,"  %s (%i)\n",file,line);
      }
      else
      {
          _snprintf(buffer,SIZE-1,"  %p (unresolved address)\n",stack[i]);
      }
      DebugOutString(buffer);
    }
    if(i>=STACKSIZE)
    {
      DebugOutString("  ...\n");
    }
  } // if output stack
    
  DebugBreak();
#else
  #ifdef DEBUG
  __assert(msg,file,line);
  #endif
#endif
}

/*------------------------------------
  mem debug
------------------------------------*/


bool _HasDebugger()
{
#ifdef WIN32
  return IsDebuggerPresent();
#else
  return false;
#endif
}

void _Sleep(const uint milliseconds) 
{
#ifdef WIN32
  Sleep( milliseconds );
#else
  struct timespec tv;
  tv.tv_sec = milliseconds/1000;
  tv.tv_nsec = (milliseconds%1000)*1000000;
  nanosleep(&tv,0);
#endif
}

#if 0
size_t _GetMemoryUsage(MemUsageMode mode)
{
#ifdef WIN32
  PROCESS_MEMORY_COUNTERS a;
  HANDLE hndl = GetCurrentProcess();
  if( GetProcessMemoryInfo( hndl, &a, sizeof(a) ) )
  {
    //printf("PeakWorkingSetSize=%u\n",a.PeakWorkingSetSize/1024/1024);
    //printf("WorkingSetSize=%u\n",a.WorkingSetSize/1024/1024);
    //printf("QuotaPeakPagedPoolUsage=%u\n",a.QuotaPeakPagedPoolUsage/1024/1024);
    //printf("QuotaPagedPoolUsage=%u\n",a.QuotaPagedPoolUsage/1024/1024);
    //printf("QuotaPeakNonPagedPoolUsage=%u\n",a.QuotaPeakNonPagedPoolUsage/1024/1024);
    //printf("QuotaNonPagedPoolUsage=%u\n",a.QuotaNonPagedPoolUsage/1024/1024);
    //printf("PagefileUsage=%u\n",a.PagefileUsage/1024/1024);
    //printf("PeakPagefileUsage=%u\n",a.PeakPagefileUsage/1024/1024);
    //printf("PrivateUsage=%u\n",a.PrivateUsage/1024/1024);
    return mode==MEMUSAGE_CURRENT ? a.WorkingSetSize : a.PeakWorkingSetSize;
  }
  else return 0;
#else
  using std::ios_base;
  using std::ifstream;
  using std::string;

  double vm_usage     = 0.0;
  double resident_set = 0.0;

  // 'file' stat seems to give the most reliable results
  //
  ifstream stat_stream("/proc/self/stat",ios_base::in);

  // dummy vars for leading entries in stat that we don't care about
  //
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;

  // the two fields we want
  //
  unsigned long vsize;
  long rss;

  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
              >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
              >> utime >> stime >> cutime >> cstime >> priority >> nice
              >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

  stat_stream.close();

  long page_size = sysconf(_SC_PAGE_SIZE);
  vm_usage     = vsize;
  resident_set = rss * page_size;
  if (mode == MEMUSAGE_CURRENT)
    return resident_set;
  else if (mode == MEMUSAGE_VM)
    return vm_usage;
  else
    return 0;
#endif
}
#endif

MemUsage GetMemoryUsage_()
{
  // based on http://www.locklessinc.com/articles/memory_usage/
  using std::ios_base;
  using std::ifstream;
  using std::string;

  MemUsage u;

  ifstream stat_stream("/proc/self/status", ios_base::in);
  stat_stream.exceptions( ifstream::failbit | ifstream::badbit | ifstream::eofbit); // exception in case of error

  string s;
  int found = 0;
  while (found != (1|2|4|8))
  {
    stat_stream >> s;
    if (s == "VmPeak:") {
      stat_stream >> u.vmem_peak;
      found |= 1;
    } else if (s == "VmSize:") {
      stat_stream >> u.vmem;
      found |= 2;
    } else if (s == "VmRSS:") {
      stat_stream >> u.rss;
      found |= 4;
    } else if (s == "VmHWM:") {
      stat_stream >> u.rss_peak;
      found |= 8;
    }
  }

  u.rss *= 1024;
  u.rss_peak *= 1024;
  u.vmem *= 1024;
  u.vmem_peak *= 1024;
  
  return u;
};


//bool OpenFileChooser( std::string &sfile )
//{
//  OPENFILENAME ofn;       // common dialog box structure
//  char szFile[4096];       // buffer for file name
//  //HWND hwnd;              // owner window
//
//  ClearMem( szFile, sizeof(szFile) );
//  CopyMem( sfile.c_str(), szFile, std::min(sizeof(szFile)-1,sfile.size()) );
//
//  // Initialize OPENFILENAME
//  ZeroMemory(&ofn, sizeof(ofn));
//  ofn.lStructSize = sizeof(ofn);
//  ofn.hwndOwner = NULL;
//  ofn.lpstrFile = szFile;
//  //
//  // Set lpstrFile[0] to '\0' so that GetOpenFileName does not 
//  // use the contents of szFile to initialize itself.
//  //
//  //ofn.lpstrFile[0] = '\0';
//  ofn.nMaxFile = sizeof(szFile);
//  ofn.lpstrFilter = "All\0*.*\0Text\0*.TXT\0";
//  ofn.nFilterIndex = 1;
//  ofn.lpstrFileTitle = NULL;
//  ofn.nMaxFileTitle = 0;
//  ofn.lpstrInitialDir = "./";
//  ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
//
//  // Display the Open dialog box. 
//
//  if (GetOpenFileName(&ofn)==TRUE) 
//  {
//    sfile.assign( szFile );
//    return true;
//  }
//  return false;
//}


void _System(const char *command) 
{
  system(command);
}


//void DebugLastError() 
//{ 
//    // Retrieve the system error message for the last-error code
//    LPVOID lpMsgBuf;
//    LPVOID lpDisplayBuf;
//    DWORD dw = GetLastError(); 
//
//    FormatMessage(
//        FORMAT_MESSAGE_ALLOCATE_BUFFER | 
//        FORMAT_MESSAGE_FROM_SYSTEM |
//        FORMAT_MESSAGE_IGNORE_INSERTS,
//        NULL,
//        dw,
//        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
//        (LPTSTR) &lpMsgBuf,
//        0, NULL );
//
//    // Display the error message and exit the process
//
//    lpDisplayBuf = (LPVOID)LocalAlloc(LMEM_ZEROINIT, 
//        (lstrlen((LPCTSTR)lpMsgBuf)+40)*sizeof(TCHAR)); 
//    StringCchPrintf((LPTSTR)lpDisplayBuf, 
//        LocalSize(lpDisplayBuf) / sizeof(TCHAR),
//        TEXT("%d: %s"), 
//        dw, lpMsgBuf); 
//    OutputDebugString((LPCTSTR)lpDisplayBuf); 
//
//    LocalFree(lpMsgBuf);
//    LocalFree(lpDisplayBuf);
//}
#endif


inline int SplitBBoxAlongAxis(const BBox3 &bb, int axis, int p, DynArray<BBox3> &res, int sharedBoundary)
{
  int l=p, r=p;
  if(sharedBoundary<0) --l;
  else if(sharedBoundary>0) ++r;
  if(bb.max[axis]<=l || bb.min[axis]>=r) {
    res.push_back(bb); return 1;
  }
  else {
    //Assert(bb.max[axis]-bb.min[axis]>=2);
    if(bb.min[axis] <= l)
      res.push_back(BBox3().Add(bb.min).Add(replace(bb.max,axis, l))); // the left box
    if(bb.max[axis] >= r)
      res.push_back(BBox3().Add(replace(bb.min,axis, r)).Add(bb.max)); // the right box
    return 2;
  }
}

inline int SplitBBoxAlongAxis(const BBox3 &bb, int axis, int p1, int p2, DynArray<BBox3> &res, int boundarySeparation)
{
  int loffset = 0, roffset = 0;
  switch(boundarySeparation)
  {
  case BBOXB:
    loffset = 1;
    roffset = -1;
    break;
  case BBOXA:
    loffset = -1;
    roffset = 1;
    break;
  }
  myAssert(p1 <= p2);
  int n = 0;
  n += SplitBBoxAlongAxis(bb, axis, p1, res, loffset); // left and middle
  BBox3 right = res.back();
  res.pop_back();
  --n;
  n += SplitBBoxAlongAxis(right, axis, p2, res, roffset); // middle and right;
  return n;
}

inline bool _OverlapsForSplit(const BBox3 &b, const BBox3 &c, int axis)
{
  for(int i=0; i<3; ++i)
  {
    if(i==axis) continue;
    if((b.max[i] <= c.min[i] || b.min[i] >= c.max[i]) && b.max[i]>b.min[i]) return false;
  }
  return true;
}

int SplitBBoxByBox(const BBox3 &bb, const BBox3 &splitbox, DynArray<BBox3> &res, bool includeUnion, int boundarySeparation)
{
  int q;
  int istart = res.size();

  DynArray<BBox3> tmp(27,ConsTags::RESERVE);
  res.push_back(bb);

  for(int axis=0; axis<3; ++axis)
  {
    q = 0;
    for(int i=istart; i<res.size(); ++i)
    {
      if(_OverlapsForSplit(res[i], splitbox, axis))
      {
        q += SplitBBoxAlongAxis(res[i], axis, splitbox.min[axis], splitbox.max[axis], tmp, boundarySeparation);
      }
      else
      {
        q += 1;
        tmp.push_back(res[i]);
      }
    }
    res.resize(istart);
    for(int i=0; i<tmp.size(); ++i)
    {
      if(axis==2 && !includeUnion && splitbox.Includes(tmp[i])) continue;
      res.push_back(tmp[i]);
    }
    tmp.remove_all();
  }
  return q;
}


