

// CC0: Time for some inner reflections
// An evolution of: https://www.shadertoy.com/view/7dKBDt
// After a tip from shane found knighty's shader: https://www.shadertoy.com/view/MsKGzw
// Knighty's shaders allows one to experiment with many cool polyhedras

// Original inspiration from:https://www.youtube.com/watch?v=qNoQXF2dKBs

// ------------------------------------------------------------------------------------
// Here are some parameters to experiment with

#define INNER_SPHERE
//#define GOT_BEER

const float poly_U        = 1.0;  // [0, inf]
const float poly_V        = 1.0;  // [0, inf]
const float poly_W        = 2.0;  // [0, inf]
const int   poly_type     = 5;    // [2, 5]

float zoom = 3.0;
// ------------------------------------------------------------------------------------

//#define PI          3.141592654
#define TAU         (2.0*PI)

#define RESOLUTION  RENDERSIZE

#define TOLERANCE       0.00001
#define MAX_RAY_LENGTH  20.0
#define MAX_RAY_MARCHES 60
#define NORM_OFF        0.001
#define ROT(a)          mat2(cos(a), sin(a), -sin(a), cos(a))
#define MAX_BOUNCES     6

// License: Unknown, author: knighty, found: https://www.shadertoy.com/view/MsKGzw
const float poly_cospin   = cos(PI/float(poly_type));
const float poly_scospin  = sqrt(0.75-poly_cospin*poly_cospin);
const vec3  poly_nc       = vec3(-0.5, -poly_cospin, poly_scospin);
const vec3  poly_pab      = vec3(0., 0., 1.);
const vec3  poly_pbc_     = vec3(poly_scospin, 0., 0.5);
const vec3  poly_pca_     = vec3(0., poly_scospin, poly_cospin);
const vec3  poly_p        = normalize((poly_U*poly_pab+poly_V*poly_pbc_+poly_W*poly_pca_));
const vec3  poly_pbc      = normalize(poly_pbc_);
const vec3  poly_pca      = normalize(poly_pca_);

const float initt = 0.125; 

// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
const vec4 hsv2rgb_K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www);
  return c.z * mix(hsv2rgb_K.xxx, clamp(p - hsv2rgb_K.xxx, 0.0, 1.0), c.y);
}
// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
//  Macro version of above to enable compile-time constants
#define HSV2RGB(c)  (c.z * mix(hsv2rgb_K.xxx, clamp(abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www) - hsv2rgb_K.xxx, 0.0, 1.0), c.y))


// License: Unknown, author: nmz (twitter: @stormoid), found: https://www.shadertoy.com/view/NdfyRM
vec3 sRGB(vec3 t) {
  return mix(1.055*pow(t, vec3(1./2.4)) - 0.055, 12.92*t, step(t, vec3(0.0031308)));
}

// License: Unknown, author: Matt Taylor (https://github.com/64), found: https://64.github.io/tonemapping/
vec3 aces_approx(vec3 v) {
  v = max(v, 0.0);
  v *= 0.6f;
  float a = 2.51f;
  float b = 0.03f;
  float c = 2.43f;
  float d = 0.59f;
  float e = 0.14f;
  return clamp((v*(a*v+b))/(v*(c*v+d)+e), 0.0f, 1.0f);
}

// License: Unknown, author: Unknown, found: don't remember
float tanh_approx(float x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  float x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}

// http://mercury.sexy/hg_sdf/
vec2 mod2(inout vec2 p, vec2 size) {
  vec2 c = floor((p + size*0.5)/size);
  p = mod(p + size*0.5,size) - size*0.5;
  return c;
}

float circle(vec2 p, float r) {
  return length(p) - r;
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float box(vec2 p, vec2 b) {
  vec2 d = abs(p)-b;
  return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}

// License: MIT, author: Inigo Quilez, found: https://www.iquilezles.org/www/articles/smin/smin.htm
float pmin(float a, float b, float k) {
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/intersectors/intersectors.htm
float rayPlane(vec3 ro, vec3 rd, vec4 p) {
  return -(dot(ro,p.xyz)+p.w)/dot(rd,p.xyz);
}

// License: Unknown, author: knighty, found: https://www.shadertoy.com/view/MsKGzw
void poly_fold(inout vec3 pos) {
  vec3 p = pos;

  for(int i = 0; i < poly_type; ++i){
    p.xy  = abs(p.xy);
    p    -= 2.*min(0., dot(p,poly_nc)) * poly_nc;
  }
  
  pos = p;
}

float poly_plane(vec3 pos) {
  float d0 = dot(pos, poly_pab);
  float d1 = dot(pos, poly_pbc);
  float d2 = dot(pos, poly_pca);
  float d = d0;
  d = max(d, d1);
  d = max(d, d2);
  return d;
}

float poly_corner(vec3 pos) {
  float d = length(pos) - .1;
  return d;
}

float dot2(vec3 p) {
  return dot(p, p);
}

float poly_edge(vec3 pos) {
  float dla = dot2(pos-min(0., pos.x)*vec3(1., 0., 0.));
  float dlb = dot2(pos-min(0., pos.y)*vec3(0., 1., 0.));
  float dlc = dot2(pos-min(0., dot(pos, poly_nc))*poly_nc);
  return sqrt(min(min(dla, dlb), dlc))-0.0125;
}

float poly_planes(vec3 pos, out vec3 pp) {
  poly_fold(pos);
  pos -= poly_p;

  pp = pos;
  return poly_plane(pos);
}

float poly_edges(vec3 pos, out vec3 pp) {
  poly_fold(pos);
  pos -= poly_p;

  pp = pos;
  return poly_edge(pos);
}


float blobs(vec2 p) {
  // Generates a grid of dots
  vec2 bp = p;
  vec2 bn = mod2(bp, vec2(3.0));

  vec2 dp = p;
  vec2 dn = mod2(dp, vec2(0.25));
  float ddots = length(dp);
  
  // Blobs
  float dblobs = 1E6;
  for (int i = 0; i < 5; ++i) {
    float dd = circle(bp-1.0*vec2(sin(smoothTimeC*0.5+float(i)), sin(float(i*i)+smoothTimeC*sqrt(0.5))), 0.1);
    dblobs = pmin(dblobs, dd, 0.35*(1.0+basshits));
  }

  float d = 1E6;
  d = min(d, ddots);
  // Smooth min between blobs and dots makes it look somewhat amoeba like
  d = pmin(d, dblobs, 0.35);
  return d;
}

vec3 skyColor(vec3 ro, vec3 rd) {
  const vec3 gcol = HSV2RGB(vec3(0.45, 0.6, 1.0));
  vec3 col = clamp(vec3(0.0025/abs(rd.y))*gcol, 0.0, 1.0);
  
  float tp0  = rayPlane(ro, rd, vec4(vec3(0.0, 1.0, 0.0), 14.0));
  float tp1  = rayPlane(ro, rd, vec4(vec3(0.0, -1.0, 0.0), 16.0));
  float tp = tp1;
  tp = max(tp0,tp1);
  if (tp > 0.0) {
    vec3 pos  = ro + tp*rd;
    const float fz = 0.25;
    const float bz = 1.0/fz;
    vec2 bpos = pos.xz/bz;
    float db = blobs(bpos)*bz;
    db = abs(db);
    vec2 pp = pos.xz*fz;
    float m = 0.5+0.25*(sin(3.0*pp.x+smoothTimeC*2.1)+sin(3.3*pp.y+smoothTimeC*2.0));
    m *= m;
    m *= m;
    pp = fract(pp+0.5)-0.5;
    float dp = pmin(abs(pp.x), abs(pp.y), 0.125);
    dp = min(dp, db);
    vec3 hsv = vec3(0.4+mix(0.15,0.0, m), tanh_approx(mix(50.0, 10.0, m)*dp), 1.0);
    vec3 pcol = 1.5*hsv2rgb(hsv)*exp(-mix(30.0, 10.0, m)*dp);
    
    float f = 1.0-tanh_approx(0.1*length(pos.xz));
    col = mix(col, pcol , f);
  }


  if (tp1 > 0.0) {
    vec3 pos  = ro + tp1*rd;
    vec2 pp = pos.xz;
    float db = box(pp, vec2(6.0, 9.0))-1.0;
    
    col += vec3(2.0)*gcol*rd.y*smoothstep(0.25, 0.0, db);
    col += vec3(0.8)*gcol*exp(-0.5*max(db, 0.0));
  }


  return col;
}

float dfExclusion(vec3 p, out vec3 pp) {
  return -poly_edges(p/zoom, pp)*zoom;
}

float shape(vec3 p) {
  
  vec3 pp;
  return poly_planes(p/zoom, pp)*zoom;
}

float df0(vec3 p) {
  float d0 = shape(p);
  float d = d0;
  return d;
}

float df1(vec3 p) {
  float d0 = -shape(p);
  float d = d0;
#if defined(INNER_SPHERE)
  float d1 = length(p) - 2.;
  d = min(d, d1);
#endif
  return d;
}

vec3 normal1(vec3 pos) {
  vec2  eps = vec2(NORM_OFF,0.0);
  vec3 nor;
  nor.x = df1(pos+eps.xyy) - df1(pos-eps.xyy);
  nor.y = df1(pos+eps.yxy) - df1(pos-eps.yxy);
  nor.z = df1(pos+eps.yyx) - df1(pos-eps.yyx);
  return normalize(nor);
}

float rayMarch1(vec3 ro, vec3 rd) {
  float t = 0.0;
  for (int i = 0; i < MAX_RAY_MARCHES; i++) {
    if (t > MAX_RAY_LENGTH) {
      t = MAX_RAY_LENGTH;    
      break;
    }
    float d = df1(ro + rd*t);
    if (d < TOLERANCE) {
      break;
    }
    t  += d;
  }
  return t;
}

vec3 normal0(vec3 pos) {
  vec2  eps = vec2(NORM_OFF,0.0);
  vec3 nor;
  nor.x = df0(pos+eps.xyy) - df0(pos-eps.xyy);
  nor.y = df0(pos+eps.yxy) - df0(pos-eps.yxy);
  nor.z = df0(pos+eps.yyx) - df0(pos-eps.yyx);
  return normalize(nor);
}

float rayMarch0(vec3 ro, vec3 rd) {
  float t = 0.0;
  for (int i = 0; i < MAX_RAY_MARCHES; i++) {
    if (t > MAX_RAY_LENGTH) {
      t = MAX_RAY_LENGTH;    
      break;
    }
    float d = df0(ro + rd*t);
    if (d < TOLERANCE) {
      break;
    }
    t  += d;
  }
  return t;
}

vec3 render1(vec3 ro, vec3 rd) {
  vec3 agg = vec3(0.0, 0.0, 0.0);
  float tagg = initt;
  vec3 ragg = vec3(1.0);

  for (int bounce = 0; bounce < MAX_BOUNCES; ++bounce) {
    float mragg = max(max(ragg.x, ragg.y), ragg.z);
    if (mragg < 0.1) break;
    float st = rayMarch1(ro, rd);
    tagg += st;
    vec3 sp = ro+rd*st;
    vec3 spp;
    float de = dfExclusion(sp, spp);
    vec3 sn = normal1(sp);
    
    float si = cos(5.0*TAU*zoom*spp.z-0.5*sp.y+smoothTimeC);
    const vec3 lcol = vec3(1.0, 1.5, 2.0)*0.8;
    float lf = mix(0.0, 1.0, smoothstep(0., 0.9, si));
    
    vec3 gcol = ragg*lcol*exp(8.0*(min(de-0.2, 0.0)));
    // Will never miss
    if (de < 0.0) {
      agg += gcol;
      ragg *= vec3(0.5, 0.6,0.8);
    } else {
      agg += gcol*lf;
      agg += ragg*lcol*1.5*lf;
      ragg = vec3(0.0);
    }
    
    rd = reflect(rd, sn);
    ro = sp+initt*rd;
    tagg += initt;
  }
#if defined(GOT_BEER)
  return agg*exp(-.5*vec3(0.3, 0.15, 0.1)*tagg);
#else  
  return agg;
#endif
}

vec3 render0(vec3 ro, vec3 rd) {
  vec3 skyCol = skyColor(ro, rd);

  vec3 col = skyCol;

  float st = rayMarch0(ro, rd);
  vec3 sp = ro+rd*st;
  vec3 sn = normal0(sp);
    vec3 spp;
  float de = dfExclusion(sp, spp);
  float ptime = mod(smoothTimeC, 30.0);
  if (st < MAX_RAY_LENGTH) {
    float sfre = 1.0+dot(rd, sn);
    sfre *= sfre;
    sfre = mix(0.1, 1.0, sfre); 
    vec3 sref   = reflect(rd, sn);
    vec3 srefr  = refract(rd, sn, 0.9);
    vec3 ssky = sfre*skyColor(sp, sref);

    if (de > 0.0) {
      col = ssky;
    } else {
      col = 0.5*sfre*ssky;
      vec3 col1 = (1.0-sfre)*render1(sp+srefr*initt, srefr);
      col += col1;
    }
    
  }

  return col;
}

vec3 effect(vec2 p) {
  p += Mirror*p*_uvc*PI*PI;
  p += _uvc*PI*Zoom;

  vec3 ro = 0.8*vec3(0.0, 4.0, 5.0);
  const vec3 la = vec3(0.0, 0.0, 0.0);
  const vec3 up = vec3(0.0, 1.0, 0.0);
  float a = 0.5*(-0.5+0.5*sin(0.123*smoothTime));
  float b = 0.1*smoothTime;
  /*
  if (_mouse.x > 0.0) {
    vec2 m = (2.0*_mouse.xy-RENDERSIZE.xy)/RENDERSIZE.y;
    // Get angle from mouse position
    a =-2.0*m.y;
    b =-2.0*m.x;
  }
  */
  ro.yz *= ROT(a);
  ro.xz *= ROT(b);
  ro.xy += Kaleidescope*(_uvc*PI);
  vec3 ww = normalize(la - ro);
  vec3 uu = normalize(cross(up, ww ));
  vec3 vv = normalize(cross(ww,uu));
  float fov = tan(TAU/6.0);
  vec3 rd = normalize(-p.x*uu + p.y*vv + fov*ww);
  rd.xy +=(_uvc*PI*rd.xy)*Flip;

  vec3 col = render0(ro, rd);
  
  return col;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  vec2 q = fragCoord/RENDERSIZE.xy;
  vec2 p = -1. + 2. * q;
  p.x *= RESOLUTION.x/RESOLUTION.y;
  vec3 col = vec3(0.0);
  col = effect(p);
  col *= smoothstep(0.0, 1.0, smoothTimeB);
  col = aces_approx(col); 
  col = sRGB(col);
  fragColor = vec4(col, 1.0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}