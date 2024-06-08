

			//******** Common Code Begins ********

// CC0: Easter tweaks
//  Easter holiday and I got some time to tinker

#define TIME        TIME
#define RESOLUTION  RENDERSIZE
#define PI          3.141592654
#define TAU         (2.0*PI)

// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
const vec4 hsv2rgb_K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www);
  return c.z * mix(hsv2rgb_K.xxx, clamp(p - hsv2rgb_K.xxx, 0.0, 1.0), c.y);
}
// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
//  Macro version of above to enable compile-time constants
#define HSV2RGB(c)  (c.z * mix(hsv2rgb_K.xxx, clamp(abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www) - hsv2rgb_K.xxx, 0.0, 1.0), c.y))

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

// License: MIT OR CC-BY-NC-4.0, author: mercury, found: https://mercury.sexy/hg_sdf/
float mod1(inout float p, float size) {
  float halfsize = size*0.5;
  float c = floor((p + halfsize)/size);
  p = mod(p + halfsize, size) - halfsize;
  return c;
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/articles/intersectors/
float rayPlane(vec3 ro, vec3 rd, vec4 p) {
  return -(dot(ro,p.xyz)+p.w)/dot(rd,p.xyz);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/articles/distfunctions/
float box(vec2 p, vec2 b) {
  vec2 d = abs(p)-b;
  return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}
  
// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float hexagon(vec2 p, float r, float rr) {
  r -= rr;
  const vec3 k = 0.5*vec3(-sqrt(3.0), 1.0, sqrt(4.0/3.0));
  p = abs(p);
  p -= 2.0*min(dot(k.xy,p),0.0)*k.xy;
  p -= vec2(clamp(p.x, -k.z*r, k.z*r), r);
  return length(p)*sign(p.y)-rr;
}




			//******** BuffA Code Begins ********

// CC0: Easter tweaks
//  Easter holiday and I got some time to tinker

#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))

#define TOLERANCE           0.0001
#define MAX_RAY_LENGTH      40.0
#define MAX_RAY_MARCHES     80
#define NORM_OFF            0.005

const vec3 sunDir       = normalize(vec3(0.,0., 1.));
const vec3 sunCol       = HSV2RGB(vec3(0.6, .95, 1E-2))*1.;
const vec3 bottomBoxCol = HSV2RGB(vec3(0.7, 0.80, 0.5))*1.;
const vec3 topBoxCol    = HSV2RGB(vec3(0.57, 0.90, 1.))*1.;
const vec3 glowCol      = HSV2RGB(vec3(0.98, 0.825, 2E-3))*1.;
  
float g_gd;
float g_gz;

float segmentz(vec3 p, float r, float l) {
  p.z = abs(p.z)-l*0.5;
  float d0 = length(p.xy)-r;
  float d1 = length(p)-r;
  return p.z > 0. ? d1: d0;
}


float df(vec3 p, float t) {
  float gz = p.z;
  p.z -= -PI*TIME;
  mat2 r = ROT(-0.1*p.z);
  vec2 p0 = p.xy;
  p0 *= r;
  vec2 s0 = sign(p0);
  p0 = abs(p0);
  p0 -= 1.;
  float d0 = length(p0) - 0.1;
  vec2 p1 = p0;
  p1 -= 0.25;
  float off = 2.*s0.x+1.*s0.y;
  float p1z = p.z+TIME*8.+4.*off;
  float n1z = mod1(p1z, 24.0);
  vec3 p1_ = vec3(p1, p1z);
  float d1 = segmentz(p1_, 0.01, 8.);
  float d = d0;

  d = min(d, d1);

  if (d1 < g_gd) {
    g_gd = d1;
    g_gz = t;
  }
  
  float d2 = hexagon(p.xy*r, 2., .5);
  d = min(d, -d2);
  return d;
}

float df(vec3 p) {
  return df(p, 1.);
}

// Back stepping technique by IQ
#define BACKSTEP
float rayMarch(vec3 ro, vec3 rd, float tinit) {
  float t = tinit;
  const float tol = TOLERANCE;
#if defined(BACKSTEP)
  vec2 dti = vec2(1e10,0.0);
#endif  
  int i = 0;
  for (i = 0; i < MAX_RAY_MARCHES; ++i) {
    float d = df(ro + rd*t,t);
#if defined(BACKSTEP)
    if (d<dti.x) { dti=vec2(d,t); }
#endif  
    if (d < TOLERANCE || t > MAX_RAY_LENGTH) {
      break;
    }
    t += d;
  }
#if defined(BACKSTEP)
  if(i==MAX_RAY_MARCHES) { t=dti.y; };
#endif  
  return t;
}


vec3 normal(vec3 pos) {
  vec2  eps = vec2(NORM_OFF,0.0);
  vec3 nor;
  nor.x = df(pos+eps.xyy) - df(pos-eps.xyy);
  nor.y = df(pos+eps.yxy) - df(pos-eps.yxy);
  nor.z = df(pos+eps.yyx) - df(pos-eps.yyx);
  return normalize(nor);
}

vec3 render0(vec3 ro, vec3 rd) {
  vec3 col = vec3(0.0);
  
  float tp0  = rayPlane(ro, rd, vec4(vec3(0.0, -1.0, 0.0), -5.0));
  float tp1  = rayPlane(ro, rd, vec4(vec3(0.0, -1.0, 0.0), 6.0));

  if (tp0 > 0.0) {
    col += bottomBoxCol*exp(-0.5*(length((ro + tp0*rd).xz)));
  }

  if (tp1 > 0.0) {
    vec3 pos  = ro + tp1*rd;
    vec2 pp = pos.xz;
    float db = box(pp, vec2(5.0, 9.0))-3.0;
    
    col += topBoxCol*rd.y*rd.y*smoothstep(0.25, 0.0, db);
    col += 0.2*topBoxCol*exp(-0.5*max(db, 0.0));
    col += 0.05*sqrt(topBoxCol)*max(-db, 0.0);
  }

  col += sunCol/(1.00035-dot(sunDir, rd));
  return clamp(col, 0., 10.); 
}

vec3 render1(vec3 ro, vec3 rd) {
  int ii;
  g_gd = 1E3;
  float t = rayMarch(ro, rd, 0.);
  float gd = g_gd;
  float gz = g_gz;

  vec3 sky = render0(ro, rd);
  vec3 col = sky;
  
  float fo = exp(-0.125*max(t-10., 0.));
  float gfo = exp(-0.125*max(gz-10., 0.));
  
  if (t < MAX_RAY_LENGTH) {
    vec3 p = ro+rd*t;
    vec3 n = normal(p);
    vec3 r = reflect(rd, n);
    vec3 rcol = render0(p, r);

    vec3 drx = dFdx(r);
    vec3 dry = dFdy(r);
    float drxy = (dot(drx,drx)+dot(dry,dry));
    

    g_gd = 1E3;
    float rt = rayMarch(p, r, 0.5);
    float rgd = g_gd;
    float rgz = g_gz;

    float fre = 1.0+dot(rd,n);
    fre *= fre;
    fre = mix(0.3, 1.0, fre);
    
    col = vec3(0.);
    col += rcol*fre;
    col += smoothstep(2E-4, 0., drxy)*4.*gfo*glowCol/max(rgd, 4E-3)*fre;
  }

  col = mix(sky, col, fo);
  col += gfo*glowCol/max(gd*abs(gd), 2E-4);

  return col;
}

vec3 effect(vec2 p) {
  const float fov = 2.;

  const vec3 up = vec3(0., 1., 0.);
  const vec3 ro   = vec3(0.0, .5, -3.0);
  const vec3 la   = vec3(0.0);

  const vec3 ww = normalize(la-ro);
  const vec3 uu = normalize(cross(up, ww));
  const vec3 vv = cross(ww, uu);
  vec3 rd = normalize(-p.x*uu + p.y*vv + fov*ww);

  vec3 col = vec3(0.0); 
  col = render1(ro, rd);
  col -= 1E-2*vec3(1.,2.,3.).zyx*(length(p)+0.25); 
  col = aces_approx(col);
  col = sqrt(col);
  return col;
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  vec2 q = fragCoord/RESOLUTION.xy;
  vec2 p = -1. + 2. * q;
  p.x *= RESOLUTION.x/RESOLUTION.y;

  vec3 col = effect(p);

  fragColor = vec4(col, 1.0);
	return fragColor; 
 } 



// CC0: Easter tweaks
//  Easter holiday and I got some time to tinker

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  vec2 isz = 1./RESOLUTION.xy;
  vec2 q = fragCoord*isz;
  fragColor = texture(BuffA, q);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}