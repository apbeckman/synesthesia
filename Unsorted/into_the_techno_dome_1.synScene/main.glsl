

// License CC0: Into the techno dome
// A simple extension of "Follow the Light" I made earlier.
// The tunnel forks and the fork is picked randomly.
// Thought it turned out nice.
// Based on: https://www.shadertoy.com/view/XsBXWt

#define PI              3.141592654
#define TAU             (2.0*PI)
#define TIME            bass_time
#define RESOLUTION      RENDERSIZE
#define ROT(a)          mat2(cos(a), sin(a), -sin(a), cos(a))
#define TOLERANCE       0.000001
#define MAX_RAY_LENGTH  36.0
#define MAX_RAY_MARCHES 100
#define NORM_OFF        0.001
#define PCOS(x)         (0.5 + 0.5*cos(x))

#define TWISTS

#if defined(TWISTS)
#define PATHA (0.75*vec2(0.1147, 0.2093))
#define PATHB (0.5*vec2(13.0, 3.0))
vec3 cam(float z)  {
    return vec3(sin(z*PATHA)*PATHB, z);
}

vec3 dcam(float z)  {
    return vec3(PATHA*PATHB*cos(PATHA*z), 1.0);
}

vec3 ddcam(float z)  {
    return vec3(-PATHA*PATHA*PATHB*sin(PATHA*z), 0.0);
}
#else
vec3 cam(float z)  {
    return vec3(0.0, 0.0, z);
}

vec3 dcam(float z)  {
    return vec3(0.0, 0.0, 1.0);
}

vec3 ddcam(float z)  {
    return vec3(0.0);
}
#endif

// License: Unknown, author: Unknown, found: don't remember
float tanh_approx(float x) {
//  return tanh(x);
  float x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/articles/smin
float pmin(float a, float b, float k) {
  float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
  return mix( b, a, h ) - k*h*(1.0-h);
}

float pabs(float a, float k) {
  return -pmin(a, -a, k);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/articles/spherefunctions
float sphered(vec3 ro, vec3 rd, vec4 sph, float dbuffer) {
    float ndbuffer = dbuffer/sph.w;
    vec3  rc = (ro - sph.xyz)/sph.w;
  
    float b = dot(rd,rc);
    float c = dot(rc,rc) - 1.0;
    float h = b*b - c;
    if( h<0.0 ) return 0.0;
    h = sqrt( h );
    float t1 = -b - h;
    float t2 = -b + h;

    if( t2<0.0 || t1>ndbuffer ) return 0.0;
    t1 = max( t1, 0.0 );
    t2 = min( t2, ndbuffer );

    float i1 = -(c*t1 + b*t1*t1 + t1*t1*t1/3.0);
    float i2 = -(c*t2 + b*t2*t2 + t2*t2*t2/3.0);
    return (i2-i1)*(3.0/4.0);
}

// "Amazing Surface" fractal
// https://www.shadertoy.com/view/XsBXWt
/*
vec4 formula(vec4 p) {
  p.xz = abs(p.xz+1.)-abs(p.xz-1.)-p.xz;
  p.y-=.25;
  p.xy*=ROT(radians(30.0));
  p=p*2.0/clamp(dot(p.xyz,p.xyz),0.24,1.0);
  return p;
}
*/
vec4 formula(vec4 p) {
  p.xz = abs(p.xz+1.)-abs(p.xz-1.)-p.xz;
  p.y-=.25;
  p.xy*=ROT(radians(30.0));
  p=p*2.0/clamp(dot(p.xyz+cos(smoothTimeC*0.0125)*0.125,p.xyz+sin(smoothTimeC*0.0125)*0.125),(0.7524-Fractal),1.0-0.1*Fractal);
  
  return p;
}

vec3  g_trap0 = vec3(0.0);

float fractal(vec3 pos) {
  vec3 tpos =pos;
  const float sz = 6.0;
  tpos.z    = abs(0.5*sz-mod(tpos.z, sz));
  vec4 p    = vec4(tpos,1.);
  
  vec3 trap0pos = vec3(-2., 0.2, -3.0);
  vec3 trap0 = vec3(1E6);
  
  for (int i=0; i < 3+(Iters); ++i) {
    p = formula(p);
    trap0 = min(trap0, abs(p.xyz-trap0pos));
  }
  g_trap0 = trap0;
  
  float fr=(length(max(vec3(0.),p.xyz-1.5))-1.0)/p.w;

  return fr;
}

float df(vec3 p) {
  // Space distortion found somewhere on shadertoy, don't remember where
  vec3 wrap = cam(p.z);
  vec3 wrapDeriv = normalize(dcam(p.z));
  p.xy -= wrap.xy;
  p -= wrapDeriv*dot(vec3(p.xy, 0), wrapDeriv)*0.5*vec3(1,1,-1);
/*
#if defined(TWISTS)
  vec3 ddcam = ddcam(p.z);
  p.xy *= ROT(-16.0*ddcam.x);
#endif
*/
#if defined(TWISTS)
  vec3 ddcam = ddcam(p.z);
  p.xy *= ROT((-16.0*(1.0+PI*Twist))*ddcam.x);
#endif

  // Splits the tunnel
  const float splitDist = 50.0;
  float mz = mod(p.z, splitDist);
  float n  = floor(p.z/splitDist);
  float h  = hash(n);
  float off = 1.75*smoothstep(15.0, 35.0, mz);

  p.x -= h>0.5 ? off : -off;
  p.x = abs(p.x);
  p.x -= 1.0+off;
  p.y = -pabs(p.y, 1.5);
  p.y -= -1.5;

  return fractal(p); 
}

float rayMarch(vec3 ro, vec3 rd, out int iter) {
  float t = 0.0;
  int i = 0;
  for (i = 0; i < MAX_RAY_MARCHES; i++) {
    float d = df(ro + rd*t);
    if (d < TOLERANCE || t > MAX_RAY_LENGTH) break;
    t += d;
  }
  iter = i;
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

vec3 render(vec3 ro, vec3 rd) {
  vec3 lightPos = cam(ro.z+20.0);
  float alpha   = 0.05*TIME;
  
  const vec3 skyCol = vec3(0.0);

  int iter    = 0;
  float t     = rayMarch(ro, rd, iter);
  vec3  trap0 = g_trap0;
  float highs =pow ((syn_HighLevel*0.5+0.5*syn_MidHighLevel+syn_Level*0.25)*(1.0+0.5*syn_Intensity), 2.);
  float pulse = smoothstep(0.1, 1.0, (sin(smoothTime)*0.)+(TAU*highs));

  //float pulse = smoothstep(0.0, 1.0, sin(TAU*TIME*0.25));
  float sr    = mix(2.0, 3.0,  pulse);
  float sd    = sphered(ro, rd, vec4(lightPos, sr), t);

  vec3 bgcol  = vec3(2.750, 2.730, 2.975).zyx;

  //vec3 bgcol  = vec3(2.0, 1.0, 0.75).zyx;
 // bgcol.xz = _rotate(bgcol.xz*0.4+trap0.xz*(sin(_uvc)/(PI*PI)), sin(_uvc.x)*pow(PI, -2));

  vec3 gcol   = mix(1.0, 1.75, pulse)*sd*sd*bgcol;

  if (t >= MAX_RAY_LENGTH) {
    return gcol;
  }

  vec3 pos    = ro + t*rd;
  vec3 nor    = normal(pos);
  vec3 refl   = reflect(rd, nor);
  float ii    = float(iter)/float(MAX_RAY_MARCHES);

  vec3 lv   = lightPos - pos;
  float ll2 = dot(lv, lv);
  float ll  = sqrt(ll2);
  vec3 ld   = lv / ll;

  float fre = abs(dot(rd, nor));
  fre *= fre;
  fre *= fre;
  float dm  = 4.0/ll2;
  float dif = pow(max(dot(nor,ld),0.0), 1.0);  
  float spe = fre*pow(max(dot(refl, ld), 0.), 10.);
  float fo  = smoothstep(0.9, 0.4, t/MAX_RAY_LENGTH);
  float ao  = 1.0-ii;

  vec3 col = vec3(0.0);
  col += pow(smoothstep(0.5, 1.0, trap0.x*0.25)*1.3, mix(6.0, 2.0, pulse+highs*0.1))*0.5*bgcol*mix(0.75, 2.25, pulse);
  col.xy = 2.0*_rotate(col.xy*0.4*trap0.xz*(sin(_uvc)/pow(PI, 2.5)), sin(_uvc.x)*pow(PI, -2));

  //col += smoothstep(0.7, 0.6, trap0.z)*smoothstep(0.4, 0.5, trap0.z)*ao*bgcol*mix(0.2, 1.4, pulse);
  col += spe*bgcol*mix(0.66, 1.75, pulse);
  col *= 1.0-sd*sd;
  col *= fo;
  //col += gcol;
  return col;
}

vec3 effect3d(vec2 p, vec2 q) {
  float z   = TIME*2.5;
  
  vec3 cam = cam(z);
  vec3 dcam = dcam(z);
  vec3 ddcam= ddcam(z);
 // p.y*=1.0-0.75* Mirror.y;
 // p.x*=1.0-0.75* Mirror.x;

  vec3 ro = cam;
  vec3 ww = normalize(dcam);
  vec3 uu = normalize(cross(vec3(0.0,1.0,0.0)+ddcam*06.0, ww ));
  vec3 vv = normalize(cross(ww,uu));
  float fov = (2.0-FOV)/tanh(TAU/6.0);
  vec3 rd = normalize(-p.x*uu + p.y*vv + fov*ww );
    rd.yz = _rotate(rd.yz, -lookXY.y*PI+_uvc.y*PI*Perspective.y*0.5);
    rd.xz = _rotate(rd.xz, lookXY.x*PI+_uvc.x*PI*Perspective.x);

  return render(ro, rd);
}

// License: Unknown, author: nmz (twitter: @stormoid), found: https://www.shadertoy.com/view/NdfyRM
float sRGB(float t) { return mix(1.055*pow(t, 1./2.4) - 0.055, 12.92*t, step(t, 0.0031308)); }
// License: Unknown, author: nmz (twitter: @stormoid), found: https://www.shadertoy.com/view/NdfyRM
vec3 sRGB(in vec3 c) { return vec3 (sRGB(c.x), sRGB(c.y), sRGB(c.z)); }

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  vec2 q = fragCoord/RESOLUTION.xy;
  vec2 p = -1. + 2. * q;
  p.x *= RESOLUTION.x/RESOLUTION.y;
  vec3 col = effect3d(p, q);
  col = sRGB(col);
  fragColor = vec4(col, 1.0);
	return fragColor; 
 } 




vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}