

// CC0: Appolloian with a twist II
//  Playing around with shadows in 2D
//  Needed a somewhat more complex distance field than boxes
//  The appolloian fractal turned out quite nice so while
//  similar to an earlier shader of mine I think it's
//  distrinctive enough to share
#define RESOLUTION  RENDERSIZE
#define TIME        TIME
#define MAX_MARCHES 120
#define TOLERANCE   0.00001
#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))
//#define PI          3.141592654
#define TAU         (2.0*PI)

const mat2 rot0 = ROT(0.0);
mat2 g_rot0 = rot0;
mat2 g_rot1 = rot0;
mat2 g_rot2 = rot0;
mat2 g_rot3 = rot0;

// License: Unknown, author: nmz (twitter: @stormoid), found: https://www.shadertoy.com/view/NdfyRM
float sRGB(float t) { return mix(1.055*pow(t, 1./2.4) - 0.1055, 2.92*t, step(t, 0.0021308))*(pow(1.0+highhits*0.5, 2.)); }
// License: Unknown, author: nmz (twitter: @stormoid), found: https://www.shadertoy.com/view/NdfyRM
vec3 sRGB(in vec3 c) { return vec3 (sRGB(c.x), sRGB(c.y), sRGB(c.z)); }

// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
const vec4 hsv2rgb_K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + hsv2rgb_K.xyz) * 7.0 - hsv2rgb_K.www);
  return c.z * mix(hsv2rgb_K.xxx, clamp(p - hsv2rgb_K.xxx, 0.0, 1.0), c.y);
}

float apolloian(vec3 p, float s, out float h) {
  float scale = 1.0;
  for(int i=0; i < 9; ++i) {
    p = -1.0 + 2.0*fract(0.5*p+0.5);
    float r2 = dot(p,p);
    float k  = s/r2;
    p       *= k;
    scale   *= k;
  }
  
  vec3 ap = abs(p/scale);  
  float d = length(ap.xy);
  d = min(d, ap.z);

  float hh = 0.0;
  if (d == ap.z){
    hh += 0.125;
  }
  h = hh;
  return d;
}

float df(vec2 p, out float h) {
  float fz = 1.0-0.0;
  float z = (Z*1.55)*fz;
   z *= Zoom;
  p /= z;
  vec3 p3 = vec3(p,0.1);
  p3.xz*=g_rot0;
  p3.yz*=g_rot3;
  float d = apolloian(p3, 1.0/fz, h);
  d *= z;
  return d;
}

float shadow(vec2 lp, vec2 ld, float mint, float maxt) {
  const float ds = 1.0-0.4;
  float t = mint;
  float nd = 1E7;
  float h;
  const float soff = 0.025;
  const float smul = .5;
  for (int i=0; i < MAX_MARCHES; ++i) {
    vec2 p = lp + ld*t;
    float d = df(p, h);
    if (d < TOLERANCE || t >= maxt) {
      float sd = 1.0-exp(-smul*max(t/maxt-soff, 0.0));
      return t >= maxt ? mix(sd, 1.0, smoothstep(0.0, 0.025, nd)) : sd;
    }
    nd = min(nd, d);
    t += ds*d;
  }
  float sd = 1.0-exp(-smul*max(t/maxt-soff, 0.0));
  return sd;
}

vec3 effect(vec2 p, vec2 q) {
  float aa = 2.0/RESOLUTION.y;
  float a = 0.25*smoothTimeC;
  float b = 0.07225*smoothTimeB;
  float c = 0.07225*smoothTime;

  g_rot0 = ROT(0.25*a); 
  g_rot1 = ROT(sqrt(0.25)*b);
  g_rot2 = ROT(sqrt(0.25)*(a+b)*0.5)*(pow(.70+highhits*0.3, 2.));
  g_rot3 = ROT(sqrt(0.25)*c);

  vec2  lightPos  = vec2(0.0, 1.0);
  lightPos        *= (g_rot1);
  lightPos.xy += PosXY.xy*0.95;

  vec2  lightDiff = lightPos - p;
  float lightD2   = dot(lightDiff,lightDiff);
  float lightLen  = sqrt(lightD2);
  vec2  lightDir  = lightDiff / lightLen;
  vec3  lightPos3 = vec3(lightPos, 0.0);
  vec3  p3        = vec3(p, -1.0);
  float lightLen3 = distance(lightPos3, p3);
  vec3  lightDir3 = normalize(lightPos3-p3);
  vec3  n3        = vec3(0.0, 0.0, 1.0);
  float diff      = max(dot(lightDir3, n3), 0.0);

  float h;
  float d   = df(p, h);
  float ss  = shadow(p,lightDir, 0.0025, lightLen);
  vec3 bcol = hsv2rgb(vec3(fract(h-0.2*length(p)+0.25*smoothTimeB), 0.666, 1.0));

  vec3 col = vec3(0.0);
  col += mix(0., 1.0, diff)*0.5*mix(0.1, 1.0, ss)/(lightLen3*lightLen3);
  col += exp(-(500.0*(1.0+abs(FOV-2.0)))*abs(d))*sqrt(bcol);
  col += exp(-40.0*max(lightLen-0.02, 0.0));
 
  return col;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  vec2 q = fragCoord/RESOLUTION.xy;
  vec2 p = -1. + 2. * q;
  //p.xy+= _uvc.xy*FOV*0.75;
  p.xy += PI*PosXY.xy;
  p.x *= (RESOLUTION.x)/(RESOLUTION.y);
 

  p.xy += _rotate(_uvc.xy*sin(p.xy*PI), smoothTimeC*0.25)*Warp;
  p.xy*=normalize(1.0+Kaleido*vec2(-_rotate(PI*_uvc*p.xy*0.5/(RENDERSIZE.xy*0.5)*RENDERSIZE,smoothTime*0.1+sin(smoothTimeC*0.1)*5.)));
  p.xy+= _uvc.xy*FOV*0.75;
  vec3 col = effect(p, q);
  //col *= mix(0.0, 1.0, smoothstep(0.0, 4.0, smoothTime)); //no fade-in
  col = sRGB(col);
  
  fragColor = vec4(col, 1.0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}