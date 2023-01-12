

// License CC0: Voronoi Glass Panes
//  Tinkering with the voronoi pattern from a Shane shader on saturday evening
//  Kind of hackish overall but that never stopped me from sharing a shader before
//  Shane shader, it's worth a look: https://www.shadertoy.com/view/Mld3Rn

#define RESOLUTION  RENDERSIZE
#define TIME        TIME
#define PI          3.141592654
#define TAU         (2.0*PI)
#define PI_2        (0.5*3.141592654)
#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))
#define DOT2(x)     dot(x,x)

const float planeDist = 1.0-0.25;

float g_hmul = 1.0;

// License: Unknown, author: nmz (twitter: @stormoid), found: https://www.shadertoy.com/view/NdfyRM
float sRGB(float t) { return mix(1.055*pow(t, 1./2.4) - 0.055, 12.92*t, step(t, 0.0031308)); }
// License: Unknown, author: nmz (twitter: @stormoid), found: https://www.shadertoy.com/view/NdfyRM
vec3 sRGB(in vec3 c) { return vec3 (sRGB(c.x), sRGB(c.y), sRGB(c.z)); }

// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
const vec4 hsv2rgb_K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www);
  return c.z * mix(hsv2rgb_K.xxx, clamp(p - hsv2rgb_K.xxx, 0.0, 1.0), c.y);
}
// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
//  Macro version of above to enable compile-time constants
#define HSV2RGB(c)  (c.z * mix(hsv2rgb_K.xxx, clamp(abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www) - hsv2rgb_K.xxx, 0.0, 1.0), c.y))

const float zoomOuter = .75;
const float zoomInner = 0.12;

float tanh_approx(float x) {
//  return tanh(x);
  float x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}

// License: Unknown, author: Unknown, found: don't remember
vec4 alphaBlend(vec4 back, vec4 front) {
  float w = front.w + back.w*(1.0-front.w);
  vec3 xyz = (front.xyz*front.w + back.xyz*back.w*(1.0-front.w))/w;
  return w > 0.0 ? vec4(xyz, w) : vec4(0.0);
}

// License: Unknown, author: Unknown, found: don't remember
vec3 alphaBlend(vec3 back, vec4 front) {
  return mix(back, front.xyz, front.w);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

vec2 hash2(vec2 p) {
  p = vec2(dot (p, vec2 (127.1, 311.7)), dot (p, vec2 (269.5, 183.3)));
  return fract (sin (p)*43758.5453123);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/articles/distfunctions2d
float hex(vec2 p, float r ) {
  const vec3 k  = 0.5*vec3(-sqrt(3.0), 1.0, sqrt(4.0/3.0));
  p = abs(p);
  p -= 2.0*min(dot(k.xy,p),0.0)*k.xy;
  p -= vec2(clamp(p.x, -k.z*r, k.z*r), r);
  return length(p)*sign(p.y);
}

vec3 offset(float z) {
  float a = z;
  vec2 p = -0.1*(vec2(cos(a), sin(a*sqrt(2.0))) + vec2(cos(a*sqrt(0.75)), sin(a*sqrt(0.5))));
  return vec3(p, z);
}

vec3 doffset(float z) {
  float eps = 0.05;
  return (offset(z + eps) - offset(z - eps))/(2.0*eps);
}

vec3 ddoffset(float z) {
  float eps = 0.05;
  return (doffset(z + eps) - doffset(z - eps))/(2.0*eps);
}

vec3 skyColor(vec3 ro, vec3 rd) {
  float ld = max(dot(rd, vec3(0.0, 0.0, 1.0)), 0.0);
  vec3 scol = HSV2RGB(vec3(0.1, 0.25, 0.9));
  return scol*tanh_approx(3.0*pow(ld, 100.0));
}

float voronoi2(vec2 p){
  vec2 g = floor(p), o; p -= g;
  
  vec3 d = vec3(1);
    
  for(int y = -1; y <= 1; y++){
    for(int x = -1; x <= 1; x++){
      o = vec2(x, y);
      o += hash2(g + o) - p;
      d.z = dot(o, o); 
      d.y = max(d.x, min(d.y, d.z));
      d.x = min(d.x, d.z); 
    }
  }
  
  return max(d.y/1.2 - d.x, 0.0)/1.2;
}

float hf2(vec2 p) {
  const float zo = zoomOuter;
  const float zi = zoomInner;
  
  p /= zo;
  p /= zi;
  
  float d = -voronoi2(p);
  d *= zi*zo;
  
  float h = 0.2*tanh_approx(3.0*smoothstep(0.0, 1.0*zo*zi, -d));
  
  return h*zo*zi;
}

float height(vec2 p) {
  return -hf2(p)*g_hmul;
}

vec3 normal(vec2 p, float eps) {
  vec2 v;
  vec2 w;
  vec2 e = vec2(0.00001, 0);
  
  vec3 n;
  n.x = height(p + e.xy) - height(p - e.xy);
  n.y = height(p + e.yx) - height(p - e.yx);
  n.z = -2.0*e.x;
  
  return normalize(n);
}

vec4 plane(vec3 pro, vec3 ro, vec3 rd, vec3 pp, vec3 off, float aa, float n_, out vec3 pnor) {
  float h0 = hash(n_);
  float h1 = fract(7793.0*h0);
  float h2 = fract(6337.0*h0);
  
  vec2 p = (pp-off*vec3(1.0, 1.0, 0.0)).xy;
  const float s = 1.0;
  vec3 lp1 = vec3(5.0,  1.0, 0.0)*vec3(s, 1.0, s)+pro;
  vec3 lp2 = vec3(-5.0, 1.0, 0.0)*vec3(s, 1.0, s)+pro;
  const float hsz = 0.2;
  float hd = hex(p.yx, hsz);

  g_hmul = smoothstep(0.0, 0.125, (hd-hsz/2.0));

  p += vec2(h0,h1)*20.0;
  p *= mix(0.5, 1.0, h2);
  float he  = height(p);
  vec3  nor = normal(p,2.0*aa);
  vec3 po   = pp;

  pnor = nor;

  vec3 ld1 = normalize(lp1 - po);
  vec3 ld2 = normalize(lp2 - po);
  
  float diff1 = max(dot(nor, ld1), 0.0);
  float diff2 = max(dot(nor, ld2), 0.0);
  diff1 = ld1.z*nor.z;;

  vec3  ref   = reflect(rd, nor);
  float ref1  = max(dot(ref, ld1), 0.0);
  float ref2  = max(dot(ref, ld2), 0.0);

  const vec3 mat   = HSV2RGB(vec3(0.55, 0.45, 0.05));
  const vec3 lcol1 = HSV2RGB(vec3(0.6, 0.25, 0.9));
  const vec3 lcol2 = HSV2RGB(vec3(0.1, 0.65, 0.9));
  
  float hf = smoothstep(0.0, 0.00012, -he);
  vec3 lpow1 = 1.0*lcol1/DOT2(ld1);
  vec3 lpow2 = 1.0*lcol2/DOT2(ld2);
  vec3 col = vec3(0.0);
  col += hf*mat*diff1*diff1*lpow1;
  col += hf*mat*diff2*diff2*lpow2;
  float spes = 20.0;
  col += pow(ref1, spes)*lcol1;
  col += pow(ref2, spes)*lcol2;

  float t = 1.0;
  t *= smoothstep(aa, -aa, -(hd-hsz/2.0)); //last term is center thiccness
  t *= mix(1.0, 0.75, hf);
  
  return vec4(col, t);
}

vec3 color(vec3 ww, vec3 uu, vec3 vv, vec3 pro, vec3 ro, vec2 p) {
  float lp = length(p);
  vec2 np = p + 1.0/RESOLUTION.xy;
  float rdd = 2.0+tanh_approx(length(0.25*p));
  
  vec3 rd = normalize(p.x*uu + p.y*vv + rdd*ww);
  vec3 nrd = normalize(np.x*uu + np.y*vv + rdd*ww);

  const int furthest = 5;
  const int fadeFrom = max(furthest-2, 0);

  const float fadeDist = planeDist*float(furthest - fadeFrom);
  float nz = floor(ro.z / planeDist);

  vec3 skyCol = skyColor(ro, rd);

  vec4 acol = vec4(0.0);
  const float cutOff = 0.98;
  bool cutOut = false;

  // Steps from nearest to furthest plane and accumulates the color 
  for (int i = 1; i <= furthest; ++i) {
    float pz = planeDist*nz + planeDist*float(i);

    float pd = (pz - ro.z)/rd.z;

    if (pd > 0.0 && acol.w < cutOff) {
      vec3 pp = ro + rd*pd;
      vec3 npp = ro + nrd*pd;

      float aa = 3.0*length(pp - npp);

      vec3 off = offset(pp.z);

      vec3 pnor = vec3(0.0);
      vec4 pcol = plane(pro, ro, rd, pp, off, aa, nz+float(i), pnor);
      
      vec3 refr = refract(rd, pnor, 1.0-0.075);
      if (pcol.w > (1.0-cutOff)&&refr != vec3(0.0)) {
        rd = refr;
      }

      float dz = pp.z-ro.z;
      const float fi = -0.; 
      float fadeIn = smoothstep(planeDist*(float(furthest)+fi), planeDist*(float(fadeFrom)-fi), dz);
      float fadeOut = smoothstep(0.0, planeDist*0.1, dz);
      pcol.w *= fadeOut*fadeIn;

      acol = alphaBlend(pcol, acol);
    } else {
      cutOut = true;
      acol.w = acol.w > cutOff ? 1.0 : acol.w;
      break;
    }

  }

  vec3 col = alphaBlend(skyCol, acol);
// To debug cutouts due to transparency  
//  col += cutOut ? vec3(1.0, -1.0, 0.0) : vec3(0.0);
  return col;
}

vec3 effect(vec2 p, vec2 q) {
  float z   = 0.33*planeDist*smoothTime;

  vec3 pro  = offset(z-1.0);
  vec3 ro   = offset(z);
  vec3 dro  = doffset(z);
  vec3 ddro = ddoffset(z);

  vec3 ww = normalize(dro);
  vec3 uu = normalize(cross(normalize(vec3(0.0,1.0,0.0)+ddro), ww));
  vec3 vv = cross(ww, uu);

  vec3 col = color(ww, uu, vv, pro, ro, p);
  
  return col;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  vec2 q = fragCoord/RENDERSIZE.xy;
  vec2 p = -1. + 2. * q;
  p.x *= RESOLUTION.x/RESOLUTION.y;
  
  vec3 col = effect(p, q);
  col *= smoothstep(0.0, 4.0, smoothTime);
  col = sRGB(col);
  
  fragColor = vec4(col, 1.0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}