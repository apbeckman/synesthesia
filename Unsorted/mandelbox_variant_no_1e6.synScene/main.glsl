

// License CC0: Mandelbox variant no 1E6
//  Another variant of Mandelbox by Evilryu: https://www.shadertoy.com/view/XdlSD4

#define TOLERANCE       0.0001
#define MAX_RAY_LENGTH  8.0
#define MAX_RAY_MARCHES 100
#define TIME            TIME
#define RESOLUTION      RENDERSIZE
#define ROT(a)          mat2(cos(a), sin(a), -sin(a), cos(a))
#define PI              3.141592654
#define TAU             (2.0*PI)

const float fixed_radius2 = 1.9;
const float min_radius2   = 0.5;
const vec3  folding_limit = vec3(1.0);
const float scale         = -2.8;
const int   max_iter      = 120;
const vec3  bone          = vec3(0.89, 0.855, 0.788);

float tanh_approx(float x) {
//  return tanh(x);
  float x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}

float pmin(float a, float b, float k) {
  float h = clamp(0.5+0.5*(b-a)/k, 0.0, 1.0);
  
  return mix(b, a, h) - k*h*(1.0-h);
}

vec3 pmin(vec3 a, vec3 b, vec3 k) {
  vec3 h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0);
  return mix(b, a, h) - k*h*(1.0-h);
}

void sphere_fold(inout vec3 z, inout float dz) {
    float r2 = dot(z, z);
    if(r2 < min_radius2) {
        float temp = (fixed_radius2 / min_radius2);
        z *= temp;
        dz *= temp;
    } else if(r2 < fixed_radius2) {
        float temp = (fixed_radius2 / r2);
        z *= temp;
        dz *= temp;
    }
}

void box_fold(float k, inout vec3 z, inout float dz) {
  vec3 zz = sign(z)*pmin(abs(z), folding_limit, vec3(k));
  z = zz * 2.0 - z;
}

float sphere(vec3 p, float t) {
  return length(p)-t;
}

float boxf( vec3 p, vec3 b, float e)
{
  p = abs(p  )-b;
  vec3 q = abs(p+e)-e;
  return min(min(
      length(max(vec3(p.x,q.y,q.z),0.0))+min(max(p.x,max(q.y,q.z)),0.0),
      length(max(vec3(q.x,p.y,q.z),0.0))+min(max(q.x,max(p.y,q.z)),0.0)),
      length(max(vec3(q.x,q.y,p.z),0.0))+min(max(q.x,max(q.y,p.z)),0.0));
}

float mb(vec3 z) {
    vec3 offset = z;
    float dr = 1.0;
    float fd = 0.0;
    const float k = 0.05;
    for(int n = 0; n < 5; ++n) {
        box_fold(k/dr, z, dr);
        sphere_fold(z, dr);
        z = scale * z + offset;
        dr = dr * abs(scale) + 1.0;        
        float r1 = sphere(z, 5.0);
        float r2 = boxf(z, vec3(5.0), 0.5);
        float r = n < 4 ? r2 : r1;
        float dd = r / abs(dr);
        if (n < 3 || dd < fd) {
          fd = dd;
        }
    }
    return fd;
}

float df(vec3 p) { 
    const float z = 0.3;
    p.y -= 1.0;
    float d1 = mb(p/z)*z;
    return d1; 
} 

float rayMarch(in vec3 ro, in vec3 rd, out int iter) {
  float t = 0.0;
  int i = 0;
  for (i = 0; i < MAX_RAY_MARCHES; i++) {
    float distance = df(ro + rd*t);
    if (distance < TOLERANCE || t > MAX_RAY_LENGTH) break;
    t += distance;
  }
  iter = i;
  return t;
}

vec3 normal(in vec3 pos) {
  vec3  eps = vec3(.0005,0.0,0.0);
  vec3 nor;
  nor.x = df(pos+eps.xyy) - df(pos-eps.xyy);
  nor.y = df(pos+eps.yxy) - df(pos-eps.yxy);
  nor.z = df(pos+eps.yyx) - df(pos-eps.yyx);
  return normalize(nor);
}

float softShadow(in vec3 pos, in vec3 ld, in float ll, float mint, float k) {
  const float minShadow = 0.25;
  float res = 1.0;
  float t = mint;
  for (int i=0; i<24; i++) {
    float distance = df(pos + ld*t);
    res = min(res, k*distance/t);
    if (ll <= t) break;
    if(res <= minShadow) break;
    t += max(mint*0.2, distance);
  }
  return clamp(res,minShadow,1.0);
}

vec3 postProcess(in vec3 col, in vec2 q)  {
  col=pow(clamp(col,0.0,1.0),vec3(1.0/2.2)); 
  col=col*0.6+0.4*col*col*(3.0-2.0*col);  // contrast
  col=mix(col, vec3(dot(col, vec3(0.33))), -0.4);  // satuation
  col*=0.5+0.5*pow(19.0*q.x*q.y*(1.0-q.x)*(1.0-q.y),0.7);  // vigneting
  return col;
}

vec3 hsv2rgb(vec3 c) {
  const vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec3 render(in vec3 ro, in vec3 rd) {
  vec3 lightPos = 2.0*vec3(1.5, 3.0, 1.0);

  vec3 skyCol = vec3(0.0);

  int iter = 0;
  float t = rayMarch(ro, rd, iter);

  float ifade = 1.0-tanh_approx(3.0*float(iter)/float(MAX_RAY_MARCHES));

  vec3 pos = ro + t*rd;    
  vec3 nor = vec3(0.0, 1.0, 0.0);
  
  vec3 color = vec3(0.0);
  
  float h = fract(TIME/30.0);
  
  if (t < MAX_RAY_LENGTH && pos.y > 0.0) {
    // Ray intersected object
    nor       = normal(pos);
    vec3 hsv  = (vec3(fract(h - 0.6 + 0.4+1.75*t), 1.0-ifade, 1.0));
    color     = hsv2rgb(hsv);
  } else if (pos.y > 0.0) {
    // Ray intersected sky
    return skyCol*ifade;
  } else {
    // Ray intersected plane
    t   = -ro.y/rd.y;
    pos = ro + t*rd;
    nor = vec3(0.0, 1.0, 0.0);
    vec2 pp = pos.xz*1.5;
    float m = 0.5+0.25*(sin(3.0*pp.x+TIME*2.1)+sin(3.3*pp.y+TIME*2.0));
    m *= m;
    m *= m;
    pp = fract(pp+0.5)-0.5;
    float dp = pmin(abs(pp.x), abs(pp.y), 0.025);
    vec3 hsv = vec3(0.4+mix(0.15,0.0, m), tanh_approx(mix(100.0, 10.0, m)*dp), 1.0);
    color = 5.5*hsv2rgb(hsv)*exp(-mix(30.0, 10.0, m)*dp);
  }

  vec3 lv   = lightPos - pos;
  float ll2 = dot(lv, lv);
  float ll  = sqrt(ll2);
  vec3 ld   = lv / ll;
  float sha = softShadow(pos, ld, ll, 0.01, 64.0);

  float dm  = min(1.0, 40.0/ll2);
  float dif = max(dot(nor,ld),0.0)*dm;
  float spe = pow(max(dot(reflect(-ld, nor), -rd), 0.), 10.);
  float l   = dif*sha;

  float lin = mix(0.2, 1.0, l);

  vec3 col = lin*color + spe*sha;

  float f = exp(-20.0*(max(t-3.0, 0.0) / MAX_RAY_LENGTH));
    
  return mix(skyCol, col , f)*ifade;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  vec2 q=fragCoord.xy/RESOLUTION.xy; 
  vec2 p = -1.0 + 2.0*q;
  p.x *= RESOLUTION.x/RESOLUTION.y;

  vec3 ro = 0.4*vec3(2.0, 0, 0.2)+vec3(0.0, 1.25, 0.0);
  float tm = mix(10.33, 11.7, 0.5+0.5*sin(TIME*TAU/60.0));
  ro.xz *= ROT(sin(tm));
  ro.yz *= ROT(sin(tm*sqrt(0.5))*0.25);

  vec3 ww = normalize(vec3(0.0, 1.0, 0.0) - ro);
  vec3 uu = normalize(cross( vec3(0.0,1.0,0.0), ww));
  vec3 vv = normalize(cross(ww,uu));
  vec3 rd = normalize( p.x*uu + p.y*vv + (2.0+0.5*tanh_approx(length(p)))*ww);

  vec3 col = render(ro, rd);
  col = clamp(col, 0.0, 1.0);
  fragColor = vec4(postProcess(col, q),1.0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}