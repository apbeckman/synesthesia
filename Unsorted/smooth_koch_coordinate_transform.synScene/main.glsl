

// License CC0: Smooth koch coordinate transform
//  Playing around with koch mappings and attempt to smooth the distance
//  function between cells

//#define TIME        smoothTime*0.35
#define TTIME       (smoothTimeC*TAU)/4.
#define RESOLUTION  RENDERSIZE
#define PI          3.141592654
#define TAU         (2.0*PI)
#define N(a)        vec2(sin(a), cos(a))
#define LESS(a,b,c) mix(a,b,step(0.,c))
#define SABS(x,k)   LESS((.5/(k))*(x)*(x)+(k)*.5,abs(x),abs(x)-(k))
#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))
#define L2(x)       dot(x, x)
#define PSIN(x)     (0.5 + 0.5*sin(x))

float tanh_approx(float x) {
//  return tanh(x);
  float x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0*(1.0+basshits*2.), 1.0*(1.0+basshits*2.));
}
vec2 mod2_1(inout vec2 p) {
  vec2 c = floor(p + 0.5);
  p = fract(p + 0.5) - 0.5;
  return c;
}

float hash(vec2 co) {
  return fract(sin(dot(co, vec2(12.9898, 58.233))) * 13758.5453);
}

float pmin(float a, float b, float k) {
  float h = clamp(0.5+0.5*(b-a)/k, 0.0, 1.0);
  
  return mix(b, a, h) - k*h*(1.0-h);
}

float pmax(float a, float b, float k) {
  return -pmin(-a, -b, k);
}

float pabs(float a, float k) {
  return pmax(a, -a, k);
}

float truchet_lw = 0.06125/(1.0+pow(syn_BassLevel*4.+syn_MidLevel*2., 2.));

const mat2[] truchet_rots = mat2[](ROT(0.0*PI/2.0), ROT(1.00*PI/2.0), ROT(2.0*PI/2.0), ROT(3.0*PI/2.0));

float circle(vec2 p, float r) {
  return length(p) - r;
}

float modMirror1(inout float p, float size) {
  float halfsize = size*0.5;
  float c = floor((p + halfsize)/size);
  p = mod(p + halfsize,size) - halfsize;
  p *= mod(c, 2.0)*2.0 - 1.0;
  return c;
}

vec2 toPolar(vec2 p) {
  return vec2(length(p), atan(p.y, p.x));
}

vec2 toRect(vec2 p) {
  return vec2(p.x*cos(p.y), p.x*sin(p.y));
}

float smoothKaleidoscope(inout vec2 p, float sm, float rep) {
  vec2 hp = p;

  vec2 hpp = toPolar(hp);
  float rn = modMirror1(hpp.y, TAU/rep);

  float sa = PI/rep - SABS(PI/rep - abs(hpp.y), sm);
  hpp.y = sign(hpp.y)*(sa);

  hp = toRect(hpp);

  p = hp;

  return rn;
}

// Attempt for classic koch mapping, thanks to art of code YT
void koch(inout vec2 p) {
  const float a = PI*2.0/3.0;
  const vec2 nn = N(PI*5.0/6.0);
  const vec2 n = N(a);

  p.y -= sqrt(0.75);
  p.x = abs(p.x);
  p.x += -1.5;
  p   -= nn*max(0.0, dot(p, nn))*2.0;
  p.x -= -1.5;

  const int levels = 4;
  
  p.x  = abs(p.x);
  p.x -= 0.5;
  p   -= n*min(0.0, dot(p, n))*2.0;


  float s = 1.0;
  for (int i = 1; i < levels; ++i) {
    s /= 3.0;
    p *= 3.0;
    p.x -= 1.5;
    
    p.x  = abs(p.x);
    p.x -= 0.5;
    p   -= n*min(0.0, dot(p, n))*2.0;
  }

  p *= s;  
}

// Attempt to smooth koch mapping
void skoch(inout vec2 p) {
  const float a       = PI*2.0/3.0;
  const vec2  nn      = N(PI*5.0/6.0);
  const vec2  n       = N(a);
  const int   levels  = 4;
  const float k       = 0.0125;

  p.xy = p.yx;
  smoothKaleidoscope(p, k, 6.0);
  p.xy = p.yx;
  p.y -= 0.875+sin(smoothTimeC*0.001);

  p.x  = pabs(p.x, k*2.0);
  p.x -= 0.5;
  p   -= n*pmin(0.0, dot(p, n), k)*2.0;

  float s = 1.0;
  for (int i = 1; i < levels; ++i) {
    s /= 3.0;
    p *= 3.0;
    p.x -= 1.5;
    
    p.x  = pabs(p.x, k/s*2.0);
    p.x -= 0.5+sin(smoothTimeC*0.001);
    p   -= n*pmin(0.0, dot(p, n), k/s)*2.0+sin(smoothTimeC*0.001);
  }

  p *= s;  
}


vec2 truchet_cell0(vec2 p, float h) {
  float d0  = circle(p-vec2(0.5), 0.5);
  float d1  = circle(p+vec2(0.5), 0.5);

  float d = 1E6;
  d = min(d, d0);
  d = min(d, d1);
  return vec2(d, 1E6); // 1E6 gives a nice looking bug, 1E4 produces a more "correct" result
}

vec2 truchet_cell1(vec2 p, float h) {
  float d0  = abs(p.x);
  float d1  = abs(p.y);
  float d2  = circle(p, mix(0.2, 0.4, h));

  float d = 1E6;
  d = min(d, d0);
  d = min(d, d1);
  d = min(d, d2);
  return vec2(d, d2+truchet_lw);
}

float truchet(vec2 p) {
  vec2 np = mod2_1(p);
  float r = hash(np);

  p *= truchet_rots[int(r*4.0)];
  float rr = fract(r*31.0);
  vec2 cd0 = truchet_cell0(p, rr);
  vec2 cd1 = truchet_cell1(p, rr);
  vec2 d0 = mix(cd0, cd1, vec2(fract(r*13.0) > 0.5));

  float d = 1E6;
  d = min(d, d0.x);
  d = abs(d) - truchet_lw;

  return d;
}

float snowFlake(vec2 p) {
  const float s = 0.2;

  vec2 kp = p;
  skoch(kp);
  kp -= 0.7345;
  kp *= ROT(1.0);
  koch(kp);
  kp += smoothTimeC*0.025;

  float d = truchet(kp/s)*s;
  return d-0.0025;
}

float df(vec2 p) {
  const float rep = 10.0;
  const float sm  = 0.05*6.0/rep;
  mat2 rot = ROT(TTIME/240.0);
  p *= rot;
  smoothKaleidoscope(p, sm, rep);
  p *= rot;
  const float ss = 0.55;
  return snowFlake(p/ss)*ss;
}

float height(vec2 p) {
  float d = df(p);
  float h = tanh_approx(smoothstep(0.02, 0.0, d)*mix(2.0, 6.0, PSIN(TTIME/60.0)));
  return -h*0.01;
}

vec3 normal(vec2 p) {
  vec2 v;
  vec2 w;
  vec2 e = vec2(4.0/RESOLUTION.y, 0);
  
  vec3 n;
  n.x = height(p + e.xy) - height(p - e.xy);
  n.y = 2.0*e.x;
  n.z = height(p + e.yx) - height(p - e.yx);
  
  return normalize(n);
}

vec3 postProcess(vec3 col, vec2 q) {
  col = clamp(col, 0.0, 1.0);
  col = pow(col, 1.0/vec3(2.2));
  col = col*0.6+0.4*col*col*(3.0-2.0*col);
  col = mix(col, vec3(dot(col, vec3(0.33))), -0.4);
  col *=0.5+0.5*pow(19.0*q.x*q.y*(1.0-q.x)*(1.0-q.y),0.7);
  return col;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  const float s = 1.0;
  const vec3 lp1 = vec3(1.0, 1.25, 1.0)*vec3(s, 1.0, s);
  const vec3 lp2 = vec3(-1.0, 1.25, 1.0)*vec3(s, 1.0, s);

  vec2 q = fragCoord/RESOLUTION.xy;
  vec2 p = -1. + 2. * q;
  p.x *= RESOLUTION.x/RESOLUTION.y;

  float aa = 2.0/RESOLUTION.y;

  vec3 col = vec3(0.0);
  float d = df(p);
  float h = height(p);
  vec3  n = normal(p);

  vec3 ro = vec3(0.0, -10.0, 0.0);
  vec3 pp = vec3(p.x, 0.0, p.y);

  vec3 po = vec3(p.x, h, p.y);
  vec3 rd = normalize(ro - po);

  vec3 ld1 = normalize(lp1 - po);
  vec3 ld2 = normalize(lp2 - po);
  
  float diff1 = max(dot(n, ld1), 0.0);
  float diff2 = max(dot(n, ld2), 0.0);

  vec3  rn    = n;
  vec3  ref   = reflect(rd, rn);
  float ref1  = max(dot(ref, ld1), 0.0);
  float ref2  = max(dot(ref, ld2), 0.0);

  vec3 lcol1 = vec3(1.5, 1.5, 2.0).zyx;
  vec3 lcol2 = vec3(2.0, 1.5, 0.75).zyx;
  vec3 lpow1 = 0.15*lcol1/L2(ld1);
  vec3 lpow2 = 0.5*lcol2/L2(ld2);
  vec3 dm = vec3(1.0)*tanh(-h*10.0+0.125);
  col += dm*diff1*diff1*lpow1;
  col += dm*diff2*diff2*lpow2;
  vec3 rm = vec3(1.0)*mix(0.25, 1.0, tanh_approx(-h*1000.0));
  col += rm*pow(ref1, 10.0)*lcol1;
  col += rm*pow(ref2, 10.0)*lcol2;
  col = postProcess(col, q);  
  
  fragColor = vec4(col, 1.0);
	return fragColor; 
 } 




vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}