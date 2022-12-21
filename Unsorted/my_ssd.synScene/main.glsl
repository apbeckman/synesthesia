

// CC0: Seven Segment Display
//  Wanted to create a seven segment display to show digits
//  I have seen versions of this on shadertoy before but gave it a go myself
//  as a fun little challenge.

#define TIME        TIME
#define RESOLUTION  RENDERSIZE

// License: Unknown, author: Unknown, found: don't remember
float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

// License: MIT OR CC-BY-NC-4.0, author: mercury, found: https://mercury.sexy/hg_sdf/
float mod1(inout float p, float size) {
  float halfsize = size*0.5;
  float c = floor((p + halfsize)/size);
  p = mod(p + halfsize, size) - halfsize;
  return c;
}

float dsegmentx(vec2 p, vec2 dim) {
  p.x = abs(p.x);
  float o = 0.5*max(dim.x-dim.y, 0.0);
  if (p.x < o) {
    return abs(p.y) - dim.y;
  }
  return length(p-vec2(o, 0.0))-dim.y;
}

vec3 digit(vec3 col, vec2 p, float aa, float n) {
  const int[16] digits = int[16](
    0x7D // 0
  , 0x50 // 1
  , 0x4F // 2
  , 0x57 // 3
  , 0x72 // 4
  , 0x37 // 5
  , 0x3F // 2
  , 0x51 // 7
  , 0x7F // 8
  , 0x77 // 9
  , 0x7B // A
  , 0x3E // B
  , 0x2D // C
  , 0x5E // D
  , 0x2F // E
  , 0x2B // F
  ); 
  const vec2 dim = vec2(0.75, 0.075);
  const float eps = 0.001;
  vec2 ap = abs(p);
  if (ap.x > (0.5+dim.y+eps)) return col;
  if (ap.y > (1.0+dim.y+eps)) return col;
  float m = mod(floor(n), 16.0);
  int digit = digits[int(m)];

  vec2 cp = (p-0.5);
  vec2 cn = round(cp);

  vec2 p0 = p;
  p0.y -= 0.5;
  p0.y = p0.y-0.5;
  float n0 = round(p0.y);
  p0.y -= n0;
  float d0 = dsegmentx(p0, dim);

  vec2 p1 = p;
  vec2 n1 = sign(p1); 
  p1 = abs(p1);
  p1 -= 0.5;
  p1 = p1.yx;
  float d1 = dsegmentx(p1, dim);
  
  vec2 p2 = p;
  p2.y = abs(p.y);
  p2.y -= 0.5;
  p2 = abs(p2);
  float d2 = dot(normalize(vec2(1.0, -1.0)), p2);

  float d = d0;
  d = min(d, d1);
  const vec3 acol = vec3(1.0, 0., 0.25);
  const vec3 icol = acol*0.1;

  float sx = 0.5*(n1.x+1.0) + (n1.y+1.0);
  float sy = -n0;
  float s  = d2 > 0.0 ? (3.0+sx) : sy;
  // Praying bit shift operations aren't TOO slow
  vec3 scol = ((digit & (1 << int(s))) == 0) ? icol : acol;  

  col = mix(col, scol, smoothstep(aa, -aa, d));
  return col;
}
vec3 digit(vec3 col, vec2 p, float n) {
  float aa = fwidth(p.y);
  return digit(col, p, aa, n);
}

vec3 effect(vec2 p, vec2 cq) {
  vec2 pp = p;
  pp *= 0.33*(-length(cq)+5.0);
  pp += 0.1*TIME;
  pp *= 10.0;

  float nx = mod1(pp.x, 1.5*3.0);
  float hx = hash(nx+123.4);

  pp.y += 8.0*(hx-0.25)*TIME;
  float ny = mod1(pp.y, 3.0);  
  float hy = hash(ny+456.7);
  float n = TIME*(hx+hy)+16.0*fract(hx+hy);

  vec3 col = vec3(0.025)*vec3(1.0, 0.0, 0.5);
  col = digit(col, pp, n+1234.5);
  col *= smoothstep(1.5, 0.5, length(cq));
  col = sqrt(col.yxz);
  return col;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  vec2 q = fragCoord/RESOLUTION.xy;
  vec2 p = -1. + 2. * q;
  vec2 pp = p;
  p.x *= RESOLUTION.x/RESOLUTION.y;

  vec3 col = effect(p, pp);  
  
  fragColor = vec4(col, 1.0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}