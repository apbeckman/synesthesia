vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


////////////////////////////////////////////////////////////////////////////////
//
// Rainbow Temple, Ungolfed. Matthew Arcus, mla, 2020.
//
// Ungolfed and rationalized version of GregRostami's Rainbow Temple:
// https://www.shadertoy.com/view/3sjcDV, which is derived from vahidk's
// Infinite Tunnel: https://www.shadertoy.com/view/tdSyDG which is derived
// from yoshin4004's original: https://twitter.com/yosshin4004/status/1251357672504360966
//
// <mouse>: look around
// f: floor
// l: change light direction
// v: vary sphere overlap
// z: move along z axis
//
// Logic is mostly the same, have dinked with the lighting a little.
//
////////////////////////////////////////////////////////////////////////////////

// Rotate vector p by angle t.
vec2 rotate(vec2 p, float t) {
  return cos(t)*p + sin(t)*vec2(-p.y,p.x);
}

//float PI = 3.14159;

vec3 transform(vec3 p) {
  if (iMouse.x > 0.0) {
    // Full range of rotation across the screen.
    float phi = (2.0*iMouse.x-RENDERSIZE.x)/RENDERSIZE.x*PI;
    float theta = (2.0*iMouse.y-RENDERSIZE.y)/RENDERSIZE.y*PI;
    p.yz = rotate(p.yz,-theta);
    p.zx = rotate(p.zx,phi);
  }
  return p;
}

const int CHAR_F = 70;
const int CHAR_L = 76;
const int CHAR_V = 86;
const int CHAR_Z = 90;
/*
bool key(int key) {
  return texelFetch(iChannel3,ivec2(key,2),0).x != 0.0;
}
*/
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  // Orthographic projection - this is the ray direction
  vec3 d = vec3((2.0*fragCoord-RENDERSIZE.xy)/RENDERSIZE.y, 1);
  d = transform(d);
      d.yz = _rotate(d.yz, lookXY.y*PI);
    d.xy = _rotate(d.xy, lookXY.x*PI);

  vec3 p = vec3(0); // Ray marching origin
    p.z = 0.2*TIME;


  vec3 q,p0 = p;
  vec3 lightdir = vec3(1,1,1);
  lightdir.xz = rotate(lightdir.xz,0.123*script_time);
  for (int i = 0; i < 120; i++) {
    float lfactor = 0.5;
    float overlap = 0.65;
    overlap += 0.1*sin(0.1*smooth_basstime);
    float floorheight = 0.45;
    float s = overlap-length(fract(p+0.5)-0.5);
    s = min(s, p.y + floorheight);
    p += lfactor*d*s;
    if (i == 80) {
      // We must have hit the surface by now.
      q = p; // Remember position
      p -= d*.01; // Back up a little to stop self-shadowing
      d = lightdir; // Light direction (for shadows)
    }
  }
  ivec3 u = ivec3(q*5e2); // u in [0..255]
  u ^= u.yzx;             // Mix up for chequer effect
  u ^= u.zxy;
  u &= 255;
  vec3 col = vec3(u)/255.0;
  col *= 0.25*(length(p-q) + 0.9);
  col += 0.041*(length(p-p0));
  fragColor = vec4(col,1);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}