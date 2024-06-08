#include "lygia/sdf.glsl"
#include "lygia/space/lookAt.glsl"
#include "lygia/space/ratio.glsl"
#include "lygia/color/space/linear2gamma.glsl"
// #include "lygia/space/center.glsl"
float rep(float p, float d) {
	return mod(p - d * .5, d) - d * .5;
}



void mo(inout vec2 p, vec2 d) {
	p.x = abs(p.x) - d.x;
	p.y = abs(p.y) - d.y;
	if(p.y > p.x);
		p = p.yx;
}

void amod(inout vec2 p, float m) {
	float a = rep(atan(p.x, p.y), m);
	p = vec2(cos(a), sin(a)) * length(p);
}
#include "hg_sdf.glsl"
// vec3 rd = 
#define RESOLUTION          RENDERSIZE

// #include "lygia/lighting/raymarch.glsl"
vec2 pixel = vec2(1.0) / RENDERSIZE;
vec2 st = gl_FragCoord.xy * pixel;
vec2 uv = ratio(st, RENDERSIZE);
vec3 cam = vec3(4.5, 4.2, 4.5) * 30.;
vec3 ta = vec3(0.0);
vec3 up = vec3(0.0, 1.0, 0.0);
mat3 viewMatrix = lookAt(cam, ta, up);
vec3 rd = viewMatrix * normalize(vec3(uv, -1.0));
#define RAYMARCH_SAMPLES    100
#define RAYMARCH_MULTISAMPLE 3
#define RAYMARCH_BACKGROUND ( vec3(0.1, 0.1, .10) + ray.y * 0.8 )
// #define RAYMARCH_CAMERA_MATRIX_FNC(ro, ta) raymarchCamera(ro, ta)
// #define RAYMARCH_RENDER_FNC(ro, rd) raymarchDefaultRender(vec3(1.0), rd)  
#define RAYMARCH_AMBIENT    vec3(0.7, 0.9, 1.0)
#define RAYMARCH_MAX_DIST 50
#include "lygia/lighting/raymarch.glsl"
// mat3 ca = RAYMARCH_CAMERA_MATRIX_FNC(ro, ta);
vec3 ones = vec3(1);
vec3 blank = vec3(0.);
vec3 axisCorrection(vec3 p) {
	p.xz = _rotate(p.xz, .75 * PI);
	p.yz = _rotate(p.yz, -.34 * PI);
	float x = p.x;
	float y = p.y;
	float z = p.z;
	p.x = x;
	p.y = z;
	p.z = y;
	return p;
}
float smin(float c, float b, float k) {
	float h = max(k - abs(c - b), 0.) / k;
	return min(c, b) - h * h * h * k * 1. / 6.;
}
vec4 raymarchMap(in vec3 pos) {
	pos = axisCorrection(pos);
	vec3 modParams = vec3(4.0, 0.0, 7.0);
	pos.xz = _rotate(pos.xz, PI*0.25);
	pos.xy = _rotate(pos.xy, TIME*0.3);
	vec3 cPos = pos;
	pR(cPos.xy, cPos.z*4.);
	// pos.y = smin(pos.y, -pos.y, 0.5);
	// pos.y = smin(pos.y+0.5, -pos.y-0.5, 0.5);
	// pos.yz = _rotate(pos.yz, TIME*0.2);
	// pos.x = smin(pos.x+.5, -pos.x-0.5, 0.5);
	// cPos.yz = _rotate(cPos.yz, 0.5*PI);
	// pos.x = mod(pos.x, 10.);
	pR(pos.xy, (pos.z*0.52)+pos.z*.1+_noise(_uv + TIME));
	pos.z += fract(smoothTimeC*0.00125)*90;
	pos.yz = _rotate(pos.yz, 0.5 * PI);
	pos.y = mod(pos.y, 9.00);
	pos.z = pModPolar(pos.xz, 6);
	cPos.z += fract(-smoothTime*0.0125)*55.;
	cPos.z = mod(cPos.z, 5.5);
	// cpos.z; 
	pos.x+= 1.5;
	vec2 uno = abs(_rotate(pos.xy * .5, -0.25 * smoothTime));
	vec2 dos = abs(_rotate(pos.xy * .5, -0.25 * smoothTime - 2.0));
	pos.x = fOpUnionStairs(uno.x, dos.x, stair_radius, stair_count);
	pos.y = fOpUnionStairs(uno.y, dos.y, stair_radius, stair_count);
	vec2 tres = abs(_rotate(pos.xy * .5, -0.5 * smoothTime - 4.0));
	pos.x -= 1.50 * modscale;
	vec4 res = vec4(1.);
	vec2 hexSize = vec2(hex_width, hex_length);
	float c = hexPrismSDF(cPos, hexSize*vec2(0.2, 1.0));
	float hex = hexPrismSDF(pos, hexSize);
	hex = fOpUnionStairs(hex, c, .25, 5.);
	vec3 hex_offset = vec3(2.0, 0.0, 0.0);
	float hexR = hexPrismSDF(pos + hex_offset, hexSize);
	float hexL = hexPrismSDF(pos - hex_offset, hexSize);
	float hex_duo = min(hexL, hexR);
	res = opUnion(res, vec4(ones, hex));
	// res = opUnion(vec4(ones, c), res);
	return res;
}

vec3 test = vec3(0.);
vec4 renderMain(void) {


	vec4 fragColor = vec4(0);
	vec3 color = vec3(0);
	vec2 pixel = vec2(1.0) / RENDERSIZE;
	vec2 st = gl_FragCoord.xy * pixel;
	vec2 uv = ratio(st, RENDERSIZE);
	color = raymarch(cam, ta, uv).rgb;
	color = linear2gamma(color);

	fragColor = vec4(color, 1.0);
	return fragColor;
}