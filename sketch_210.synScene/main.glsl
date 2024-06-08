#include "lygia/sdf.glsl"
#include "lygia/space/lookAt.glsl"
#include "lygia/space/ratio.glsl"
#include "lygia/color/space/linear2gamma.glsl"
// #include "lygia/space/center.glsl"

#include "hg_sdf.glsl"
// vec3 rd = 
#define RESOLUTION          RENDERSIZE

// #include "lygia/lighting/raymarch.glsl"
vec2 pixel = vec2(1.0) / RENDERSIZE;
vec2 st = gl_FragCoord.xy * pixel;
vec2 uv = ratio(st, RENDERSIZE);
vec3 cam = vec3(4.5, 4.2, 4.5) * 20.;
vec3 ta = vec3(0.0);
vec3 up = vec3(0.0, 1.0, 0.0);
mat3 viewMatrix = lookAt(cam, ta, up);
vec3 rd = viewMatrix * normalize(vec3(uv, -1.0));
#define RAYMARCH_SAMPLES    100
#define RAYMARCH_MULTISAMPLE 4
#define RAYMARCH_BACKGROUND ( vec3(0.1, 0.1, .10) + ray.y * 0.8 )
#define RAYMARCH_CAMERA_MATRIX_FNC(ro, ta) raymarchCamera(ro, ta)
#define RAYMARCH_RENDER_FNC(ro, rd) raymarchDefaultRender(vec3(1.0), rd)  
#define RAYMARCH_AMBIENT    vec3(0.7, 0.9, 1.0)
#define RAYMARCH_MAX_DIST 55
#include "lygia/lighting/raymarch.glsl"
// mat3 ca = RAYMARCH_CAMERA_MATRIX_FNC(ro, ta);
vec3 ones = vec3(1);
vec3 blank = vec3(0.);
vec3 axisCorrection(vec3 p){
	p.xz = _rotate(p.xz, .75*PI);
	p.yz = _rotate(p.yz, -.34*PI);
	float x = p.x;
	float y = p.y;
	float z = p.z;
	p.x = x;
	p.y = z;
	p.z = y;
	return p;
}
vec3 camera(vec3 ro, vec2 uv, vec3 ta) {
	uv = mix(uv, uv + _uvc * PI, fov);
	vec3 fwd = normalize(ta - ro);
	vec3 left = cross(vec3(0, 1, 0), fwd);
	vec3 up = cross(fwd, left);
	return normalize(fwd + uv.x * left + up * uv.y);
}
vec4 raymarchMap(in vec3 pos){
	pos = axisCorrection(pos);
	pos.z = mod(pos.z, 4.5);
	pos.yz = _rotate(pos.yz, 0.5*PI);
	vec4 res = vec4(1.);
	vec2 hexSize = vec2(hex_width, 3.5);
	vec3 hex_offset = vec3(4.0, 0.0, 0.0);
	float hex = hexPrismSDF(pos + hex_offset, hexSize);
	hex = min(hex, hexPrismSDF(pos - hex_offset, hexSize));
	res = opUnion(res, vec4(ones, hex));
	return res; 

}
vec3 test = vec3(0.);
vec4 renderMain(void)
{
	float dt = (bass_time) * 3.;
	vec3 ro = vec3(0, 0,  dt);
	vec3 ta = vec3(0, 0, dt);
	vec3 rd;
	vec4 fragColor = vec4(0);
	vec2 uv = _xy;
	vec3 color = vec3(0);
	vec2 pixel = vec2(1.0) / RENDERSIZE;
	vec2 st = gl_FragCoord.xy * pixel;
	// cam.z -= TIME;
	rd = camera(cam, uv, ta);
	uv = ratio(st, RENDERSIZE);

	color = raymarch(rd,ta, uv).rgb;
	color = linear2gamma(color);

	fragColor = vec4(color, 1.0);
	return fragColor;
}