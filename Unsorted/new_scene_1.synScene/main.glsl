
#include "hg_sdf.glsl"
#include "lygia/sdf.glsl"
#include "sqTile.glsl"
#include "../space/mirrorTile.glsl"
#include "../space/depth2viewZ.glsl"
// #include "../space/screen2viewPosition.glsl"
#include "../sample.glsl"
#include "../math/lengthSq.glsl"
#include "lygia/draw/fill.glsl"
#include "lygia/space/lookAt.glsl"
#include "lygia/space/ratio.glsl"
#include "lygia/space/xyz2equirect.glsl"
#include "lygia/math/inverse.glsl"
#include "lygia/sdf/opExtrude.glsl"
#include "lygia/math/const.glsl"
#include "lygia/math/mod2.glsl"
#define CAMERA_POSITION     cam

// #include "lygia/space/center.glsl"
// #include "lygia/space/screen2viewPosition.glsl"
#define RESOLUTION          RENDERSIZE
#define RAYMARCH_SAMPLES    90
#define RAYMARCH_MULTISAMPLE 3
#define RAYMARCH_BACKGROUND ( vec3(0.1, 0.1, .10) + ray.y * 0.8 )
// #define RAYMARCH_CAMERA_MATRIX_FNC(cam, ta)
// #define RAYMARCH_RENDER_FNC(ro, rd)
#define RAYMARCH_AMBIENT    vec3(0.7, 0.9, 1.0)
#define RAYMARCH_MAX_DIST 40
#define M 60


// #define SAMPLE2DCUBE_CELL_SIZE 64.0
// #define SAMPLE2DCUBE_CELLS_PER_SIDE 8.0
#define SAMPLE2DCUBE_FNC(TEX, UV) sampleBicubic(TEX, UV, vec2(512.0))
// #define SAMPLE2DCUBE_FNC(TEX, UV) sampleSmooth(TEX, UV, vec2(512.0))
// #define SAMPLE2DCUBE_FNC(TEX, UV) sampleNearest(TEX, UV, vec2(512.0))

#include "lygia/math/cubic.glsl"
#include "lygia/math/quartic.glsl"
#include "lygia/math/quintic.glsl"
// #include "lygia/sdf.glsl"
#include "lygia/space/ratio.glsl"
#include "lygia/space/scale.glsl"
#include "lygia/sample/bicubic.glsl"
#include "lygia/sample/smooth.glsl"
#include "lygia/sample/nearest.glsl"
#include "lygia/sample/3DSdf.glsl"
#include "lygia/math/saturate.glsl"
#include "lygia/lighting/raymarch.glsl"
#include "lygia/color/space/linear2gamma.glsl"

// vec3 cam = vec3(vec2(1., -1.), -100);


float smin(float c, float b, float k) {
	float h = max(k - abs(c - b), 0.) / k;
	return min(c, b) - h * h * h * k * 1. / 6.;
}
vec2 smin2(vec2 a, vec2 b, float k){
	vec2 mixed;
	mixed.x = smin(a.x, b.x, k);
	mixed.y = smin(a.y, b.y, k);
	return mixed;
}
vec3 smin3(vec3 a, vec3 b, float k){
	vec3 mixed;
	mixed.x = smin(a.x, b.x, k);
	mixed.y = smin(a.y, b.y, k);
	mixed.z = smin(a.z, b.z, k);
	return mixed;
}
mat3 wiggleMatrix(vec3 axis, float angle) {
//  from 'hex array' 
    float s = sin(angle);
    float c = cos(angle);
    float oc = 2.0 - c;
    
    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}
//3D Rotation
vec3 rot3d(in vec3 x, in vec3 r) {
	x.xy = _rotate(x.xy, r.z);
	x.yz = _rotate(x.yz, r.x);
	x.zx = _rotate(x.zx, r.y);
	return x;
}
#include "lygia/lighting/raymarch.glsl"

float checkBoard(vec2 uv, vec2 _scale) {
	uv.xy += fract(TIME*0.001);
	for(int i = 0; i <= 2; i++) {
		uv += _noise(abs(_uv) + smoothTimeC * 0.01);
	}
	uv = floor(fract(uv * _scale) * 12.0);
	return min(1.0, uv.x + uv.y) - (uv.x * uv.y);
}

vec4 raymarchMap(in vec3 pos) {
	vec3 ones, box_pos, hex_pos, hex1_pos, hex2_pos, rotate, gem_pos, repeat;
	vec2 hex_size, st;
	st = (_xy * 1/RENDERSIZE);
	st = ratio(st, RENDERSIZE);
	vec4 t = mirrorTile(st*4.0);
	hex_size = vec2(0.5, 5.0);

	ones = vec3(1.0);
	box_pos = pos;
	box_pos += 4.;
	box_pos.y += 1.;
	gem_pos = box_pos;
	vec3 box2_pos = box_pos;
	vec3 box3_pos = box_pos;
	box_pos = smin3(box_pos, -box_pos, 0.125);
	// box_pos.x = smin(box_pos.x, -box_pos.x, 0.25);
	box_pos.xz = _rotate(box_pos.xz, smoothTimeC*0.05);
	// box_pos.y = smin(box_pos.y, -box_pos.y, 0.25);
	box_pos.yz = _rotate(box_pos.yz, smoothTime*0.05);
	// box_pos.z = smin(box_pos.z, -box_pos.z, 0.25);
	// box_pos.x = smin(box_pos.x, -box_pos.x, 0.25);
	box_pos = smin3(box_pos*0.9, -box_pos, 0.125);
	box_pos.xy = _rotate(box_pos.xy, smoothTimeC*0.05);
	box_pos = smin3(box_pos, -box_pos, 0.125);
	gem_pos.xz = _rotate(gem_pos.xz, smoothTimeC*0.035);
	gem_pos = abs(gem_pos);
	gem_pos.yz = _rotate(gem_pos.yz, smoothTime*0.035);
	gem_pos = abs(gem_pos);
	gem_pos.xy = _rotate(gem_pos.xy, smoothTimeC*0.035);
	gem_pos = abs(gem_pos);
	box2_pos.xz = _rotate(box2_pos.xz, smoothTime*0.05);
	box2_pos.yz = _rotate(box2_pos.yz, smoothTimeC	*0.05);
	box2_pos.xy = _rotate(box2_pos.xy, smoothTime*0.05);
	box3_pos.yz = _rotate(box3_pos.yz, smoothTime*0.03);
	box3_pos.xy = _rotate(box3_pos.xy, smoothTimeC*0.03);
	box3_pos.xz = _rotate(box3_pos.xz, smoothTimeC*0.03);
	hex_pos = pos;
	rotate = vec3(xy1.x, xy1.y, 0);
	hex_pos.xz = _rotate(hex_pos.xz, .75*PI);
	hex_pos.yz = _rotate(hex_pos.yz, -.37*PI);
	pos.xz = _rotate(pos.xz, .5*PI);
	pos.yz = _rotate(pos.yz, .5*PI);
	hex_pos.z -= 0.75;
	hex_pos = rot3d(hex_pos, vec3(.25, 0., 0.));
	pR(hex_pos.xz, hex_pos.y*0.0125);
	hex_pos.y -= fract(smoothTimeC*0.02)*60;
	hex_pos.y = mod(hex_pos.y, 6.0);
	// hex_pos.z = wiggleMatrix(vec3(0.0, 0.0, 2.0), 1.0 * sin(hex_pos.z*2.5));
	// / hex_pos += vec3(10,0,0);
	hex_pos = opRepeat(hex_pos, repeat);
	// hex_pos.xy = _rotate(hex_pos.xy, sin(TIME*0.2)*0.2002);
	hex_pos.z += 1.25;
	hex_pos.x = fOpUnionChamfer(hex_pos.x, -hex_pos.x, 0.2);
	hex_pos.xz = _rotate(hex_pos.xz, 0.5);
	hex_pos.z = fOpUnionChamfer(hex_pos.z, -hex_pos.z, 0.2);
	hex_pos.xz = _rotate(hex_pos.xz, TIME*0.1);
	pR(hex_pos.xz, hex_pos.y*-0.025);
	hex_pos.x = fOpUnionChamfer(hex_pos.x, -hex_pos.x, 0.2);
	hex_pos.xz = _rotate(hex_pos.xz, -TIME*0.2);
	hex_pos.z = fOpUnionChamfer(hex_pos.z, -hex_pos.z, 0.2);
	hex_pos.xz = _rotate(hex_pos.xz, -0.1);
	// hex_pos.xz = _rotate(hex_pos.xz, (TIME*0.2)*0.2002);
	hex_pos.z = pModPolar(hex_pos.xz,4);
	hex_pos -= vec3(1.0);
	// hex_pos.z = pModPolar(hex_pos.xz, 10.)-5.;
	hex_pos = opRepeat(hex_pos, vec3(3.5, 0., 3.5));
	// hex_pos.x -= vec3(-1., 0., 0.);
	hex_pos.xy = abs(hex_pos.xy);
	// hex_pos = rot3d(hex_pos, vec3(0., 0. ,0.+sin(TIME)*0.5 ));
	// hex2_pos = hex_pos;
	hex_pos.x -= 1.90;
	hex2_pos.x += 02.0;
	// hex_pos = opRepeat(hex_pos, repeat);
	
	// hex2_pos.y = mod(hex2_pos.y + TIME, 2.5);
	// hex1_pos.xy = _rotate(hex1_pos.xy, .5);
	// hex2_pos.xy = _rotate(hex2_pos.xy, -.5);
	// box_pos = rot3d(box_pos, rotate);
	// box_pos.x = mod(box_pos.x, 1.);
	// hex_pos.x = abs(1.0-hex_pos.x);
	// hex_pos = rot3d(hex_pos, hex_rotate);
	vec3 bfh = vec3(0.75, 0.75, 0.75);
	vec3 bf2h = vec3(0.5, 0.5, 0.5);
	vec4 res = vec4(1.);
	float hexP = hexPrismSDF(hex_pos, hex_size);
	float hexP2 = hexPrismSDF(hex2_pos, hex_size);
	// hexP = opRepeat(hex1_pos, vec3(1.0, 0.0, 0.0));
	float dodecahedron = dodecahedronSDF(gem_pos, 0.5);
	float boxframe = boxFrameSDF(box_pos, bfh*1.2, .1);
	
	
	float boxframe2 = boxFrameSDF(box2_pos, bf2h, .05);
	
	// float boxframe3 = boxFrameSDF(box3_pos*0.45, bf2h*01.2, .05);
	// float check = checkBoard(vec2(_noise(pos.xz)), vec2(5.0));
	// res.z = wiggleMatrix(vec3(0.0, 1.0, 0.0), 1.0*sin(res.z*0.5)); 
	boxframe = fOpUnionStairs(boxframe, boxframe2, .1, 4.);
	// boxframe = fOpUnionStairs(boxframe, boxframe3, .4, 4.);
	boxframe = fOpUnionColumns(boxframe, dodecahedron, .12, 4.);
	float final = fOpUnionRound(boxframe, hexP, 0.1);
	// barrel distortion 
	// final = barrel(vec2(0.5, 0.5) , final);
	res = opUnion(res, vec4(ones, final));
	// res = opUnion(res, vec4(vec3(ones),hexP)); 
	// res = opUnion(res, vec4(vec3(ones),hexP2));
	// res = opUnion(res, vec4(ones, boxframe2));
	// res = opUnion(res, vec4(ones, boxframe3));
	// res = opUnion(res, vec4(ones, dodecahedron));
	// res = opUnion(res, vec4(1.0, 1.0, 1.0, sample3DSdf(u_tex0, c_pos)));
	// res = opUnion(res, planeSDF,)
	// res = opUnion(res, vec4(vec3(0.5 + check * 0.5), planeSDF(pos + 1.0)), 0.1);

	return res;
}

vec4 renderMain() {
	vec3 color;
	vec3 cam = vec3(4.5, 4.2, 4.5) * 10.;
	vec3 eye = cam;
	vec3 up = vec3(0.0, 1.0, 0.0);
	vec3 target = vec3(0.0, 0.0, 0.0);
	vec3 forward = vec3(0.0, 0.0, -1.0);
	// vec3 forward = normalize(vec3(eye - target));
	mat3 rmc = raymarchCamera(vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 0.0), up);
	// mat3 look = lookAt(forward, up);
	vec4 fragColor = vec4(0.);
	vec2 fragCoord = _xy;
	float d2z = depth2viewZ(01., 0., 1.+TIME);
	// cam.z = d2z;
	// cam = (look[1]);
	color = vec3(0.0);
	vec2 pixel = vec2(1.0) / RENDERSIZE;
	vec2 st = gl_FragCoord.xy * pixel;
	vec2 uv = ratio(st, RENDERSIZE);
	vec3 s2v = cart2polar( vec3(_uvc, 1.0));
	vec3 rd = vec3(2.*fragCoord - RENDERSIZE.xy, RENDERSIZE.y);
    // rd = normalize(vec3(rd.xy, sqrt(max(rd.z*rd.z - dot(rd.xy, rd.xy)*.2, 0.))));
	vec2 mo = _mouse.xy * pixel;
	float time = 32.0 + TIME * 1.5;
	
	// cam *= rmc;
	// vec3 cam = vec3(rd);
	
	// cam.z + TIME;
	// cam = normalize(cam);
	color = raymarch(cam, uv).rgb;
	color = linear2gamma(color);

	fragColor = vec4(color, 1.0);
	return fragColor;
}
    // vec3 n = normal(p);
    // vec3 l = normalize(lightPos - p);
    // float diff = clamp(dot(n, l), 0.0, 1.0);
    // vec4 fragColor = vec4(vec3(diff), 1.0);
