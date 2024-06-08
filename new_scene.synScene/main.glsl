#include "lygia/sdf.glsl"
#include "lygia/space/lookAt.glsl"
#include "lygia/space/ratio.glsl"
#include "lygia/color/space/linear2gamma.glsl"
// #include "lygia/space/center.glsl"
#include "lygia/sample.glsl"
// #include "lygia/map.glsl"

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
	// uv = mix(uv, uv + _uvc * PI, fov);
	vec3 fwd = normalize(ta - ro);
	vec3 left = cross(vec3(0, 1, 0), fwd);
	vec3 up = cross(fwd, left);
	return normalize(fwd + uv.x * left + up * uv.y);
}

vec3 tri(in vec3 x){return abs(x-floor(x)-.5);} // Triangle function.

float surfFunc(in vec3 p){
    
	return dot(tri(p*0.5 + tri(p*0.25).yzx), vec3(0.666));
}
float map(vec3 p){
    float sf = surfFunc(p - vec3(0, cos(p.z/3.)*.15, 0));

//     
    // Round tunnel.
    // For a round tunnel, use the Euclidean distance: length(tun.y)
    vec2 tun = (p.xy - (p.z))*vec2(0.5, 0.7071);
    float n = 1.- length(tun) + (0.5 - sf);
    return min(n, p.y + 1.0);  

    
/*
    // Rounded square tunnel using Minkowski distance: pow(pow(abs(tun.x), n), pow(abs(tun.y), n), 1/n)
    vec2 tun = abs(p.xy - path(p.z))*vec2(0.5, 0.7071);
    tun = pow(tun, vec2(4.));
    float n =1.-pow(tun.x + tun.y, 1.0/4.) + (0.5 - sf);
    return min(n, p.y + FH);
*/
 
}
vec4 raymarchMap(in vec3 pos){
	// pos = axisCorrection(pos);
	// pos.z = mod(pos.z, 4.5);
	pos.yz = _rotate(pos.yz, 0.5*PI);
	// pos.z= pModPolar(pos.xz,6);
	// pos.x -= 1.0;
		// pos.y += 0.9;
	vec4 res = vec4(1.);
	vec2 Size = vec2(.5, 1.);
	// vec3 offset = vec3(-.850, 0.0, 0.0);
	float ics = fOctahedron(pos, Size.x, Size.y);
	// ics = min(ics, fOctahedron(pos - offset, Size.x, Size.y));
	// vec3 n = sampleNormal();
	// res = castRay(pos, n);
	res = opUnion(res, vec4(ones, ics));
	return res; 
}
vec3 test = vec3(0.);
vec4 renderMainImage(void)
{
	vec4 fragColor = vec4(0);
	vec2 fragCoord = _xy;
	float a = (bass_time) * 3.;
	vec2 uv = (fragCoord - RENDERSIZE.xy*0.5)/RENDERSIZE.y;
	/* 
	)--]--[-----from Shane: "Abstract Corridor"-------]--[--(
	)--]--[-- https://www.shadertoy.com/view/MlXSWX --]--[--(	
	*/
	vec3 camPos = vec3(0.0, 0.0, a*5.); // Camera position, doubling as the ray origin.
	vec3 lookAt = camPos + vec3(0.0, 0.1, 0.5+TIME);  // "Look At" position.
 
    // Light positioning. One is a little behind the camera, and the other is further down the tunnel.
 	vec3 light_pos = camPos + vec3(0.0, 0.125, -0.125);// Put it a bit in front of the camera.
	vec3 light_pos2 = camPos + vec3(0.0, 0.0, 6.0);// Put it a bit in front of the camera.

	// Using the Z-value to perturb the XY-plane.
	// Sending the camera, "look at," and two light vectors down the tunnel. The "path" function is 
	// synchronized with the distance function.
	// lookAt.xy += path(lookAt.z);
	// camPos.xy += path(camPos.z);
    // Using the above to produce the unit ray-direction vector.
    float FOV = PI/3.; // FOV - Field of view.
    vec3 forward = normalize(lookAt-camPos);
    vec3 right = normalize(vec3(forward.z, 0., -forward.x )); 
    vec3 up = cross(forward, right);

	// rd - Ray direction.
    vec3 rd = normalize(forward + FOV*uv.x*right + FOV*uv.y*up);
	float t = 0.0, dt;
	for(int i=0; i<128; i++){
		dt = map(camPos + rd*t);
		if(dt<0.005 || t>150.){ break; } 
		t += dt*0.75;
	}
	
	vec3 color = vec3(0);
	vec2 pixel = vec2(1.0) / RENDERSIZE;
	vec2 st = gl_FragCoord.xy * pixel;
	uv = ratio(st, RENDERSIZE);
	color = raymarch(cam, uv).rgb;
	color = linear2gamma(color);
	fragColor = vec4(color, 1.0);
	return fragColor;
}
vec4 renderMain() {
	if(PASSINDEX == 0) {
		return renderMainImage();
	}
	// if(PASSINDEX == 1) {
		// return renderMainImage();
	// }
}