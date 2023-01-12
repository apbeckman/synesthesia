


////////////////////////////////////////////////////////////////////
// DiamondVision  by mojovideotech
//
// based on :
// glslsandbox.com\/e#46411.0
//
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0
////////////////////////////////////////////////////////////////////

#ifdef GL_ES
precision highp float;
#endif
#define		eConst  2.7182818284590452353602874713527
#define    	c30    	cos(30.)     	// cos 30
#define    	twpi   	PI*2.0       // two pi, 2*pi
#define 	pi   	PI 		// pi
#define		phi		2.61803398875
#define		piphi	pi * phi		// pi*phi
#define 	sqpi 	pow(pi, 0.5)		// square root of pi
//#define 	phi   	phi		// golden ratio
#define 	cucuphi	pow(pow(phi, 1.0/3.0), 1.0/3.0)		// cube root of cube root of phi
#define		epi		eConst/PI		// e/pi
#define 	time2	bass_time*2.0+50.0
#define 	time3	smoothTimeC*0.0125+50.0


vec2 rotate(in vec2 r, in float o) { return vec2(cos(o)*r.x + sin(o)*r.y, -sin(o)*r.x + cos(o)*r.y); }

float torus(in vec3 pos, in vec2 tor) { 
	vec2 qt = abs(vec2(max(abs(pos.x), abs(pos.z))-tor.x, pos.y));
	return max(qt.x, qt.y)-tor.y;
}

float trap(in vec3 tp) {
	return abs(min(torus(tp, vec2(epi, 0.125)), max(abs(tp.z)-0.125, abs(tp.x)-0.125)))-0.0025;
}

float map(in vec3 pm) {
	float c = dot(abs((pm.yz)),vec2(0.5))-0.05;
	vec3 m = abs(1.0-mod(pm,2.0));
	m.yz = rotate(m.yz+test2, sqrt(time3));
	float e = 9999.9999, f = 1.0;
	for (float i = 0.0; i < 4.0+its; i++) {
		m.xz = rotate(m.xz, radians(i*(1./3.)+time3));
		m.zy = rotate(m.yz, radians((i+i)*(2./3.)+time3*phi));
		m = abs(1.0-mod(m+i/(3.0+test+(0.5*sin(spin_time*0.5))*2.),2.0));
		m *=abs(sqrt(m)*sqpi);
		f *= 0.5;
		e = min(e, trap(m) * f);
	}
	return max(e, -c);
}

vec3 hsv(in float h, in float s, in float v) {
	return mix(vec3(1.0), clamp((abs(fract(h + vec3(3.0, 2.0, 1.0) / 3.0) * 6.0 - 3.0) - 1.0), 0.0 , 1.0), s) * v;
}

vec3 intersect(in vec3 rayOrigin, in vec3 rayDir) {
	float d = 1.0, it = 0.0;
	vec3 p = rayOrigin, col = vec3(0.0);
	float md = phi+sin(time2*0.25)*0.25;
	for (int i = 0; i < 50; i++) {		
		if (d < 0.000099999) continue;
		d = map(p);
		p += d*rayDir; 
		md = min(md, d);
		it++;
	}
	if (d < 0.0001) {
		float x = (it/49.0);
		float y = (d-0.01)/0.01/(49.0);
		float z = (0.01-d)/0.01/49.0;
		float q = 1.0-x-y*2.+z;
		col = hsv(q*0.82+0.5, 1.0-q*epi, q);
	} 
		col += hsv(d, 1.0,2.0)*md*(28.0);
	return col;
}
/*
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    
    // Unit direction ray vector: Note the absence of a divide term. I came across
    // this via a comment Shadertoy user "coyote" made. I'm pretty happy with this.
    vec3 rd = vec3(2.*fragCoord - RENDERSIZE.xy, RENDERSIZE.y);
    
    // Barrel distortion;
    rd = normalize(vec3(rd.xy, sqrt(max(rd.z*rd.z*FOVmod - dot(rd.xy, rd.xy)*.2, 0.))));
    
    // Rotating the ray with Fabrice's cost cuttting matrix. I'm still pretty happy with this also. :)
    //vec2 m = sin(vec2(1.57079632, 0) + TIME/8.);
    //rd.xy = rd.xy*mat2(m.xy, -m.y, m.x);
    //rd.xz = rd.xz*mat2(m.xy, -m.y, m.x);
    rd.yz = _rotate(rd.yz, lookXY.y*PI);
    rd.xy = _rotate(rd.xy, -1.0*lookXY.x*PI);
    rd.xz = _rotate(rd.xz, lookZ*PI);
    // Ray origin: Sending it along the Z-axis.
    //vec3 ro = vec3(0, 0, TIME*rate*0.25);
    vec3 ro = vec3(0, 0, ((smoothTime) * rate));
*/

vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 ps = -1.0 + 2.0 * _xy.xy / RENDERSIZE.xy;
	ps.x *= RENDERSIZE.x / RENDERSIZE.y;
	vec3 up = vec3(0, -1, 0);
	vec3 cd = vec3(1, 0, 0);
	vec3 co = vec3(( bass_time*0.125+50.), 0, 0);
	vec3 uw = normalize(cross(up, co));
	vec3 vw = normalize(cross(cd, uw));
	//vec3 rd = normalize(uw * ps.x + vw * ps.y + cd*(1.0-length(ps)*phi));
    vec3 rd = vec3(2.*fragCoord - RENDERSIZE.xy, RENDERSIZE.y);

	rd = normalize(vec3(rd.xy*(1.0+(Mirror*(-1+_uvc*PI))), sqrt(max(rd.z*rd.z - dot(rd.xy, rd.xy)*.2, 0.)*(FOVmod-Flip*0.5))*FOVmod));
//    rd = normalize(vec3(rd.xy*(1.0+(Mirror*(-1+_uvc*PI))), (rd.z - length(rd.xy-_uvc*PI*Fisheye)*.125)));

	rd.yz = _rotate(rd.yz, lookXY.y*PI+_uvc.y*Flip*PI);
    rd.xy = _rotate(rd.xy, lookXY.x*PI);
    rd.xz = _rotate(rd.xz, Rotate*PI);

	out_FragColor = vec4(vec3(intersect(co, rd)),1.0);

return out_FragColor; 
 } 
 