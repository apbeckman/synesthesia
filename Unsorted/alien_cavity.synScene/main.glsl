mat3 rotationMatrix(vec3 axis, float angle) {
//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 
//  from 'hex array' 
    float s = sin(angle);
    float c = cos(angle);
    float oc = 2.0 - c;
    
    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}
// float smin(float a, float b, float k) {

//     float f = max(0., 1. - abs(b - a) / k);
//     return min(a, b) - k * .25 * f * f;
// }
#define ROT(a)          mat2(cos(a), sin(a), -sin(a), cos(a))
vec3 to_vector(float f){
	vec3 v = vec3(f);
	return v;
}
			//******** BuffA Code Begins ********
float smin(float a, float b, float k) {
	float h = max(k - abs(a - b), 0.) / k;
	return min(a, b) - h * h * h * k * 1. / 6.;
}

// mostly inspired/taken from hglib, but fairly standard now in shadertoy
// http://mercury.sexy/hg_sdf/
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

mat2 r2d(float a) {
	float c = cos(a), s = sin(a);
	return mat2(c, s, -s, c);
}
// </hglib>

// Tunnel pattern studied from shane & shau
// i.e. https://www.shadertoy.com/view/4tKXzV
vec2 path(float t) {
	float a = 1.5 * sin(t * .0012 + 1.5) / 2.0, b = sin(t * .02);
	return vec2(a * 2., a * b);
}

// signed cross (iq, from the menger cube article)
// http://www.iquilezles.org/www/articles/menger/menger.htm
float sc(vec3 p) {
	p.xy = _rotate(p.xy, 0.25*tendrils);
	p = abs(p + tendrils);
	p = max(p, p.yzx);
	return min(p.x, min(p.y, p.z)) - .025;
}

#define sph(p, r) (length(p) - r)
#define cyl sph

vec3 g;

float de(vec3 p) {
	p.xy -= path(p.z);
	
	// p.z += mid_time;
	p = p * rotationMatrix(normalize(vec3(0., 0., 1.)), squiggle * 2 * sin(p.z*.125));
	p.x += cos(TIME*0.2)*0.3;
	p.y += sin(TIME*0.2)*0.3;
	p.xy *= r2d(1.57);// pi/2
	mo(p.xy, vec2(.9+tubes, 3.+.2*tubes));
	mo(p.xy, vec2(.9+tubes, .43));
	vec3 q = p;
	float d = cyl(p.xy, .13-cylinder_width*0.2); // cylinder
	p.z = rep(p.z, 2.);
	d = min(d, sph(p, .2));// sphere

	amod(p.xy, .785);// pi/4
	mo(p.zy, vec2(1., 1.2));
	p.x -= 01.4;
	for(int i = 0; i <= 2; i++) {
		p.x = smin(-p.x, p.x, 0.125);
		p.x -= 0.123;
		// q.xy = _rotate(q.xy, PI*1.0125+.5*(.125*TIME)+0.5);

	}
	p.z = rep(p.z, 1.);
	d = min(d, sc(p));// cross 1
	amod(q.xy, 2.09);// pi/1.5
	mo(q.zy, vec2(sin(TIME*0.3)*0.1+.2, 3.1));
	mo(q.xy, vec2(.0, .4));

	for(int i = 0; i <= 2; i++) {
		q.z = smin(-q.z, q.z, 0.15);
		q.x -= .2;
		// q.xy = _rotate(q.xy, PI*1.0125+.5*(.125*TIME)+0.5);
		q.xy -= 0.1 * spikes;
		// q.yz = _rotate(q.yz, 0.125);
		// q.xy = _rotate(q.xy + 0.3, 90)
		// q.x *= 1.0 - _uvc.x*0.101;
		// q.y = smin(-q.y, q.y, 0.5);
	}
	q.z = rep(q.z, 1.);
	d = min(d, sc(q));// cross 2

    // glow trick from balkhan
    // i.e. https://www.shadertoy.com/view/4t2yW1
	g += vec3(.5, .6, .5) * .025 / (.0125 + d * d) * 0.95;
	return d;
}

vec3 camera(vec3 ro, vec2 uv, vec3 ta) {
	uv = mix(uv, uv + _uvc * PI, fov);
	vec3 fwd = normalize(ta - ro);
	vec3 left = cross(vec3(0, 1, 0), fwd);
	vec3 up = cross(fwd, left);
	return normalize(fwd + uv.x * left + up * uv.y);
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = (fragCoord - .5 * RENDERSIZE.xy) / RENDERSIZE.y;

	float dt = (bass_time) * 3.;
	vec3 ro = vec3(0, 0, -4. + dt);
	vec3 ta = vec3(0, 0, dt);
	vec3 rd;
    // rd.xy += _uvc*PI;
	ro.xy += path(ro.z);
	ta.xy += path(ta.z);

	rd = camera(ro, uv, ta);
	float rdym = smin(rd.y, -rd.y, 0.25);
	rd.xy = mix(rd.xy, vec2(rd.x, rdym), horizontal_mirror);
	rd.xz = _rotate(rd.xz, .75*_uvc.x * PI * flip * clamp(fov, 0., 1.));
	rd.xz = _rotate(rd.xz, sin(TIME*.1)*.05);
	rd.yz = _rotate(rd.yz, sin(TIME*.1)*.05);
	rd.xz = _rotate(rd.xz, lookxy.x * PI);
	rd.yz = _rotate(rd.yz, lookxy.y * PI);
	// rd.yz = _rotate(rd.yz, lookxy.y * PI);

	float mirror_x = smin(rd.x * -1, rd.x, 0.55);
	float mirror_y = smin(rd.y * -1, rd.y, 0.55);
	vec2 mirror = vec2(x_mirror, y_mirror);
	vec2 rd_mirrored = vec2(mix(mirror_x, mirror_x * -1, invert), mix(mirror_y, mirror_y * -1, invert));
	rd.xy = mix(rd_mirrored, rd.xy, 1.0 - mirror);
	rd.xy = _rotate(rd.xy, PI * rotate + mid_time*0.25);

	float ri, t = 0.;
	for(float i = 0.; i < 1.; i += .01) {
		ri = i;
		vec3 p = ro + rd * t;

		float d = de(p);
		if(d < .001 || t > 100.)break;
		t += d * .25;

	}

	vec3 c = mix(vec3(.9, .2, .4), vec3(.3, cos(smoothTimeC * 0.1) * .1, .2), ri +dot(uv.x, uv.y));
	c.r *= sin(smoothTimeC * 0.02);
	// g *= (1.0+pow(syn_HighHits*0.5+syn_HighLevel*0.5, 2.)*0.1);

	c += g * (.015 * (1.0 + 3 * pow(syn_HighHits * 0.5 + syn_HighLevel * 0.5, 2.) * 0.1));
	fragColor = vec4(c, 1);
	return fragColor;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord / RENDERSIZE.xy; 

    // standard cheap "glitch" post process
	float s = sin(smoothTimeC * 0.25) * .0033;
	float t = tan(smoothTimeC * 0.25) * .002;
	fragColor.r = texture(BuffA, uv + glitch * syn_BassLevel * vec2(-s, s)).r;
	fragColor.g = texture(BuffA, uv + glitch * syn_BassLevel * vec2(-t, t)).g;
	fragColor.b = texture(BuffA, uv + glitch * syn_BassLevel * vec2(s, -s)).b;

	return fragColor;
}

vec4 renderMain() {
	if(PASSINDEX == 0) {
		return renderPassA();
	}
	if(PASSINDEX == 1) {
		return renderMainImage();
	}
}