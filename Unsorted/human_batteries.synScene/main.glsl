

// 'Human Batteries' dean_the_coder (Twitter: @deanthecoder)
// https://www.shadertoy.com/view/wlycRR
//
// "The human body generates more bio electricity than a
//  120 volt battery and over 25000 BTUs of body heat.
//  Combined with a form of fusion the machines had found
//  all the energy they would ever need."
//   - The Matrix (1999)
//
// Tricks to get the performance:
//   - There's only one 'pod'.
//	 The others are all cloned by folding space and
//	 mirroring the axis.
//   - There's only one column.
//	 Again, mirroring the axis x/z axis creates more.
//   - Only two lighting bolts.
//	 Axis mirroring again...
//
// Thanks to Evvvvil, Flopine, Nusan, BigWings, Iq, Shane
// and a bunch of others for sharing their knowledge!

// License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License

#define MIN_DIST		 .0015
#define MAX_DIST		 80.0
#define MAX_STEPS		 100.0
#define SHADOW_STEPS	 20.0
#define MAX_SHADOW_DIST  10.0
#define BODY_STEPS	     30.0
#define MAX_BODY_DIST    1.7

float T; // Global time.
vec2 g = vec2(0); // Glow for lightning and pods.

struct Hit {
	float d; // SDF distance.
	int id;  // Material ID.
	vec3 uv; // Ray position.
};

// Thnx Dave_Hoskins - https://www.shadertoy.com/view/4djSRW
#define HASH  p = fract(p * .1031); p *= p + 3.3456; return fract(p * (p + p));

vec2 hash22(vec2 p) { HASH }
vec4 hash44(vec4 p) { HASH }

float n31(vec3 p) {
	// Thanks Shane - https://www.shadertoy.com/view/lstGRB
	const vec3 s = vec3(7, 157, 113);
	vec3 ip = floor(p);
	p = fract(p);
	p = p * p * (3. - 2. * p);

	vec4 h = vec4(0, s.yz, s.y + s.z) + dot(ip, s);
	h = mix(hash44(h), hash44(h + s.x), p.x);

	h.xy = mix(h.xz, h.yw, p.y);
	return mix(h.x, h.y, p.z);
}

float n21(vec2 p) { return n31(vec3(p, 1)); }

float smin(float a, float b, float k) {
	float h = clamp(.5 + .5 * (b - a) / k, 0., 1.);
	return mix(b, a, h) - k * h * (1. - h);
}

void minH(inout Hit a, Hit b) {
	if (b.d < a.d) a = b;
}

mat2 rot(float a) {
	float c = cos(a), s = sin(a);
	return mat2(c, s, -s, c);
}

vec2 opModPolar(vec2 p, float n)
{
	float angle = 3.141 / n,
		  a = mod(atan(p.y, p.x), 2. * angle) - angle;
	return length(p) * vec2(cos(a), sin(a));
}

float sdBox(vec3 p, vec3 b) {
	vec3 q = abs(p) - b;
	return length(max(q, 0.)) + min(max(q.x, max(q.y, q.z)), 0.);
}

float sdCapsule(vec3 p, float h, float r) {
	p.x -= clamp(p.x, 0., h);
	return length(p) - r;
}

float sdTorus(vec3 p, vec2 t) {
  return length(vec2(length(p.xz) - t.x, p.y)) - t.y;
}

vec3 getRayDir(vec3 ro, vec3 lookAt, vec2 uv) {
	vec3 f = normalize(lookAt - ro),
		 r = normalize(cross(vec3(0, 1, 0), f));
	return normalize(f + r * uv.x + cross(f, r) * uv.y);
}

Hit sdPod(vec3 p) {
	// Pod.
	float d = exp(-p.x * .6) * .3 - .07 * p.x * sin(4.7 * p.x) + .06, d2;
	Hit h = Hit(sdCapsule(p, 1.5, d), 3, p);

	// The end bit.
	d2 = sdBox(p - vec3(1.6, 0, 0), vec3(.05)) - .03;

	// Feeding pipes.
	d2 = min(min(d2, sdTorus(p + vec3(.6, 0, 0), vec2(.8, .05))),
			 sdTorus(p.xzy + vec3(1.75, 0, 0), vec2(2, .1)));

	// Pod 'cage'.
	p.yz = abs(p.yz);
	p.y -= d * .5;
	p.z -= d;
	d = sdCapsule(p, 1.65, .03 + abs(sin(p.x * 30.)) * .005);
	minH(h, Hit(min(d2, d) * .9, 2, p));

	g.y += .0001 / (.2 + d * d);

	return h;
}

float sdBody(vec3 p) {
	// Map world to pod point.
	p.xz = opModPolar(abs(p.xz) - 20., 60.);
	p.y = mod(p.y, 3.) - 1.5;
	p.x -= 12.1;

	// Head
	float d = length(p) - .07;

	// Torso
	p.z = abs(p.z);
	p.x -= .1;
	p.xy *= rot(-.2);
	d = smin(d, sdCapsule(p, .3, p.x * .09 + .02), .1);

	// Legs
	p.x -= .35;
	p.z -= .06;
	p.xy *= rot(.3);
	p.xz *= rot(.15);
	return smin(d, sdCapsule(p, .5, .01), .2);
}

float sdBolts(vec3 p, float i) {
	p = mix(p, p.zyx, step(12., p.x));

	float d, t = T + sign(p.z) + i,
		  r = n21(vec2(t, i * .2)) - .5;

	p.x += 10.;
	p.y += mod(70. - t * 10. * i, 70.) - 18.
		   + n21(vec2(p.x, t * 15.2 + i)) * 2.;

	p.z = abs(p.z) - 24. - r * 5.;

	p.xy *= rot(r);
	p.xz *= rot((i - 1.) * .5);

	d = sdCapsule(p, 20., .01);
	g.x += .01 / (.01 + d * d);

	return d * .6;
}

// Map the scene using SDF functions.
Hit map(vec3 p) {
	// Lightning.
	// We only render two bolts - The others
	// are added using axis reflection.
	Hit h = Hit(min(sdBolts(p, 1.), sdBolts(p, .4)), 6, p);

	// Main cylinders.
	// Only one created - We use axis reflection again.
	p.xz = abs(p.xz) - 20.;
	minH(h, Hit(length(p.xz) - 12. + sin(p.y * 2.09 + 4.7) * .2, 1, p));

	// Pods.
	// Again, only one is ever rendered!
	p.xz = opModPolar(p.xz, 60.);
	p.x -= 11.8;
	p.y = mod(p.y, 3.) - 1.5;
	minH(h, sdPod(p));

	return h;
}

vec3 calcN(vec3 p, float t) {
	float h = t * .5;
	vec3 n = vec3(0);
	for (int i = min(int(FRAMECOUNT), 0); i < 4; i++) {
		vec3 e = .005773 * (2. * vec3((((i + 3) >> 1) & 1), (i >> 1) & 1, i & 1) - 1.);
		n += e * map(p + e * h).d;
	}

	return normalize(n);
}

float calcShadow(vec3 p, vec3 ld) {
	// Thanks iq.
	float s = 1., t = .1;
	for (float i = 0.; i < SHADOW_STEPS; i++)
	{
		float h = map(p + ld * t).d;
		s = min(s, 15. * h / t);
		t += h;
		if (s < .01 || t > MAX_SHADOW_DIST) break;
	}

	return clamp(s, 0., 1.);
}

// March through the pod - Basically the same as a shadow.
float bodyTint(vec3 p, vec3 rd) {
	float s = 1.;
	const float stp = MAX_BODY_DIST / BODY_STEPS;
	for (float t = 0.; t < MAX_BODY_DIST; t += stp)
	{
		float h = sdBody(p + rd * t);
		s = min(s, 20. * h / t);
		if (s < .01) break;
	}

	return clamp(s, 0., 1.);
}

// Quick ambient occlusion.
float ao(vec3 p, vec3 n, float h) { return map(p + h * n).d / h; }

// Sub-surface scattering. (Thanks Evvvvil)
float sss(vec3 p, vec3 ld, float h) { return smoothstep(0.0, 1.0, map(p + ld * h).d / h); }

/**********************************************************************************/

vec3 vignette(vec3 c, vec2 fc) {
	vec2 q = fc.xy / RENDERSIZE.xy;
	c *= .5 + .5 * pow(16. * q.x * q.y * (1. - q.x) * (1. - q.y), .4);
	return c;
}

vec3 lights(vec3 p, vec3 rd, float d, Hit h) {
	vec3 ld = normalize(vec3(6, 3, -10) - p),
		 n = calcN(p, d), c;
	float ss = 0., // Sub-surface scatter
		  sp = 10.; // Specular multiplier.

	if (h.id == 1) {
		// Column.
		// Small ridges added using a quick bump map.
		c = vec3(.01);
		n.y += (abs(sin(h.uv.y * 31.)) - .5) * .4;
		n = normalize(n);
	} else if (h.id == 2) {
		// Feeder Pipes.
		c = vec3(.02);
		sp = 50.; // Slightly increased specular to add 'shine'.
	} else if (h.id == 3) {
		// Pod pink.
		c = vec3(1, .32, .27) * smoothstep(1.55, 1.3, h.uv.x);
		c *= .2 + .8 * smoothstep(0., .2, h.uv.x);

		// Sub-surface scattering.
		ss = sss(p, ld, .45);
	} else c = vec3(1);

	// Ambient occlusion.
	float ao = dot(vec2(ao(p, n, .2), ao(p, n, .5)), vec2(.1, .2)),

	// Primary light.
	l1 = max(0., .1 + .9 * dot(ld, n))
		 * (.3 + .7 * calcShadow(p, ld)) // ...with shadow.
		 * (.3 + .7 * ao), // ...and _some_ AO.

	// Secondary(/bounce) light.
	l2 = max(0., .1 + .9 * dot(ld * vec3(-1, 0, -1), n)) * .3,

	// Fresnel
	fre = smoothstep(.7, 1., 1. + dot(rd, n)) * .5,

	// Specular.
	spe = smoothstep(0., 1., pow(max(0., dot(rd, reflect(ld, n))), 30.)) * sp;

	// Combine into final color.
	return mix((l1 + (l2 + spe) * ao + ss) * c * vec3(2, 1.6, 1.4),
			   vec3(.01),
			   fre);
}

vec3 march(vec3 ro, vec3 rd) {
	// Raymarch.
	vec3 p, c = vec3(0);

	float d = .01;
	Hit h;
	for (float i = 0.; i < MAX_STEPS; i++) {
		p = ro + rd * d;
		h = map(p);

		if (abs(h.d) < MIN_DIST || d > MAX_DIST)
			break;

		d += h.d; // No hit, so keep marching.
	}

	// Stash the glow before calculating normals/etc alters it.
	vec2 gg = g;

	if (d < MAX_DIST) {
		c = lights(p, rd, d, h);

		// If we hit a pod, do a mini-march through it
		// to build a tint of the 'bodies'.
		if (h.id == 3)
			c *= .4 + .6 * bodyTint(p - rd, rd);
	}

	c *= exp(d * d * -.001);

	// Raymarch to get fog depth.
#define FOG_STEPS      8.0
	float maxD = min(d, MAX_DIST),
		  d2 = 0.01;
	for (float i = 0.; i < FOG_STEPS; i++) {
		p = ro + rd * d2;
		float fd = 0.7 * abs(p.y + 15.0 + 10.0 * n31(vec3(p.xz * 0.1, T * 0.1)));
		if (abs(fd) < MIN_DIST || d2 > maxD)
			break;
		d2 += fd;
	}

	// Mix in the fog color.
	d = smoothstep(0., 1., (min(d, MAX_DIST) - d2) / 30.);
	c = mix(c,
			vec3(.25, .3, .4),
			d * d // Base fog.
			* n21(p.xz * .06) // Low frequency patches.
			* (.2 + .8 * n31(vec3(p.xz + T, T) * .2)) // Finer details.
			);

	return c // Base color.
		   + gg.x * vec3(.4, .6, 1) // Lighting bolt glow.
		   + gg.y * vec3(1, .32, .27); // Faint pod glow.
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fc = _xy;

	T = mod(smoothTime, 40.) + 1.3;

	float dim = 1. - abs(cos(clamp(T, -1., 1.) * 1.57)),
		  t = .5 + .5 * cos(.157 * T);
	t = mix(smoothstep(0., 1., t), t, t);
	vec3 ro = vec3(-6.5, 25. - t * 20., 21);
	ro.xz += t * 25.;

	vec2 uv = (fc - .5 * RENDERSIZE.xy) / RENDERSIZE.y;
	vec3 col = march(ro, getRayDir(ro, vec3(0, -t * 10., 0), uv));

	// Output to screen.
	fragColor = vec4(vignette(pow(col * dim, vec3(.45)), fc), 0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}