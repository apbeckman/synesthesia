

////////////////////////////////////////////////////////////////////
// KaliCircuitsExplorer  by mojovideotech
//
// based on :
// shadertoy.com/XlX3Rj  by Kali
//
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0
////////////////////////////////////////////////////////////////////

#ifdef GL_ES
precision highp float;
#endif
mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }


// IQ's vec2 to float hash.
float hash21(vec2 p){  return fract(sin(dot(p, vec2(27.619, 57.583)))*43758.5453); }


// vec2 to vec2 hash.
vec2 hash22B(vec2 p) { 

    // Faster, but doesn't disperse things quite as nicely. However, when framerate
    // is an issue, and it often is, this is a good one to use. Basically, it's a tweaked 
    // amalgamation I put together, based on a couple of other random algorithms I've 
    // seen around... so use it with caution, because I make a tonne of mistakes. :)
    float n = sin(dot(p, vec2(1, 113)));
    p = fract(vec2(262144, 32768)*n)*2. - 1.; 
    return sin(p*6.2831853 + smoothTimeC*0.125);   
}

// vec2 to vec2 hash.
vec2 hash22C(vec2 p) { 

    // Faster, but doesn't disperse things quite as nicely. However, when framerate
    // is an issue, and it often is, this is a good one to use. Basically, it's a tweaked 
    // amalgamation I put together, based on a couple of other random algorithms I've 
    // seen around... so use it with caution, because I make a tonne of mistakes. :)
    float n = sin(dot(p, vec2(289, 41)));
    return fract(vec2(262144, 32768)*n)*2. - 1.;
    
    // Animated.
    p = fract(vec2(262144, 32768)*n)*2. - 1.; 
    return sin(p*6.2831853 + smoothTime*0.0125); 
}


// Based on IQ's gradient noise formula.
float n2D3G( in vec2 p ){
   
    // Cell ID and local coordinates.
    vec2 i = floor(p); p -= i;
    
    // Four corner samples.
    vec4 v;
    v.x = dot(hash22C(i), p);
    v.y = dot(hash22C(i + vec2(1, 0)), p - vec2(1, 0));
    v.z = dot(hash22C(i + vec2(0, 1)), p - vec2(0, 1));
    v.w = dot(hash22C(i + 1.), p - 1.);

    // Cubic interpolation.
    p = p*p*(3. - 2.*p);
    
    // Bilinear interpolation -- Along X, along Y, then mix.
    return mix(mix(v.x, v.y, p.x), mix(v.z, v.w, p.x), p.y);
    
}

// Two layers of noise.
float fBm(vec2 p){ return n2D3G(p)*.57 + n2D3G(p*2.)*.28 + n2D3G(p*4.)*.15; }


#define 	pi   	PI 	// pi
	
float S = (101.0*(1.0+highhits) + glow) * (intensity);

vec3 color = vec3(0.0);

void formula(vec2 z, float t) 
{
	S *= (1.0+0.05*highhits);
	
	float M = 0.0;
	float o, ot2, ot=ot2=1000.0;
	float K = floor(loops/4.0)+floor(5.0 * zoom);
	for (int i=0; i<11; i++) {
		z = abs(z) / clamp(dot(z, z), 0.1, 0.5) - t;
		float l = length(z);
		o = min(max(abs(min(z.x, z.y)), -l + 0.25), abs(l - 0.25));
		ot = min(ot, o);
		ot2 = min(l * 0.1, ot2);
		M = max(M, float(i) * (1.0 - abs(sign(ot - o))));
		if (K <= 0.0) break;
		K -= 1.0;
	}
	M += 1.0;
	float w = (intensity) * M+highhits*0.00125;
	float circ = pow(max(0.0, w - ot2) / w, 6.0);
	S += max(pow(max(0.0, w - ot) / w, 0.25), circ);
	vec3 col = normalize(0.1 + vec4(01.45, 0.75, M * 0.1, 1.0).rgb);

	color += col * (0.4 + mod(M / 9.0 - t * pulse + ot2 * 2.0, 1.0));
	color += vec3(1.0, 0.7, 0.3) * circ * (10.0 - M) * (3.0+highhits*0.25);
}


vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);

	float R = 0.0;
	float N = (smoothTimeC*0.075) * 0.01 * rate;
	float T = 2.0 * rate;
	if (N > 6.0 * rate) { 
		R += 1.0;
		N -= (R * 8.0 * rate);
	}
	if (N < 6.0 * rate) T += N;
	else  T = 8.0 * rate - N;
	float Z = (1.05-zoom);
	vec2 pos = _xy.xy / RENDERSIZE.xy - 0.5;
		//pos *= (1.0+Warp*(((_uvc.xy))-1.));

	pos.x *= RENDERSIZE.x/RENDERSIZE.y;
	vec2 uv = pos ;
	uv -= _uvc*Zoom2*PI;

	//uv *= -1.0 + Warp * (1.0 - (_uvc*PI-uv));
	uv -= center;
	uv = _rotate(uv, Rotate);

	float sph = length(uv)*0.1; 
	sph = sqrt(1.0 - sph * sph) * 2.0 ;
	float a = T * pi;
	float b = a + T*2.;
	float c = cos(a) + sin(b);

	uv *= mat2(cos(b), sin(b), -sin(b), cos(b));

	//uv *= 1.0+ Warp*(_uvc*fBm(_rotate(_uvc*uv+(fBm(vec2(_rotate(_uvc, smoothTime*0.05) ))), (smoothTime*0.0125))-1.));
	uv *= mat2(cos(a),-sin(a), sin(a),cos(a));
	uv -= vec2(sin(c), cos(c)) / pi;

	uv *= Z;

	float pix = 0.5 / RENDERSIZE.x * Z / sph;
	float dof = (zoom * focus) + (T * 0.25);
	float L = floor(loops);
	for (int aa=0; aa<24; aa++) {
		vec2 aauv = floor(vec2(float(aa) / 6.0, mod(float(aa), 6.0)+_uvc*PI));
		formula(uv + aauv * pix * dof, T);
		if (L <= 0.0) break;
		L -= 1.0;
	}
	S /= floor(loops); 
	color /= floor(loops);
	vec3 colo = mix(vec3(0.125), color, S) * (1.0 - length(pos))*min(1.,abs(.5+mod(smoothTimeB+.5,1.))*10.);	
	//colo*= 1.0 +highhits*0.5;
	colo *=vec3(1.2, 1.1, 1.0);

	out_FragColor = sqrt(max(vec4(colo, 1.0), 0.0) -0.2);

return out_FragColor; 
 } 
