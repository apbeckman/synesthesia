

////////////////////////////////////////////////////////////
// UltimateKaliCircuits  by mojovideotech
//
// based on :
// shadertoy/XlX3Rj  by Kali
//
// Creative Commons Attribution-NonCommercial-ShareAlike 3.0
////////////////////////////////////////////////////////////



#define 	pisq  	9.869604401089359	// pi squared, pi^2
#define 	twpi  	6.283185307179586  	// two pi, 2*pi
#define 	phicu  	4.236067977499791 	// phi cubed, phi^3
#define 	cuphi  	1.173984996705329   // cube root of phi
#define 	rctwpi  0.159154943091895	// reciprocal of twpi, 1/twpi      
#define		r36		0.027777777777778


vec3 color = vec3(0.0,0.5,0.0);
float S = glow;
float T = rate*smoothTime*0.0025;

vec2 hash22(vec2 p) { return fract(vec2(21.0, 97.0)*sin(dot(p, vec2(92.0, 61.0)))); }

void formula(vec2 z, float f) 
{
	float m = 0.0; 
	float o, ot2, ot=ot2=10000.0;
	for (int i=0; i<11; i++) {
		z = abs(z)/clamp(dot(z,z), 0.1, 0.5)-f;
		float l = length(z);
		o = min(max(abs(min(z.x, z.y)),-l+0.25), abs(l-0.25));
		ot = min(ot, o);
		ot2 = min(l*0.1, ot2);
		m = max(m, float(i)*(1.0-abs(sign(ot-o))));
	}
	m += 1.0;
	float w = intensity*m*2.0;
	float circ = pow(max(0.0,w-ot2)/w,6.0);
	S += max(pow(max(0.0,w-ot)/w,0.25),circ);
	vec3 col = vec3(hash22(z),f);
	col += highhits*1.6;
    color += col*(0.4+mod(m/9.0-(smoothTimeC*0.00095)*trace+ot2*2.0, 1.0)*1.6);
	color += vec3(1.0, 0.7, 0.3)*circ*(10.0-m)*3.0
			 *smoothstep(0.0, 0.5, vec3(f, _uv));
}

vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);

	vec2 pos = 2.0 * _xy.xy - RENDERSIZE.xy;
	pos /= max(RENDERSIZE.x,RENDERSIZE.y);
	vec2 uv = pos-center;
	uv *= 4.0-zoom;
	float a = T + mod(T, floor(runtime))*cuphi;
	float b = a*phicu;
	uv *= mat2(cos(b),sin(b),-sin(b),cos(b));
	uv += vec2(sin(a),cos(a*cuphi))*pisq;
	uv *= offset;
	float pix = cuphi/RENDERSIZE.x*offset;
	float c = 1.5+mod(floor(T), 16.0)*0.125;
	for (int aa=0; aa<36; aa++) {
		vec2 aauv = floor(vec2(float(aa)*rctwpi, mod(float(aa), twpi)));
		formula(uv+aauv*pix, c);
	}
	S *= r36, color *= r36;
	vec3 colo = mix(vec3(0.025), color, S)*(1.5-length(pos)); 	
	colo *= vec3(1.2, 1.1, 1.0);
	out_FragColor = vec4(colo, 1.0);

return out_FragColor; 
 } 
