//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


/*
	PsyVR test
*/

#define R RENDERSIZE.xy
#define T TIME * .5

// toggle for psychedelic madness
//#define ENABLE_COLOR_CYCLE 0

// FabriceNeyret2 
#define hue(v)  (.5 + cos(6.3 * (v) + vec4(0, 23, 21, 0)))

int id = -1;

mat2 rotate(float a) {
	float c = cos(a),
		s = sin(a);
	return mat2(c, s, -s, c);
}

float random(in vec2 st) {
    return fract(sin(dot(st.xy, vec2(12.9898, 78.233))) * 43758.1);
}

float noise(vec2 p) {
	vec2 i = ceil(p);
    vec2 f = fract(p);
    vec2 u = f * f * (3. - 2. * f);
   	float a = random(i);
    float b = random(i + vec2(1., 0.));
    float c = random(i + vec2(0., 1.));
    float d = random(i + vec2(1., 1.));
    return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

float fbm(in vec2 p) { 
	float s = .0;
	float m = .0;
	float a = .5;	
	for(int i = 0; i < 8; i++) {
		s += a * noise(p);
		m += a;
		a *= .35;
		p *= 2.;
	}
	return s / m;
}

vec3 renderFractal(vec2 uv) {

    vec3 color = vec3(0.);
    vec2 p = uv;
	
    // per channel iters
    float t = T;
    for (int c = 0; c < 2; c++) {
    
        t += .1; // time offset per channel
        
		float l = 0.;
        float s = 1.;
        for (int i = 0; i < 8; i++) {
            // from Kali's fractal iteration
            p = abs(p) / dot(p, p);
            p -= s;
            p *= rotate(t * .5);
            s *= .8;
            l += (s  * .08) / length(p)*syn_HighLevel;
        }
        color[c] += l;
    
    }

	return color;

}

float map(vec3 p) {
	
    float m = 1000.;
    
    vec3 q = p;
    float k = fbm(q.xz + fbm(q.xz + smoothTimeC *0.125));
   	
    q.y += 1.1;
    float d = dot(q, vec3(0., 1., 0.)) + k;
	d = min(5. - d, d);
    if (d < m) { 
        m = d;
        id = 1;
    }
    
    q = p;
    q.xz = mod(q.xz + 2., 4.) - 2.;
    d = min(d, length(q.xz) - .5);
    if (d < m) { 
        m = d;
        id = 2;
    }
    
    return m;
}

vec3 render(vec3 ro, vec3 rd) {

    vec3 col = vec3(0.);
	vec3 p;
    
	float t = 0.;
	for (int i = 0; i < 256; i++) {
		p = ro + rd * t;
		float d = map(p);
		if (d < .001 || t > 50.) break;
		t += .5 * d;
#if ENABLE_COLOR_CYCLE 
        col += .012 * hue(d * .5 + T * .8).rgb;
#else
        col += .0075 * hue(d).rgb;
#endif
	}
    col /= 1.5;
    
    vec3 tex =  renderFractal(fract(.1 * p.xz) - .5);
    if (id == 1) col += tex / (1. + t * t * .5);
    if (id == 2) col += abs(.1 / sin(2. * abs(pow(p.y+_uv.y, 2.)) - T*3.)) * vec3(0., 1., 1.)*(1.0+pow(syn_HighLevel*0.5+syn_MidHighLevel*0.5, 1.0+syn_Intensity));
    
	return col;

}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 I = _xy;


    vec2 uv = (2. * I - R)
        / R.y;
	vec3 col = vec3(0.);
	
	vec3 ro = vec3(2., 1., bass_time * 2.);
	vec3 rd = vec3(uv, 1.);
	
    vec3 pc = render(ro, rd);
    
    fragColor = vec4(pc, 1.);
    return fragColor;
}
/*
void mainVR(out vec4 fragColor, vec2 I, vec3 roVR, vec3 rdVR) {        
    vec3 col = render(roVR + vec3(2., 1., T * 2.), rdVR);
	fragColor = vec4(col, 1.);
	return fragColor; 
 } 
*/

vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}