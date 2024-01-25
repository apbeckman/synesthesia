

			//******** BuffA Code Begins ********

#define G(ic,x) texture(ic, x)
//#define iC0 BuffD
//#define iC1 BuffA
#define A BuffA
#define B BuffB
#define C BuffC
#define D BuffD
#define o0 1.0*(1.0+0.1*syn_Intensity*syn_BassLevel)
#define o1 3.0*(1.0+0.1*syn_Intensity)
#define stddev 3.5

float gaussian(float x, float s) {
    return exp(-x*x/(s*s));
}

vec4 gaussian(vec4 x, float s) {
    return exp(-x*x/(s*s));
}

vec2 wrap(vec2 x) {
    return mod(mod(x, 1.0) + 1.0, 1.0);
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    fragCoord+= _uvc*(1.+syn_Intensity)*Zoom;
    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    
    vec2 ix = vec2(1.0 / RENDERSIZE.x, 0.0);
    vec2 iy = vec2(0.0, 1.0 / RENDERSIZE.y);
    vec4 i0 = vec4(-4.0, -3.0, -2.0, -1.0);
    vec4 i1 = vec4(1.0, 2.0, 3.0, 4.0);
    vec4 g0 = gaussian(i0, stddev);
    vec4 g1 = g0.wzyx;
    float g = gaussian(0.0, stddev);
    float sum = (2.0 * dot(g0, vec4(1.0)) + g);
    g0 /= sum;
    g1 /= sum;
    g /= sum;

    // 2 complete blur passes
    vec4 leftX0  = g0 * vec4(G(D, wrap(uv + o0 * ix * i0.x)).x, G(D, wrap(uv + o0 * ix * i0.y)).x, G(D, wrap(uv + o0 * ix * i0.z)).x, G(D, wrap(uv + o0 * ix * i0.w)).x);
    vec4 rightX0 = g1 * vec4(G(D, wrap(uv + o0 * ix * i1.x)).x, G(D, wrap(uv + o0 * ix * i1.y)).x, G(D, wrap(uv + o0 * ix * i1.z)).x, G(D, wrap(uv + o0 * ix * i1.w)).x); 
    float centerX0 = g * G(D, uv).x;
    float sumX0 = centerX0 + dot(leftX0, vec4(1.0)) + dot(rightX0, vec4(1.0));

    vec4 leftY0  = g0 * vec4(G(A, wrap(uv + o0 * iy * i0.x)).x, G(A, wrap(uv + o0 * iy * i0.y)).x, G(A, wrap(uv + o0 * iy * i0.z)).x, G(A, wrap(uv + o0 * iy * i0.w)).x);
    vec4 rightY0 = g1 * vec4(G(A, wrap(uv + o0 * iy * i1.x)).x, G(A, wrap(uv + o0 * iy * i1.y)).x, G(A, wrap(uv + o0 * iy * i1.z)).x, G(A, wrap(uv + o0 * iy * i1.w)).x); 
    float centerY0 = g * G(A, uv).x;
    float sumY0 = centerY0 + dot(leftY0, vec4(1.0)) + dot(rightY0, vec4(1.0));

    vec4 leftX1  = g0 * vec4(G(A, wrap(uv + o1 * ix * i0.x)).y, G(A, wrap(uv + o1 * ix * i0.y)).y, G(A, wrap(uv + o1 * ix * i0.z)).y, G(A, wrap(uv + o1 * ix * i0.w)).y);
    vec4 rightX1 = g1 * vec4(G(A, wrap(uv + o1 * ix * i1.x)).y, G(A, wrap(uv + o1 * ix * i1.y)).y, G(A, wrap(uv + o1 * ix * i1.z)).y, G(A, wrap(uv + o1 * ix * i1.w)).y); 
    float centerX1 = g * G(A, uv).y;
    float sumX1 = centerX1 + dot(leftX1, vec4(1.0)) + dot(rightX1, vec4(1.0));

    vec4 leftY1  = g0 * vec4(G(A, wrap(uv + o1 * iy * i0.x)).z, G(A, wrap(uv + o1 * iy * i0.y)).z, G(A, wrap(uv + o1 * iy * i0.z)).z, G(A, wrap(uv + o1 * iy * i0.w)).z);
    vec4 rightY1 = g1 * vec4(G(A, wrap(uv + o1 * iy * i1.x)).z, G(A, wrap(uv + o1 * iy * i1.y)).z, G(A, wrap(uv + o1 * iy * i1.z)).z, G(A, wrap(uv + o1 * iy * i1.w)).z); 
    float centerY1 = g * G(A, uv).z;
    float sumY1 = centerY1 + dot(leftY1, vec4(1.0)) + dot(rightY1, vec4(1.0));
    
    fragColor = vec4(sumX0, sumY0, sumX1, sumY1);
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

#define G(ic,x) texture(ic, x)
//#define iC0 BuffA
//#define iC1 BuffB
#define o01 o0*(9.)
#define o11 o1*(9.)
//#define stddev 2.5
/*
float gaussian(float x, float s) {
    return exp(-x*x/(s*s));
}

vec4 gaussian(vec4 x, float s) {
    return exp(-x*x/(s*s));
}

vec2 wrap(vec2 x) {
    return mod(mod(x, 1.0) + 1.0, 1.0);
}
*/
vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    
    vec2 ix = vec2(1.0 / RENDERSIZE.x, 0.0);
    vec2 iy = vec2(0.0, 1.0 / RENDERSIZE.y);
    vec4 i0 = vec4(-4.0, -3.0, -2.0, -1.0);
    vec4 i1 = vec4(1.0, 2.0, 3.0, 4.0);
    vec4 g0 = gaussian(i0, stddev);
    vec4 g1 = g0.wzyx;
    float g = gaussian(0.0, stddev);
    float sum = 2.0 * dot(g0, vec4(1.0)) + g;
    g0 /= sum;
    g1 /= sum;
    g /= sum;

    // 2 complete blur passes
    vec4 leftX0  = g0 * vec4(G(A, wrap(uv + o01 * ix * i0.x)).w, G(A, wrap(uv + o01 * ix * i0.y)).w, G(A, wrap(uv + o01 * ix * i0.z)).w, G(A, wrap(uv + o01 * ix * i0.w)).w);
    vec4 rightX0 = g1 * vec4(G(A, wrap(uv + o01 * ix * i1.x)).w, G(A, wrap(uv + o01 * ix * i1.y)).w, G(A, wrap(uv + o01 * ix * i1.z)).w, G(A, wrap(uv + o01 * ix * i1.w)).w); 
    float centerX0 = g * G(A, uv).w;
    float sumX0 = centerX0 + dot(leftX0, vec4(1.0)) + dot(rightX0, vec4(1.0));

    vec4 leftY0  = g0 * vec4(G(B, wrap(uv + o01 * iy * i0.x)).x, G(B, wrap(uv + o01 * iy * i0.y)).x, G(B, wrap(uv + o01 * iy * i0.z)).x, G(B, wrap(uv + o01 * iy * i0.w)).x);
    vec4 rightY0 = g1 * vec4(G(B, wrap(uv + o01 * iy * i1.x)).x, G(B, wrap(uv + o01 * iy * i1.y)).x, G(B, wrap(uv + o01 * iy * i1.z)).x, G(B, wrap(uv + o01 * iy * i1.w)).x); 
    float centerY0 = g * G(B, uv).x;
    float sumY0 = centerY0 + dot(leftY0, vec4(1.0)) + dot(rightY0, vec4(1.0));

    vec4 leftX1  = g0 * vec4(G(B, wrap(uv + o11 * ix * i0.x)).y, G(B, wrap(uv + o11 * ix * i0.y)).y, G(B, wrap(uv + o11 * ix * i0.z)).y, G(B, wrap(uv + o11 * ix * i0.w)).y);
    vec4 rightX1 = g1 * vec4(G(B, wrap(uv + o11 * ix * i1.x)).y, G(B, wrap(uv + o11 * ix * i1.y)).y, G(B, wrap(uv + o11 * ix * i1.z)).y, G(B, wrap(uv + o11 * ix * i1.w)).y); 
    float centerX1 = g * G(B, uv).y;
    float sumX1 = centerX1 + dot(leftX1, vec4(1.0)) + dot(rightX1, vec4(1.0));

    vec4 leftY1  = g0 * vec4(G(B, wrap(uv + o11 * iy * i0.x)).z, G(B, wrap(uv + o11 * iy * i0.y)).z, G(B, wrap(uv + o11 * iy * i0.z)).z, G(B, wrap(uv + o11 * iy * i0.w)).z);
    vec4 rightY1 = g1 * vec4(G(B, wrap(uv + o11 * iy * i1.x)).z, G(B, wrap(uv + o11 * iy * i1.y)).z, G(B, wrap(uv + o11 * iy * i1.z)).z, G(B, wrap(uv + o11 * iy * i1.w)).z); 
    float centerY1 = g * G(B, uv).z;
    float sumY1 = centerY1 + dot(leftY1, vec4(1.0)) + dot(rightY1, vec4(1.0));
    
    fragColor = vec4(sumX0, sumY0, sumX1, sumY1);
	return fragColor; 
 } 


			//******** BuffC Code Begins ********

#define G(ic,x) texture(ic, x)
//#define iC0 BuffB
//#define iC1 BuffC
#define o02 o01*(9.)
#define o12 o02*(9.)
//#define stddev 2.5
/*
float gaussian(float x, float s) {
    return exp(-x*x/(s*s));
}

vec4 gaussian(vec4 x, float s) {
    return exp(-x*x/(s*s));
}

vec2 wrap(vec2 x) {
    return mod(mod(x, 1.0) + 1.0, 1.0);
}
*/
vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    
    vec2 ix = vec2(1.0 / RENDERSIZE.x, 0.0);
    vec2 iy = vec2(0.0, 1.0 / RENDERSIZE.y);
    vec4 i0 = vec4(-4.0, -3.0, -2.0, -1.0);
    vec4 i1 = vec4(1.0, 2.0, 3.0, 4.0);
    vec4 g0 = gaussian(i0, stddev);
    vec4 g1 = g0.wzyx;
    float g = gaussian(0.0, stddev);
    float sum = 2.0 * dot(g0, vec4(1.0)) + g;
    g0 /= sum;
    g1 /= sum;
    g /= sum;

    // 2 complete blur passes
    vec4 leftX0  = g0 * vec4(G(B, wrap(uv + o02 * ix * i0.x)).w, G(B, wrap(uv + o02 * ix * i0.y)).w, G(B, wrap(uv + o02 * ix * i0.z)).w, G(B, wrap(uv + o02 * ix * i0.w)).w);
    vec4 rightX0 = g1 * vec4(G(B, wrap(uv + o02 * ix * i1.x)).w, G(B, wrap(uv + o02 * ix * i1.y)).w, G(B, wrap(uv + o02 * ix * i1.z)).w, G(B, wrap(uv + o02 * ix * i1.w)).w); 
    float centerX0 = g * G(B, uv).w;
    float sumX0 = centerX0 + dot(leftX0, vec4(1.0)) + dot(rightX0, vec4(1.0));

    vec4 leftY0  = g0 * vec4(G(C, wrap(uv + o02 * iy * i0.x)).x, G(C, wrap(uv + o02 * iy * i0.y)).x, G(C, wrap(uv + o02 * iy * i0.z)).x, G(C, wrap(uv + o02 * iy * i0.w)).x);
    vec4 rightY0 = g1 * vec4(G(C, wrap(uv + o02 * iy * i1.x)).x, G(C, wrap(uv + o02 * iy * i1.y)).x, G(C, wrap(uv + o02 * iy * i1.z)).x, G(C, wrap(uv + o02 * iy * i1.w)).x); 
    float centerY0 = g * G(C, uv).x;
    float sumY0 = centerY0 + dot(leftY0, vec4(1.0)) + dot(rightY0, vec4(1.0));

    vec4 leftX1  = g0 * vec4(G(C, wrap(uv + o12 * ix * i0.x)).y, G(C, wrap(uv + o12 * ix * i0.y)).y, G(C, wrap(uv + o12 * ix * i0.z)).y, G(C, wrap(uv + o12 * ix * i0.w)).y);
    vec4 rightX1 = g1 * vec4(G(C, wrap(uv + o12 * ix * i1.x)).y, G(C, wrap(uv + o12 * ix * i1.y)).y, G(C, wrap(uv + o12 * ix * i1.z)).y, G(C, wrap(uv + o12 * ix * i1.w)).y); 
    float centerX1 = g * G(C, uv).y;
    float sumX1 = centerX1 + dot(leftX1, vec4(1.0)) + dot(rightX1, vec4(1.0));

    vec4 leftY1  = g0 * vec4(G(C, wrap(uv + o12 * iy * i0.x)).z, G(C, wrap(uv + o12 * iy * i0.y)).z, G(C, wrap(uv + o12 * iy * i0.z)).z, G(C, wrap(uv + o12 * iy * i0.w)).z);
    vec4 rightY1 = g1 * vec4(G(C, wrap(uv + o12 * iy * i1.x)).z, G(C, wrap(uv + o12 * iy * i1.y)).z, G(C, wrap(uv + o12 * iy * i1.z)).z, G(C, wrap(uv + o12 * iy * i1.w)).z); 
    float centerY1 = g * G(C, uv).z;
    float sumY1 = centerY1 + dot(leftY1, vec4(1.0)) + dot(rightY1, vec4(1.0));
    
    fragColor = vec4(sumX0, sumY0, sumX1, sumY1);
	return fragColor; 
 } 


			//******** BuffD Code Begins ********

#define Pr 0.1
#define Pg 0.187
#define Pb 01.814
#define saturation 0.995
#define darkening 0.0075
#define BIAS(x,b) (x / ((1./b - 2.)*(1.-x))+1. )
#define rate 0.0035// +sin(smoothTime*0.001)
#define blur 0.00812

float hash( vec2 p ) {
	float h = dot(p,vec2(127.1,311.7));	
    return fract(sin(h)*43758.5453123);
}

vec4 renderPassD() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec2 texel = 1.0 / RENDERSIZE.xy;
    vec4 b0 = vec4(texture(BuffA, uv).yw, texture(BuffB, uv).yw);
    vec2 b1 = texture(BuffC, uv).yw;
    
    const float _K0 = -20.0/6.0; // center weight
    const float _K1 = 4.0/8.0; // edge-neighbors
    const float _K2 = 1.0/6.0; // vertex-neighbors
    
    // 3x3 neighborhood coordinates
    float step_x = texel.x;
    float step_y = texel.y;
    vec2 n  = vec2(0.0, step_y);
    vec2 ne = vec2(step_x, step_y);
    vec2 e  = vec2(step_x, 0.0);
    vec2 se = vec2(step_x, -step_y);
    vec2 s  = vec2(0.0, -step_y);
    vec2 sw = vec2(-step_x, -step_y);
    vec2 w  = vec2(-step_x, 0.0);
    vec2 nw = vec2(-step_x, step_y);

    vec4 is =    texture(BuffD, uv);
    vec4 is_n =  texture(BuffD, uv+n);
    vec4 is_e =  texture(BuffD, uv+e);
    vec4 is_s =  texture(BuffD, uv+s);
    vec4 is_w =  texture(BuffD, uv+w);
    vec4 is_nw = texture(BuffD, uv+nw);
    vec4 is_sw = texture(BuffD, uv+sw);
    vec4 is_ne = texture(BuffD, uv+ne);
    vec4 is_se = texture(BuffD, uv+se);

    // laplacian of all components
    vec4 lapl  = _K0*is + _K1*(is_n + is_e + is_w + is_s) + _K2*(is_nw + is_sw + is_ne + is_se)*(0.990+highhits*0.2);
    
    vec3 weights[6]; 
    weights[0] = vec3(1.0, -1.0, 1.0); 
    weights[1] = vec3(2.0, 2.0, 1.0); 
    weights[2] = vec3(3.0, -2.0, -4.0); 
    weights[3] = vec3(4.0, 3.0, 6.0); 
    weights[4] = vec3(5.0, 5.0, 3.0); 
    weights[5] = vec3(6.0, 3.0, -2.0);

    // difference of gaussians
    float dogs[6];
    dogs[0] = is.x - b0.x;
    dogs[1] = b0.x - b0.y;
    dogs[2] = b0.y - b0.z;
    dogs[3] = b0.z - b0.w;
    dogs[4] = b0.w - b1.x;
    dogs[5] = b1.x - b1.y;
    
    float lowest_variation = 1000000.0;
    vec3 diff = vec3(0.0);
    for(int i = 0; i < 5; i++) {
        float variation = abs(dogs[i]);
        if( variation < lowest_variation )
        {
            lowest_variation = variation;
            diff = sign(dogs[i]) * weights[i];
        }
    }
    
    vec4 p = vec4(sqrt(is.x*is.x*Pr + is.y*is.y*Pg + is.z*is.z*Pb));
//    vec4 desaturated = vec4(p) + (is - vec4(p)) * saturation;
      vec4 desaturated = vec4(p) + (is - vec4(p)) * saturation;
      
    vec4 eps = vec4(0.1);
    
    // initialize with noise
    if(FRAMECOUNT<=10) {
        fragColor = vec4(hash(uv));
    } else {
        
        if(distance(fragCoord.xy, RENDERSIZE.xy) < 15.0) {
            fragColor = (vec4(1.0) - eps) * is + eps;    
            fragColor -= distance(fragCoord, _mouse.xy)*fragColor;
        } else {
            fragColor = clamp(desaturated + rate *(1.0+basshits)* vec4(diff, 0.0) + blur * lapl - darkening, -1.0, 1.0);
        }
    }
    

	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord.xy / RENDERSIZE.xy;
	fragColor = 0.5 + 0.5 * texture(BuffD, uv);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderPassC();
	}
	if(PASSINDEX == 3){
		return renderPassD();
	}
	if(PASSINDEX == 4){
		return renderMainImage();
	}
}