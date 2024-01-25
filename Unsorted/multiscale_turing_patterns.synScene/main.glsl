

			//******** BuffA Code Begins ********
vec4 mediaEdges = texture(media_pass_fx, _uv);
float media_lum = sin(PI*0.75*length(mediaEdges));
float dist(vec2 p0, vec2 pf) {
    return sqrt((pf.x - p0.x) * (pf.x - p0.x) + (pf.y - p0.y) * (pf.y - p0.y));
}
float d2 = dist(RENDERSIZE.xy * 0.5, _xy.xy) * 0.004;

#define G(ic,x) texture(ic, x)
#define iC0 BuffD
#define iC1 BuffA
#define iC2 BuffC
#define iC3 BuffB
#define o0 1.0
#define o1 3.0

#define stddev 2.5+0.25*syn_MidLevel 

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
    
    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec2 ix = vec2(1.0 / RENDERSIZE.x, 0.0);
    //ix += syn_BassLevel*vec2(_noise(ix));
    vec2 iy = vec2(0.0, 1.0 / RENDERSIZE.y);
    //iy += syn_BassLevel*vec2(_noise(iy));

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
    vec4 leftX0  = g0 * vec4(G(iC0, wrap(uv + o0 * ix * i0.x)).x, G(iC0, wrap(uv + o0 * ix * i0.y)).x, G(iC0, wrap(uv + o0 * ix * i0.z)).x, G(iC0, wrap(uv + o0 * ix * i0.w)).x);
    vec4 rightX0 = g1 * vec4(G(iC0, wrap(uv + o0 * ix * i1.x)).x, G(iC0, wrap(uv + o0 * ix * i1.y)).x, G(iC0, wrap(uv + o0 * ix * i1.z)).x, G(iC0, wrap(uv + o0 * ix * i1.w)).x); 
    float centerX0 = g * G(iC0, uv).x;
    float sumX0 = centerX0 + dot(leftX0, vec4(1.0)) + dot(rightX0, vec4(1.0));

    vec4 leftY0  = g0 * vec4(G(iC1, wrap(uv + o0 * iy * i0.x)).x, G(iC1, wrap(uv + o0 * iy * i0.y)).x, G(iC1, wrap(uv + o0 * iy * i0.z)).x, G(iC1, wrap(uv + o0 * iy * i0.w)).x);
    vec4 rightY0 = g1 * vec4(G(iC1, wrap(uv + o0 * iy * i1.x)).x, G(iC1, wrap(uv + o0 * iy * i1.y)).x, G(iC1, wrap(uv + o0 * iy * i1.z)).x, G(iC1, wrap(uv + o0 * iy * i1.w)).x); 
    float centerY0 = g * G(iC1, uv).x;
    float sumY0 = centerY0 + dot(leftY0, vec4(1.0)) + dot(rightY0, vec4(1.0));

    vec4 leftX1  = g0 * vec4(G(iC1, wrap(uv + o1 * ix * i0.x)).y, G(iC1, wrap(uv + o1 * ix * i0.y)).y, G(iC1, wrap(uv + o1 * ix * i0.z)).y, G(iC1, wrap(uv + o1 * ix * i0.w)).y);
    vec4 rightX1 = g1 * vec4(G(iC1, wrap(uv + o1 * ix * i1.x)).y, G(iC1, wrap(uv + o1 * ix * i1.y)).y, G(iC1, wrap(uv + o1 * ix * i1.z)).y, G(iC1, wrap(uv + o1 * ix * i1.w)).y); 
    float centerX1 = g * G(iC1, uv).y;
    float sumX1 = centerX1 + dot(leftX1, vec4(1.0)) + dot(rightX1, vec4(1.0));

    vec4 leftY1  = g0 * vec4(G(iC1, wrap(uv + o1 * iy * i0.x)).z, G(iC1, wrap(uv + o1 * iy * i0.y)).z, G(iC1, wrap(uv + o1 * iy * i0.z)).z, G(iC1, wrap(uv + o1 * iy * i0.w)).z);
    vec4 rightY1 = g1 * vec4(G(iC1, wrap(uv + o1 * iy * i1.x)).z, G(iC1, wrap(uv + o1 * iy * i1.y)).z, G(iC1, wrap(uv + o1 * iy * i1.z)).z, G(iC1, wrap(uv + o1 * iy * i1.w)).z); 
    float centerY1 = g * G(iC1, uv).z;
    float sumY1 = centerY1 + dot(leftY1, vec4(1.0)) + dot(rightY1, vec4(1.0));
    
    fragColor = vec4(sumX0, sumY0, sumX1, sumY1);

	return fragColor; 
 

} 


			//******** BuffB Code Begins ********

#define G(ic,x) texture(ic, x)
//#define iC0 BuffA
//#define iC1 BuffB
#define o2 9.0
#define o3 27.0
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
    vec4 leftX0  = g0 * vec4(G(BuffA, wrap(uv + o2 * ix * i0.x)).w, G(BuffA, wrap(uv + o2 * ix * i0.y)).w, G(BuffA, wrap(uv + o2 * ix * i0.z)).w, G(BuffA, wrap(uv + o2 * ix * i0.w)).w);
    vec4 rightX0 = g1 * vec4(G(BuffA, wrap(uv + o2 * ix * i1.x)).w, G(BuffA, wrap(uv + o2 * ix * i1.y)).w, G(BuffA, wrap(uv + o2 * ix * i1.z)).w, G(BuffA, wrap(uv + o2 * ix * i1.w)).w); 
    float centerX0 = g * G(BuffA, uv).w;
    float sumX0 = centerX0 + dot(leftX0, vec4(1.0)) + dot(rightX0, vec4(1.0));

    vec4 leftY0  = g0 * vec4(G(BuffB, wrap(uv + o2 * iy * i0.x)).x, G(BuffB, wrap(uv + o2 * iy * i0.y)).x, G(BuffB, wrap(uv + o2 * iy * i0.z)).x, G(BuffB, wrap(uv + o2 * iy * i0.w)).x);
    vec4 rightY0 = g1 * vec4(G(BuffB, wrap(uv + o2 * iy * i1.x)).x, G(BuffB, wrap(uv + o2 * iy * i1.y)).x, G(BuffB, wrap(uv + o2 * iy * i1.z)).x, G(BuffB, wrap(uv + o2 * iy * i1.w)).x); 
    float centerY0 = g * G(BuffB, uv).x;
    float sumY0 = centerY0 + dot(leftY0, vec4(1.0)) + dot(rightY0, vec4(1.0));

    vec4 leftX1  = g0 * vec4(G(BuffB, wrap(uv + o3 * ix * i0.x)).y, G(BuffB, wrap(uv + o3 * ix * i0.y)).y, G(BuffB, wrap(uv + o3 * ix * i0.z)).y, G(BuffB, wrap(uv + o3 * ix * i0.w)).y);
    vec4 rightX1 = g1 * vec4(G(BuffB, wrap(uv + o3 * ix * i1.x)).y, G(BuffB, wrap(uv + o3 * ix * i1.y)).y, G(BuffB, wrap(uv + o3 * ix * i1.z)).y, G(BuffB, wrap(uv + o3 * ix * i1.w)).y); 
    float centerX1 = g * G(BuffB, uv).y;
    float sumX1 = centerX1 + dot(leftX1, vec4(1.0)) + dot(rightX1, vec4(1.0));

    vec4 leftY1  = g0 * vec4(G(BuffB, wrap(uv + o3 * iy * i0.x)).z, G(BuffB, wrap(uv + o3 * iy * i0.y)).z, G(BuffB, wrap(uv + o3 * iy * i0.z)).z, G(BuffB, wrap(uv + o3 * iy * i0.w)).z);
    vec4 rightY1 = g1 * vec4(G(BuffB, wrap(uv + o3 * iy * i1.x)).z, G(BuffB, wrap(uv + o3 * iy * i1.y)).z, G(BuffB, wrap(uv + o3 * iy * i1.z)).z, G(BuffB, wrap(uv + o3 * iy * i1.w)).z); 
    float centerY1 = g * G(BuffB, uv).z;
    float sumY1 = centerY1 + dot(leftY1, vec4(1.0)) + dot(rightY1, vec4(1.0));
    
    fragColor = vec4(sumX0, sumY0, sumX1, sumY1);
    fragColor += _edgeDetect(BuffB)*edgeMix*syn_MidLevel*0.25;


	return fragColor; 
 } 


			//******** BuffC Code Begins ********

#define G(ic,x) texture(ic, x)
//#define iC0 BuffB
//#define iC1 BuffC
#define o4 81.0
#define o5 243.0
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
    vec4 leftX0  = g0 * vec4(G(BuffB, wrap(uv + o4 * ix * i0.x)).w, G(BuffB, wrap(uv + o4 * ix * i0.y)).w, G(BuffB, wrap(uv + o4 * ix * i0.z)).w, G(BuffB, wrap(uv + o4 * ix * i0.w)).w);
    vec4 rightX0 = g1 * vec4(G(BuffB, wrap(uv + o4 * ix * i1.x)).w, G(BuffB, wrap(uv + o4 * ix * i1.y)).w, G(BuffB, wrap(uv + o4 * ix * i1.z)).w, G(BuffB, wrap(uv + o4 * ix * i1.w)).w); 
    float centerX0 = g * G(BuffB, uv).w;
    float sumX0 = centerX0 + dot(leftX0, vec4(1.0)) + dot(rightX0, vec4(1.0));

    vec4 leftY0  = g0 * vec4(G(BuffC, wrap(uv + o4 * iy * i0.x)).x, G(BuffC, wrap(uv + o4 * iy * i0.y)).x, G(BuffC, wrap(uv + o4 * iy * i0.z)).x, G(BuffC, wrap(uv + o4 * iy * i0.w)).x);
    vec4 rightY0 = g1 * vec4(G(BuffC, wrap(uv + o4 * iy * i1.x)).x, G(BuffC, wrap(uv + o4 * iy * i1.y)).x, G(BuffC, wrap(uv + o4 * iy * i1.z)).x, G(BuffC, wrap(uv + o4 * iy * i1.w)).x); 
    float centerY0 = g * G(BuffC, uv).x;
    float sumY0 = centerY0 + dot(leftY0, vec4(1.0)) + dot(rightY0, vec4(1.0));

    vec4 leftX1  = g0 * vec4(G(BuffC, wrap(uv + o5 * ix * i0.x)).y, G(BuffC, wrap(uv + o5 * ix * i0.y)).y, G(BuffC, wrap(uv + o5 * ix * i0.z)).y, G(BuffC, wrap(uv + o5 * ix * i0.w)).y);
    vec4 rightX1 = g1 * vec4(G(BuffC, wrap(uv + o5 * ix * i1.x)).y, G(BuffC, wrap(uv + o5 * ix * i1.y)).y, G(BuffC, wrap(uv + o5 * ix * i1.z)).y, G(BuffC, wrap(uv + o5 * ix * i1.w)).y); 
    float centerX1 = g * G(BuffC, uv).y;
    float sumX1 = centerX1 + dot(leftX1, vec4(1.0)) + dot(rightX1, vec4(1.0));

    vec4 leftY1  = g0 * vec4(G(BuffC, wrap(uv + o5 * iy * i0.x)).z, G(BuffC, wrap(uv + o5 * iy * i0.y)).z, G(BuffC, wrap(uv + o5 * iy * i0.z)).z, G(BuffC, wrap(uv + o5 * iy * i0.w)).z);
    vec4 rightY1 = g1 * vec4(G(BuffC, wrap(uv + o5 * iy * i1.x)).z, G(BuffC, wrap(uv + o5 * iy * i1.y)).z, G(BuffC, wrap(uv + o5 * iy * i1.z)).z, G(BuffC, wrap(uv + o5 * iy * i1.w)).z); 
    float centerY1 = g * G(BuffC, uv).z;
    float sumY1 = centerY1 + dot(leftY1, vec4(1.0)) + dot(rightY1, vec4(1.0));
    
    fragColor = vec4(sumX0, sumY0, sumX1, sumY1)*(1.+0.0125*syn_Level);
    //fragColor *=1.0 - _edgeDetectSobel(syn_UserImage, _uv)*0.5*syn_BassLevel;
    //fragColor += _hBlurTiltShift(BuffC)*blurMix*syn_MidLevel*0.75;
    //fragColor += _edgeDetect(BuffC)*edgeMix*syn_MidLevel*0.5;

	return fragColor; 
 } 


			//******** BuffD Code Begins ********

#define Pr 0.299  //0.299
#define Pg 0.587 //.587
#define Pb 0.114 //0.114
#define saturation 0.9925
#define darkening 0.0099 //0.01
#define rate 0.005*(1.0 + 0.25*syn_BassLevel)// 0.005
#define blur 0.03*(1.0+blurMix*0.4) //0.02

float hash( vec2 p ) {
	float h = dot(p,vec2(127.1,311.7));	
    return fract(sin(h)*43758.5453123);
}

vec4 renderPassD() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    float a = .003*sin(.125*smoothTimeC-2.*length(fragCoord-0.5*RENDERSIZE.xy)/RENDERSIZE.y)*Spin;
    fragCoord *= mat2(cos(a),-sin(a),sin(a),cos(a));
    fragCoord += d2 * fisheye * + _uvc;
    fragCoord -= _uvc*(Zoom*0.875*(1+syn_BassLevel*1.25+syn_MidLevel*0.75));
    fragCoord += _uvc*(Stretch*0.75*(1+syn_BassLevel+syn_MidLevel*0.5));
    fragCoord += (moveXY*(1+syn_BassLevel+syn_MidLevel*0.5)); \
    fragCoord += vec2(_noise(sin(TIME*0.3+2)), _noise(cos(TIME*0.3-1)));
    fragCoord -= _noise(_uv+_uvc)*(syn_BassPresence*0.35+0.35*syn_Intensity+syn_BassLevel*0.3)*Warp;
  //  fragCoord += 
    fragCoord += vec2(dot(cos(_uv.x+_uvc.x),sin(_uv.y+_uvc.y)))*Pull;
    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    vec2 texel = 1.0 / RENDERSIZE.xy;
    
    vec4 b0 = vec4(texture(BuffA, uv).yw+Noise*_noise(cos(syn_Presence*0.5+syn_BassLevel*0.5)), texture(BuffB, uv).yw+Noise*_noise(sin(syn_Presence*0.5+syn_BassLevel*0.5)));
   // b0 += 10* _luminance(_textureMedia(uv));

    vec2 b1 = texture(BuffC, uv).yw;

// b1 += vec2(_noise(b1));
    const float _K0 = -20.0/6.0; // center weight
    const float _K1 = 4.0/6.0; // edge-neighbors
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
    vec4 lapl  = _K0*is + _K1*(is_n + is_e + is_w + is_s) + _K2*(is_nw + is_sw + is_ne + is_se);
    lapl*=1.0 + syn_MidLevel;


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
    
    float lowest_variation = 10000.0;
    vec3 diff = vec3(0.0);
    for(int i = 0; i < 6; i++) {
        float variation = abs(dogs[i]);
        if( variation < lowest_variation )
        {
            lowest_variation = variation;
            diff = sign(dogs[i]) * weights[i];

        }
    }
    
    vec4 p = vec4(sqrt(is.x*is.x*Pr + is.y*is.y*Pg + is.z*is.z*Pb));
    p += mediaEdges;
    vec4 desaturated = vec4(p) + (is - vec4(p)) * saturation;
    
    vec4 eps = vec4(0.1);
    
    // initialize with noise
    if(FRAMECOUNT<10) {
        fragColor = vec4(hash(uv));
    } else {
        if(distance(fragCoord.xy, _mouse.xy) < 40.0*paint_size && _mouse.z > 0) {
            fragColor = (vec4(1.0) - eps) * is + eps;    
        } else {
            fragColor = clamp(desaturated + rate * vec4(diff, 0.0) + blur * lapl - darkening, -1.0, 1.0);

          //  fragColor *=1.0 - _edgeDetectSobel(syn_UserImage, _uv)*0.085*syn_MidLevel;

        }
    }
    fragColor += length(_edgeDetect(BuffD)*edgeMix*syn_MidLevel*0.03);

	return fragColor; 
 } 


vec4 mediaPass() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec4 media = _loadMedia();
    return media*media_impact;
    //fragColor *=1.0 - _edgeDetectSobel(syn_UserImage, _uv)*0.5*syn_MidLevel;

	return fragColor; 
 } 
vec4 mediaPassFX() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec4 edges = mix(mediaPass(), _edgeDetectSobel(media_pass), edge_mix);
    return edges;
    //fragColor *=1.0 - _edgeDetectSobel(syn_UserImage, _uv)*0.5*syn_MidLevel;

	return fragColor; 
 } 
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    fragCoord = mix(fragCoord, abs(fragCoord), vec2(mirror_x, mirror_y));

	vec2 uv = fragCoord.xy / RENDERSIZE.xy;
	fragColor = 0.5 + 0.5 * texture(BuffD, uv);
    //fragColor *=1.0 - _edgeDetectSobel(syn_UserImage, _uv)*0.5*syn_MidLevel;
    fragColor.xy = abs(fragColor.xy);
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
		return mediaPass();
	}
	if(PASSINDEX == 5){
		return mediaPassFX();
	}
	if(PASSINDEX == 6){
		return renderMainImage();
	}
}