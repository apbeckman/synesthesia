

			//******** Common Code Begins ********

// For normalized fixed-point packing
//#define unpack(x) unpackSnorm2x16(floatBitsToUint(x))
//#define pack(x) uintBitsToFloat(packSnorm2x16(x))
#extension GL_ARB_shading_language_packing : enable     
#define pack(d) uintBitsToFloat(packHalf2x16(d))
#define unpack(d) unpackHalf2x16(floatBitsToUint(d))
float growthFactor = pow((syn_BassLevel*0.5)+(syn_MidLevel*0.25)+(syn_Level*0.125), 2.0);


// contrast
#define SIGMOID_CONTRAST 20

vec3 contrast(vec3 x) {
	return 1.0 / (1.0 + exp(-SIGMOID_CONTRAST * (x - 0.5)));    
}

vec3 normz(vec3 x) {
	return x == vec3(0) ? vec3(0) : normalize(x);
}

// Begin IQ's simplex noise:

// The MIT License
// Copyright © 2013 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

vec2 hash( vec2 p ) // replace this by something better
{
	p = vec2( dot(p,vec2(127.1,311.7)),
			  dot(p,vec2(269.5,183.3)) );

	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}

float noise( in vec2 p )
{
    const float K1 = 0.366025404; // (sqrt(3)-1)/2;
    const float K2 = 0.211324865; // (3-sqrt(3))/6;

	vec2 i = floor( p + (p.x+p.y)*K1 );
	
    vec2 a = p - i + (i.x+i.y)*K2;
    vec2 o = step(a.yx,a.xy);    
    vec2 b = a - o + K2;
	vec2 c = a - 1.0 + 2.0*K2;

    vec3 h = max( 0.5-vec3(dot(a,a), dot(b,b), dot(c,c) ), 0.0 );

	vec3 n = h*h*h*h*vec3( dot(a,hash(i+0.0)), dot(b,hash(i+o)), dot(c,hash(i+1.0)));

    return dot( n, vec3(70.0) );
	
}

// End IQ's simplex noise


// GGX from noby's Goo shader https://www.shadertoy.com/view/lllBDM
float G1V(float dnv, float k){
    return 1.0/(dnv*(1.0-k)+k);
}

float ggx(vec3 n, vec3 v, vec3 l, float rough, float f0){
    float alpha = rough*rough;
    vec3 h = normalize(v+l);
    float dnl = clamp(dot(n,l), 0.0, 1.0);
    float dnv = clamp(dot(n,v), 0.0, 1.0);
    float dnh = clamp(dot(n,h), 0.0, 1.0);
    float dlh = clamp(dot(l,h), 0.0, 1.0);
    float f, d, vis;
    float asqr = alpha*alpha;
    const float pi = 3.14159;
    float den = dnh*dnh*(asqr-1.0)+1.0;
    d = asqr/(pi * den * den);
    dlh = pow(1.0-dlh, 5.0);
    f = f0 + (1.0-f0)*dlh;
    float k = alpha/1.0;
    vis = G1V(dnl, k)*G1V(dnv, k);
    float spec = dnl * d * f * vis;
    return spec;
}

			//******** BuffA Code Begins ********

#define STEPS 50  // advection steps

#define ts 0.25*Amp    // advection curl
#define cs -2.0   // curl scale
#define ls 0.05*Amp   // laplacian scale
#define ps -2.0*pow(Divergence, 01.25)   // laplacian of divergence scale
#define ds -0.4/(Divergence)   // divergence scale
#define dp -0.03  // divergence update scale
#define pl 0.3    // divergence smoothing
#define amp 1.0*(0.999+growthFactor*0.001)   // self-amplification
#define upd 0.4   // update smoothing

#define _D 0.6    // diagonal weight

#define _K0 -20.0/6.0 // laplacian center weight
#define _K1 4.0/6.0   // laplacian edge-neighbors
#define _K2 1.0/6.0   // laplacian vertex-neighbors

#define _G0 0.25*(weight)      // gaussian center weight
#define _G1 0.125     // gaussian edge-neighbors
#define _G2 0.0625    // gaussian vertex-neighbors

bool reset() {
    return 0.0 > 0.5;
}

vec2 normz(vec2 x) {
	return x == vec2(0.0) ? vec2(0.0) : normalize(x);
}

#define T(d) texture(BuffA, fract(aUv+d)).xyz

vec3 advect(vec2 ab, vec2 vUv, vec2 texel, out float curl, out float div, out vec3 lapl, out vec3 blur) {
    
    vec2 aUv = vUv - ab * texel;
    vec4 t = vec4(texel, -texel.y, 0.0);

    vec3 uv =    T( t.ww); vec3 uv_n =  T( t.wy); vec3 uv_e =  T( t.xw);
    vec3 uv_s =  T( t.wz); vec3 uv_w =  T(-t.xw); vec3 uv_nw = T(-t.xz);
    vec3 uv_sw = T(-t.xy); vec3 uv_ne = T( t.xy); vec3 uv_se = T( t.xz);
    
    curl = uv_n.x - uv_s.x - uv_e.y + uv_w.y + _D * (uv_nw.x + uv_nw.y + uv_ne.x - uv_ne.y + uv_sw.y - uv_sw.x - uv_se.y - uv_se.x);
    div  = uv_s.y - uv_n.y - uv_e.x + uv_w.x + _D * (uv_nw.x - uv_nw.y - uv_ne.x - uv_ne.y + uv_sw.x + uv_sw.y + uv_se.y - uv_se.x);
    lapl = _K0*uv + _K1*(uv_n + uv_e + uv_w + uv_s) + _K2*(uv_nw + uv_sw + uv_ne + uv_se);
    blur = _G0*uv + _G1*(uv_n + uv_e + uv_w + uv_s) + _G2*(uv_nw + uv_sw + uv_ne + uv_se);
    
    return uv;
}

vec2 rot(vec2 v, float th) {
	return vec2(dot(v, vec2(cos(th), -sin(th))), dot(v, vec2(sin(th), cos(th)))); 
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    vec2 vUv = fragCoord.xy / RENDERSIZE.xy;
    vec2 texel = 1. / RENDERSIZE.xy;
    
    vec3 lapl, blur;
    float curl, div;
    
    vec3 uv = advect(vec2(0), vUv, texel, curl, div, lapl, blur)*(0.9995+syn_MidLevel*0.00025+syn_BassLevel*0.00025);

    float sp = ps * lapl.z;
    float sc = (cs * curl)*(0.99925+syn_BassLevel*0.0175);
	float sd = uv.z + dp * div + pl * lapl.z;
    vec2 norm = normz(uv.xy);

    vec2 off = uv.xy;
    vec2 offd = off;
    vec3 ab = vec3(0);

    for(int i = 0; i < STEPS; i++) {
        //advect(off, vUv, texel, curl, div, lapl, blur);
        vec4 samp = texture(BuffB, fract(vUv - off * texel));
        curl = samp.w;
        blur = samp.xyz;
        offd = rot(offd,ts*curl);
        off += offd;
    	ab += blur / float(STEPS);  
    }
    
    vec2 tab = (amp * ab.xy*(0.9925+growthFactor*0.0725) + (ls * lapl.xy + norm * sp)*0.99999+syn_Level*0.00001 + (uv.xy * ds * sd)*(0.999125+syn_MidLevel*0.000895));    
    vec2 rab = rot(tab,sc)*impulse;
    
    vec3 abd = mix(vec3(rab,sd), uv, upd+syn_Level*0.025);
    
    if (_mouse.z > 0.0) {
    	vec2 d = (fragCoord.xy - _mouse.xy) / RENDERSIZE.xy;
        vec2 m = (Size*0.5) * normz(d) * exp(-length(d) / (0.02/(Size)));
       //m *= Size;
        abd.xy += m;
        uv.xy += m;
    }
    
    // initialize with noise
    if(uv == vec3(0) || reset()) {
        vec3 rnd = vec3(noise(8.0 * vUv + 1.1), noise(8.0 * vUv + 2.2), noise(8.0 * vUv + 3.3));
        fragColor = vec4(rnd, 0);
    } else {
        abd.z = clamp(abd.z, -1.0, 1.0);
        abd.xy = clamp(length(abd.xy) > 1.0 ? normz(abd.xy) : abd.xy, -1.0, 1.0)*(0.99+normalize(growthFactor)*0.01);
        fragColor = vec4(abd, 0.0);
    }

	return fragColor; 
 } 


			//******** BuffB Code Begins ********

// This computes the laplacian of the input
/*
#define _G0 0.25      // gaussian center weight
#define _G1 0.125     // gaussian edge-neighbors
#define _G2 0.0625    // gaussian vertex-neighbors
*/
#define _D 0.6    // diagonal weight

#define TB(d) texture(BuffA, fract(vUv+d)).xyz

vec4 curl_gaussian(vec2 fragCoord, vec2 RENDERSIZE) {
    vec2 texel = 1.0 / RENDERSIZE;
    vec2 vUv = fragCoord * texel;
    vec4 t = vec4(texel, -texel.y, 0.0);

    vec3 uv =    TB( t.ww); vec3 uv_n =  TB( t.wy); vec3 uv_e =  TB( t.xw);
    vec3 uv_s =  TB( t.wz); vec3 uv_w =  TB(-t.xw); vec3 uv_nw = TB(-t.xz);
    vec3 uv_sw = TB(-t.xy); vec3 uv_ne = TB( t.xy); vec3 uv_se = TB( t.xz);
    
    float curl = uv_n.x - uv_s.x - uv_e.y + uv_w.y + _D * (uv_nw.x + uv_nw.y + uv_ne.x - uv_ne.y + uv_sw.y - uv_sw.x - uv_se.y - uv_se.x);
    vec3 gaussian = _G0*uv + _G1*(uv_n + uv_e + uv_w + uv_s) + _G2*(uv_nw + uv_sw + uv_ne + uv_se);
    
    return vec4(gaussian, curl);
}

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    fragColor = curl_gaussian(fragCoord, RENDERSIZE.xy);
	return fragColor; 
 } 


			//******** BuffC Code Begins ********

bool resetC() {
    return 0.0 > 0.5;
}

vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord / RENDERSIZE.xy;
    vec2 texel = 1.0 / RENDERSIZE.xy;
    
/*
U:
 1.9000778e-006   2.6484959e-003  -1.0391083e-001
-2.0540590e-003  -4.4125709e-003  -3.4489894e-001
-1.3499238e-002  -6.8390049e-002  -4.8179728e-001
-5.1257182e-002  -2.5511910e-001  -3.1508410e-001
-1.5871959e-001  -5.5398879e-001   1.1805352e-001
-4.7194022e-001  -1.2936001e-001   1.1920299e-001
-7.0606907e-001   4.6167731e-001  -1.4625093e-001
-4.7194022e-001  -1.2936001e-001   1.1920299e-001
-1.5871959e-001  -5.5398879e-001   1.1805352e-001
-5.1257182e-002  -2.5511910e-001  -3.1508410e-001
-1.3499238e-002  -6.8390049e-002  -4.8179728e-001
-2.0540590e-003  -4.4125709e-003  -3.4489894e-001
 1.9000778e-006   2.6484959e-003  -1.0391083e-001

V:
 6.2750203e-003  5.2398670e-002  3.7276962e-001  
-1.6414278e-003  4.2486224e-002  5.4995743e-001
-4.3498466e-002 -1.0892533e-001  2.4023362e-001
-1.3135171e-001 -3.3953004e-001 -7.8265086e-004
-3.0484343e-001 -5.0984393e-001  1.8311873e-002
-6.2280256e-001  3.2918550e-001 -2.3270335e-002
-5.8224388e-076  3.2916004e-064 -1.0109051e-055
 6.2280256e-001 -3.2918550e-001  2.3270335e-002
 3.0484343e-001  5.0984393e-001 -1.8311873e-002
 1.3135171e-001  3.3953004e-001  7.8265088e-004  
 4.3498466e-002  1.0892533e-001 -2.4023362e-001
 1.6414278e-003 -4.2486224e-002 -5.4995743e-001
-6.2750203e-003 -5.2398670e-002 -3.7276962e-001

diag(S):

  5.2045614e-001
  4.5787111e-002
  5.3607463e-003
  1.3379961e-003
  2.4305267e-005
  6.5520767e-008
  2.7241624e-013
  1.8098574e-013
  4.2556708e-014
  3.2104951e-014
  5.6655005e-016
  2.1958056e-018
  1.0637097e-030

*/

    float p_y3[13] = float[](-1.0391083e-001, -3.4489894e-001, -4.8179728e-001, -3.1508410e-001,  1.1805352e-001,  1.1920299e-001, -1.4625093e-001,  1.1920299e-001,  1.1805352e-001, -3.1508410e-001, -4.8179728e-001, -3.4489894e-001, -1.0391083e-001);
 	float p_y2[13] = float[](2.6484959e-003, -4.4125709e-003, -6.8390049e-002, -2.5511910e-001, -5.5398879e-001, -1.2936001e-001, 4.6167731e-001, -1.2936001e-001, -5.5398879e-001, -2.5511910e-001, -6.8390049e-002, -4.4125709e-003, 2.6484959e-003);
    float p_y1[13] = float[](1.9000778e-006, -2.0540590e-003, -1.3499238e-002, -5.1257182e-002, -1.5871959e-001, -4.7194022e-001, -7.0606907e-001, -4.7194022e-001, -1.5871959e-001, -5.1257182e-002, -1.3499238e-002, -2.0540590e-003,  1.9000778e-006);

    float p_x3[13] = float[](3.7276962e-001,  5.4995743e-001,  2.4023362e-001, -7.8265086e-004,  1.8311873e-002, -2.3270335e-002, -1.0109051e-055,  2.3270335e-002, -1.8311873e-002,  7.8265088e-004, -2.4023362e-001, -5.4995743e-001, -3.7276962e-001);
    float p_x2[13] = float[](5.2398670e-002,  4.2486224e-002, -1.0892533e-001, -3.3953004e-001, -5.0984393e-001,  3.2918550e-001,  0.0, -3.2918550e-001,  5.0984393e-001,  3.3953004e-001,  1.0892533e-001, -4.2486224e-002, -5.2398670e-002);
    float p_x1[13] = float[](6.2750203e-003, -1.6414278e-003, -4.3498466e-002, -1.3135171e-001, -3.0484343e-001, -6.2280256e-001, 0.0, 6.2280256e-001, 3.0484343e-001, 1.3135171e-001, 4.3498466e-002, 1.6414278e-003, -6.2750203e-003);
        
    float s_i[3] = float[](  5.2045614e-001, 4.5787111e-002, 5.3607463e-003);
    
    float g_x[13] = float[](1.8154960e-002, 5.1439053e-002, 1.1757498e-001, 2.2045309e-001, 3.4292702e-001, 4.4580513e-001, 
         4.8633287e-001, 4.4580513e-001, 3.4292702e-001, 2.2045309e-001, 1.1757498e-001, 5.1439053e-002, 1.8154960e-002);  

    #define RANGE 6
    
    #define Po(m,n) texture(BuffA, fract(uv + texel * vec2(m,n)))
    
    vec2 P1 = vec2(0);
    vec2 P2 = vec2(0);
    vec2 P3 = vec2(0);
    float G = 0.0;
    float Gw = 0.0;
    for (int i = -RANGE; i <= RANGE; i++) {
        int index = RANGE + i;

        vec2 t = Po(i,0).xy;
        float g = texture(BuffD, fract(uv + texel * vec2(i,0))).x;
        
        P1 += vec2(p_x1[index], p_y1[index]) * t;
        P2 += vec2(p_x2[index], p_y2[index]) * t;
        P3 += vec2(p_x3[index], p_y3[index]) * t;
        
        Gw += g_x[index];
        G  += g_x[index] * g;
    }
    
    G /= Gw;
    
    if(reset()) {
        fragColor = vec4(0);
    } else {
        fragColor = vec4(pack(P1),pack(P2),pack(P3), G);
    }

	return fragColor; 
 } 


			//******** BuffD Code Begins ********

bool resetD() {
    return 0.0 > 0.5;
}

vec4 renderPassD() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord / RENDERSIZE.xy;
    vec2 texel = 1.0 / RENDERSIZE.xy;
    
    float p_y3[13] = float[](-1.0391083e-001, -3.4489894e-001, -4.8179728e-001, -3.1508410e-001,  1.1805352e-001,  1.1920299e-001, -1.4625093e-001,  1.1920299e-001,  1.1805352e-001, -3.1508410e-001, -4.8179728e-001, -3.4489894e-001, -1.0391083e-001);
 	float p_y2[13] = float[](2.6484959e-003, -4.4125709e-003, -6.8390049e-002, -2.5511910e-001, -5.5398879e-001, -1.2936001e-001, 4.6167731e-001, -1.2936001e-001, -5.5398879e-001, -2.5511910e-001, -6.8390049e-002, -4.4125709e-003, 2.6484959e-003);
    float p_y1[13] = float[](1.9000778e-006, -2.0540590e-003, -1.3499238e-002, -5.1257182e-002, -1.5871959e-001, -4.7194022e-001, -7.0606907e-001, -4.7194022e-001, -1.5871959e-001, -5.1257182e-002, -1.3499238e-002, -2.0540590e-003,  1.9000778e-006);

    float p_x3[13] = float[](3.7276962e-001,  5.4995743e-001,  2.4023362e-001, -7.8265086e-004,  1.8311873e-002, -2.3270335e-002, -1.0109051e-055,  2.3270335e-002, -1.8311873e-002,  7.8265088e-004, -2.4023362e-001, -5.4995743e-001, -3.7276962e-001);
    float p_x2[13] = float[](5.2398670e-002,  4.2486224e-002, -1.0892533e-001, -3.3953004e-001, -5.0984393e-001,  3.2918550e-001,  0.0, -3.2918550e-001,  5.0984393e-001,  3.3953004e-001,  1.0892533e-001, -4.2486224e-002, -5.2398670e-002);
    float p_x1[13] = float[](6.2750203e-003, -1.6414278e-003, -4.3498466e-002, -1.3135171e-001, -3.0484343e-001, -6.2280256e-001, 0.0, 6.2280256e-001, 3.0484343e-001, 1.3135171e-001, 4.3498466e-002, 1.6414278e-003, -6.2750203e-003);
        
    float s_i[3] = float[](  5.2045614e-001, 4.5787111e-002, 5.3607463e-003);
    
    float g_x[13] = float[](1.8154960e-002, 5.1439053e-002, 1.1757498e-001, 2.2045309e-001, 3.4292702e-001, 4.4580513e-001, 
         4.8633287e-001, 4.4580513e-001, 3.4292702e-001, 2.2045309e-001, 1.1757498e-001, 5.1439053e-002, 1.8154960e-002);  

    //#define RANGE 6
    
    #define Pc(m,n) texture(BuffC, fract(uv + texel * vec2(m,n)))
    
    vec2 P = vec2(0);
    float G = 0.0;
    float Gw = 0.0;
    for (int i = -RANGE; i <= RANGE; i++) {
        int index = RANGE + i;
        
        vec4 tx = Pc(0,i);
        vec2 t1 = unpack(tx.x);
        vec2 t2 = unpack(tx.y);
        vec2 t3 = unpack(tx.z);

        float g = tx.w;
        
        P += s_i[0] * vec2(p_x1[index], p_y1[index]).yx * t1;
        P += s_i[1] * vec2(p_x2[index], p_y2[index]).yx * t2;
        P += s_i[2] * vec2(p_x3[index], p_y3[index]).yx * t3;
        Gw += g_x[index];
        G  += g_x[index] * g;
    }
    
    G /= Gw;

    if(resetD()) {
        fragColor = vec4(0);
    } else {
        fragColor = vec4(-(P.x + P.y) + G);
    }

	return fragColor; 
 } 


/* 
	Created by Cornus Ammonis (2021)
	Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

    Fork of "Perceptual Depth Poisson Solver" by cornusammonis. https://shadertoy.com/view/MlByW3
    2020-10-10 22:58:33

	This is an update of my previously-posted "Perceptual Depth Poisson Solver" shader.
	As in the previous technique, a Poisson solver kernel is generated by recursively 
    solving multiple Jacobi iterations. See this shader (https://www.shadertoy.com/view/MdSczK)
    by Robert Schuetze for another example of a multi-step kernel Poisson solver method. 
    
    This approach sacrifices some accuracy for greatly improved performance by precomputing
    separable Poisson solver kernels using singular value decomposition.
    
    Here, the Laplacian calculation step and Poisson solver step are combined together by 
    convolution, resulting in a single step with two kernels, one each for the x and y 
    components of the input vector field. The combined kernels are then decomposed into separable 
    kernels using singular value decomposition. Here, the 3 largest singular values are taken, 
    and the sum of the 3 resulting separable convolutions is used to approximate the 2 full 
    Poisson kernels. 
    
    The result of the first pass (Buffer C) is packed to half floats in order to store each of the 6
    values (3 each for x and y components) from the first convolution pass in a single buffer.
    
    I have also made some general performance improvements to the original shader by offloading
    some computation in the dynamical system in Buffer A to a second pass in Buffer B. These 
    changes are ancillary to the Poisson-SVD method demonstrated here.

	Comment out "#define POISSON" below to render using the original vector map without using the
    Poisson solver.
*/

// displacement (for texturing)
#define DISP 0.02

// bump mapping scale
#define BUMP 1.5*(1.0+basshits*0.5)

// mip level
#define MIP 0.0

// comment to use the original vector field without running through the Poisson solver
#define POISSON

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 texel = 1. / RENDERSIZE.xy;
    vec2 uv = fragCoord.xy / RENDERSIZE.xy;

    vec2 n  = vec2(0.0, texel.y);
    vec2 e  = vec2(texel.x, 0.0);
    vec2 s  = vec2(0.0, -texel.y);
    vec2 w  = vec2(-texel.x, 0.0);

    #ifdef POISSON
        float d   = texture(BuffD, uv).x;
        float d_n = texture(BuffD, fract(uv+n)).x;
        float d_e = texture(BuffD, fract(uv+e)).x;
        float d_s = texture(BuffD, fract(uv+s)).x;
        float d_w = texture(BuffD, fract(uv+w)).x; 

        float d_ne = texture(BuffD, fract(uv+n+e)).x;
        float d_se = texture(BuffD, fract(uv+s+e)).x;
        float d_sw = texture(BuffD, fract(uv+s+w)).x;
        float d_nw = texture(BuffD, fract(uv+n+w)).x; 

        float dxn[3];
        float dyn[3];

        dyn[0] = d_nw - d_sw;
        dyn[1] = d_n  - d_s; 
        dyn[2] = d_ne - d_se;

        dxn[0] = d_ne - d_nw; 
        dxn[1] = d_e  - d_w; 
        dxn[2] = d_se - d_sw; 
    #else
        vec2 d   = texture(BuffA, uv).xy;
        vec2 d_n = texture(BuffA, fract(uv+n)).xy;
        vec2 d_e = texture(BuffA, fract(uv+e)).xy;
        vec2 d_s = texture(BuffA, fract(uv+s)).xy;
        vec2 d_w = texture(BuffA, fract(uv+w)).xy; 

        vec2 d_ne = texture(BuffA, fract(uv+n+e)).xy;
        vec2 d_se = texture(BuffA, fract(uv+s+e)).xy;
        vec2 d_sw = texture(BuffA, fract(uv+s+w)).xy;
        vec2 d_nw = texture(BuffA, fract(uv+n+w)).xy; 

        float dxn[3];
        float dyn[3];

        dyn[0] = d_n.y;
        dyn[1] = d.y; 
        dyn[2] = d_s.y;

        dxn[0] = d_e.x; 
        dxn[1] = d.x; 
        dxn[2] = d_w.x; 
    #endif
    
    #define I(d_x,d_y) texture(image3, fract(vec2(0.5) + DISP * vec2(d_x,d_y)), MIP).xyz

    vec3 i   = I(dxn[0],dyn[0]);
    vec3 i_n = I(dxn[1],dyn[1]);
    vec3 i_e = I(dxn[2],dyn[2]);
    vec3 i_s = I(dxn[1],dyn[2]);
    vec3 i_w = I(dxn[2],dyn[0]);
    
    vec3 ib = 0.4 * i + 0.15 * (i_n+i_e+i_s+i_w);

    vec3 ld = normz(vec3(0.5+0.5*vec2(cos(smoothTimeB/4.0), sin(smoothTimeB/4.0)) - uv, -1.));
    
    float spec = 0.0;    
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            vec3 dxy = normalize(vec3(-BUMP*vec2(dxn[i], dyn[j]), -1.0));
            spec += ggx(dxy, vec3(0,0,-1), ld, 0.4, 0.1) / 9.0;
        }
    }

    // end bumpmapping section

    vec3 tc = 0.9*contrast(0.9*ib);

    fragColor = vec4((tc + vec3(0.9, 0.85, 0.8)*spec),1.0);
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