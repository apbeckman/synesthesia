vec3 colorPaletteChooser(float colReg, float var){
  vec3 paletteCol = vec3(0.5);
  if (colReg < 1.0){
    paletteCol = _palette(var, vec3(0.500, 0.500, 0.520), vec3(0.500, 0.500, 0.500), vec3(0.780, 0.765, 0.750), vec3(0.360, 0.570, 0.680));
  } else if (colReg < 2.0){
    paletteCol = _palette(var, vec3(1.025, 0.361, 0.703), vec3(0.900, 0.506, 0.724), vec3(0.720, 1.005, 0.075), vec3(1.000, 0.950, 0.590));
  } else if (colReg < 3.0){
    paletteCol = _palette(1.0-var,vec3(0.540, -0.500, 0.340), vec3(0.630, 0.816, 0.680), vec3(1.410, 0.855, 0.570), vec3(0.370, 0.620, 0.630));
  } else if (colReg < 4.0){
    paletteCol = _palette(var,vec3(0.500, 0.500, 0.500), vec3(0.500, 0.500, 0.500), vec3(0.500, 0.825, 0.750), vec3(0.500, 0.500, 0.500));
  } else if (colReg < 5.0){
    paletteCol = _palette(1.0-var,vec3(0.000, 0.580, 0.453), vec3(0.848, 0.703, 0.110), vec3(0.700, 0.175, 0.542), vec3(0.000, 0.182, 0.915));
  } else if (colReg < 6.0){
    paletteCol = _palette(var,vec3(0.411, -0.070, 0.500), vec3(0.350, 0.880, 0.471), vec3(1.398, 1.363, 0.678), vec3(0.652, 0.250, 0.809)).gbr;
  } else if (colReg < 7.0){
    paletteCol = _palette(1.0-var,vec3(1.040, 0.180, 0.260), vec3(0.053, 0.775, 0.330), vec3(0.142, 0.523, 0.800), vec3(0.242, 0.887, 0.000));
  }
  return paletteCol;
}

			//******** BuffA Code Begins ********

float sHighHits = syn_HighHits*audio_reactivity;
float sBassHits = syn_BassHits*audio_reactivity;

// ---------------------------------------------
float or = 9.0*scale*scale;         // outer gaussian std dev
float ir = 3.0*scale*scale;          // inner gaussian std dev
float b1 = 0.19*birth_lower-        0.18*regime_toggle  ;                // birth1
float b2 = 0.222*birth_upper+       0.15*regime_toggle*(1.0-sBassHits*3.5);// birth2
float s1 = 0.267*survival_lower-                (1.0-sBassHits);                // survival1
float s2 = 0.445*survival_upper+    regime_toggle*0.2+  sBassHits*0.2;                  // survival2
float dt = (0.2*timestep+sHighHits*0.3);          // timestep
float alpha_n = 0.017*flow_resistance;   // sigmoid width for outer fullness
float alpha_m = 0.112*space_filling;   // sigmoid width for inner fullness
// ---------------------------------------------
// float or = 18.0; 
// float ir = 6.0; 
// float b1 = 0.2; 
// float b2 = 0.215; 
// float s1 = 0.25; 
// float s2 = 0.5; 
// float dt = 0.2; 
// float alpha_n = 0.02; 
// float alpha_m = 0.11;

bool reset() {
    return ((FRAMECOUNT < 2)||(reset_sim>0.5));
}

// the logistic function is used as a smooth step function
float sigma1(float x,float a,float alpha) 
{ 
    return 1.0 / ( 1.0 + exp( -(x-a)*4.0/alpha ) );
}

float sigma2(float x,float a,float b,float alpha)
{
    return sigma1(x,a,alpha) 
        * ( 1.0-sigma1(x,b,alpha) );
}

float sigma_m(float x,float y,float m,float alpha)
{
    return x * ( 1.0-sigma1(m,0.5,alpha) ) 
        + y * sigma1(m,0.5,alpha);
}

// the transition function
// (n = outer fullness, m = inner fullness)
float s(float n,float m)
{
    return sigma2( n, sigma_m(b1,s1,m,alpha_m), 
        sigma_m(b2,s2,m,alpha_m), alpha_n );
}

#define T(d) texture(BuffA, fract(uv+d)).x

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 tx = 1.0 / RENDERSIZE.xy;
    vec2 uv = fragCoord.xy * tx;
    const float _K0 = -20.0/6.0; // center weight
    const float _K1 = 4.0/6.0;   // edge-neighbors
    const float _K2 = 1.0/6.0;   // vertex-neighbors
    
    
	// We can optionally add a laplacian to the update rule
	// to decrease the appearance of aliasing, but we also
	// introduce subtle anisotropy by doing so.
        // vec4 t = vec4(tx, -tx.y, 0.0);
        // float u =    T( t.ww); float u_n =  T( t.wy); float u_e =  T( t.xw);
        // float u_s =  T( t.wz); float u_w =  T(-t.xw); float u_nw = T(-t.xz);
        // float u_sw = T(-t.xy); float u_ne = T( t.xy); float u_se = T( t.xz);
        // float lapl  = _K0*u + _K1*(u_n + u_e + u_w + u_s) + _K2*(u_nw + u_sw + u_ne + u_se);
       

    vec4 current = texture(BuffA, uv);
    vec2 fullness = texture(BuffC, uv).xy;
    
    float delta =  2.0 * s( fullness.x, fullness.y ) - 1.0;
    float new = clamp( current.x + dt * delta, 0.0, 1.0 );
    
    // if(paint_on > 0.0) {
    //     // from chronos' SmoothLife shader https://www.shadertoy.com/view/XtdSDn
    //     float dst = length(fragCoord.xy - paint_xy*RENDERSIZE);
    //     if(dst <= or) {
    //     	new = step((ir+1.5), dst) * (1.0 - step(or, dst));
    //     }
    // }
    
    vec4 init = texture(image47, uv*1.4);
    if(FRAMECOUNT < 10 || reset() || (init != vec4(0) && current.w == 0.0)) {
    	fragColor = vec4(2.0*init.xyz,1.0);    
    } else {
    	fragColor = vec4(new, fullness, current.w);
    }

    fragColor -= cycling_rings*fragColor*sin(length(_uvc)*10.0+TIME*2.0)*0.1;

    float lum = dot(_loadUserImage().rgb,vec3(1.0))/3.0;
    fragColor.r = max(fragColor.r, pow(lum,contrast*contrast));
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

#define SQRT_2_PI 2.50662827463

// ---------------------------------------------
// const float or2 = 18.0;         // outer gaussian std dev
// const float ir2 = 6.0;          // inner gaussian std dev
const int   oc = 50;           // sample cutoff
// ---------------------------------------------


vec2 gaussian1(float i, vec2 a, vec2 d) {
     return a * exp( -(i*i) / d );
}

vec2 gaussian1d1(sampler2D sam, vec2 sigma, vec2 uv, vec2 tx) {
    vec2 a = vec2(1.0 / (sigma * SQRT_2_PI));
    vec2 d = vec2(2.0 * sigma * sigma);
    vec2 acc = vec2(0.0);
    vec2 sum = vec2(0.0);
    
    // centermost term
    acc += a * texture(sam, uv).x;
    sum += a;

    // sum up remaining terms symmetrically
    for (int i = 1; i <= oc; i++) {
        float fi = float(i);
        vec2 g = gaussian1(fi, a, d);
        vec2 posL = fract(uv - tx * fi);
        vec2 posR = fract(uv + tx * fi);
        acc += g * (texture(sam, posL).x + texture(sam, posR).x);
        sum += 2.0 * g;
    }

    return acc / sum;
}

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 tx = 1.0 / RENDERSIZE.xy;
    vec2 uv = fragCoord.xy * tx;
        uv += xy_flow*1;

    tx = (mod(float(FRAMECOUNT),2.0) < 1.0) ? vec2(tx.x,0) : vec2(0,tx.y);
    vec2 x_pass = gaussian1d1(BuffA, vec2(or, ir), uv, tx);
    fragColor = vec4(x_pass,0,0);
	return fragColor; 
 } 


			//******** BuffC Code Begins ********

#define SQRT_2_PI 2.50662827463

// ---------------------------------------------
// const float or3 = 18.0;         // outer gaussian std dev
// const float ir3 = 6.0;          // inner gaussian std dev
// ---------------------------------------------


vec2 gaussian(float i, vec2 a, vec2 d) {
     return a * exp( -(i*i) / d );
}

vec2 gaussian1d(sampler2D sam, vec2 sigma, vec2 uv, vec2 tx) {
    vec2 a = vec2(1.0 / (sigma * SQRT_2_PI));
    vec2 d = vec2(2.0 * sigma * sigma);
    vec2 acc = vec2(0.0);
    vec2 sum = vec2(0.0);
    
    // centermost term
    acc += a * texture(sam, uv).x;
    sum += a;

    // sum up remaining terms symmetrically
    for (int i = 1; i <= oc; i++) {
        float fi = float(i);
        vec2 g = gaussian(fi, a, d);
        vec2 posL = fract(uv - tx * fi);
        vec2 posR = fract(uv + tx * fi);
        acc += g * (texture(sam, posL).xy + texture(sam, posR).xy);
        sum += 2.0 * g;
    }

    return acc / sum;
}

vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 tx = 1.0 / RENDERSIZE.xy;
    vec2 uv = fragCoord.xy * tx;
            uv += xy_flow*1;

    tx = (mod(float(FRAMECOUNT),2.0) < 1.0) ? vec2(0,tx.y) : vec2(tx.x,0);
    vec2 y_pass = gaussian1d(BuffB, vec2(or, ir), uv, tx);
    fragColor = vec4(y_pass,0,0);
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord.xy / RENDERSIZE.xy;
	vec4 col = texture(BuffA, uv);
    
    float semiBlur = col.b;
    float blur = col.g;
    float hard = col.r;

    float colSel = mix(mix(blur, semiBlur, syn_Presence*audio_reactivity), hard, syn_HighHits*audio_reactivity);
    colSel = mix(colSel, (col.b+col.g+col.r)*(0.33+syn_BassLevel*0.3), max(no_blur, 1.0-audio_reactivity));
    // vec3 palCol = vec3(0.0,0.0,blur)*syn_BassLevel+colSel*_palette(colSel, vec3(1.439, -0.012, 0.314), vec3(0.746, 0.702, 0.276), vec3(0.680, 0.321, 1.190), vec3(0.407, 0.661, 0.759));
    vec3 palCol = vec3(0.0,0.0,blur)*syn_BassLevel+colSel*colorPaletteChooser(color_regime, colSel);

    vec3 normCol = col.x*vec3(1.0) + col.y*vec3(1,0.5,0) + col.z*vec3(0,0.5,1);
    fragColor = vec4(mix(palCol, normCol, greyscale), 1.0);
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
		return renderMainImage();
	}
}