

			//******** BuffA Code Begins ********

#define RAYMARCH_MAX_STEPS 100
#define RAYMARCH_MAX_DIST 100.
#define RAYMARCH_SURFACE_DIST .001
#define SPEED 1.0
#define ROTSPEED 0.25
#define CAMDIST 4.0
#define color vec3(1.0, 0.2, 0.2)//vec3(0.1, 1.0, 0.4)
#define time mod(TIME, 60.0)
#define PI 3.14159265359
#define HALF_PI 1.57079632679
#define QUARTER_PI 0.78539816339

vec3 randomVector(vec3 p) {
	vec3 a = fract(p.xyz*vec3(123.34, 234.34, 345.65));
    a += dot(a, a+34.45);
    vec3 v = fract(vec3(a.x*a.y, a.y*a.z, a.z*a.x));
    return v;
}

vec3 lightPos1 = vec3(0.0, 0.0, 2.9457);
vec3 lightPos2 = vec3(0.0, 0.0, 4.8569);
vec3 lightPos3 = vec3(0.0, 0.0, 6.7142);
float light1 = 0.0;
float light2 = 0.0;
float light3 = 0.0;

vec3 lightPosGlobal = vec3(0.0, 0.0, 30.0);
float lightGlobal = 0.0;

float cameraMovement = 0.0;

float beamsShine = 0.0;
float tmpo = 0.0;

mat2 rotate(float a) {
    float s=sin(a), c=cos(a);
    return mat2(c, -s, s, c);
}

// This creates repetitions of something along the way
vec3 repeat(vec3 p, vec3 s) {
	return (fract(p/s - 0.5) - 0.5) * s;
}

vec2 repeat(vec2 p, vec2 s) {
	return (fract(p/s - 0.5) - 0.5) * s;
}

float repeat(float p, float s) {
	return (fract(p/s - 0.5) - 0.5) * s;
}

float sdCylinder(vec2 p, float r) {
    return length(p) - r;
}

vec3 tunnel(vec3 p) {
	vec3 off = vec3(0.0);
  	float dd = (p.z * 0.01) ;
	dd = floor(dd*1.0) + smoothstep(0.0, 1.0, smoothstep(0.0, 1.0, fract(dd*1.0)));
    dd *= 20.1;
    dd += time*0.1;
		// dd *= 0.1;
  	off.x += sin(dd) * 6.0;
  	off.y = sin(dd * 0.7) * 6.0;
  	return off;
}

vec3 navigate(vec3 p) {
  p += tunnel(p);
  p.xy *= rotate((p.z * ROTSPEED) + (time*ROTSPEED));
  p.y -= 0.3;
  return p;
}
float Hash21(vec2 p) {
    p = fract(p*vec2(123.34,233.53));
    p += dot(p, p+23.234);
    return fract(p.x*p.y);
}

float sdBox(vec3 p, vec3 s) {
    p = abs(p)-s;
	return length(max(p, 0.))+min(max(p.x, max(p.y, p.z)), 0.);
}

float SDG(vec3 p, float scale, float thickness, float bias, float lx, float ly) {
  vec3 p2 = p;
  p2 *= scale;
  float ls = max(lx, ly)*1.6;
  float gyroid = abs(dot(sin(p2*lx), cos(p2.zxy*ly))-bias) / (scale*ls) - thickness;
  return gyroid;
}

float smin(float a, float b, float h) {
	float k = clamp((a-b) / h * 0.5 + 0.5, 0.0, 1.0);
	return mix(a, b, k) - k * (1.0 - k) * h;
}

float getDist(vec3 p) {
    vec3 p2 = p;
    p2 = navigate(p2);

    float box = sdBox(p2, vec3(0.0));

    float lz = fract((p.z/100.0) * 0.02);
    float t = time*0.5;
    float lx = 1.25 + ((sin((lz+t) * 0.2576) * 0.5) + 0.5) * 0.25;
    float ly = 1.25 + ((cos((lz+t) * 0.1987) * 0.5) + 0.5) * 0.25;

    float g1 = SDG(p2, 0.543, 0.1, 1.4, lx, ly);

    float g2 = SDG(p2, 10.756, 0.03, 0.3, 1.0, 2.0);
    float g3 = SDG(p2, 20.765, 0.03, 0.3, 1.0, 1.0);
    float g4 = SDG(p2, 40.765, 0.03, 0.3, 1.0, 1.0);
    float g5 = SDG(p2, 60.765, 0.03, 0.3, 1.0, 1.0);
    float g6 = SDG(p2, 120.765, 0.03, 0.3, 1.0, 1.0);

    g1 -= g2*0.4;
    g1 -= g3*0.3;
    g1 -= g4*0.2;
    g1 -= g5*0.2;
    g1 += g6*0.1;

    float d = g1*0.6;//max(box, g1*0.4);
	
    // light bulbs
 	float dl1 = length(lightPos1 - p) - 0.1;
	light1 += 0.02/(0.01+dl1);
	d = smin(d, dl1, 0.5*2.0);
 	float dl2 = length(lightPos2 - p) - 0.1;
	light2 += 0.02/(0.01+dl2);
	d = smin(d, dl2, 0.5*2.0);
 	float dl3 = length(lightPos3 - p) - 0.1;
	light3 += 0.02/(0.01+dl3);
	d = smin(d, dl3, 0.5*2.0);

	// Beams
	vec3 p4 = p2;
	p4.xy *= rotate(p4.z * 0.05);
	p4.z = repeat(p4.z, 10.0);
	p4.x += sin(p4.y*0.3 + p2.z * 0.08 + time * 0.5) * 2.0;
	float beams = sdCylinder(p4.xz, 0.5);
	tmpo = beams;
	beams -= g2*0.4;
	beams -= g3*0.3;
    beams -= g4*0.2;
    beams -= g5*0.2;
    beams += g6*0.1;
	beams *= 0.8;
	d = smin(d, beams, 1.0);

  return d;
}

float rayMarch(vec3 ro, vec3 rd) {
	float dO = RAYMARCH_SURFACE_DIST;

    for (int i = 0; i < RAYMARCH_MAX_STEPS; i++) {
			vec3 p = ro + rd * dO;
      float dS = getDist(p);
      dO += dS;
      if (dO > RAYMARCH_MAX_DIST || dS < RAYMARCH_SURFACE_DIST) break;
    }

    return dO;
}

vec3 getNormal(vec3 p) {
	float d = getDist(p);
    vec2 e = vec2(0.01, 0.0);

    vec3 n = d - vec3(
        getDist(p - e.xyy),
        getDist(p - e.yxy),
        getDist(p - e.yyx));

    return normalize(n);
}

vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 f = normalize(l-p),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f,r),
        c = f*z,
        i = c + uv.x*r + uv.y*u,
        d = normalize(i);
    return d;
}

vec3 farColor(vec3 rd) {
  vec3 col = vec3(0.0);
  float y = rd.y * 0.5 + 0.5;
  return color * 0.1;
}

float getLight(vec3 p, vec3 lightPos) {
    vec3 l = normalize(lightPos - p);
    vec3 n = getNormal(p);

    float dif = dot(n, l);

    return dif;
}

vec3 calcLight(vec3 lp) {
  vec3 lp2 = randomVector(lp);
  float lpr = lightPos2.z;
  float lpRotation = lp2.z + time * 0.5;
  lp2.x = sin(lpRotation);
  lp2.y = cos(lpRotation);
  lp2.z = -CAMDIST+2.0 + (sin((lpr + time)*0.25) * 0.5 + 0.5) * 4.0;
  lp2.z += cameraMovement;
  lp2 -= tunnel(lp2);
  return lp2;
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    vec2 uv = (gl_FragCoord.xy-.5 * RENDERSIZE.xy) / RENDERSIZE.y;

    vec3 col = vec3(0.0);


    //


    // camera
    float t = time * SPEED;
    cameraMovement = t;
    vec3 ro = vec3(0.0, 0.0, -CAMDIST);
    ro.z += cameraMovement;
    ro -= tunnel(ro);
    vec3 ta = vec3(0.0);
    ta.z += cameraMovement;
    ta -= tunnel(ta);

    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
    vec3 vv = normalize(cross(uu, ww));

    float fov = 1.0;
    vec3 rd = normalize(uv.x * uu + uv.y * vv + ww * fov);
    
    
    //
    

    // lights 
    float lightRotation = time * 0.5;
    lightPos1.x = sin(lightRotation);
    lightPos1.y = cos(lightRotation);
    lightPos1.z = -CAMDIST+2.0 + (sin(time*0.25) * 0.5 + 0.5) * 8.0;
    lightPos1.z += cameraMovement;
    lightPos1 -= tunnel(lightPos1);

    lightPos2 = calcLight(lightPos2);
    lightPos3 = calcLight(lightPos3);

    float d = rayMarch(ro, rd);

    if(d<RAYMARCH_MAX_DIST) {
        vec3 p = ro + rd * d;
        vec3 n = getNormal(p);

        // calc lights
        float dif1 = clamp(getLight(p, lightPos1), 0.0, 1.0);
        col = vec3(dif1)*normalize(pow(color, vec3(0.585)));
        float dif2 = clamp(getLight(p, lightPos2), 0.0, 1.0);
        col = vec3(dif2)*normalize(pow(color, vec3(0.585)));
        float dif3 = clamp(getLight(p, lightPos3), 0.0, 1.0);
        col = vec3(dif3)*normalize(pow(color, vec3(0.585)));

        p = navigate(p);

        float g2 = SDG(p, 10.756, 0.03, 0.3, 1.0, 1.0);
        col *= smoothstep(-0.1, 0.06, g2);

        float cw = -0.04 + smoothstep(0.0, -0.5, n.y) * 0.08;
        float sh = smoothstep(cw, -0.03, g2);

        float pz = pow(p.z, 10.0/p.z);
        float g3 = SDG(p, 5.756, 0.03, 0.0, 1.0, 1.0);
        float g4 = SDG(p, 4.756, 0.03, 0.0, 1.0, 1.0);

        sh *= g3 * g4 * 0.8 + 0.2 * smoothstep(0.2, 0.0, n.y);

        col += sh * normalize(pow(color, vec3(0.5)))*20.0;

    }

    float mindist = 3.0;
    float maxdist = 8.0;
    col = mix(col, farColor(rd), smoothstep(mindist, maxdist, d));

    vec3 l1 = light1 * color*0.2;
    col += l1;
    //
    vec3 l2 = light2 * color*0.2;
    col += l2;
    //
    vec3 l3 = light3* color*0.2;
    col += l3;


    col = 1.0 - exp(-col * 3.0);
    col = pow(col, vec3(1.6));


    fragColor = vec4(col, d);
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

//  Blur  effect
//  Edited  from  https://www.shadertoy.com/view/XdfGDH
precision mediump float;

float  normpdf(in float  x,  in float  sigma)  {
    return  0.39894 * exp(-0.5 * x * x / (sigma * sigma)) / sigma;
}

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


	vec2  uv = fragCoord.xy / RENDERSIZE.xy;

    vec3  c = texture(BuffA, uv).rgb;

    vec2  center = vec2(0.5, 0.5);
    center = vec2(0.5, 0.5);

    float  d = smoothstep(0.3, 1.0, 0.1 + distance(center, uv));

    //  grain  effect
    float  strength = 8.0;
    float  x = (uv.x + 4.0) * (uv.y + 4.0) * (TIME * 10.0);
    vec3  grain = vec3(mod((mod(x, 13.0) + 1.0) * (mod(x, 123.0) + 1.0), 0.01) - 0.005) * strength;

    const int  mSize = 9;
    const int  kSize = (mSize - 1) / 2;
    float  kernel[mSize];
    vec3  final_colour = vec3(0.0);

    //create  the  1-D  kernel
    float  Z = 0.0;
    float sigma = 0.1 + texture(BuffA, uv).a * 7.0;
    for (int  j = 0; j <= kSize; ++j) {
        kernel[kSize + j] = kernel[kSize - j] = normpdf(float(j), sigma);
    }

    //get  the  normalization  factor  (as  the  gaussian  has  been  clamped)
    for (int  j = 0; j < mSize; ++j) {
        Z += kernel[j];
    }

    //read  out  the  texels
    for (int  i = -kSize; i <= kSize; ++i) {
        for(int  j = -kSize; j <= kSize; ++j) {
            final_colour += kernel[kSize + j] * kernel[kSize + i] * texture(BuffA, (gl_FragCoord.xy + vec2(float(i), float(j))) / RENDERSIZE.xy).rgb;
        }
    }

    vec3  c_step_1 = final_colour / (Z * Z);

    float  nd = 1.0 - d;
    vec3 c_step_2 = clamp(c_step_1 * nd, 0.0, 1.0);

    // I don't like the image too clean
    c_step_2 += grain * 1.0;

    fragColor = vec4(c_step_2, 1.0);
	return fragColor; 
 } 



// Fork of "GEHENNA" by kesson. https://shadertoy.com/view/3stBDr
// 2023-01-11 04:35:41

/* --- GEHENNA --- */

// Copyright 2020 - Giovanni Muzio
// https://kesson.io
//
//
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
//

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord / RENDERSIZE.xy;
    vec3 c = texture(BuffB, uv).rgb;
    fragColor = vec4(c, 1.0);
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
		return renderMainImage();
	}
}