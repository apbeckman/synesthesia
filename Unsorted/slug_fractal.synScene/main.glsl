

			//******** Common Code Begins ********

#define pi 3.14159

#define thc(a,b) tanh(a*cos(b))/tanh(a)
#define ths(a,b) tanh(a*sin(b))/tanh(a)
#define sabs(x) sqrt(x*x+1e-2)
//#define sabs(x, k) sqrt(x*x+k)

#define rot(a) mat2(cos(a), -sin(a), sin(a), cos(a))

float cc(float a, float b) {
    float f = thc(a, b);
    return sign(f) * pow(abs(f), 0.25);
}

float cs(float a, float b) {
    float f = ths(a, b);
    return sign(f) * pow(abs(f), 0.25);
}

vec3 pal(in float t, in vec3 d) {
    return 0.5 + 0.5 * cos(2. * pi * (0.5 * t + d));
}

vec3 pal(in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d) {
    return a + b * cos(2. * pi * (c * t + d));
}

float h21(vec2 a) {
    return fract(sin(dot(a.xy, vec2(12.9898, 78.233))) * 43758.5453123);
}

float mlength(vec2 uv) {
    return max(abs(uv.x), abs(uv.y));
}

float mlength(vec3 uv) {
    return max(max(abs(uv.x), abs(uv.y)), abs(uv.z));
}

float sfloor(float a, float b) {
    return floor(b) + 0.5 + 0.5 * tanh(a * (fract(b) - 0.5)) / tanh(0.5 * a);
}

// From iq, k = 0.12 is good
float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0., 1.);
    return mix(b, a, h) - k * h * (1. - h);
}

float smax(float a, float b, float k) {
    float h = clamp(0.5 - 0.5 * (b - a) / k, 0., 1.);
    return mix(b, a, h) + k * h * (1. - h); 
}

#define MAX_STEPS 400
#define MAX_DIST 100.
#define SURF_DIST .001

// RayMarching from TheArtOfCode

// From BlackleMori
#define FK(k) floatBitsToInt(k*k/7.)^floatBitsToInt(k)
float hash(float a, float b) {
    int x = FK(a), y = FK(b);
    return float((x*x+y)*(y*y-x)-x)/2.14e9;
}

vec3 erot(vec3 p, vec3 ax, float ro) {
  return mix(dot(ax, p)*ax, p, cos(ro)) + cross(ax,p)*sin(ro);
}

vec3 face(vec3 p) {
     vec3 a = abs(p);
     return step(a.yzx, a.xyz)*step(a.zxy, a.xyz)*sign(p);
}

float sdBox(vec3 p, vec3 s) {
    p = abs(p)-s;
	return length(max(p, 0.))+min(max(p.x, max(p.y, p.z)), 0.);
}

float sdBox(in vec2 p, in vec2 b){
    vec2 d = abs(p)-b;
    return length(max(d,0.))+min(max(d.x,d.y),0.);
}

vec3 GetRayOrigin() {
    vec2 m = _mouse.xy/RENDERSIZE.xy;
    float r = 4.;
    float a = 0.0525 * smoothTime;
    vec3 ro = vec3(r * cos(a), 0, r * sin(a));
    //ro.yz *= rot(-m.y*3.14+1.);
    //ro.xz *= rot(-m.x*6.2831);
    return ro;
}

float GetDist(vec3 p) {
    float sd = length(p) - 0.5;
    for(int j=0; j<16; j++){
        float io = 2. * pi * float(j) / 16.;
        vec3 q =normalize(p.yzx);
        //p.xy *= rot(p.z + 0.25 * TIME);
        p.xz = abs(p.zx * rot(0.125 * smoothTimeC + io) - 0.25 * abs(p.z-p.x))-0.0;
        p.yz = abs(p.yz * rot(0.125 * smoothTimeC + io) + .25 * abs(p.z - p.y))-0.01;
      // p.y += 0.01;
      // p.x -= 0.01;
      // p.z -= 0.01;
    }
    return max(sd, 0.1 * abs(p.x + p.y + p.z) - 0.01);
   // return max(0., length(p.yxz)/10. - 0.01);
}

float RayMarch(vec3 ro, vec3 rd, float z) {
	
    float dO=0.;
    float s = sign(z);
    for(int i=0; i<MAX_STEPS; i++) {
    	vec3 p = ro + rd*dO;
        float dS = GetDist(p);
        if (s != sign(dS)) { z *= 0.5; s = sign(dS); }
        if(abs(dS)<SURF_DIST || dO>MAX_DIST) break;
        dO += dS*z; 
    }
    
    return min(dO, MAX_DIST);
}

vec3 GetNormal(vec3 p) {
	float d = GetDist(p);
    vec2 e = vec2(.001, 0);
    
    vec3 n = d - vec3(
        GetDist(p-e.xyy),
        GetDist(p-e.yxy),
        GetDist(p-e.yyx));
    
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

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-.5*RENDERSIZE.xy)/RENDERSIZE.y;
	
    vec3 ro = GetRayOrigin();
    
    vec3 rd = GetRayDir(uv, ro, vec3(0), 1. + 5. * tanh(0.25 * smoothTime));
    vec3 col = vec3(0);
   
    float d = RayMarch(ro, rd, 1.);
    
    vec3 p = ro + rd * d;
    float IOR = 1.05;
    if(d<MAX_DIST) {        
        vec3 n = GetNormal(p);
        vec3 r = reflect(rd, n);

        vec3 pIn = p - 4. * SURF_DIST * n;
        vec3 rdIn = refract(rd, n, 1./IOR);
        float dIn = RayMarch(pIn, rdIn, -1.);
        
        vec3 pExit = pIn + dIn * rdIn;
        vec3 nExit = -GetNormal(pExit); // *-1.; ?

        float dif = dot(n, normalize(p))*.5+.5;
        float dif2 = dot(n, normalize(abs(r-p)))*.5+.5;
        col = vec3(dif);
       // col *= clamp(0.25 + p.y, 0., 1.);
        float fres = pow(1. + dot(rd, n), 4.);
        float fog = 1.-exp(-length(p));

        float spec = pow(dif, 32.);
        float spec2 = pow(dif2, 1.);
        
        col *= 0.5 + 0.5 * n;
        col *= pal(n.y * 0.5 + 0.1 * smoothTimeB, vec3(0,1,2)/2.);
        col *= 1. + 0.75 * n.y;
        
        
        col += spec;
        //col = vec3(spec);
        col = pal(dif + p.y * 0.5 + 0.5 * smoothTimeB, 0.25 * vec3(2,1,1));
        col += spec;
        col += p.y -0.0;
       
        //col = sqrt(col);
        float th = .6 + .4 * thc(8., 12. * log(length(p)) - 1.5 * pi * TIME);
       
       col *= th *exp(-12. * abs(fract(log(length(p))-1. * TIME)-0.5)) * abs(n.x + n.y + n.z);
       //col *= exp(-4. * length(p)) * abs(n.y);
        col *= abs(n.x + n.y + n.z) ;
    }
    
    
    col = pow(col, vec3(.4545));	// gamma correction
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}