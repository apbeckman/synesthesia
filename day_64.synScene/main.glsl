vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** Common Code Begins ********
float smoothTime = (smooth_basstime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.95;
float smoothTimeB = (smooth_hightime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.95;
float smoothTimeC = (smooth_midtime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.95;

// Repeat space along one axis. Use like this to repeat along the x axis:
// <float cell = pMod1(p.x,5);> - using the return value is optional.
float pMod1(inout float p, float size) {
	float halfsize = size*0.5;
	float c = floor((p + halfsize)/size);
	p = mod(p + halfsize, size) - halfsize;
	return c;
}
// Shortcut for 45-degrees rotation
void pR45(inout vec2 p) {
	p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}
float fOpDifferenceColumns(float a, float b, float r, float n) {
	a = -a;
	float m = min(a, b);
	//avoid the expensive computation where not needed (produces discontinuity though)
	if ((a < r) && (b < r)) {
		vec2 p = vec2(a, b);
		float columnradius = r*sqrt(2.)/n/2.0;
		columnradius = r*sqrt(2.)/((n-1.)*2.+sqrt(2.));

		pR45(p);
		p.y += columnradius;
		p.x -= sqrt(2.)/2.*r;
		p.x += -columnradius*sqrt(2.)/2.;

		if (mod(n,2.) == 1.) {
			p.y += columnradius;
		}
		pMod1(p.y,columnradius*2.);

		float result = -length(p) + columnradius;
		result = max(result, p.x);
		result = min(result, a);
		return -min(result, b);
	} else {
		return -m;
	}
}

float fOpIntersectionColumns(float a, float b, float r, float n) {
	return fOpDifferenceColumns(a,-b,r, n);
}

// The "Stairs" flavour produces n-1 steps of a staircase:
// much less stupid version by paniq
float fOpUnionStairs(float a, float b, float r, float n) {
	float s = r/n;
	float u = b-r;
	return min(min(a,b), 0.5 * (u + a + abs ((mod (u - a + s, 2. * s)) - s)));
}

// We can just call Union since stairs are symmetric.
float fOpIntersectionStairs(float a, float b, float r, float n) {
	return -fOpUnionStairs(-a, -b, r, n);
}

float fOpDifferenceStairs(float a, float b, float r, float n) {
	return -fOpUnionStairs(-a, b, r, n);
}


			//******** BuffA Code Begins ********

// Fork of "Day 62" by jeyko. https://shadertoy.com/view/wtySRh
// 2020-02-21 08:51:00


vec3 glow = vec3(0);
vec3 glowG = vec3(0);
vec3 glowB = vec3(0);
#define pmod(p, x) mod(p, x) - 0.5*x
#define pi acos(-1.)
#define tau (2.*pi)
//#define mx (TIME*(0.8 )+ sin(TIME*2.)*0.4 + 20.*iMouse.x/RENDERSIZE.x)
//#define mx (TIME*(0.6 )+ sin(TIME*2.)*0.3 + 20.*iMouse.x/RENDERSIZE.x)

#define mx (   (smoothTime )+ sin(smoothTime*1.)*0.2 + 20.*iMouse.x/RENDERSIZE.x)

vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
    uv *= 0.8;
	vec3 dir = normalize(lookAt - ro);
	vec3 right = normalize(cross(vec3(0,1,0), dir));
	vec3 up = normalize(cross(dir, right));
    return normalize(dir + right*uv.x + up*uv.y);
}


float sdBox(vec3 p, vec3 s){
	p = abs(p) - s;
    return max(p.x, max(p.y, p.z));
}
#define pal(a,b,c,d,e) (a + b*sin(c*d + e))
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
#define modD  (1.)
#define dmin(a,b) a.x < b.x ? a : b

float zid = 0.;
vec3 pp = vec3(0);


vec2 map(vec3 p){
	vec2 d = vec2(10e6);
    
    p.xz = pmod(p.xz, modD);
    
    vec3 q = p;
    
    #define ceilH 0.5
    float ceils = p.y + ceilH;
    ceils = min(ceils, -p.y + ceilH);
    d.x = min(d.x, ceils );
    #define colW 0.0328
    //float cols = length(p.xz) - 0.06;
    float cols = max(abs(p.x) - colW, abs(p.z) - colW);
    #define floors 8.
    
    d = dmin(d, vec2(fOpUnionStairs(d.x, cols, 2.*ceilH/2., floors), 0.));
    
    // outer
    p = q;
    for(int i = 0; i < 4; i++){
    	p = abs(p);
        p.xz *= rot(0.125*pi);
        //q.x += 0.03;
    }
    p = abs(p);
    //p.xy *= rot(0.4);
    //p.y -= 0.2;
    q = p;
    p = abs(p);
    p.z -= 0.25;
    p.xz *= rot(pi*0.25);
    p.yz = abs(p.yz) - 0.02;
	float outer = max(p.y,p.z);  
    
    //q.xz *= rot(0.1);
    q.xz *= rot(0.25*pi);
    q = abs(q);
    q.xz *= rot(0.25*pi);
    q = abs(q);
    q -= 0.02;
    pp = q;
    d = dmin(d, vec2(outer, 1.));
    //d.x = min(d.x, outer);
    
    q = abs(q);
	outer = max(q.y,q.x);
    
    outer = max(outer, -q.z + 0.23) + 0.001;
    
    d = dmin(d, vec2(outer, 2.));
    //d.x = min(d.x, outer);
	//outer = max(outer,max(q.y,q.x));
    //glowB += exp(-outer*20.)*1.;
    float gBsc = 0.05/(0.01 + outer*outer*200.);
    glowB += gBsc*pal(0.0, 0.9, vec3(4.,2.,3.), vec3(1.,3.,1.6),5.9)*1.3;
    
    glowB += pow(abs(sin(p.x*4.+ TIME) ), 200.)*gBsc*pal(0.9, 0.9, vec3(4.,9.,3.), vec3(1.,3.,2.6),5.9)*2.;
    
    
    glowG += exp(-d.x*200. )*pal(0.88, 0.2, vec3(4.,2.,3.), vec3(1.,2.3,1.6),5.9)*0.7;

    d.x *= 0.7;
    glow += exp(-d.x*150. );
    return d;
}

vec2 march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	vec2 d = vec2(10e6);
	p = ro;
    hit = false; t = 0.;
    
    for(int i = 0; i < 250  ; i++){
    	d = map(p);
        if(d.x < 0.0005){
        	hit = true;
            break;
        }
        t += d.x;
        p = ro + rd*t;
    }
    
    return d;
}
vec3 getNormal(vec3 p){
	vec2 u = vec2(0.001,0.);
    return normalize(map(p).x - vec3(
    	map(p - u.xyy).x,
    	map(p - u.yxy).x,
    	map(p - u.yyx).x
    ));
}

vec3 text(vec2 t, vec3 p){
	vec3 o = vec3(0);
	
    float d = 10e6;
    t = pmod(t,1./16.);
    
    t *= 16.;
    float yid = (floor( (p.y + 1.)*16. ) );
    //t *= rot(0.25*pi);
    //t *= rot(0.5*pi);
    float W = 0.02;
    
    float modd = 0.15;
    //t = abs(mod(t,modd)/modd - 0.5);
    #define lmod(d, x) (mod(d,x)/x - 0.5)
    float sqD = max(abs(t.x), abs(t.y));
    sqD += TIME*0.2 + yid*0.04;
    float sqid = floor(sqD/modd);
    sqD = lmod(sqD, modd);
    
    d = min(d, sqD);
    //d = min(d, length(t.x) - W );
    //d = min(d, length(t.y) - W );

    
    o +=  pal(0.6, vec3(1.,0.7,0.6)*0.5, vec3(8.4 ,4.19,7.4 - yid*0.2), vec3(3.,7.,3.),-1. + smoothTimeC + sqid*0.5 + p.z + t.x*2.5 - t.y*1.5);
    o *= step(sin(sqid*40.), -0.3);

    
    float aa = 20.;
    //sqD -= modd*0.25;
    sqD -= 0.5;
    sqD = abs(sqD*1.);
    o -= exp(-sqD*aa)*4.;
    sqD -= 1.;
    sqD = abs(sqD*1.);
    o -= exp(-sqD*aa)*4.;
    //o += smoothstep(0.001,0., d);
    
    return o;
}
// Tri-Planar blending function. Based on an old Nvidia writeup:
// GPU Gems 3 - Ryan Geiss: https://developer.nvidia.com/gpugems/GPUGems3/gpugems3_ch01.html
vec3 tex3D(  in vec3 p, in vec3 n ){
   
    n = max((abs(n) - .2)*7., .001);
    n /= (n.x + n.y + n.z );  
    
	vec3 q = (text(p.yz, p)*n.x + text(p.zx,p)*n.y + text(p.xy,p)*n.z).xyz;
    
    return q;
}
vec4 nint(float t){
	vec4 a = texture(image30, vec2(floor(t)*0.02));
	vec4 b = texture(image30, vec2(floor(t+ 1.)*0.02));
    return mix(a,b, pow(smoothstep(0.,1.,fract(t)),1.));
}
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    uv *= 1. - length(uv)*0.2;
    vec3 col = vec3(0);
	
    vec3 ro = vec3(0);
    ro.y -= 0. + sin(TIME*0.25)*0.4;
    ro.z += mx;
    vec3 lookAt = vec3(0,0,ro.z + 2.);
    
    vec3 rd = getRd(ro, lookAt, uv);
   
    float r = nint(TIME/tau).x;
    rd.xy *= rot(sin(nint(0.5*TIME/tau).x*tau)*0.2);
    
    rd.xz *= rot(sin(nint(0.5*TIME/tau).y
                     *tau)*0.3);
    
    bool hit; float t; vec3 p;
    vec2 d = march(ro, rd, p, t, hit);
    vec3 n = getNormal(p);

    float fogD = smoothstep(0.,1.,t*0.18);
    if(hit){
        //col += pal(0.5, 0.4, vec3(3.,1.,1.), vec3(1.,1.,2.),2.6 - did*2.3 )*0.09;
        //col += pal(0.5, 0.6, vec3(3.,1.1,1.), vec3(1.,1.,2.),0.4 - did*0.6 )*0.09;
        //p.y -= 0.03;
        p.x -= 0.03;
        p.z -= 0.03;
    } 
    
    if(d.y == 1.){

        vec3 tt = tex3D(pp*0.1,n)*1.;
        
		col += 0.1*length(tt)* pal(0.7, vec3(1.,0.7,0.6)*0.5, vec3(8.4 ,4.19,7.4), vec3(4.,7.,3.),-1. + TIME);; 
    } 
    if(d.y == 2.){
        
        col = clamp(col, 0., 1.);
        col -= col;
		col += 0.2*pal(0.3, vec3(1.,1.,0.6)*0.5, vec3(8.4 ,4.19,7.4), vec3(4.,7.,3.),-1. + TIME);; 
    }
    col -= glowG*0.1;
    col += glow*0.03;
    col += glowB*0.03;
    if(d.y == 0.){
    	col += tex3D(p,n)*1.;
    }
    col = clamp(col, 0., 1.);

    col = mix(col,vec3(0.2,0.014,0.1)*0.2,fogD);
    //fogD = smoothstep(0.,1.,t*0.98);
    col += glow*0.003*fogD*vec3(0.2,0.064,0.1);
    //col += pal(0.5, 0.5, vec3(4.,2.,1.), vec3(1.,1.,1.),2.9 + pid.x + did*2.);
    col = pow(col, vec3(0.45));
    
    col = mix(col, vec3(0), dot(uv*0.55,uv*0.55)*2.);
    //col = mix(col, smoothstep(0.,1., col), 0.6);
//    col *= 1. - pow(abs(uv.x) - 0.25,6. )*50.;
//    col *= 1. - pow(abs(uv.y) - 0.3,3. )*50.;
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 



vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    
    
    //float m = pow(abs(sin(p.z*0.03)),10.);

    // Radial blur
    float steps = 10.;
    float scale = 0.00 + pow(length(uv - 0.5)*1.2,5.)*0.4;
    //float chromAb = smoothstep(0.,1.,pow(length(uv - 0.5), 0.3))*1.1;
    float chromAb = pow(length(uv - 0.5),1.)*3.1;
    vec2 offs = vec2(0);
    vec4 radial = vec4(0);
    for(float i = 0.; i < steps; i++){
    
        scale *= 0.91;
        vec2 target = uv + offs;
        offs -= normalize(uvn)*scale/steps;
    	radial.r += texture(BuffA, target + chromAb*1./RENDERSIZE.xy).x;
    	radial.g += texture(BuffA, target).y;
    	radial.b += texture(BuffA, target - chromAb*1./RENDERSIZE.xy).z;
    }
    radial /= steps;
    
    
    fragColor = radial*1.; 
    
    fragColor *= 1.2;
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.2);
    
    fragColor = max(fragColor, 0.);
    
    fragColor.xyz = pow(fragColor.xyz, vec3(1.5,1. + sin(TIME)*0.2,1. - cos(TIME)*0.4));

    //fragColor = pow(fragColor, vec4(0.4545 + dot(uvn,uvn)*1.7));
    fragColor *= 1. - dot(uvn,uvn)*1.8;
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}