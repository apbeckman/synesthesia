vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


// "Corona Virus" 
// by Martijn Steinrucken aka BigWings/CountFrolic - 2020
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Music: Xerxes - early-morning-crystals
//
// This effect depicts my idea of what a virus could
// look like taking a huge artistic license. I started making what
// I imagine to be a lipid bilayer and then realized.. a virus
// doesn't have one! So then I just figured I'd make it look 'mean'
// 
// At first I tried using sphere coordinates but they distort too 
// much and have poles, so I ended up with projected cubemap
// coordinates instead. I think the function WorldToCube below
// is extremely usefull if you want to bend stuff around a sphere
// without too much distortion.
//
// As usual, the code could be a lot better and cleaner but I figure
// that by the time its all clean, I've lost interest and the world
// has moved on. Better ship it while its hot ;)
//
// uncomment the MODEL define to see once particle by itself
// you can change the amount of particles by changing the
// FILLED_CELLS define


//#define MODEL
#define FILLED_CELLS .3

#define PI 3.1415
#define TAU 6.2831

#define MAX_STEPS 400
#define MAX_DIST 40.
#define SURF_DIST .01



mat2 Rot(float a) {
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

float smin( float a, float b, float k ) {
    float h = clamp( 0.5+0.5*(b-a)/k, 0., 1. );
    return mix( b, a, h ) - k*h*(1.0-h);
}

float sdCapsule(vec3 p, vec3 a, vec3 b, float r) {
	vec3 ab = b-a;
    vec3 ap = p-a;
    
    float t = dot(ab, ap) / dot(ab, ab);
    t = clamp(t, 0., 1.);
    
    vec3 c = a + t*ab;
    
    return length(p-c)-r;
}

float N31(vec3 p) {
	vec3 a = fract(vec3(p) * vec3(213.897, 653.453, 253.098));
    a += dot(a, a + 79.76);
    return fract(a.x * a.y * a.z);
}

float N21(vec2 p) {
	vec3 a = fract(vec3(p.xyx) * vec3(213.897, 653.453, 253.098));
    a += dot(a, a.yzx + 79.76);
    return fract((a.x + a.y) * a.z);
}

vec3 SphereCoord(vec3 p) {
	float x = atan(p.x, p.z);
    float y = atan(length(p.xz), p.y);
    
    return vec3(x/TAU, length(p), 2.*y/TAU);
}

// returns cubemap coordinates
// xy = uv coords for face of cube, z = cube index (-3,-2,-1, 1, 2, 3)
vec3 WorldToCube(vec3 p) {
	vec3 ap = abs(p);
    vec3 sp = sign(p);
    float m = max(ap.x, max(ap.y, ap.z));
    vec3 st;
    if(m==ap.x)
        st = vec3(p.zy, 1.*sp.x);
    else if(m==ap.y)
        st = vec3(p.zx, 2.*sp.y);
    else
        st = vec3(p.xy, 3.*sp.z);
    
    st.xy /= m;
    
    // mattz' distortion correction
    st.xy *= (1.45109572583 - 0.451095725826*abs(st.xy));
       
    return st;
}

float Lipid(vec3 p, float twist, float scale) {
    vec3 n = sin(p*20.)*.2;
    p *= scale;
    
	p.xz*=Rot(p.y*.3*twist);
    p.x = abs(p.x);
    
    float d = length(p+n)-2.;
    
    float y = p.y*.025;
    float r = .05*scale;
    float s = length(p.xz-vec2(1.5,0))-r+max(.4,p.y);
    d = smin(d, s*.9,.4);
    
    return d/scale;
}

float sdTentacle(vec3 p) {
    float offs = sin(p.x*50.)*sin(p.y*30.)*sin(p.z*20.);
    
    p.x += sin(p.y*10.+script_high_time)*.02;
    p.y *= .2;
    
    float d = sdCapsule(p, vec3(0,0.1,0), vec3(0,.8,0), .04);
    
    p.xz = abs(p.xz);
    
    d = min(d, sdCapsule(p, vec3(0,.8,0), vec3(.1,.9,.1), .01));
    d += offs*.01;
    
    return d;
}


float Particle(vec3 p, float scale, float amount) {  
    float t = (script_time*0.75+script_bass_time*0.25);
 
    vec3 st = WorldToCube(p);
    vec3 cPos = vec3(st.x, length(p), st.y);
    vec3 tPos = cPos;
    
    vec3 size = vec3(.05);
  	
    cPos.xz *= scale;
    vec2 uv = fract(cPos.xz)-.5;
    vec2 id = floor(cPos.xz);
    
    uv = fract(cPos.xz)-.5;
    id = floor(cPos.xz);
    
    
    float n = N21(id);
    
    t = (t+st.z+n*123.32)*1.3;
    float wobble = sin(t)+sin(1.3*t)*.4;
    wobble /= 1.4;
    
    wobble *= wobble*wobble;
    
    wobble = wobble*amount/scale;
    vec3 ccPos = vec3(uv.x, cPos.y, uv.y);
    vec3 sPos = vec3(0, 3.5+wobble, .0);
    
    vec3 pos = ccPos-sPos;
    
    pos.y *= scale/2.;
   
    float r = 16./scale;
    r/=sPos.y; // account for height
    float d = length(pos)-r;
    d = Lipid(pos, n, 10.)/scale;
    
    d = min(d, length(p)-.2*scale);	// inside blocker
    
    
    float tent = sdTentacle(tPos);
    d = min(d, tent);
    
    return d;
}

float dCell(vec3 p, float size) {
    p = abs(p);
    float d = max(p.x, max(p.y, p.z));
    
    return max(0., size - d);
}

float GetDist(vec3 p) {
	float t = TIME;
    
    float scale=8.;
    
    #ifndef MODEL
    p.z += t;
    vec3 id = floor(p/10.);
    p = mod(p, vec3(10))-5.;
    float n = N21(id.xz);
    p.xz *= Rot(t*.2*(n-.5));
    p.yz *= Rot(t*.2*(N21(id.zx)-.5));
    scale = mix(4., 16., N21(id.xz));//mod(id.x+id.y+id.z, 8.)*2.;
    
    n = N31(id);
    if(n>FILLED_CELLS) {			// skip certain cells
        return dCell(p, 5.)+.1;
    }
    #endif
    
   
    
    p += sin(p.x+t)*.1+sin(p.y*p.z+t)*.05;
    
   
    float surf = sin(scale+t*.2)*.5+.5;
    surf *= surf;
    surf *= 4.;
    surf += 2.;
    float d = Particle(p, scale, surf);
    
    p.xz *= Rot(.78+t*.08);
    p.zy *= Rot(.5);
    
    d = smin(d, Particle(p, scale, surf), .02);
    
    
    return d;
}

float RayMarch(vec3 ro, vec3 rd) {
	float dO=0.;
    float cone = .0005;
    for(int i=0; i<MAX_STEPS; i++) {
    	vec3 p = ro + rd*dO;
        
        float dS = GetDist(p);
        dO += dS;
        if(dO>MAX_DIST || abs(dS)<SURF_DIST+dO*cone) break;
    }
    
    return dO;
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


vec3 R(vec2 uv, vec3 p, vec3 l, vec3 up, float z) {
    vec3 f = normalize(l-p),
        r = normalize(cross(up, f)),
        u = cross(f,r),
        c = p+f*z,
        i = c + uv.x*r + uv.y*u,
        d = normalize(i-p);
    return d;
}



vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-.5*RENDERSIZE.xy)/RENDERSIZE.y;
	vec2 m = (iMouse.xy-.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    if(m.x<-0.45&&m.y<-.45) m += .5;
    
    float t = TIME;
    
    vec3 col = vec3(0);
    
    #ifdef MODEL
    vec3 ro = vec3(0, 0, -7);
    //ro.y += sin(t*.1)*3.;
    ro.yz *= Rot(-m.y*2.);
    ro.xz *= Rot(TIME*.0-m.x*6.2831);
    vec3 rd = R(uv, ro, vec3(0,0,0), vec3(0,1,0), .5);
    //ro += 3.;
    #else
    vec3 ro = vec3(0, 0, -1);
    //ro.y += sin(t*.1)*3.;
    ro.yz *= Rot(-m.y*2.);
    ro.xz *= Rot(TIME*.0-m.x*6.2831);
    
    vec3 up = vec3(0,1,0);
    up.xy *= Rot(sin(t*.1));
    vec3 rd = R(uv, ro, vec3(0,0,0), up, .5);
    
    ro.x += 5.;
    ro.xy *= Rot(t*.1);
    ro.xy -= 5.;
    #endif
    
    float d = RayMarch(ro, rd);
    
    float bg = rd.y*.5+.3;
    float poleDist = length(rd.xz);
    float poleMask = smoothstep(.5, 0., poleDist);
    bg += sign(rd.y)*poleMask;
    
    float a = atan(rd.x, rd.z);
    bg += (sin(a*5.+t+rd.y*2.)+sin(a*7.-t+rd.y*2.))*.2;
    float rays = (sin(a*5.+t*2.+rd.y*2.)*sin(a*37.-t+rd.y*2.))*.5+.5;
    bg *= mix(1., rays, .25*poleDist*(sin(t*.1)*.5+.5));//*poleDist*poleDist*.25;
    col += bg;
    
    if(d<MAX_DIST) {
    	vec3 p = ro + rd * d;
   
    	vec3 n = GetNormal(p);
       
        #ifndef MODEL
        p = mod(p, vec3(10))-5.;
        #endif
        
        float s = dot(n, normalize(p))-.4;
        float f = -dot(rd, n);
        
        col += dot(n,-rd)*.5+.5;
    	//col += (1.-f*f)*s*1.5;
        
        col *= 0.;
        float r = 3.7;
        float ao = smoothstep(r*.8, r, length(p));
        col += (n.y*.5+.5)*ao*2.;
        //col *= 2.;
        col *= smoothstep(-1., 6., p.y);
        
        //col += n*.5+.5;
    }
    
    col = mix(col, vec3(bg), smoothstep(0., 40., d));
    
    //col *= vec3(1., .9, .8);
    //col = 1.-col;
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}