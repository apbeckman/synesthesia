

#define MAX_STEPS 120
#define MAX_DIST 120.
#define SURF_DIST .000125
#define Rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))
#define antialiasing(n) n/min(RENDERSIZE.y,RENDERSIZE.x)
#define S(d,b) smoothstep(antialiasing(1.0),b,d)
#define B(p,s) max(abs(p).x-s.x,abs(p).y-s.y)
#define Tri(p,s,a) max(-dot(p,vec2(cos(-a),sin(-a))),max(dot(p,vec2(cos(a),sin(a))),max(abs(p).x-s.x,abs(p).y-s.y)))
#define DF(a,b) length(a) * cos( mod( atan(a.y,a.x)+6.28/(b*8.0), 6.28/((b*8.0)*0.5))+(b-1.)*6.28/(b*8.0) + vec2(0,11) )
#define ZERO (min(int(FRAMECOUNT),0))

float Hash21(vec2 p) {
    p = fract(p*vec2(234.56,789.34));
    p+=dot(p,p+34.56);
    return fract(p.x+p.y);
}

// thx iq! https://iquilezles.org/articles/distfunctions/
float sdTorus( vec3 p, vec2 t )
{
    vec2 q = vec2(length(p.xy)-t.x,p.z);
    return length(q)-t.y;
}

// thx iq! https://iquilezles.org/articles/distfunctions/
float sdBox( vec3 p, vec3 b )
{
    vec3 q = abs(p) - b;
    return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

// thx iq! https://iquilezles.org/articles/distfunctions/
// tweaked as the center aligned horizontal capsule. 
float sdHorizontalCapsule( vec3 p, float w, float r )
{
  p.x-= clamp( p.x, -w*0.5, w*0.5 );
  return length( p ) - r;
}

float speaker(vec3 p){
    vec3 prevP = p;
    float d = sdBox(p, vec3(0.45,0.95,0.34))-0.03;
    float d2 = length(p-vec3(0.,-0.2,-0.53))-0.38;
    
    d = max(-d2,d);
    
    d2 = sdTorus(p-vec3(0.,-0.2,-0.36),vec2(0.36-syn_HighHits*0.0125,0.03));
    d = min(d,d2);
    
    d2 = sdTorus(p-vec3(0.,-0.2,-0.32),vec2(0.32-syn_HighHits*0.0125,0.025));
    d = min(d,d2);
    d2 = length(p-vec3(0.,-0.25,-0.08))-0.12;
    d = min(d,d2);
    
    d2 = sdHorizontalCapsule(p-vec3(0.,-0.75,-0.36),0.6,0.06);
    d = max(-d2,d);
    
    d2 = length(p-vec3(0.,0.55,-0.36))-0.2;
    d = max(-d2,d);
    
    d2 = sdTorus(p-vec3(0.,0.55,-0.36),vec2(0.2-syn_HighHits*0.0125,0.03));
    d = min(d,d2);
    
    p.z-=-0.36;
    p.x = abs(p.x)-0.4;
    p.y = abs(p.y)-0.9;
    d2 = length(p)-0.03;
    d = min(d,d2);
    
    return d;
}

float speaker2(vec3 p){
    vec3 prevP = p;
    float d = sdBox(p, vec3(0.95,0.45,0.34))-0.03;
    float d2 = sdBox(p-vec3(0.,0.,-0.35), vec3(0.9,0.4,0.01))-0.03;
    d = max(-d2,d);
    
    p.x = abs(p.x);
    d2 = length(p-vec3(0.4,0.,-0.5))-0.36;
    d = max(-d2,d);
    
    d2 = sdTorus(p-vec3(0.4,0.,-0.3),vec2(0.34,0.03));
    d = min(d,d2);
    
    d2 = sdTorus(p-vec3(0.4,0.,-0.29),vec2(0.3,0.025));
    d = min(d,d2);
    d2 = length(p-vec3(0.45,0.,-0.08))-0.1;
    d = min(d,d2);
    
    p.z-=-0.3;
    p.x = abs(p.x)-0.86;
    p.y = abs(p.y)-0.36;
    d2 = length(p)-0.03;
    d = min(d,d2);    
    
    return d;
}

float speaker3(vec3 p){
    vec3 prevP = p;
    float d = sdBox(p, vec3(0.95,0.95,0.34))-0.03;
    
    float d2 = length(p-vec3(0.0,0.,-0.68))-0.66;
    d = max(-d2,d);
    
    d2 = sdTorus(p-vec3(0.0,0.,-0.35),vec2(0.64-syn_BassHits*0.035,0.05));
    d = min(d,d2);
    
    d2 = sdTorus(p-vec3(0.0,0.,-0.33),vec2(0.6-syn_BassHits*0.05,0.045));
    d = min(d,d2);


    d2 = length(p-vec3(0.0,0.,0.1))-0.2;
    d = min(d,d2);
    
    d2 = sdTorus(p-vec3(0.0,0.,-0.3),vec2(0.56-syn_BassHits*0.025,0.035));
    d = min(d,d2);    
    
    d2 = sdTorus(p-vec3(0.0,0.,-0.24),vec2(0.52-syn_BassHits*0.025,0.035));
    d = min(d,d2);        
        
    d2 = sdTorus(p-vec3(0.0,0.,-0.19),vec2(0.47-syn_BassHits*0.025,0.035));
    d = min(d,d2);  
    
    d2 = abs(length(p.xy)-0.73)-0.07;
    d = min(d,max((abs(p.z)-0.38),d2));
    
    p.z-=-0.37;
    p.x = abs(p.x)-0.86;
    p.y = abs(p.y)-0.86;
    d2 = length(p)-0.03;
    d = min(d,d2);    
    
    p = prevP;
    p.z-=-0.37;
    p.xy = DF(p.xy,3.0);
    p.xy -= vec2(0.52);
    d2 = length(p)-0.03;
    d = min(d,d2); 
    
    p = prevP;
    p.xy*=Rot(radians((smoothTimeC*accentRotate)*120.));
    p.z-=-0.37;
    p.y=abs(p.y)-0.93;
    d2 = Tri(p.xy,vec2(0.08),radians(45.));
    d = min(d,max((abs(p.z)-0.02),d2));    
    
    p = prevP;
    p.xy*=Rot(radians(90.+(smoothTimeC*accentRotate)*120.));
    p.z-=-0.37;
    p.y=abs(p.y)-0.93;
    d2 = Tri(p.xy,vec2(0.08),radians(45.));
    d = min(d,max((abs(p.z)-0.02),d2));      
    
    return d;
}

float changeSpeakers(vec3 p, float start, float speed){
    vec3 prevP = p;
    float endTime = 3.;
    float t = smoothTimeC*0.1;
    float scenes[3] = float[](0.,1.,2.);
    for(int i = 0; i<scenes.length(); i++){
        scenes[i] = mod(scenes[i]+start,endTime);
    }
    
    float scene = scenes[int(mod(t,endTime))];
    
    float d = 10.;
    if(scene<1.) {
        p.x=abs(p.x)-0.5;
        d = speaker(p);
    } else if (scene >= 1. && scene<2.){
        p.y=abs(p.y)-0.5;
        d = speaker2(p);
    } else {
        d = speaker3(p);
    }
    
    return d;
}

vec2 GetDist(vec3 p) {
    vec3 prevP = p;
    
    p.y -=smoothTime*0.5;
    vec2 id = floor(p.xy*0.5);
    p.z-=3.;
    p.xy = mod(p.xy,2.0)-1.0;

    id*=.5;
    float rand = Hash21(id);
    
    float d = 10.;
    p.z-=rand*0.3;
    if(rand<0.3) {
        d = changeSpeakers(p,1.,0.5+rand);
    } else if(rand>=0.3 && rand<0.7) {
        d = speaker3(p);
    } else {
        p.x=abs(p.x)-0.5;
        d = speaker(p);
    }
    
    return vec2(d,0);
}

vec2 RayMarch(vec3 ro, vec3 rd, float side, int stepnum) {
    vec2 dO = vec2(0.0);
    
    for(int i=0; i<stepnum; i++) {
        vec3 p = ro + rd*dO.x;
        vec2 dS = GetDist(p);
        dO.x += dS.x*side;
        dO.y = dS.y;
        
        if(dO.x>MAX_DIST || abs(dS.x)<SURF_DIST) break;
    }
    
    return dO;
}

vec3 GetNormal(vec3 p) {
    float d = GetDist(p).x;
    vec2 e = vec2(.0001, 0);
    
    vec3 n = d - vec3(
        GetDist(p-e.xyy).x,
        GetDist(p-e.yxy).x,
        GetDist(p-e.yyx).x);
    
    return normalize(n);
}

vec3 R(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 f = normalize(l-p),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f,r),
        c = p+f*z,
        i = c + uv.x*r + uv.y*u,
        d = normalize(i-p);
    return d;
}

// https://www.shadertoy.com/view/3lsSzf
float calcOcclusion( in vec3 pos, in vec3 nor )
{
    float occ = 0.0;
    float sca = 1.0;
    for( int i=ZERO; i<4; i++ )
    {
        float h = 0.01 + 0.15*float(i)/4.0;
        vec3 opos = pos + h*nor;
        float d = GetDist( opos ).x;
        occ += (h-d)*sca;
        sca *= 0.95;
    }
    return clamp( 1.0 - 2.0*occ, 0.0, 1.0 );
}

vec3 diffuseMaterial(vec3 n, vec3 rd, vec3 p, vec3 col) {
    float occ = calcOcclusion(p,n);
    vec3 diffCol = vec3(0.0);
    vec3 lightDir = normalize(vec3(1+cos(smoothTimeB*0.25)*2.,10+sin(smoothTimeB*0.25)*10,-10));
    float diff = clamp(dot(n,lightDir),0.0,1.0);
    float shadow = step(RayMarch(p+n*0.3,lightDir,1.0, 15).x,0.9);
    float skyDiff = clamp(0.5+0.5*dot(n,vec3(0,1,0)),0.0,1.0);
    float bounceDiff = clamp(0.5+0.5*dot(n,vec3(0,-1,0)),0.0,1.0);
    diffCol = col*vec3(-0.5)*diff*shadow*occ;
    diffCol += col*vec3(1.0,1.0,0.9)*skyDiff*occ;
    diffCol += col*vec3(0.3,0.3,0.3)*bounceDiff*occ;
    diffCol += col*pow(max(dot(rd, reflect(lightDir, n)), 0.0), 60.)*occ; // spec
        
    return diffCol;
}

vec3 materials(int mat, vec3 n, vec3 rd, vec3 p, vec3 col){
    col = diffuseMaterial(n,rd,p,vec3(1.3));
    return col;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-.5*RENDERSIZE.xy)/RENDERSIZE.y;
    vec2 prevUV = uv;
    vec2 m =  LookXY.xy;
    
    vec3 ro = vec3(0, 0, -1.5);
    //if(_mouse.z>0.){
        ro.yz *= Rot(m.y*PI+1.);
        ro.y = max(-0.9,ro.y);
        ro.xz *= Rot(-m.x*2*PI);
    //} else {
        /*float scene = mod(TIME,15.);
        float rotY = -10.;
        float rotX = 0.;
        if(scene>=5. && scene<10.){
            rotY = 0.;
            rotX = -30.;
        } else if(scene>=10.){
            rotY = 0.;
            rotX = 30.;
        }
        
        ro.yz *= Rot(radians(rotY));
        ro.xz *= Rot(radians(rotX));
    }
    */
    vec3 rd = R(uv, ro, vec3(0,0.0,0), 1.0);
    vec2 d = RayMarch(ro, rd, 1.,MAX_STEPS);
    vec3 col = vec3(.0);
    
    if(d.x<MAX_DIST) {
        vec3 p = ro + rd * d.x;
        vec3 n = GetNormal(p);
        int mat = int(d.y);
        col = materials(mat,n,rd,p,col);
    }
    
    // gamma correction
    col = pow( col, vec3(0.9545) );    

    fragColor = vec4(sqrt(col),1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}