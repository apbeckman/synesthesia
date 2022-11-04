vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


float t;
vec3 map;

float anim1(float x, float sm){
  float xmd = mod(x,2.) - .5;
  return smoothstep(-sm,sm,xmd) - smoothstep(-sm,sm,xmd - 1.);
}

float h(vec2 p){ return fract(sin(dot(p.xy,vec2(12.9898,78.233)))*43758.585); }
float hm(vec2 p, float v){ return fract(sin(dot(p.xy,vec2(12.9898,78.233)))*43758.585 + v); }

float c(vec2 p, float s) { return length(p) - s; }
float sq(vec2 p, float s, float r) { return length(max(abs(p) - s,0.)) - r; }
float ushape(vec2 p, float t, float sel) {
    p.x += sel;
    float d = max(abs(min(c(p - vec2(.5),.5),c(p + vec2(.5),.5))) - t,sq(p,.5,.001));
    p.x -= 1.;
    d = min(d,max(abs(min(c(p - vec2(.5,-.5),.5),c(p + vec2(.5,-.5),.5))) - t,sq(p,.5,.001)));
    return d - .05;
}
float pattern(vec2 p, float sz) {
    float pz = p.y;
    p *= sz;
    vec2 pf = floor(p); 
    p = fract(p) - .5;
    float rep = Repetitions;
    float d = ushape(p,rep*.5,step(h(pf),.75));
    int lim = int(anim1(t*.5 + pz*.25,.5)*11.);
    for(int i = 0; i < 20; i++) {
        if(i == lim) break;
        d = abs(d - rep) - rep;
    }
    return d/sz;
}

mat2 rot(float a) {
    float c = cos(a), s = sin(a);
    return mat2(c,s,-s,c);
}

float df(vec3 p) {
    p.z += t;
    p.y += .5;
    p.xy *= rot(cos(p.z));
    p.y -= .5;
    map = p;
    return max(-1.,-(length(p.xy) - 1.));
}

#define EPSI .001
vec3 normal(vec3 p){
    vec2 u = vec2(0.,EPSI); float d = df(p);
    return normalize(vec3(df(p + u.yxx),df(p + u.xyx),df(p + u.xxy)) - d);
}

#define MIN_DIST 0.
#define MAX_DIST 100.
#define MAX_STEPS 110
#define LIM .00001
float cptI = 0.;
vec3 rm(vec3 c, vec3 r) {
    vec3 color = vec3(1.0);
    vec3 p = c + r*MIN_DIST;
    for(int i = 0; i < MAX_STEPS; i++) {
        float d = df(p);
        if(d < LIM) {
            vec3 n = normal(p);
            color = n*.75+.5;
            vec3 dph = float(i)*.025*vec3(1.);
            vec2 pr = vec2(atan(map.x,map.y),map.z);
            float d1 = pattern(pr/vec2(3.14159265359*.25,1.) + vec2(t*2.,2.*t),1.);
            float aa = 1./RENDERSIZE.x;
            color = mix(color+dph,dph*(n*.5 + .5),smoothstep(-aa,aa,d1));
            color = mix(color, color + map*.5+.5, .0001);
            cptI = float(i);
            return color;
        }
        if(distance(c,p) > MAX_DIST) return color;
        p += d*r;
    }
    return color;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 st = fragCoord.xy/RENDERSIZE.xy - .5;
    st.x *= RENDERSIZE.x/RENDERSIZE.y;
    
    float tScale = .5 +  .0001*(cos(smooth_hightime*.5)*.5+.5);
    t = (TIME*tScale+(smooth_basstime*0.25))*0.5;
    vec2 xyshift = vec2(cos(t*2.)*.53,sin(t*2.)*.53);
    t -= pow(length((st + xyshift)*.5),1.3)*(.5+.5*cos(t*.5));
    vec3 c = vec3(xyshift,(-2.));
    vec3 r = normalize(vec3(st,.35));
    r.xy *= rot(-t*1.2);
   // r.xz *= rot(3.14*anim1(smooth_basstime*.05,.05));
    
    #define MBLUR 5
    vec3 color = vec3(0.);
    for(int i = 0; i < MBLUR; i++) {
        t += (.51*(1.-(cos(smooth_hightime*1.5)*.5+.5)))/60.;
        color += rm(c,r);
    }
    color /= float(MBLUR);
    fragColor = vec4(color + pow(length(st)*1.,4.),1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}