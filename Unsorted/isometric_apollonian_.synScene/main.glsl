

/**

    Isometric Apollonian
    @byt3_m3chanic | 7/18/2021

    I really like @Flopine's and others isometric stuff.
    So my attempt at isometric using an Apollonian formula
    and some linear timing.

*/
#define R   RENDERSIZE
#define M   _mouse
#define T   TIME
#define PI  3.14159265359
#define PI2 6.28318530718

#define MAX_DIST    100.
#define MIN_DIST    1e-4

float time,ga1,ga2,ga3,ga4,tmod;
mat2 turn,r45,rx,ry;
mat2 rot(float a){ return mat2(cos(a),sin(a),-sin(a),cos(a)); }
float lsp(float begin, float end, float t) { return clamp((t - begin) / (end - begin), 0.0, 1.0); }

//@iq - all sdf's
float box(vec3 p, vec3 s) 
{
    p = abs(p)-s;
    return length(max(p, 0.))+min(max(p.x, max(p.y, p.z)), 0.);
}
float leaf(vec2 p, float r, float d)
{
    p = abs(p);
    float b = sqrt(r*r-d*d);
    return ((p.y-b)*d>p.x*b) ? length(p-vec2(0.0,b))
                             : length(p-vec2(-d,0.0))-r;
}

vec3 map (vec3 p) 
{
	float scale =  1.45;
    float py= (p.y+3.75)*.025;
    
    p.yz*=r45;
    p.xz*=r45;
    p.xz*=rx;

    vec3 q = p;
	float orb =1e7;
    float ph = 1.35;
    
    p.y-=ga1*2.;

	for( int i=0; i<2;i++ ) {
		p = -1.+2.*fract(.5*p+.5);
		float r2 = dot(p,p);  
        float k = ph/r2;
		p *= k;
        scale *= k;
        orb = min(orb,r2);
	}

    float d = length(p)-1.015;
    
    float thx = .02*scale;
    float tubes = length((p.xz))-thx;
    float td=d;
    d=min(tubes,d);
    d/=scale;

    float bx = box(q,vec3(1.275,6.,1.275));
    d=max(d,bx);

    float adr1 = .75*floor((2.15*p.y+1.15)*4.);
    float adr2 = .25*floor(1.5*length(p.xz)-p.y+.5);
    float adr=mix(adr2,adr1,clamp(py,.0,1.));
    if(tubes<td)adr=5.5*clamp(py,.0,1.);
	return vec3(d ,adr, orb);
}

vec3 normal(vec3 p, float t)
{
    float e = MIN_DIST*t;
    vec2 h =vec2(1,-1)*.5773;
    vec3 n = h.xyy * map(p+h.xyy*e).x+
             h.yyx * map(p+h.yyx*e).x+
             h.yxy * map(p+h.yxy*e).x+
             h.xxx * map(p+h.xxx*e).x;
    return normalize(n);
}

vec3 hue(float t)
{ 
    const vec3 c = vec3(0.122,0.663,1.000);
    return .65 + .45*cos(13.+PI2*t*(c+vec3(0.959,0.970,0.989))); 
}

void topLayer(inout vec3 C, vec2 uv, float alpha) 
{
    float px = fwidth(uv.x);
    float d = length(uv)-.5;
    d=abs(d)-.02;
    d=smoothstep(px,-px,d);

    vec2 mv = uv;
    mv.y+=.25*sin(mv.x*2.15+time);
    float bx = length(mv.y)-.1;
    bx=abs(abs(abs(abs(bx)-.01)-.01)-.01)-.005;
    bx=smoothstep(px,-px,bx);

    float e = length(uv)-.525;
    e=smoothstep(-px,px,e);

    C =mix(C,C+.5,min(min(bx,alpha),e));
    C =mix(C,vec3(1),d);   

    uv*=rot(T*35.*PI/180.);
    
    vec2 luv = abs(uv)-vec2(.14);
    luv*=r45;
    float fl = leaf(luv,.28,.125);

    luv = abs(uv*r45)-vec2(.14);
    luv*=r45;
    fl = min(leaf(luv,.28,.125),fl);
    
    fl=abs(abs(abs(fl)-.01)-.01)-.005;
    fl=smoothstep(px,-px,fl);
    C =mix(C,C+.5,fl);  
    C = clamp(C,vec3(0.),vec3(1.));
}

vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 F = _xy;
   
    // precal       
    time = T;
    
    r45=rot(-0.78539816339);
    turn=rot(T*5.*PI/180.);

    tmod = mod(time, 10.);
    float t1 = lsp(0.0, 4.0, tmod);
    float t2 = lsp(4.0, 10.0, tmod);
    ga1 = (t1)+floor(time*.1);
    ga2 = (t2)+floor(time*.1);
    rx=rot(-t2*1.57079632679);
    
    float t3 = lsp(1.0, 4.0, tmod);
    float t4 = lsp(6.0, 9.0, tmod);
    ga3 = (t3-t4);
    
    //
    vec2 uv = (2.*F.xy-R.xy)/max(R.x,R.y);
    
    //@Flopine's isometric setup - tweaked
    //https://www.shadertoy.com/view/NtXSWS
    vec3 ro = vec3(uv*2.25,-10.);
    vec3 rd = vec3(0.,0.,1.);
    //
    
    //background
    vec3 C = mix(vec3(0.737,0.918,0.459),vec3(0.180,0.643,0.761),clamp((uv.y+1.)*.75,0.,1.));
    vec2 f=fract((uv*2.75)-T*.2)-.5;
    if(f.x*f.y>0.) C=vec3(.03);
    
    vec3 p;
    float d=0.,m,orb,alpha=1.;
    // ray marcher
    for(int i=0;i<150;i++)
    {
        p=ro+rd*d;
        vec3 ray = map(p);
        d += ray.x;
        m  = ray.y;
        orb= ray.z;
        if(abs(ray.x)<MIN_DIST ||d>MAX_DIST)break;
    }
    
    if(d<MAX_DIST)
    {
        alpha=0.;
        vec3 n = normal(p,d);
        vec3 lpos =vec3(-1.05, -.25, -6.25);
        vec3 l = normalize(lpos-p);

        float diff = clamp(dot(n,l),0.,1.);
        
        float shdw = .95;
        for( float t=.001; t < 16.; )
        {
            float h = map(p + l*t).x;
            if( h<MIN_DIST ) { shdw = 0.; break; }
            shdw = min(shdw, 24.*h/t);
            t += h * .25;
            if( shdw<MIN_DIST || t>25. ) break;
        }

        diff = mix(diff,diff*shdw,.75);

        vec3 view = normalize(p - ro);
        vec3 ret = reflect(normalize(lpos), n);
        float spec = .5 * pow(max(dot(view, ret), 0.), 20.);

        vec3 h = hue(m);

        C = diff*h+spec;
    } 
     
    topLayer(C,uv,alpha);
    C=clamp(C,vec3(.03),vec3(1.));
    C=pow(C, vec3(.4545));
    O = vec4(C,1.0);
	return O; 
 } 





vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}