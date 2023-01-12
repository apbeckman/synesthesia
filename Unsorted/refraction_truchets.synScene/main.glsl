

/**
    License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
    
    Refraction Truchets [015] (mouseable)
    11/12/21 @byt3_m3chanic
    
    An update to one of my earlier shaders, where I could barely do
    2 bounce reflections - now sporting the refraction and some 
    smooth movement.
    
    I'd like to figure out chromatic aberration for the refraction
    or maybe I should just get into ray tracing?!

    original https://www.shadertoy.com/view/sdXSW4

*/

#define R 		    RENDERSIZE
#define T 		    TIME
#define M 		    _mouse

#define PI          3.141592653
#define PI2         6.283185307

float hash21(vec2 a) { return fract(sin(dot(a,vec2(21.23,41.232)))*4123.2323); }
mat2 rot(float a) { return mat2(cos(a),sin(a),-sin(a),cos(a)); }
float lsp(float begin, float end, float t) { return clamp((t - begin) / (end - begin), .0, 1.); }
float eoc(float t) { return (t = t - 1.) * t * t + 1.; }

//@iq sdf's
float box( vec3 p, vec3 b ){
    vec3 q = abs(p) - b;
    return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}
float cyl(vec3 p, float h, float r ){
    vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(h,r);
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}
float cap( vec3 p, vec2 hr ){
    p.y -= clamp( p.y, 0.0, hr.x );
    return length( p ) - hr.y;
}
float torus( vec3 p, vec2 t) {
    vec2 q = vec2(length(p.xy)-t.x,p.z);   
    return length(q)-t.y;
}
//@iq
float smin( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); 
}

vec2 mod2(inout vec3 p, float scale){
    float hlf = scale*.5;
    vec2 id=floor((p.xy+hlf)/scale);
    p.xy = mod(p.xy+hlf,scale)-hlf;
    return id;
}

vec3 hit,hitPoint;
vec2 id,tileId;
mat2 turn,trs;
float zIndex=0.,zid=0.;
float tmod=0.,ga2=0.,ga5=0.;

// sdf scene
const float size = 4.;
const float hlf = size*.5;
const float hzf = hlf*.5;
const float sp  = 15.;
const float hsp = sp*.5;

vec2 map(vec3 p) {

    vec3 ip = p-vec3(0,0,ga2*sp);
    p.z+=4.25;
    p.z-=ga5*sp;

    zid = floor((ip.z+hsp)/sp)-.5;
    p.xy*=rot(zid*.198);

    p.z = mod(p.z+hsp,sp)-hsp;

    vec2 res = vec2(1e5,0.);

    float d = 1e5, t = 1e5;
    vec3 q = p;
    vec3 s = p-vec3(hlf,hlf,-hlf);
    vec2 cid = mod2(q,size);
    vec2 sid = mod2(s,size);
    float ht = hash21(cid+zid); 

    // build parts
    float thx = (.08+.055*sin(p.y*1.25) ) *size;
    if(ht>.5) q.x *= -1.;

    vec2 t2 = vec2(length(q.xy-hlf),length(q.xy+hlf));
    vec2 q2 = t2.x<t2.y?vec2(q.xy-hlf):vec2(q.xy+hlf);
    vec3 q3 = vec3(q2,q.z);
    
    float ti = torus(q3,vec2(hlf,thx));

    // truchet
    if(ti<res.x) {
        res = vec2(ti*.8,2.);
        hit = q;
        id = cid;
    }

    vec2 bfm = vec2(hlf*.825,hlf+.1);

    vec3 nq = cid!=vec2(0) ? q-vec3(0,0,.5) : q;
    float b = box(nq,bfm.xxy);

    b = min(box(q,bfm.yxx),b);
    b = min(box(q,bfm.xyx),b);

    float c = box(q,vec3(hlf)*.975);
    float di = max(c,-b);

    float sp = cap(s.xzy,vec2(size+.25,.675));
    di=smin(di,sp,.125);

    // box
    if(di<res.x) {
        res = vec2(di,3.);
        hit = q;
    }

    return res;
}

//Tetrahedron technique
//https://iquilezles.org/articles/normalsSDF
vec3 normal(vec3 p, float t, float mindist) {
    float e = mindist*t;
    vec2 h = vec2(1.,-1.)*.5773;
    return normalize( h.xyy*map( p + h.xyy*e ).x + 
					  h.yyx*map( p + h.yyx*e ).x + 
					  h.yxy*map( p + h.yxy*e ).x + 
					  h.xxx*map( p + h.xxx*e ).x );
}
   
vec3 hue(float t){ 
    return .65 + .35*cos(232.+PI2*t*(vec3(0.989+_uvc.x,0.977+_uvc.y,0.989)*vec3(0.925,0.722,0.435))); 
}

vec3 shade(vec3 p, vec3 rd, float d, float m, inout vec3 n, inout float fresnel, inout vec3 h) {

    n = normal(p,d,1.);

    vec3 l = normalize(vec3(5.,30,10)-p);
    float diff = clamp(dot(n,l),.1,.8);
    //@Shane specular highlights!
    float spec = pow(max(dot(reflect(l, n), rd ), .1), 32.);
        
    fresnel = pow(clamp(1.+dot(rd, n), 0., 1.), 9.);
    fresnel = mix(.0, .95, fresnel);
    // doin a second fill light
    l = normalize(vec3(-2.,020,-25)-p);
    float diff2 = clamp(dot(n,l),.0,.8)*.65;
    float spec2 = pow(max(dot(reflect(l, n), rd ), .2), 64.)*.65;

    h = vec3(.075);

    if(m==2.) {   
        hitPoint/=size;
        float dir = mod(tileId.x + tileId.y,2.) * 2. - 1.; 
        vec2 uv = hitPoint.xy-sign(hitPoint.x+hitPoint.y+.0001)*.5;

        float angle = atan(uv.x, uv.y);
        float a = sin( dir * angle * 14. + smoothTimeC * .5);
        a = abs(a)-(.145+0.35*sin(smoothTimeC*0.125)+0.25);a = abs(a)-.135;
        vec3 nz = hue(floor(a*3.5)+zIndex*2.);
        h = mix(nz, vec3(.00001), smoothstep(.1, .11, a));  
    }

    return h*diff*vec3(0.933,0.984,0.984)+spec*diff2+spec2;
}

vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 F = _xy;


    trs = rot(PI*4.5);
    turn = rot(-.1);
    // movment and timing
    float tt = bass_time*2.;
    tmod = mod(tt, 10.);
    float t9 = lsp(0.0, 9.0, tmod);
    t9 = eoc(t9);
    t9 = t9*t9*t9;  
    ga5 = t9;
    ga2 = t9+floor(tt*.1);
    
    vec3 C=vec3(.0);
    vec2 uv = (2.*F.xy-R.xy)/max(R.x,R.y);

    vec3 ro = vec3(0,0,1.),
         rd = normalize(vec3(uv,-1));
    
    float x = M.z<0.||M.xy == vec2(0) ? 0. :  (M.y/R.y * .25 - .125) * PI;
    float y = M.z<0.||M.xy == vec2(0) ? 0. : -(M.x/R.x * .5 - .25) * PI;

    mat2 rx = rot(x),ry = rot(y);
    ro.yz *= rx; ro.xz *= ry;
    rd.yz *= rx; rd.xz *= ry;

    vec3  p = ro + rd * .1;
    float atten = .95,k = 1.;
    float bt =2.;// cheat to break for reflection bounce
    float ct =6.;// cheat to break for refraction bounce
    // loop inspired/adapted from @blackle's 
    // marcher https://www.shadertoy.com/view/flsGDH
    for(int i=0;i<134;i++)
    {
        vec2 ray = map(p);
        float d = ray.x*.9;
        float m = ray.y;

        p += rd * d *k;
        
        if (d*d < 1e-8) {
            hitPoint=hit; tileId=id, zIndex= zid;
            
            vec3 n=vec3(0);
            float fresnel= 0.;
            vec3 h=vec3(0);

            C+=shade(p,rd,d,ray.y,n,fresnel,h)*atten;
            if(bt<1.||ct<1.)break;

            p += rd*.002;
            k = sign(map(p).x);

            vec3 rr = vec3(0);
            // use h value and material to decide to reflect or refract
            if(m==2.&&(h.r>.175||h.g>.175||h.b>.175)) {
                atten *=.25;
                rd=reflect(-rd,n);
                p+=n*.02;
                bt--;
            }else{
                ct--;
                atten *= .6;
                rr = refract(rd,n,.95);
                rd=mix(rd,rr,.95-fresnel);
            }
        }  
        if(distance(p,rd)>55.) { break; }
    }

    // Output to screen
    C = sqrt(smoothstep(0.,1.,C));
    O = vec4(C,1.);
	return O; 
 } 


//end


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}