

/** 
    License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
    
    🧊 Cube Stacks - its really just one layer packed with waves. 
    4/3/22 | @byt3_m3chanic
    
    Mostly a homage to @Shane's techniques for multi-tap grids and cube mapping.
    Was just a fun exercise in playing with cube maps
    
*/

#define R           RENDERSIZE
#define T           TIME
#define M           _mouse

#define PI          3.14159265359
#define PI2         6.28318530718

#define MIN_DIST    .00001
#define MAX_DIST    100.

#define ZERO (min(FRAMECOUNT,0))
// AA 1 = Off/ 2+ = On
#define AA 1

mat2 rot(float a) { 
    return mat2(cos(a),sin(a),-sin(a),cos(a)); 
}
float hash21(vec2 a){ 
    return fract(sin(dot(a, vec2(27.609, 57.583)))*43768.5453); 
}
//@iq ALL THE SDFS
float box(vec3 p, vec3 b){
    vec3 d = abs(p)-b;
    return length(max(d,0.))+min(max(d.x,max(d.y,d.z)),.0);
}
float box(vec2 p, vec2 b){
    vec2 d = abs(p)-b;
    return length(max(d,0.))+min(max(d.x,d.y),.0);
}

//globals s = saved / g = global
vec3 s_hit,g_hit;
vec2 s_id,g_id;
float travelspeed = 0.,s_hash,g_hash;

//scales for all things
const float scale = 1./.345;
const float scale_h = scale*.5;
const vec2 s = vec2(scale)*2.;

//positions and rotations
const vec2 pos = vec2(.5,-.5);
const vec2[4] ps4 = vec2[4](pos.yx,pos.xx,pos.xy,pos.yy);
const vec4 rotations = vec4(0.,1.570,3.141,4.712);

// block map
// 4 tap technique @Shane
vec2 map(vec3 q3){
    vec2 res = vec2(1e5,0.);

    q3.xz += vec2(travelspeed,-travelspeed);
    
    vec2 p,
         ip,
         id = vec2(0),
         ct = vec2(0);

    float t=1e5, y=1e5;

    for(int i =0; i<4; i++){
        ct = ps4[i]/2. -  ps4[0]/2.;
        p = q3.xz - ct*s;
        ip = floor(p/s) + .5;
        p -= (ip)*s;
        vec2 idi = (ip + ct)*s;

        float px = sin(idi.x*2.)*.4;
        float py = cos(idi.y*2.)*.4;
        float ah = -floor((px+py)*scale)*scale;
  
        vec3 q = vec3(p.x,q3.y+9.-ah,p.y);

        // just a box
        float b = box(q,vec3(scale*.4645))-.125;
        if(b<t) {
            t = b;
            s_id = idi;
            s_hit = q;
            s_hash= px+py;
        }
  
    }
    
	if(t<res.x) res = vec2(t,2.);

    return res;
}

vec2 marcher(vec3 ro, vec3 rd, inout vec3 p, inout bool hit, int steps) {
    hit = false; float d=0., m = 0.;
    for(int i=0;i<steps;i++) {
        p = ro + rd * d;
        vec2 t = map(p);
        if(abs(t.x)<d*MIN_DIST) hit = true;
        d += i<32? t.x*.35:t.x;
        m  = t.y;
        if(d>MAX_DIST) break;
    } 
    return vec2(d,m);
}

vec3 getFace(int face, vec3 p, inout vec3 mask) {
    vec2 uv = p.xy;
    float px = fwidth(uv.x)*.5;
    uv = uv/scale_h;

    vec3 color = vec3(.5);

    
    vec2 grid_uv = fract(uv)-.5;
    vec2 grid_id = floor(uv);
    
    float rnd = hash21(grid_id+g_id);
    if(rnd>.5) grid_uv.x*=-1.;

    vec2 d2 = vec2(length(grid_uv-.5),length(grid_uv+.5));
    vec2 gx = d2.x<d2.y? vec2(grid_uv-.5):vec2(grid_uv+.5);

    float fm = box(uv,vec2(scale_h*.5))-.075;
    float fx = fm;
    fm=abs(fm)-.025;

    fm=smoothstep(px,-px,fm);
    fx=smoothstep(px,-px,fx);
    color=mix(color,vec3(.05),fx);

    float circle = length(gx)-.5;
    circle=abs(circle)-.125;

    mask=vec3(circle,fm,fx);
    circle=smoothstep(px,-px,circle);
    color = mix(color, vec3(.15),circle);

    return vec3(1.-color);
}

// based on bmp mapping from @Passion
// https://www.shadertoy.com/view/ld3yDn
vec3 doBumpMap( vec3 p, vec3 n, float bf, float per, int face, inout vec3 mask){
    vec2 e = vec2(per*MIN_DIST, 0);  
    mat3 m = mat3( 
        getFace(face, p - e.xyy,mask), 
        getFace(face, p - e.yxy,mask), 
        getFace(face, p - e.yyx,mask)
    );
    vec3 g = vec3(.299,.587,.114) * m; 
    g = (g-dot(getFace(face, p,mask),vec3(.299,.587,.114)))/e.x; 
    g -= n*dot(n,g);  
    return normalize(n+g*bf);
}

vec3 normal(vec3 p, float t) {
    t*=MIN_DIST;
    float d = map(p).x;
    vec2 e = vec2(t,0);
    vec3 n = d - vec3(
        map(p-e.xyy).x,
        map(p-e.yxy).x,
        map(p-e.yyx).x
    );
    return normalize(n);
}

//@Shane AO
float calcAO(in vec3 p, in vec3 n){
    float sca = 4.5, occ = 0.;
    for( int i = 0; i<5; i++ ){
        float hr = float(i + 1)*.0029/.27; 
        float d = map(p + n*hr).x;
        occ += (hr - d)*sca;
    }
    return clamp(1. - occ, 0., 1.);
}

const vec3 FC = vec3(0.9,0.2,01.752);

vec3 hue(float t){ return .54 + .45*cos(PI2*t*(vec3(.95,.97,.98)+vec3(0.471,0.878,0.063))); }

vec3 render(inout vec3 ro, inout vec3 rd, inout vec3 ref, int bnc, inout float d) {
        
    vec3 RC=vec3(0);
    vec3 p = ro;
    float m = 0.;
    bool hit = false;

    vec2 ray = marcher(ro,rd,p, hit, 85);
    d = ray.x;
    m = ray.y;
    g_hit = s_hit;
    g_id = s_id;
    g_hash=s_hash;
    
    if(hit)
    {
        vec3 n = normal(p,d);
        vec3 lpos =  vec3(-18.*sin(T*.1),11,8.*cos(T*.1));
        vec3 l = normalize(lpos);
        
        float rhash= hash21(g_id.yx);
 
        vec3 h = vec3(.45);
        vec3 clr = hue((rhash*3.44)*.46);
        
        vec3 hp = g_hit;
        vec3 tn = n;
        
        vec3 cuv;
        int face;
        
        if(m==2.){
            tn = n;
            //@Shane https://www.shadertoy.com/view/3sVBDd
            //finding the face of a cube using normal
            vec3 aN = abs(tn);
            ivec3 idF = ivec3(tn.x<-.25? 0 : 5, tn.y<-.25? 1 : 4, tn.z<-.25? 2 : 3);
            face = aN.x>.5? idF.x : aN.y>.5? idF.y : idF.z;

            // set coords
            if(face==0) cuv = hp.yzx;
            if(face==1) cuv = hp.zxy;
            if(face==2) cuv = hp.xyz;
            if(face==3) cuv = hp.xyz;
            if(face==4) cuv = hp.xzy;
            if(face==5) cuv = hp.yzx;

            // get bump map surface
            vec3 mask = vec3(0);

            n=doBumpMap(cuv,n,.5,d,face,mask);
            vec3 txt = vec3(.7);
            vec3 dtr = texture(image6,(cuv.xy*.5)-.5).rgb*vec3(.815,.812,.839);
            vec3 otr = texture(image1,cuv.xy*.145).rgb;
            vec3 bse = rhash>.8?clr:dtr;
            vec3 sqr = vec3(.001);
            vec3 pth = rhash>.8?clr*dtr:vec3(0.902);
            
            float mx =abs(mask.x+.05)-.05;
            mask.x=smoothstep(.011,.01,mask.x);
            mx=smoothstep(.011,.01,mx);
            
            h=mix(otr,bse,mask.z>1.5?1.:0.);
            h=mix(h,pth,mask.x>.5?1.:0.);
            h=mix(h,sqr,min(1.-mask.x,mask.y)>.5?.8:0.);
     
            ref =min(1.-mask.x,mask.z)>.5? h:vec3(.0001);
        }
        
        float diff = clamp(dot(n,l),.25,.95);
        float ao = calcAO(p,n);
        float shdw = 1.;
        //@Shane - softshadows
        for( float t=.01;t<10.; ) {
            float h = map(p + l*t).x;
            if(h<MIN_DIST) { shdw = 0.; break; }
            shdw = min(shdw, 12.*h/t);
            t += h;
            if(shdw<MIN_DIST||t>12.)break;
        }
        diff = mix(diff,diff*shdw,.75);
        
        float spec = .65 * pow(max(dot(normalize(p-ro),reflect(l,n)),.2),24.);

        RC = (h*diff) + min(spec,shdw)*ao;

        ro = p+n*.05;
        rd = reflect(rd,n);

    } else {
        RC = vec3(.05);
    } 

    return RC;
}

vec3 render( vec2 uv, in vec2 F ) {

    vec3 ro = vec3(0,0,10);
    vec3 rd = normalize(vec3(uv, -1.));

    // mouse //
    float x = M.xy==vec2(0)||M.z<0. ? .07+.02*sin(T*.1) : -(M.y/R.y*.25-.125)*PI;
    float y = M.xy==vec2(0)||M.z<0. ? 0. : -(M.x/R.x*2.-1.)*PI;
    
    mat2 rx =rot(-.68+x);
    mat2 ry =rot(.58+y);
    rd.xz = _rotate(rd.xz, LookXY.x*PI);
    //ro.xz = _rotate(ro.xz, LookXY.x*PI);
    //ro.yz = _rotate(ro.yz, LookXY.y*PI);
    rd.yz = _rotate(rd.yz, LookXY.y*PI);
    //ro.zy*=rx;rd.zy*=rx;
//    /ro.xz*=ry;rd.xz*=ry;
    ro.y += 1.0;
    vec3 C=vec3(0),RC=vec3(0),ref=vec3(0),fill=vec3(1),p = ro;
    C.xy += _uvc*0.1+0.025*_rotate(_uvc, smoothTimeC);
    float d = 0., fA = 0.;
    
    // reflection loop + stop unroll
    for(int i = 0; i < 2 + min(FRAMECOUNT, 0); i++){
        RC = render(ro,rd,ref,2-i,d);
        C = mix(C,clamp(C+RC,vec3(0),vec3(1)),fill);
        fill *= ref; 
        if(i==0)fA=d;
    }
    // fog dusting
    C = mix(C,FC, 1.-exp(-.000000125*fA*fA*fA));
    return C;
}


// AA from @iq https://www.shadertoy.com/view/3lsSzf
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 F = _xy;

    travelspeed = T;
    
    vec3 C = vec3(0.0);
#if AA>1
    for( int m=ZERO; m<AA; m++ )
    for( int n=ZERO; n<AA; n++ )
    {
        // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 uv = (-R.xy + 2.0*(F+o))/max(R.x,R.y);
#else    
        vec2 uv = (-R.xy + 2.0*F)/max(R.x,R.y);
#endif

        vec3 color = render(uv,F);
        // compress        
        color = 1.35*color/(1.0+color);
        // gamma
        color = pow( color, vec3(0.4545) );

        C += color;
#if AA>1
    }
    C /= float(AA*AA);
#endif
    // Output to screen
    O = vec4(C,1.0);
	return O; 
 } 

//end

vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}