vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


/**
        UV Mapping a Truchet Tile Set
        @byt3_m3chanic 8/16/21
        
        
        Thank you @Fabrice for the knowledge and math
        Started as an experiment - how can I do this.
        https://www.shadertoy.com/view/sdtGRn
        
        And finally ended up here, it's pretty tricky as
        you have to get the closest arc and use that in
        the mapping formula.

*/


#define R           RENDERSIZE
#define T           (smooth_basstime*0.5)
#define M           iMouse

#define PI          3.14159265359
#define PI2         6.28318530718

#define MAX_DIST    20.00
#define MIN_DIST    0.001
#define SCALE       0.7500

//utils
float hash21(vec2 p){ return fract(sin(dot(p,vec2(26.34,45.32)))*4324.23); }
mat2 rot(float a){ return mat2(cos(a),sin(a),-sin(a),cos(a)); }
//globals
vec3 hit,hitP1,sid,id;
float speed,sdir,hitD,chx,checker;
mat2 t90;
float box( vec3 p, vec3 b ){
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}
//@iq torus sdf
float torus( vec3 p, vec2 t ) {
  vec2 q = vec2(length(p.xy)-t.x,p.z);
  return length(q)-t.y;
}
//make tile piece
float truchet(vec3 p,vec3 x, vec2 r) {
    return min(torus(p-x,r),torus(p+x,r));
}
//const
const float size = 1./SCALE;
const float hlf = size/2.;
const float shorten = 1.26;   
//domain rep
vec3 drep(inout vec3 p) {
    vec3 id=floor((p+hlf)/size);
    p = mod(p+hlf,size)-hlf;
    return id;
}
vec2 drep(inout vec2 p) {
    vec2 id=floor((p+hlf)/size);
    p = mod(p+hlf,size)-hlf;
    return id;
}
vec2 map(vec3 q3){
    vec2 res = vec2(100.,0.);

    float k = 5.0/dot(q3,q3); 
    q3 *= k;

    q3.z += speed;

    vec3 qm = q3;
    vec3 qd = q3+hlf;
    qd.xz*=t90;
    vec3 qid=drep(qm);
    vec3 did=drep(qd);
    
    float ht = hash21(qid.xy+qid.z);
    float hy = hash21(did.xz+did.y);
    
    float chk1 = mod(qid.y + qid.x,2.) * 2. - 1.;
    float chk2 = mod(did.y + did.x,2.) * 2. - 1.;

    // truchet build parts
    float thx = .115;
    float thz = .200;

    if(ht>.5) qm.x *= -1.;
    if(hy>.5) qd.x *= -1.;

    float t = truchet(qm,vec3(hlf,hlf,.0),vec2(hlf,thx));
    if(t<res.x) {
        sid = qid;
        hit = qm;
        chx = chk1;
        sdir = ht>.5 ? -1. : 1.;
        res = vec2(t,2.);
    }

    float d = truchet(qd,vec3(hlf,hlf,.0),vec2(hlf,thz));
    if(d<res.x) {

        sid = did;
        hit = qd;
        chx = chk2;
        sdir = hy>.5 ? -1. : 1.;
        res = vec2(d,1.);
    }

    float mul = 1.0/k;
    res.x = res.x * mul / shorten;
    
    return res;
}

// Tetrahedron technique @iq
// https://www.iquilezles.org/www/articles/normalsSDF
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
    vec3 d = vec3(0.220,0.961,0.875);
    return .375 + .375*cos(PI2*t*(vec3(.985,.98,.95)+d)); 
}

float gear(vec2 p, float radius)
{
    //length of cog
    float sp = floor(radius*PI2)*2.;
    float gs = length(p.xy)-radius;
    float at = atan(p.y,p.x);
    //gear teeth
    float gw = abs(sin(at*sp)*.15);
    gs +=smoothstep(.05,.5,gw);

    gs=max(gs, -(length(p.xy)-(radius*.45)) );
    return gs;
}

vec4 FC= vec4(.001,.001,.001,0.);
vec3 lpos = vec3(-hlf,hlf,3.85);

vec4 render(inout vec3 ro, inout vec3 rd, inout vec3 ref, bool last, inout float d) {

    vec3 C = vec3(0);
    vec3 p = ro;
    float m = 0.;
    
    for(int i=0;i<150;i++) {
        p = ro + rd * d;
        vec2 ray = map(p);
        if(abs(ray.x)<MIN_DIST*d||d>MAX_DIST)break;
        d += i<64? ray.x*.35: ray.x;
        m  = ray.y;
    } 
    // save globals
    hitP1 = hit;
    id = sid;
    hitD = sdir;
    checker=chx;
    
    float alpha = 0.;
    if(d<MAX_DIST) {
    
        vec3 p = ro + rd * d;
        vec3 n = normal(p,d);
        vec3 l = normalize(lpos-p);
        
        vec3 h = vec3(.05);

        float diff = clamp(dot(n,l),0.,1.);
        float bounce = clamp(dot(n,vec3(0.,-1.,0.)), 0.,1.);
        
        float fresnel = pow(clamp(1.+dot(rd, n), 0., 1.), 5.);
        fresnel = mix(.01, .7, fresnel);

        float shdw = 1.0;
        for( float t=.01; t < 12.; ){
            float h = map(p + l*t).x;
            if( h<MIN_DIST ) { shdw = 0.; break; }
            shdw = min(shdw, 24.*h/t);
            t += h;
            if( shdw<MIN_DIST || t>32. ) break;
        }
        
        diff += bounce;
        diff = mix(diff,diff*shdw,.65);
  
        
        vec3 view = normalize(p - ro);
        vec3 ret = reflect(normalize(lpos), n);
        float spec =  0.5 * pow(max(dot(view, ret), 0.), (m==2.||m==4.)?24.:64.);

        // uv mapping stuff here
        
        if(m==2.) {
           vec3 hp = hitP1*hitD;
            
            // get closest arc
            vec2 d3 = vec2(length(hp-hlf), length(hp+hlf));
            vec3 g3 = d3.x<d3.y? vec3(hp-hlf) : vec3(hp+hlf);
            
            //angle for the tube
            float angle = atan(g3.y,g3.x)/PI2;
            //angle for the arc/truchet path
            float gz =  atan( hp.z,  length(g3.yx)-hlf ) / PI2;
            //make the uv
            vec2 uv = vec2(angle,gz);

            if(hitD<1. ^^ checker>0.) uv.y*=-1.;
            // make truchet design down
            float px  = .0125;
            
            vec2 scaler = vec2(28.,6.);
            vec2 grid = fract(uv.xy*scaler)-.5;
            vec2 cid   = floor(uv.xy*scaler);

            if(hash21(cid)<.5) grid.x*=-1.;
            
            h = vec3(.6); 

            vec2 d2 = vec2(length(grid-.5), length(grid+.5));
            vec2 gx = d2.x<d2.y? vec2(grid-.5) : vec2(grid+.5);

            float circle = length(gx)-.5;
            float center =smoothstep(.03-px,px,abs(abs(abs(circle)-.2)-.1)-.025);
            h = mix(vec3(0),h, center);
           
            ref = vec3(clamp(1.-center,0.,1.))-fresnel;
        }

        if(m==1.) {
            vec3 hp = hitP1*hitD;
            
            // get closest arc
            vec2 d3 = vec2(length(hp-hlf), length(hp+hlf));
            vec3 g3 = d3.x<d3.y? vec3(hp-hlf) : vec3(hp+hlf);
    
            //angle for the tube
            float angle = atan(g3.y,g3.x)/PI2;
            //angle for the arc/truchet path
            float gz =  atan( hp.z,  length(g3.yx)-hlf ) / PI2;
            //make the uv
            vec2 uv = vec2(angle,gz);

            if(hitD<1. ^^ checker>0.) uv.y*=-1.;
            // make truchet design down
            float px  = .0125;
           
            vec2 scaler = vec2(28.,10.);
            vec2 grid = fract(uv.xy*scaler)-.5;
            vec2 cid   = floor(uv.xy*scaler);
            float hs = hash21(cid);
            if(hs<.5) grid.x*=-1.;

            vec2 d2 = vec2(length(grid-.5), length(grid+.5));
            vec2 gx = d2.x<d2.y? vec2(grid-.5) : vec2(grid+.5);

            float circle = length(gx)-.5;
            float center =smoothstep(.03-px,px,abs(circle)-.15);
            h = mix(hue(length(p.zy*.3)*3.) ,hue(length(p.zx*.5)*2.), center);
            
            //debug with london texture 
            //h=texture(iChannel0,uv*vec2(8.,4.)).rgb;
            
            // super freakout animation
            float chk = mod(cid.y + cid.x,2.) * 2. - 1.;
            vec2 arc = grid-sign(grid.x+grid.y+.001)*.5;
            float angle2 = atan(arc.x, arc.y);
            float width = .2;
            float d = length(arc);
            // coord checker
            float tm = T*.25;
            vec2 tuv = vec2(
                fract(1.*chk*angle2/1.57+tm),
                (d-(.5-width))/(2.*width)*2.
            );
            tuv.y-=.5;
            vec2 tid = vec2(
                floor(1.*chk*angle2/1.57+tm),
                floor(d-(.5-width))/(2.*width)
            );

            tuv.xy*=vec2(2.,.5);
            tuv.x=mod(tuv.x+.5,1.)-.5;
            // float ddt = length(tuv.xy-vec2(0,.25))-.25;
            
            // new gear spin thing
            // comment out to change back to dots - and use above.
            vec2 gvec = tuv.xy-vec2(0,.25);
            float dir = (chk>0.^^ hs>.5) ? -1.:1.;
            gvec*= rot( (script_high_time*0.75) * dir);
            float ddt = gear(gvec,.45);
            // end
            
            ddt = smoothstep(-px,px,min(ddt,center));
            h = mix(h,vec3(0.),ddt);
            
            
            ref = vec3(clamp(1.-center,0.,1.))-fresnel;
        }

        C = diff*h+spec;
        if(last) C = mix(FC.rgb,C,  exp(-.05*d*d*d));
    
        ro = p+n*.002;
        rd = reflect(rd,n);
    } else {
        C = FC.rgb;
    }
     
    return vec4(C,alpha);
}

vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 F = _xy;

    // precal
    t90 = rot(90.*PI/180.);
    speed = T*.225;
    //
    vec2 uv = (2.*F.xy-R.xy)/max(R.x,R.y);
    vec3 ro = vec3(0,0,3.5);
    vec3 rd = normalize(vec3(uv,-1));
    //mouse rotation
    float x = -.305;
    float y = -.750;
    
    if(M.z>0.){
        x = -(M.y/R.y * .5 - .25) * PI;
        y = -(M.x/R.x * .5 - .25) * PI;
    }
    
    mat2 rx = rot(x);
    mat2 ry = rot(y);
    
    ro.yz *= rx;
    rd.yz *= rx;
    ro.xz *= ry;
    rd.xz *= ry;
    
    // pre render
    vec3 C = vec3(0);
    vec3 ref=vec3(0);
    vec3 fil=vec3(1.);
    
    float d =0.;
    float bounces = 2.;
    for(float i=0.; i<bounces; i++) {
        vec4 pass = render(ro, rd, ref, i==bounces-1., d);
        C += pass.rgb*fil;
        fil*=ref;
        if(i==0.) FC = vec4(FC.rgb,exp(-.145*d*d*d));
    }

    //mixdown
    C = mix(C,FC.rgb,1.-FC.w);
    C = clamp(C,vec3(0),vec3(1));
    C = pow(C, vec3(.4545));
    O = vec4(C,1.0);
	return O; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}