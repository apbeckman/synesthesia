

			//******** BuffA Code Begins ********

/**
    License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
    
    Playing with some abstract forms and shapes
    Standard truchet patten in 3d (grid) and
    then warped nothing too exciting
    
    8/14/21 @byt3_m3chanic
*/

#define R           RENDERSIZE
#define T           TIME
#define M           _mouse

#define PI          3.14159265359
#define PI2         6.28318530718

#define MAX_DIST    20.00
#define MIN_DIST    0.001
#define SCALE       0.6500

//utils
float hash21(vec2 p){ return fract(sin(dot(p,vec2(26.34,45.32)))*4324.23); }
mat2 rot(float a){ return mat2(cos(a),sin(a),-sin(a),cos(a)); }

//globals
vec3 hit,hitP1,sid,id;
float glow,speed;
mat2 t90;

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
const float dbl = 1./size;
const float shorten = 1.26;  

//domain rep
vec3 drep(inout vec3 p) {
    vec3 id=floor((p+hlf)/size);
    p = mod(p+hlf,size)-hlf;
    return id;
}

vec2 map(vec3 q3, float gw){
    vec2 res = vec2(100.,0.);

    float k = 5.0/dot(q3,q3); 
    q3 *= k;

    q3.z += speed*.4;

    vec3 qm = q3;
    vec3 qd = q3+hlf;
    qd.zy*=t90;
    vec3 qid=drep(qm);
    vec3 did=drep(qd);
    
    float ht = hash21(qid.xy+qid.z);
    float hy = hash21(did.xz+did.y);
    
    // truchet build parts
    float thx = (.0750+.0475*sin(T+(q3.x+qid.z)*5.25) ) *size;
    float thz = (.0115+.0325*sin(T*3.+(q3.x+did.z)*6.5) ) *size;
    float thd = (.0515+.0325*sin(T*3.+(q3.x+did.z)*6.5) ) *size;
    if(ht>.5) qm.x *= -1.;
    if(hy>.5) qd.x *= -1.;
    
    // ring movement
    float dir = mod(did.x+did.y,2.)<.5? -1. : 1.;
    if(mod(did.z,2.)<1.)dir*=-1.;
    vec2 d2 = vec2(length(qd-hlf), length(qd+hlf));  
    
    vec2 pp = d2.x<d2.y? vec2(qd - hlf) : vec2(qd + hlf);
    pp *= rot(speed*dir);
    
    float a = atan(pp.y, pp.x);
    float amt = 6.;
    a = (floor(a/PI2*amt) + .5)/amt;

    vec2 qr = rot(-a*PI2)*pp; 
    qr.x -= hlf;
    vec3 npos = vec3(qr.x,qr.y,qd.z);

    //truchets1
    float t = truchet(qm,vec3(hlf,hlf,.0),vec2(hlf,thx));
    if(t<res.x) {
        sid = qid;
        hit = qm;
        res = vec2(t,2.);
    }
    //truchets2
    float d = truchet(qd,vec3(hlf,hlf,.0),vec2(hlf,thz));
    if(d<res.x) {
        sid = did;
        hit = qd;
        res = vec2(d,1.);
    }
    //rings
    float f = truchet(npos.xzy,vec3(.0,0,0),vec2(thd,.025));
    f=max(f,-(d-.01));
    if(f<res.x) {
        sid = did;
        hit = qd;
        res = vec2(f,3.);
    }
    //glows
    if (gw==1.) glow += smoothstep(.1,.25,.002/(.0165+f*f)); 

    float mul = 1.0/k;
    res.x = res.x * mul / shorten;
    
    return res;
}

// Tetrahedron technique @iq
// https://iquilezles.org/articles/normalsSDF
vec3 normal(vec3 p, float t)
{
    float e = MIN_DIST*t;
    vec2 h =vec2(1,-1)*.5773;
    vec3 n = h.xyy * map(p+h.xyy*e,0.).x+
             h.yyx * map(p+h.yyx*e,0.).x+
             h.yxy * map(p+h.yxy*e,0.).x+
             h.xxx * map(p+h.xxx*e,0.).x;
    return normalize(n);
}

vec3 FC= vec3(.001,.001,.001);
vec3 gcolor = vec3(0.145,0.659,0.914);
vec3 lpos = vec3(.0,.001,3.85);

vec3 render(inout vec3 ro, inout vec3 rd, inout vec3 ref, inout float d) {
    vec3 C = vec3(0);
    vec3 p = ro;
    float m = 0.;
    
    for(int i=0;i<120;i++) {
        p = ro + rd * d;
        vec2 ray = map(p,1.);
        if(abs(ray.x)<MIN_DIST*d||d>MAX_DIST)break;
        d += i<32? ray.x*.35: ray.x*.95;
        m  = ray.y;
    } 
    hitP1 = hit;
    id = sid;
    
    float alpha = 0.;
    if(d<MAX_DIST) {
    
        vec3 p = ro + rd * d;
        vec3 n = normal(p,d);
        vec3 l = normalize(lpos-p);
        
        vec3 h = vec3(.5);

        float diff = clamp(dot(n,l),0.,1.);
        float fresnel = pow(clamp(1.+dot(rd, n), 0., 1.), 12.);
        fresnel = mix(.01, 1., fresnel);

        float shdw = 1.0;
        for( float t=.01; t < 12.; ){
            float h = map(p + l*t,0.).x;
            if( h<MIN_DIST ) { shdw = 0.; break; }
            shdw = min(shdw, 24.*h/t);
            t += h;
            if( shdw<MIN_DIST || t>32. ) break;
        }
        diff = mix(diff,diff*shdw,.65);
        
        vec3 view = normalize(p - ro);
        vec3 ret = reflect(normalize(lpos), n);
        float spec =  0.5 * pow(max(dot(view, ret), 0.), 24.);

        if(m==1.) {
            vec3 hp = hitP1/(1./SCALE);
            float dir = mod(id.x + id.y,2.) * 2. - 1.;  
            vec2 uv = hp.xy-sign(hp.x+hp.y+.001)*.5;
            float angle = atan(uv.x, uv.y);
            float a = sin(dir*angle*4.+T*1.25);
            a = abs(abs(abs(a)-.5)-.24)-.124;
            h = mix(gcolor, gcolor*.2, smoothstep(.01, .02, a));   
            ref = h-fresnel;
        }

        if(m==2.) {
            vec3 hp = hitP1/(1./SCALE);
            float dir = mod(id.x + id.y,2.) * 2. - 1.;  
            if(mod(id.z,2.)<1.)dir*=-1.;
            vec2 uv = hp.xy-sign(hp.x+hp.y+.001)*.5;
            float angle = atan(uv.x, uv.y);
            float a = sin(dir*angle*8.+T*0.5);
            a = abs(abs(abs(abs(a)-.75)-.5)-.25)-.06;
            h = mix( vec3(.01), vec3(0.400,0.761,0.941) , smoothstep(.02, .01, a)); //vec3(0.400,0.761,0.941)
            ref = h-fresnel;
        }
        if(m==3.) {
            h = vec3(0.906,0.757,0.894); 
            ref = vec3(.0);
        }
        
        C = diff*h+spec;
        C = mix(FC.rgb,C,  exp(-.375*d*d*d));
    
        ro = p+n*.002;
        rd = reflect(rd,n);
        
    } else {
        C = FC.rgb;
    }
     
    return vec3(C);
}

void topLayer(inout vec3 C, vec2 uv, float alpha) 
{
    float px = fwidth(uv.x);
    float md = mod(T*.1,2.);
    float zw = md<1.? 2.: 1.25*sin(fract(T*.1));
    float d = length(uv)-zw;
    d=abs(d)-.002;
    d=smoothstep(px,-px,d);
    C =mix(C,gcolor,d);   
}

vec4 renderPassA() {
	vec4 O = vec4(0.0);
	vec2 F = _xy;

    // precal
    t90 = rot(90.*PI/180.);
    speed = T*.35;
    FC = (texture(BuffA,F.xy/R.xy).rgb)*.9;
    gcolor = vec3(0.145,0.659,0.914);
    //
    
    vec2 uv = (2.*F.xy-R.xy)/max(R.x,R.y);
    vec3 ro = vec3(0,0,3.25);
    vec3 rd = normalize(vec3(uv,-1));

    float x = M.xy == vec2(0) ? -.005 : -(M.y/R.y * .25 - .125) * PI;
    float y = M.xy == vec2(0) ? 0.015 : -(M.x/R.x * .25 - .125) * PI;
    
    mat2 rx = rot(x);
    mat2 ry = rot(y);
    
    ro.yz *= rx;
    rd.yz *= rx;
    ro.xz *= ry;
    rd.xz *= ry;
    
    vec3 C = vec3(0);
    vec3 ref=vec3(0);
    
    float glowMask =0.;
    float d =0.;

    C = render(ro, rd, ref, d).rgb;


    vec3 gcolor = vec3(0.145,0.659,0.914);

    glowMask = clamp(glow,.0,1.);
    C = mix(C,gcolor*glow,glowMask);
    
    float vw = .5+.5*sin(uv.x*.075+T*.3);
    float fade = clamp((d*.045)+vw,0.,1.);
    fade=abs(abs(fade)-.005)-.001;
    C = mix(C,gcolor,smoothstep(.11,.1,fade)); 
    
    topLayer(C,uv,d);
    C = clamp(C,vec3(0),vec3(1));
    O = vec4(C,1.0);
	return O; 
 } 


/**
    License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
    
    Abstract Truchet Inversion 
    8/14/21 @byt3_m3chanic

    (mouseable)
    
    suggest music to play (since soundcloud isnt working in shader)
    https://soundcloud.com/relapserecords/07-survive-low-fog
*/
float offset[3] = float[]( 0.01, 1.75, 4.95 );
float weight[3] = float[]( 0.2, 0.35, 0.075 );

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
  
	vec2 uv = fragCoord.xy/RENDERSIZE.xy;

    vec3 C = texture(BuffA, uv).rgb;
    vec3 R = texture(BuffA, uv).rgb * weight[0];
    for (int i=1; i<3; i++) {
      R += texture(BuffA, uv + vec2(offset[i])/RENDERSIZE.xy, 0.0).rgb * weight[i];
      R += texture(BuffA, uv - vec2(offset[i])/RENDERSIZE.xy, 0.0).rgb * weight[i];
    }
    
    // mask for effect and mixdown 
    float dt = distance(uv.xy,vec2(.5))*.32;
    dt = smoothstep(0.82,.965,1.-dt);
    vec3 Color = mix(R,C,dt);
    
    // output
    Color=pow(Color, vec3(.4545));
    fragColor = vec4(Color,1.);

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