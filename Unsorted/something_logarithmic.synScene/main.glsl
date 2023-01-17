

/**
    License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
        
    Something Logarithmic
    @byt3_m3chanic 08/24/21
    
    More ray marching a log spherial mapped animated
    truchet tile system. More learning shader and 
    something to keep in my toolbox for later use.
    
    fun to play with tile size (sz) and density
    
    found post online on Log Spherical Warping 
    https://www.osar.fr/notes/logspherical/
    
*/

#define R			RENDERSIZE
#define T			TIME
#define M			_mouse

#define PI2			6.28318530718
#define PI			3.14159265358

#define MAX_DIST 	85.
#define MIN_DIST	.00001

float hash21(vec2 p){ return fract(sin(dot(p,vec2(26.34,45.32)))*4324.23); }
mat2 rot(float a){ return mat2(cos(a),sin(a),-sin(a),cos(a)); }

float box(vec3 p, vec3 b) {
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

// constants 
float sz = 2.,hl = sz*.5;
vec2 boxSize = vec2(sz*.475,.12);

float density = (Truchet)*10.;
float dshalf = density/2.;

//globals and stuff
vec3 hit, hitPoint;
vec2 sid, cellId;
float shorten = 1., lpscale, movement, _movement,  trackspeed, time;
mat2 turn;

vec2 map(vec3 q){
    vec2 res = vec2(1e5,0.);

    vec3 p = q;
    p.xz*=turn;

    // log-spherical map
    float r = length(p);
    float mul = r/lpscale;
    p = vec3(log(r), acos(p.y / r ), atan(p.z, p.x));
    p *= lpscale;
    p -= vec3(_movement,dshalf,hl);
    
    vec2 id = floor((p.xz+hl)/sz);
    p.xz = mod(p.xz+hl,sz)-hl;

    float hs = hash21(id);
    float dir = mod(id.y + id.x,2.) <.5 ? -1. : 1.;

    if(hs>.5) p.x *= -1.;
    // get closest point and make vector for 1/4 torus
    vec2 d2 = vec2(length(p.xz - hl), length(p.xz + hl));  
    vec2 pp = d2.x<d2.y ? vec2(p.xz - hl) : vec2(p.xz + hl);
    float pth = abs(min(d2.x, d2.y) - hl);

    float tr = length(vec2(pth, p.y+.15))-.25;
    float bx = box(p,boxSize.xyx)-.0125;

    bx = max(bx,-tr);
    if(bx<res.x) {
        res = vec2(bx,3.);
        sid = id;
        hit = p;
    }
    
    //balls
    pp *= rot(trackspeed*dir);
     
    float amt = 4.;
    float dbl = 8.;
    
    float a = atan(pp.y, pp.x);
    // make id
    float ai = floor(dir*a/PI*amt);
    a = (floor(a/PI2*dbl) + .5)/dbl;
    vec2 qr = rot(-a*PI2)*pp; 
    qr.x -= hl;

    vec3 bq = vec3(qr.x, p.y+.125, qr.y);
   
    float sph = length(bq)-.1;
    
    if(sph<res.x) {
        res = vec2(sph,4.);
        sid = vec2(ai,dir);
        hit = bq;
    }

    res.x *= mul/shorten;
    return res;
} 
// @iq https://iquilezles.org/articles/normalsSDF
vec3 normal(vec3 p, float t){
    float e = t*MIN_DIST;
    vec2 h = vec2(1.,-1.)*.5773;
    return normalize( h.xyy*map( p + h.xyy*e).x + 
					  h.yyx*map( p + h.yyx*e).x + 
					  h.yxy*map( p + h.yxy*e).x + 
	 				  h.xxx*map( p + h.xxx*e).x );
}
//@iq https://iquilezles.org/articles/palettes
const vec3 c = vec3(0.959,0.970,0.989),
           d = vec3(0.067,0.812,0.910);
vec3 hue(float t){ 
    return .45 + .45*cos(PI2*t*(c+d) ); 
}

vec4 FC= vec4(0.000,0.000,0.000,0.);
vec3 render(inout vec3 ro, inout vec3 rd, inout vec3 ref, inout float d) {
    
    vec3 C = vec3(0);
    vec3 p = ro;
    float m = 0.;
    
    // marcher
    for(int i=0;i<175;i++) {
        p = ro + rd * d;
        vec2 ray = map(p);
        if(abs(ray.x)<MIN_DIST*d||d>MAX_DIST)break;
        d += i<128? ray.x*.5: ray.x * .75;
        m  = ray.y;
    } 
    
    cellId = sid;
    hitPoint = hit;
    
    float alpha = 0.;
    if(d<MAX_DIST) {
    
        vec3 p = ro + rd * d;
        vec3 n = normal(p,d);
        vec3 lpos = vec3(2.0,5.0,3.85);
        vec3 l = normalize(lpos-p);
        
        vec3 h = vec3(.5);

        float diff = clamp(dot(n,l),.03,1.);
        float fresnel = pow(clamp(1.+dot(rd, n), 0., 1.), 9.);
        fresnel = mix(.01, .9, fresnel);

        float shdw = 1.;
        float t=.0;
        for( float i=.01; i < 32.;i++ ){
            float h = map(p + l*t).x;
            if( h<MIN_DIST ) { shdw = 0.; break; }
            shdw = min(shdw, 24.*h/t);
            t += h;
            if( shdw<MIN_DIST || t>32. ) break;
        }
        diff = mix(diff,diff*shdw,.65);
        
        vec3 view = normalize(p - ro);
        vec3 ret = reflect(normalize(lpos), n);
        float spec =  0.5 * pow(max(dot(view, ret), 0.), 24.);
        
        h = vec3(1);

        if(m==3.) {
            vec3 uv = hitPoint;
            float px  = fwidth(uv.x*1.75);
            
            vec2 id = cellId;
            vec2 grid = uv.xz;
            float hs = hash21(id);
            float chk = mod(id.y + id.x,2.) * 2. - 1.;

            vec2 d2 = vec2(length(grid-hl), length(grid+hl));
            vec2 gx = d2.x<d2.y? vec2(grid-hl) : vec2(grid+hl);
            float pth = abs(min(d2.x,d2.y)-hl);
            vec2 vuv = vec2(pth, uv.z);
            
            float back = length(gx)-hl;
            back=(chk>0.^^ hs>.5) ? smoothstep(-px,px,back) : smoothstep(px,-px,back);
            
            vec2 pid = floor(grid*4.);
            vec2 puv = fract(grid*4.)-.5;
            float fs = hash21(pid);
            if(fs>.5)puv.x*=-1.;
            vec2 d5 = vec2(length(puv-.5), length(puv+.5));
            vec2 kx = d5.x<d5.y? vec2(puv-.5) : vec2(puv+.5);
            
            float ptrn = length(kx)-.5;
            ptrn = smoothstep(px,-px,abs(abs(abs(abs(ptrn)-.1)-.1)-.1)-.05);
            vec3 c2 = hue((10.+cellId.x*.05) );
            vec3 c3 = hue((1. -cellId.x*.075));
            //c2.xy *= PI*_uvc*0.5;
            vec3 c4 = hue((cellId.x)*.05);
            
            h = mix(c3,mix(c3, c2*PI,ptrn),back);
            
            float circle4;
            float circle2 = length(gx)-hl;
            float circle3 = smoothstep(px,-px,abs(circle2)-.13);
            circle4 = smoothstep(px,-px,abs(abs(abs(circle2)-.2)-.065)-.012);
            circle2 = smoothstep(px,-px,abs(abs(circle2)-.125)-.15);
            
            h=mix(h,vec3(.4),circle2);
            
            h=mix(h,c4,circle3);
            h=mix(h,vec3(0),circle4);
 
            ref = mix(vec3(0),h-fresnel,circle2);
        }

        if(m==4.) {
            h=mod(cellId.x+25.,2.)==0.?vec3(.03):vec3(.9);
            ref = h-fresnel;
        }
        
        C = diff*h+spec;
        
        ro = p+n*.001;
        rd = reflect(rd,n);
        
    } else {
        C = FC.rgb;
    }
    return C;
}

vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 F = _xy;

    // precal
    time = TIME;
  
    lpscale = (floor((Truchet)*density)/(PI));
    turn = rot(spin_time*5.*PI/180.);

    trackspeed = .5*mid_time;
    movement = .95*smoothTimeC*lpscale*.125;
    _movement = bass_time*lpscale*.25;
    
    vec2 uv = (2.*(F.xy)-R.xy)/max(R.x,R.y);
    uv.xy += _uvc*PI*FOV;
    vec3 ro = vec3(0,0,9.);
    vec3 rd = normalize(vec3(uv,-1));
    
    float y = M.xy == vec2(0) ? 0. : (M.x/R.x * 1. - .5) * PI;
    
    mat2 rx = rot(.8+.2*(2.0*Look));
    mat2 ry = rot(y);
    
    ro.yz *= rx;ro.xz *= ry;
    rd.yz *= rx;rd.xz *= ry;
    
    vec3 C = vec3(0);
    vec3 ref=vec3(0);
    vec3 fill=vec3(1.);
    
    float d =0.;
    //@BigWIngs - reflection loop
    for(float i=0.; i<2.; i++) {
        vec3 pass = render(ro, rd, ref, d);
        C += pass*fill;
        fill*=ref;
        if(i==0.) FC = vec4(FC.rgb,exp(-.00015*d*d*d));
    }

    C = mix(C,FC.rgb,1.-FC.w);
    C = clamp(C,vec3(MIN_DIST),vec3(1));
    C = pow(C, vec3(.4545));
    O = vec4(C,1.0);
	return O; 
 } 

// end

vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}