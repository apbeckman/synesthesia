

			//******** BuffA Code Begins ********

/**
    Analogous Oceans | some abstract marching
    @byt3_m3chanic | 07/10/21
*/

#define R   RENDERSIZE
#define M   _mouse
#define T   TIME
#define PI  3.14159265359
#define PI2 6.28318530718

#define MAX_DIST    100.
#define MIN_DIST    .0001

float hash21(vec2 p){ return fract(sin(dot(p,vec2(26.34,45.32)))*4324.23); }
mat2 rot(float a){ return mat2(cos(a),sin(a),-sin(a),cos(a)); }

float noise (in vec2 _st) {
    vec2 i = floor(_st);
    vec2 f = fract(_st);

    // Four corners in 2D of a tile
    float a = hash21(i);
    float b = hash21(i + vec2(1.0, 0.0));
    float c = hash21(i + vec2(0.0, 1.0));
    float d = hash21(i + vec2(1.0, 1.0));

    vec2 u = f * f * (3. - 2.0 * f);
    return mix(a, b, u.x) +
            (c - a)* u.y * (1.0 - u.x) +
            (d - b) * u.x * u.y;
}

//@iq
float box(vec3 p, vec3 s) {
    p = abs(p)-s;
    return length(max(p, 0.))+min(max(p.x, max(p.y, p.z)), 0.)-.025;
}

vec3 hit=vec3(0),hitPoint,gid,sid;
vec3 speed = vec3(0);
float glow =0.0,wtime;
const float size = .75;
const float hlf = size/2.;
const float dbl = size*2.;
          
vec2 map(in vec3 p, float sg) {
    vec2 res = vec2(1e5,0.);
    p += speed;

    float id,nz=0.;
    vec3 q;
    for(int i = 0; i<2; i++)
    {
        float cnt = i<1 ? size : dbl;
        q = vec3(p.x-cnt,p.yz);
        id = floor(q.x/dbl) + .5;
        q.x -= (id)*dbl;
        float qf = (id)*dbl + cnt;

        nz = noise(vec2(qf*.22+wtime,q.z*.22)) * 2.;
        q.y += nz+nz*sin(q.z*.55);
        q.y += nz+nz*cos(q.x*.35);
        float d = box(q+vec3(0.,3.25,0.),vec3(.125,.25,50.));
        d=max(d,-(p.y+6.4));
        
        if(d<res.x){
            res = vec2(d,1.);
            hitPoint = q+vec3(qf,qf,float(i));
            gid = vec3(qf,nz,float(i));
            if(sg==1. && mod(qf,6.)==0.) glow += .001/(.0005+d*d);
        }
        
    }

    float fl = p.y+6.5;
    if(fl<res.x){
        res=vec2(fl,2.);
        hitPoint=p;
    }
    return res;
}
//Tetrahedron technique
//https://iquilezles.org/articles/normalsSDF
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

vec2 marcher(vec3 ro, vec3 rd, int maxsteps, float sg){
	float d = 0.;
    float m = 0.;
    for(int i=0;i<maxsteps;i++){
    	vec2 ray = map(ro + rd * d, sg);
        if(ray.x<MIN_DIST*d||d>MAX_DIST) break;
        d += i<42?ray.x*.4:ray.x;
        m  = ray.y;
    }
	return vec2(d,m);
}

// Tri-Planar blending function. GPU Gems 3 - Ryan Geiss:
vec3 tex3D(sampler2D t, in vec3 p, in vec3 n ){
    n = max(abs(n), MIN_DIST);
    n /= dot(n, vec3(1));
	vec3 tx = texture(t, p.yz).xyz;
    vec3 ty = texture(t, p.zx).xyz;
    vec3 tz = texture(t, p.xy).xyz;
    return (tx*tx*n.x + ty*ty*n.y + tz*tz*n.z);
}

vec4 FC = vec4(0.019,0.019,0.019,0.);

vec3 hue(float t){ 
    const vec3 c = vec3(0.757,0.110,0.992);
    return .65 + .45*cos(13.+PI2*t*(c*vec3(0.959,0.970,0.989))); 
}

vec4 render(inout vec3 ro, inout vec3 rd, inout vec3 ref, bool last, inout float d, vec2 uv) {

    vec3 C = vec3(0);
    vec2 ray = marcher(ro,rd,150, 1.);
    hitPoint = hit;
    sid=gid;
    d = ray.x;
    float m = ray.y;
    float alpha = 0.;
    
    if(d<MAX_DIST)
    {
        vec3 p = ro + rd * d;
        vec3 n = normal(p,d);
        vec3 lpos =vec3(1.,11,-3.);
        vec3 l = normalize(lpos-p);
        
        vec3 h = vec3(.5);
        vec3 hp = hitPoint;

        float diff = clamp(dot(n,l),0.,1.);
        float fresnel = pow(clamp(1.+dot(rd, n), 0., 1.), 9.);
        fresnel = mix(.01, .7, fresnel);

        float shdw = 1.0;
        for( float t=.01; t < 10.; )
        {
            float h = map(p + l*t,0.).x;
            if( h<MIN_DIST ) { shdw = 0.; break; }
            shdw = min(shdw, 18.*h/t);
            t += h;
            if( shdw<MIN_DIST || t>32. ) break;
        }

        diff = mix(diff,diff*shdw,.75);

        vec3 view = normalize(p - ro);
        vec3 ret = reflect(normalize(lpos), n);
        float spec = 0.3 * pow(max(dot(view, ret), 0.), 20.);
        
        vec3 clr = hue(sid.x*.1+sid.y)*.35;
        // materials
        if(m==1.){
            h= mod(sid.x,6.)==0. ? clr : tex3D(image3,hp*.1,n).rgb*clr;
            
            ref = h-fresnel;
            C = (diff*h)+spec;
        }
        if(m==2.){
            h=texture(image3,hp.xz*.1).rgb;
            ref = vec3(1.)-fresnel;
            C = (diff*h)+spec;
        }

        C = mix(FC.rgb,C,exp(-.00015*d*d*d));
        ro = p+n*.01;
        rd = reflect(rd,n);
        C+=clamp(glow,0.,.95)*clr;
    } 
    return vec4(C,alpha);
}

vec4 renderPassA() {
	vec4 O = vec4(0.0);
	vec2 F = _xy;
   
    float T2 = T*2.;
    float rrt = T2*PI/180.;
    wtime=T*.1;
    speed = vec3(T2,0,5.*sin(T*.25));
    
    vec2 uv = (2.*F.xy-R.xy)/max(R.x,R.y);
    vec3 ro = vec3(0,.25,1.5);
    vec3 rd = normalize(vec3(uv,-1));
    
    //camera
    mat2 rx = rot(.75);
    mat2 ry = rot(rrt);
    
    ro.yz *= rx;
    rd.yz *= rx;
    ro.xz *= ry;
    rd.xz *= ry;
    
    // reflection loop (@BigWings)
    vec3 C = vec3(0);
    vec3 ref=vec3(0), fil=vec3(1);
    float d =0.;
    float numBounces = 2.;
    for(float i=0.; i<numBounces; i++) {
        vec4 pass = render(ro, rd, ref, i==numBounces-1., d, uv);
        C += pass.rgb*fil;
        fil*=ref;
        if(i==0.) FC = vec4(FC.rgb,exp(-.000015*d*d*d));
    }
    C = mix(C,FC.rgb,1.-FC.w);
    C = clamp(C,vec3(.01),vec3(1));
    //C=pow(C, vec3(.4545));
    O = vec4(C,1.0);
	return O; 
 } 





/**
    Analogous Oceans | some abstract marching
    @byt3_m3chanic | 07/10/21
*/

float offset[3] = float[]( 0.0, .7215, 1.73 );
float weight[3] = float[]( 0.2, 0.35, 0.0735 );

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
    float dt = distance(uv.xy,vec2(.5))*.35;
    dt = smoothstep(.80,.92,1.-dt);

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