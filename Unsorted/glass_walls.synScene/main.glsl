

/** 
    License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License

    Beveled Glass 5/23/22 @byt3_m3chanic
    
    Just some glass with a refractive surface - animated textured wall behind.
    AA makes it pretty - but can be slow on lower end machines.
    
    big thank you to, and a little bit of help from:
    @blackle, @tdhooper, @iq, @Drakyen
    
*/

#define R 		RENDERSIZE
#define T 		TIME
#define M 		_mouse

#define PI2     PI*2.0
//#define PI      3.14159265359

// AA Setting - comment out to disable
#define AA 2

float time;
mat2 rot(float a) { return mat2(cos(a),sin(a),-sin(a),cos(a)); }
float hash21(vec2 a) { return fract(sin(dot(a,vec2(21.23,41.232)))*43758.5453); }

//@iq https://iquilezles.org/articles/palettes
vec3 hue(float t){ return .45+.4*cos( PI2*t*vec3(.95,.97,.88)*vec3(0.110,0.584,0.949) ); }


float box( vec3 p, vec3 b ){
    vec3 q = abs(p) - b;
    return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float box( vec2 p, vec2 b){
    vec2 d = abs(p)-b;
    return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}

//global
vec3 hit=vec3(0),hitPoint=vec3(0);
float radius = 1.25;
     
vec2 scaler =vec2(6.*Scale);
float pattern(vec3 pt){
    vec3 hp = pt;
    vec2 uv   = hp.xy;

    vec2 grid = fract(uv.xy*scaler)-.5;
    vec2 id   = floor(uv.xy*scaler);
    float hs = hash21(id);

    if(hs>.5) grid.x*=-1.;
    float chk = mod(id.y + id.x,2.) * 2. - 1.;

    vec2 d2 = vec2(length(grid-.5), length(grid+.5));
    vec2 gx = d2.x<d2.y? vec2(grid-.5) : vec2(grid+.5);

    float circle = length(gx)-.5;
    circle=abs(circle)-.075;

    if(hs>.82){
     circle = abs(length(gx.x)-.5)-.075;
     circle = min(circle,abs(length(gx.y)-.5)-.075);
    }
    
    float px = 45.5/R.x;
    circle=abs(abs(circle)-.08)-.01;
    circle=smoothstep(-px,.175-px,circle);

    float h = mix(.0, .1,circle);
    return clamp(h*.1,.0,.1);
}

mat2 spin;

vec2 map(vec3 p){
    vec2 res = vec2(1e5,0.);
    p.xy*=spin;
    
    vec3 q = p-vec3(.2*sin(smoothTime*.05),.2*cos(smoothTime*.0481),1);
    vec3 r = p+vec3(0,0,7);

    float hmp = pattern(q);
    q.z+=.085*sin(smoothTimeC*.25+q.x*3.7);

    float glass = abs(q.z)-hmp;

    if(glass<res.x) {
        res = vec2(glass,2.);
        hit=q;
    }    

    float background = box(r,vec3(12.,8.,.1));
    if(background<res.x) {
        res = vec2(background,3.);
        hit=r;
    }    

    return res;
}

//Tetrahedron technique
//https://iquilezles.org/articles/normalsSDF
vec3 normal(vec3 p, float t, float mindist) {
    float e = mindist*t;
    vec2 h = vec2(1.0,-1.0)*0.5773;
    return normalize( h.xyy*map( p + h.xyy*e ).x + 
					  h.yyx*map( p + h.yyx*e ).x + 
					  h.yxy*map( p + h.yxy*e ).x + 
					  h.xxx*map( p + h.xxx*e ).x );
}

vec3 shade(vec3 p, vec3 rd, float d, float m, inout vec3 n) {
    n = normal(p,d,1.);
    vec3 lpos = vec3(.1,9,7);
    vec3 l = normalize(lpos-p);

    float diff = clamp(dot(n,l),0.,1.);
    vec3 h = vec3(.0005);
    
    if(m==3.) {
        vec2 uv = hitPoint.xy*.075;
        uv.x+=smoothTime*.0125;
        
        float scale = 10.;
        vec2 id = floor(uv*scale);
        uv = fract(uv*scale)-.5;

        float hs = hash21(id);
        vec3 hue = 0.5 + 0.5*cos(hs+vec3(0,1,2)-id.xyx*.2);
         
        if(hs>.55) {
            scale*=.75;
            id   = floor(uv*2.);
            uv = fract(uv*2.)-.5;
            hue *= 0.5 + 0.5*cos(id.x+vec3(0,1,2)+id.y);
        }
    
        h =vec3(0);

        float px = 1.5*scale/R.x;
        float d = box(uv,vec2((fract(hs*3.1)>.9)?.3:.325))-(.015*scale);
        if(fract(hs*3.1)>.9) d=abs(d)-.05;
        d = smoothstep(px,-px,d);

        h = mix (h,hue,d);
    }
    
    
    return h*pow(h, vec3(1.25))*2.5;
}

vec3 renderAll( vec2 uv){    

    vec3 C=vec3(.0);
    vec3 ro = vec3(0,0,1.45),
         rd = normalize(vec3(uv,-1));

    float x = M.xy == vec2(0) ? .0 : -(M.y/R.y * .05 - .025) * PI;
    float y = M.xy == vec2(0) ? .0 :  (M.x/R.x * .05 - .025) * PI;

    mat2 rx = rot(x), ry = rot(y);
    ro.yz *= rx;ro.xz *= ry;
    rd.yz *= rx;rd.xz *= ry;
    
    vec3  p = ro + rd * .1;
    float atten = 1.;
    float k = 1.;
    float iv= 1.;
    float bnc = 5.;
    float alpha = 1.;
    
    // loop inspired/adapted from @blackle's 
    // marcher https://www.shadertoy.com/view/flsGDH
    for(int i=0;i<110;i++)
    {
        vec2 ray = map(p);
        vec3 n=vec3(0);
        float fresnel=0.;
        float d = i<32?ray.x*.85:ray.x;
        float m = ray.y;
        p += rd * d *k;
        
        if (d*d < 1e-7) {
            bnc--;
            hitPoint=hit;

            alpha*=d;
            
            C+=shade(p,rd,d,m,n)*atten;
            if(m==3.)break;
            atten *= .7;
            p += rd*.001;
            k = sign(map(p).x);
            
            vec3 rf=refract(rd,n,iv > 0. ? 1./1.2 : 1.1);
            iv *= -1.;
            if(length(rf) == 0.) rf = reflect(rd,n);
            rd=rf;
            p+=-n*.001;
   
            if(bnc<1.) break;
        }  
        if(distance(p,rd)>100.) { break; }
       
    }

   return clamp(C,vec3(0),vec3(1));
}

float vmul(vec2 v) {return v.x * v.y;}
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
 

    vec3 col = vec3(.00); 
    
    float mTime = TIME;
    time = mTime;    
    spin = rot(.5*sin(T*.2));
    vec2 o = vec2(0);

    // AA and motion blur from iq https://www.shadertoy.com/view/3lsSzf
    // set AA above renderFull
    #ifdef AA
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
    	// pixel coordinates
    	o = vec2(float(m),float(n)) / float(AA) - 0.5;
    	// time coordinate (motion blurred, shutter=0.5)
    	float d = 0.5*vmul(sin(mod(fragCoord.xy * vec2(147,131), vec2(PI * 2.))));
    	time = mTime - 0.1*(1.0/24.0)*(float(m*AA+n)+d)/float(AA*AA-1);
    #endif
        //time = mod(time, 1.);
    	vec2 p = (-RENDERSIZE.xy + 2. * (fragCoord + o)) / RENDERSIZE.x;
    	col += renderAll(p);
        
    #ifdef AA
    }
    col /= float(AA*AA);
    #endif

    col = pow( col, vec3(0.4545) );
    fragColor = vec4(col, 0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}