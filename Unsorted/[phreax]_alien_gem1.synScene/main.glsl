

/* Creative Commons Licence Attribution-NonCommercial-ShareAlike 
   phreax 2022
*/

// /#define PI 3.141592
#define TAU (2.*PI)
#define SIN(x) (sin(x)*.5+.5)
#define PHI 1.618033988749895
#define BUMP_EPS 0.004
#define FAR 100.


float tt, g_mat;
float closest = FAR;


mat2 rot(float a) { return mat2(cos(a), -sin(a), sin(a), cos(a)); }

float fadeInOut(float t1, float t2, float fadeTime) {
    return smoothstep(t1, t1+fadeTime, smoothTimeB)-smoothstep(t2-fadeTime, t2, smoothTimeB);
}

// zucconis spectral palette https://www.alanzucconi.com/2017/07/15/improving-the-rainbow-2/
vec3 bump3y (vec3 x, vec3 yoffset)
{
    vec3 y = 1. - x * x;
    y = clamp((y-yoffset), vec3(0), vec3(1));
    return y;
}

const highp float NOISE_GRANULARITY = 0.5/255.0;

highp float random(highp vec2 coords) {
   return fract(sin(dot(coords.xy, vec2(12.9898,78.233))) * 43758.5453);
}

vec3 spectral_zucconi6(float x) {
    x = fract(x);
    const vec3 c1 = vec3(3.54585104, 2.93225262, 2.41593945);
    const vec3 x1 = vec3(0.69549072, 0.49228336, 0.27699880);
    const vec3 y1 = vec3(0.02312639, 0.15225084, 0.52607955);
    const vec3 c2 = vec3(3.90307140, 3.21182957, 3.96587128);
    const vec3 x2 = vec3(0.11748627, 0.86755042, 0.66077860);
    const vec3 y2 = vec3(0.84897130, 0.88445281, 0.73949448);
    return bump3y(c1 * (x - x1), y1) + bump3y(c2 * (x - x2), y2) ;
}

// Amazing fractal texture from jarble https://www.shadertoy.com/view/csl3zl
vec2 triangle_wave(vec2 a){
    
    vec2 a2 = vec2(1.,0.25),
    
    a1 = a-a2;
    
    return abs(fract((a1)*(a2.x+a2.y))-.5);
}

float getDepth(sampler2D sampler, vec2 uv) {
    uv = clamp(uv, vec2(.0), vec2(1.));
    vec2 eps = vec2(0.004, 0); 
    float vc = texture(sampler, uv).r;
    float vu = texture(sampler, uv-eps.yx).r;
    float vb = texture(sampler, uv+eps.yx).r;
    float vr = texture(sampler, uv+eps.xy).r;
    float vl = texture(sampler, uv-eps.xy).r;
    
    return clamp((vc + vu + vb + vr + vl)/5., 0., .4);
    
}

float rect( vec2 p, vec2 b, float r ) {
    vec2 d = abs(p) - (b - r);
    return length(max(d, 0.)) + min(max(d.x, d.y), 0.) - r;
}


vec3 kalei(vec3 p) {
    float iter = 3.;


    float s= 1.;
    for(float i=0.; i< iter; i++) {
        p = abs(p);  
        p.xz *= rot(i/iter*PI+.2*tt);
        p -= .1;
         
    }

    return clamp(p, -1e5, 1e5);

}


// n-fold symmetry by mla
vec2 foldSym(vec2 p, float N) {

    float t = atan(p.x,-p.y);
    t = mod(t+PI/N,2.0*PI/N)-PI/N;
    return length(p)*vec2(cos(t),sin(t));
}

vec3 transform(vec3 p) {

    p *= 1.5;
    
    p.xy = foldSym(p.xy, 6.);
    p.xy = abs(p.xy)-0.25*SIN(.3*tt);

    p = kalei(p);
    
    
    p.yz *= rot(PI*.5);
    p.xz *= rot(.2*tt);

    float r = 1.;
    vec2 cp = vec2(length(p.xz)-r, p.y);
    

    float rev = 2.;
    float a = atan(p.z, p.x);
    
 
    cp *= rot(rev*a+.3*tt);
    cp= abs(cp) - mix(.4, 1., SIN(tt));
    cp= mix(cp, abs(cp) - .2, smoothstep(10., 11., bass_time));
    cp *= rot(.5*tt);


    return vec3(cp, p.z);
}

vec3 transform2(vec3 p) {

    p *= .8;

    p.xy = foldSym(p.xy, 3.); 
    p.xy = abs(p.xy)-0.3*SIN(-.3*tt);  
    p.yz *= rot(PI*.25);
    

    float r = 1.;
    vec2 cp = vec2(length(p.xz)-r, p.y);
 

    float rev = 1.;
    float a = atan(p.z, p.x);
    
 
    cp *= rot(rev*a+.3*tt);
    cp= abs(cp) - mix(.4, 1., SIN(tt));
    // cp=abs(cp) - .5;
    cp *= rot(-.5*tt);

    return vec3(cp, p.z);
}

float smin(float a, float b, float k) {
  float h = clamp((a-b)/k * .5 + .5, 0.0, 1.0);
  return mix(a, b, h) - h*(1.-h)*k;
}

float map(vec3 p) {   
  
    vec3 bp = p;
    float edge = 0.05;
   

    vec2 cp = transform(p).xy;
    vec2 cp2 = transform2(p).xy;
   
    float dr = rect(cp.xy, vec2(.3, .3), edge);
    float dr2 = rect(cp2.xy, vec2(.08), 0.02);
   
    float d = smin(dr, dr2,.8);
    
    g_mat = dr < dr2 ? 0. : 1.;
    
    return .5*d;
}


vec3 getNormal(vec3 p) {

    vec2 eps = vec2(0.001, 0.0);
    return normalize(vec3(map(p + eps.xyy) - map(p - eps.xyy),
                          map(p + eps.yxy) - map(p - eps.yxy),
                          map(p + eps.yyx) - map(p - eps.yyx)
                         )
                     );
}

float bumpSurf3D( in vec3 p){

    p.z += .3*tt;
    p = abs(mod(p*4., 2.*0.125)-0.0125);
    
    float x = min(p.x,min(p.z, p.y))/0.03125;

    return clamp(x, 0., 1.);


}

// Standard function-based bump mapping function (from Shane)
vec3 doBumpMap(in vec3 p, in vec3 nor, float bumpfactor){
    
    const float eps = BUMP_EPS;
    float ref = bumpSurf3D(p);                 
    vec3 grad = vec3( bumpSurf3D(vec3(p.x-eps, p.y, p.z))-ref,
                      bumpSurf3D(vec3(p.x, p.y-eps, p.z))-ref,
                      bumpSurf3D(vec3(p.x, p.y, p.z-eps))-ref )/eps;                     
          
    grad -= nor*dot(nor, grad);          
                      
    return normalize( nor + bumpfactor*grad );
	
}
// from iq
float softshadow( in vec3 ro, in vec3 rd, float mint, float maxt, float k )
{
    float res = 1.0;
    float ph = 1e20;
    for( float t=mint; t<maxt; )
    {
        float h = map(ro + rd*t);
        if( h<0.001 )
            return 0.0;
        float y = h*h/(2.0*ph);
        float d = sqrt(h*h-y*y);
        res = min( res, k*d/max(0.0,t-y) );
        ph = h;
        t += h;
    }
    return res;
}


vec2 raymarch(vec3 ro, vec3 rd, float steps) {

    float mat = 0.,
          t   = 0.,
          d   = 0.;
    vec3 p = ro;
    
    
    for(float i=.0; i<steps; i++) {
    
        d = map(p);
        mat = g_mat;  // save global material
        
        closest = min(closest, d/t);
        if(abs(d) < 0.0001 || t > FAR) break;
 
        t += d;
        p += rd*d;
       
    }
    
    return vec2(t, mat);
}

float n21(vec2 p) {
      return fract(sin(dot(p, vec2(524.423,123.34)))*3228324.345);
}

float noise(vec2 n) {
    const vec2 d = vec2(0., 1.0);
    vec2 b = floor(n);
    vec2 f = mix(vec2(0.0), vec2(1.0), fract(n));
    return mix(mix(n21(b), n21(b + d.yx), f.x), mix(n21(b + d.xy), n21(b + d.yy), f.x), f.y);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


   	vec2 uv = (fragCoord - .5*RENDERSIZE.xy)/RENDERSIZE.y;

    tt = .5*bass_time;
    vec3 ro = vec3(uv*6.,-4.),
          rd = vec3(0,0,1.),
          lp = vec3(3., 0., -2),
          lp2 = vec3(-3., 0., -2);

    vec3 col;
    

       
    float mat = 0.,
          t   = 0.,
          d   = 0.;
 

    vec2 e = vec2(0.0035, -0.0035);
     
    // background color
    vec3 c1 = vec3(0.165,0.208,0.298);
    vec3 c2 = vec3(0.180,0.337,0.337);
    
    // light color
    vec3 lc1 = vec3(0.573,0.424,0.976);
    vec3 lc2 = vec3(0.573,0.922,0.969);
    
    
    // currently only one pass
    for(float i = 0.; i < 1.; i++) {
        float steps = i > 0. ? 50. : 250.;
        vec2 rm = raymarch(ro, rd, steps);
        mat = rm.y;
        
        vec3 p = ro + rm.x*rd;
        
        vec3 n = normalize( e.xyy*map(p+e.xyy) + e.yyx*map(p+e.yyx) +
                                e.yxy*map(p+e.yxy) + e.xxx*map(p+e.xxx));


        
        if(mat == 1.) {
        
            vec3 pt = transform2(p);
            vec3 np = transform2(n);
            n = doBumpMap(pt, n, mix(0., .001, length(p.xy)-.5));
        }
    
     
        if(rm.x < 50.) {
        
            vec3 l = normalize(lp-p);
            vec3 l2 = normalize(lp2-p);
            float dif = max(dot(n, l), .0);
            float dif2 = max(dot(n, l2), .0);
            float spe = pow(max(dot(reflect(-rd, n), -l), .0),40.);

            float sss = smoothstep(0., 1., map(p + l * .4)) / .4;
            vec3 n2 = n;
            n2.xy += noise(p.xy) * .5 - .025;
            n2 = normalize(n2);
            float height = atan(n2.y, n2.x);

            vec3 iri = spectral_zucconi6(height*1.11)*smoothstep(.8, .2, abs(n2.z))-.02;
       
          
            col += dif*lc1*.5 + .5*dif2*lc2 + .5*iri ;
           
          
            if(mat > 0.) {
                n += 111.4*texture(image3, n.xy*200.).rgb;
                rd = reflect(rd, n);
                
               // rd.yz *= rot(PI*.8);
         
                //vec3 refl = texture(image3, rd).rgb;
                
               // refl *= mix(vec3(1), spectral_zucconi6(n.x*n.y*3.), .3); // reflect rainbows too
               // col = mix(col, refl, .6);
    
              
            } 
            if(mat == 1.)
              col = mix(col, col*spectral_zucconi6(p.z*.2+length(p.xy*.2)*.2+.05*tt), .95);
            if(mat == 0.) 
              col = mix(col, col*(mix(.6, .0, length(p*.5)-.3+.3*sin(tt))+vec3(1.000,0.506,0.239)*spectral_zucconi6(length(p.xy*.4)+.6+.1*tt)), 1.);

        } else {
            col =  mix(c1-.2, c2, (.3-pow(dot(uv, uv), .8)))*.0; // background
            // outer glow (inspired from https://www.shadertoy.com/view/ldB3Rz)
            float f = 1.0 - clamp(closest * 0.5, 0.0, 1.0);		

            float glowAmount = 0.0;
       
            glowAmount +=  pow(f, 400.0) * (0.08+.1*SIN(tt));
            vec3 glowColor = spectral_zucconi6(length(p.xy*.2)+.6);
            col += glowColor * glowAmount;
        } 
    
    }

    
    col += mix(-NOISE_GRANULARITY, NOISE_GRANULARITY, random(uv)); // dithering
    
    col *= 1.8;
    col *= mix(.2, 1., (1.3-pow(dot(uv, uv), .5))); // vignette
    col = pow(col, vec3(.7)); // gamma
    

    
    fragColor = vec4(col, 1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}