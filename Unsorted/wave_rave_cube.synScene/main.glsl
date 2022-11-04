

			//******** Common Code Begins ********


// Display the individual hexagon grid cells.
//#define SHOW_HEX_GRID
    
// Flat shading. No lighting.
//#define FLAT_SHADING
 
// Apply hatching
#define HATCHING

// Tile scale.
#define GSCALE vec2(1./5.)*vec2(.5, .8660254)


//  Vertices and edge midpoints: Clockwise from the bottom left. -- Basically, the ones 
// above rotated anticlockwise. :)
vec2[6] vID = vec2[6](vec2(-.5, -1./3.)/vec2(.5, 1), vec2(-.5, 1./3.)/vec2(.5, 1), vec2(0, 2./3.)/vec2(.5, 1), 
                      vec2(.5, 1./3.)/vec2(.5, 1), vec2(.5, -1./3.)/vec2(.5, 1), vec2(0, -2./3.)/vec2(.5, 1));
//vec2[6] eID = vec2[6](vec2(-.5, 0)/vec2(.5, 1), vec2(-.25, .5)/vec2(.5, 1), vec2(.25, .5)/vec2(.5, 1), vec2(.5, 0)/vec2(.5, 1), 
                      //vec2(.25, -.5)/vec2(.5, 1), vec2(-.25, -.5)/vec2(.5, 1));

 
// Standard 2D rotation formula.
mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }


// IQ's vec2 to float hash.
float hash21(vec2 p){ 
    
    p = mod(p, 256.);
    // An annoying, but necessary, hack for systems with less sin
    // function accuracy. If anyone knows a way around it, feel 
    // free to let me know.
    //p = floor(p*1048576.)/1048576.;
    return fract(sin(dot(p, vec2(27.649, 57.583)))*43758.5453); 
}

/*
// UE4 random function: I like this because it incorporates a modulo
// 128 wrap, so in theory, things shouldn't blow up with increasing input.
// Also, in theory, you could tweak the figures by hand to get a really
// scrambled output... When I'm feeling less lazy, I might do that.
float hash21(vec2 p) {

    p = fract(p/128.)*128. - vec2(59.340627, 73.465623);
    return fract(dot(p.xyx*p.xyy, vec3(20.390625, 80.703127, 12.4281203)));
    
}
*/


/*
// This is IQ's WebGL 2 hash formula: I've modified it slightly to 
// take in the normal decimal floats that we're used to passing. It 
// works here, I think, but I'd consult the experts before using it.
//
// I remember reading through a detailed explanation of the C++ hash 
// we all used to use many years ago (which the following would be
// similar to), but have long since forgotten why it works. By the 
// way Nimitz and Dave Hoskins have formulae on Shadertoy that's worth
// looking up.
//
// Integer Hash - III - IQ
// https://www.shadertoy.com/view/4tXyWN
float hash21(vec2 p){
    
    uvec2 q = uvec2(ivec2(p*1023.));
	q = 1103515245U*((q>>1U)^q.yx);
    uint n = 1103515245U*(q.x^(q.y>>3U));
    return float(n)/float(0xffffffffU);
}
*/


// IQ's distance to a regular pentagon, without trigonometric functions. 
// Other distances here:
// https://iquilezles.org/articles/distfunctions2d
//
#define NV2 4
//
float sdPoly4(in vec2 p, in vec2[NV2] v){

    int num = v.length();
    float d = dot(p - v[0],p - v[0]);
    float s = 1.0;
    for( int i = 0, j = num - 1; i < num; j = i, i++){
    
        // distance
        vec2 e = v[j] - v[i];
        vec2 w =    p - v[i];
        vec2 b = w - e*clamp(dot(w, e)/dot(e, e), 0., 1. );
        d = min( d, dot(b,b) );

        // winding number from http://geomalgorithms.com/a03-_inclusion.html
        bvec3 cond = bvec3( p.y>=v[i].y, p.y<v[j].y, e.x*w.y>e.y*w.x );
        if( all(cond) || all(not(cond)) ) s*=-1.0;  
    }
    
    return s*sqrt(d);
}


// Determines which side of a line a pixel is on. Zero is the threshold.
float line(vec2 p, vec2 a, vec2 b){
     return ((b.x - a.x)*(p.y - a.y) - (b.y - a.y)*(p.x - a.x));
}

/*
// IQ's unsigned box formula.
float sBoxS(in vec2 p, in vec2 b, in float sf){

  vec2 d = abs(p) - b + sf;
  return min(max(d.x, d.y), 0.) + length(max(d, 0.)) - sf;
}
*/
// IQ's standard box function.
float sBox(in vec2 p, in vec2 b){
   
    vec2 d = abs(p) - b;
    return min(max(d.x, d.y), 0.) + length(max(d, 0.));
}

// This will draw a box (no caps) of width "ew" from point "a "to "b". I hacked
// it together pretty quickly. It seems to work, but I'm pretty sure it could be
// improved on. In fact, if anyone would like to do that, I'd be grateful. :)
float lBox(vec2 p, vec2 a, vec2 b, float ew){
    
    float ang = atan(b.y - a.y, b.x - a.x);
    p = rot2(ang)*(p - mix(a, b, .5));

    vec2 l = vec2(length(b - a), ew);
    return sBox(p, (l + ew)/2.) ;
}

 


// Fork of "Box Singularity" by Tater. https://shadertoy.com/view/7dVGDd
// 2021-10-20 06:09:34

#define MDIST 350.0
#define STEPS 200.0
#define pi 3.1415926535
#define rot(a) mat2(cos(a),sin(a),-sin(a),cos(a))
#define pmod(p,x) (mod(p,x)-0.5*(x))


//iq palette
vec3 pal( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d ){
    return a + b*cos(2.*pi*(c*t+d));
}
float h21 (vec2 a) {
    return fract(sin(dot(a.xy,vec2(12.9898,78.233)))*43758.5453123);
}
float h11 (float a) {
    return fract(sin((a)*12.9898)*43758.5453123);
}
float box(vec3 p, vec3 b){
    vec3 d = abs(p)-b;
    return max(d.x,max(d.y,d.z));
}
//iq box sdf
float ebox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float swave(float x, float a){
    return (sin(x*pi/3.-pi/2.)/sqrt(a*a+sin(x*pi/3.-pi/2.)*sin(x*pi/3.-pi/2.))+1./sqrt(a*a+1.))*0.5;
}
vec3 rdg = vec3(0);
//No Cell bounds SDF of the blocks
float nsdf = 0.;
bool rnsdf = false;
vec2 blocks(vec3 p, vec3 scl, vec3 rd){
    float t = smoothTimeC*0.25;
    
    bvec3 isEdge = bvec3(true);
    // the domain of the fractal being generated
    // will be modified in the iteration part
    vec3 dMin = vec3(-0.5) * scl;
    vec3 dMax = vec3(0.5) * scl;
    vec3 dMini = dMin;
    vec3 dMaxi = dMax;
    
    float id = 0.;
    float seed = floor(t/8.)+0.2;
    
    float ITERS = 5.;
    
    vec3 dim = dMax - dMin;
    //Big thanks for @0b5vr for cleaner version of subdiv algo
    for (float i = 0.; i < ITERS; i++) {

        vec3 divHash = vec3(
            0.49,
            0.501,
            0.511
        );
        vec3 divide = divHash * dim + dMin;
        // update the box domain
        dMax = mix( dMax, divide, step( p, divide ));
        dMin = mix( divide, dMin, step( p, divide ));

        float pad = 0.01;
        if(dMaxi.x>dMax.x+pad&&dMini.x<dMin.x-pad)isEdge.x=false;
        if(dMaxi.y>dMax.y+pad&&dMini.y<dMin.y-pad)isEdge.y=false;
        if(dMaxi.z>dMax.z+pad&&dMini.z<dMin.z-pad)isEdge.z=false;
        
        // id will be used for coloring and hash seeding
        
        vec3 diff = mix( -divide, divide, step( p, divide));
        id = length(diff + 10.0);
    
        // recalculate the dimension
        dim = dMax - dMin;
    }
    float volume = dim.x*dim.y*dim.z;
    vec3 center = (dMin + dMax)/2.0;
    float b = 0.;
    if(any(isEdge)) {
    
        float expand = 1.05+0.45*(sin(length(center+100.0)*0.35+t*1.0*pi/3.)*0.5+0.5);
        //expand = 1.5; 
        if(isEdge.x){
        dim.x+=abs(center.x*expand-center.x)*2.0;
        center.x*=expand;
        }
        if(isEdge.y){
        dim.y+=abs(center.y*expand-center.y)*2.0;
        center.y*=expand;
        }
        if(isEdge.z){
        dim.z+=abs(center.z*expand-center.z)*2.0;
        center.z*=expand;
        }
        //id = 1.;
    }
    
    vec3 edgeAxis = mix(dMin, dMax, step(0.0, rd));
    vec3 dAxis = abs(p - edgeAxis) / (abs(rd) + 1E-4);
    float dEdge = min(dAxis.x,min(dAxis.y,dAxis.z));
    b= dEdge;
    vec3 d = abs(center);
    dim-=0.4;
    float a = ebox(p-center,dim*0.5)-0.2;


    nsdf = a;
    a = min(a, b);

    
    id = h11(id)*1000.0;

    return vec2(a,id);
}

vec3 map(vec3 p){
    float t = smoothTimeB*0.25;

    vec3 po = p;
    vec2 a = vec2(1);

    vec3 scl = vec3(15,15,15);
    vec3 rd2 = rdg;
    p.yz*=rot(t*0.5*pi/3.);
    rd2.yz*=rot(t*0.5*pi/3.);
    p.xy*=rot(t*0.5*pi/3.);
    rd2.xy*=rot(t*0.5*pi/3.);
    a = blocks(p,scl,rd2)+0.01;
    
   
    a.x = max(box(p,vec3(scl*1.3)),a.x);
    

    return vec3(a,nsdf);
}
vec3 norm(vec3 p){
    vec2 e = vec2(0.01,0.);
    return normalize(map(p).x-vec3(
    map(p-e.xyy).x,
    map(p-e.yxy).x,
    map(p-e.yyx).x));
}
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    float t = smoothTime*0.125;
    vec3 col = vec3(0);
    vec3 ro = vec3(0,11.5,-20)*1.8; //+- for closeness
    if(_mouse.z>0.){
    ro.yz*=rot(2.0*(_mouse.y/RENDERSIZE.y-0.5));
    ro.zx*=rot(-7.0*(_mouse.x/RENDERSIZE.x-0.5));
    }
   // else ro.xz*=rot(t*0.3);
    ro.xz*=rot(-pi/4.);
    vec3 lk = vec3(0,0.,0);
    vec3 f = normalize(lk-ro);
    vec3 r = normalize(cross(vec3(0,1,0),f));
    vec3 rd = normalize(f*(0.99)+uv.x*r+uv.y*cross(f,r));    
    rdg = rd;
    vec3 p = ro;
    float dO = 0.;
    vec3 d = vec3(0);
    bool hit = false;
    for(float i = 0.; i<STEPS; i++){
        p = ro+rd*dO;
        d = map(p);
        dO+=d.x*0.99;
        if(abs(d.x)<0.0001||i==STEPS-1.){
            hit = true;
            break;
        }
        if(d.x>MDIST){
            dO=MDIST;
            break;
        }
    }
    if(hit){
        vec3 ld = normalize(vec3(0.5,0.9,-0.9));
        vec3 n = norm(p);
        vec3 r = reflect(rd,n);
        vec3 e = vec3(0.5);
        vec3 al = pal(fract(d.y+t*1.0/12.0)*0.6-0.2,e*1.5,e,e*2.0,vec3(0,0.33,0.66));
        col = al;
        
        //lighting EQs from @blackle
        float spec = length(sin(r*5.)*.5+.5)/sqrt(3.);
        float fres = 1.-abs(dot(rd,n))*0.9;
        
        float diff = length(sin(n*2.)*.5+.7)/sqrt(3.);
        
        #define AO(a,n,p) smoothstep(-a,a,map(p+n*a).z)
        float ao = AO(0.3,n,p)*AO(.5,n,p)*AO(.9,n,p);
        col = al*diff+pow(spec,5.0)*fres*diff*vec3(1.000,0.878,0.792);
        col*=pow(ao,0.4);
        
    }
    //col = pow(col,vec3(0.99));
    vec3 bg = mix(vec3(0.345,0.780,0.988),vec3(0.424,0.059,0.925),length(uv));
    col = mix(col,bg,dO/MDIST);
    fragColor = vec4(col,1.0);
    return fragColor;
}
/*
#define AA 2.0
#define ZERO min(0.0,TIME)
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    float px = 1.0/AA;
    vec4 col = vec4(0);
    
    if(AA==1.0) {render(col,fragCoord); fragColor = col; return;}
    
    for(float i = ZERO; i <AA; i++){
        for(float j = ZERO; j <AA; j++){
            vec4 col2;
            vec2 coord = vec2(fragCoord.x+px*i,fragCoord.y+px*j);
            render(col2,coord);
            col.rgb+=col2.rgb;
            rdg = vec3(0);
        }
    }
    col/=AA*AA;
    fragColor = vec4(col);
	return fragColor; 
 } 

*/


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}