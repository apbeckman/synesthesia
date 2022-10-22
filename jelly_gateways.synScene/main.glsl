

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

 


#define MDIST 100.0
#define STEPS 128.0
#define pi 3.1415926535
#define rot(a) mat2(cos(a),sin(a),-sin(a),cos(a))
#define pmod(p,x) (mod(p,x)-0.5*(x))

vec3 rdg = vec3(0);
vec3 hsv(vec3 c){
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}
//iq box sdf
float ebox(vec3 p, vec3 b){
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}
float ebox(vec2 p, vec2 b){
  vec2 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,q.y),0.0);
}

float lim(float p, float s, float lima, float limb){
    return p-s*clamp(round(p/s),lima,limb);
}
float idlim(float p, float s, float lima, float limb){
    return clamp(round(p/s),lima,limb);
}

float dibox(vec3 p,vec3 b,vec3 rd){
    vec3 dir = sign(rd)*b;   
    vec3 rc = (dir-p)/rd;
    float dc = min(rc.y,rc.z)+0.01;
    return dc;
}
float easeOutBounce(float x) {
    float n1 = 7.5625;
    float d1 = 2.75;
    if (x < 1. / d1) {
        return n1 * x * x;
    } 
    else if (x < 2. / d1) {
        return n1 * (x -= 1.5 / d1) * x + 0.75;
    } 
    else if (x < 2.5 / d1) {
        return n1 * (x -= 2.25 / d1) * x + 0.9375;
    } 
    else {
        return n1 * (x -= 2.625 / d1) * x + 0.984375;
    }
}
vec3 map(vec3 p){
    float t = -smoothTime*0.2;
    vec3 rd2 = rdg;
    vec2 a = vec2(1);
    vec2 b = vec2(2);
    p.xz*=rot(t*0.3*pi/3.);
    rd2.xz*=rot(t*0.3*pi/3.);
    //p.xz*=rot(pi/4.);
    //rd2.xz*=rot(pi/4.); 
    vec3 po = p;
    float dsz = 0.45;
    float m = 2.42-dsz;
    float bs = 1.-dsz*0.5;
    
    //VERTIAL TRANSLATION
    p.y+=t*m;
    
    //VERTIAL REP
    float id1 = floor(p.y/m);
    p.y = pmod(p.y,m);
    
    //ROTATE EACH LAYER
    p.xz*=rot(id1*pi/2.);
    rd2.xz*=rot(id1*pi/2.);

    vec3 p2 = p; //dibox p1
    
    //Auxillary boxes positions
    vec3 p3 = p;
    vec3 rd3 = rd2;
     
    p3.xz*=rot(pi/2.);
    rd3.xz*=rot(pi/2.);
    vec3 p4 = p3; 
    
    
    //HORIZONTAL REP
    p2.z = pmod(p2.z-m*0.5,m);
    p4.z = pmod(p4.z-m*0.5,m);
    
    float cnt = 100.;
    float id2 = idlim(p.z,m,-cnt,cnt);
    float id3 = idlim(p3.z,m,-cnt,cnt);
    p.z = lim(p.z,m,-cnt,cnt);
    p3.z = lim(p3.z,m,-cnt,cnt);
    
    
    //CLOSING ANIMATION 
    float close = max((id1-t)*1.,-2.);
    float close2 = clamp(max((id1-t-0.3)*1.,-2.)*1.4,0.,1.);
    close+=id2*0.025;
    close = clamp(close*1.4,0.,1.);
    close = easeOutBounce(close);
    //close = 1.0-easeOutBounce(1.-close);

    
    
    //CLOSING OFFSET
    p.x = abs(p.x)-34.5*0.5-0.25*7.;
    p.x-=close*34.5*0.52-0.055;
    
    p3.x = abs(p3.x)-36.5;

    p.x-=((id1-t)*0.55)*close*2.4;
    p3.x-=((id1-t)*0.55)*close2*2.4;
    //WAVEY
    p.x+=(sin(id1+id2-t*6.0)*0.18+4.)*close*2.4;
    p3.x+=(sin(id1+id3-t*6.0)*0.18+4.)*smoothstep(0.,1.,close2)*2.4;
    
    
    //BOX SDF
    a = vec2(ebox(p,vec3(7.5*2.5,bs,bs))-0.2,id2);
    
    //AUXILLARY BOX
    b = vec2(ebox(p3,vec3(7.5*2.5,bs,bs))-0.2,id3);
    
    a=(a.x<b.x)?a:b;
    //ARTIFACT REMOVAL
    float c = dibox(p2,vec3(1,1,1)*m*0.5,rd2)+.1;
    //ARTIFACT REMOVAL 2
    c = min(c,dibox(p4,vec3(1,1,1)*m*0.5,rd3)+.1);
    

    float nsdf = a.x;
    
    a.x = min(a.x,c); //Combine artifact removal
    a.y = id1;
    return vec3(a,nsdf);
}
vec3 norm(vec3 p){
    vec2 e = vec2(0.005,0);
    return normalize(map(p).x-vec3(
    map(p-e.xyy).x,
    map(p-e.yxy).x,
    map(p-e.yyx).x));
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    vec3 col = vec3(0);
    vec3 ro = vec3(0,13,-5)*1.5;
    if(_mouse.z>0.){
    ro.yz*=rot(1.0*(_mouse.y/RENDERSIZE.y-0.2));
    ro.zx*=rot(-7.0*(_mouse.x/RENDERSIZE.x-0.5));
    }
    vec3 lk = vec3(0,0,0);
    vec3 f = normalize(lk-ro);
    vec3 r = normalize(cross(vec3(0,1,0),f));
    vec3 rd = normalize(f*(0.5)+uv.x*r+uv.y*cross(f,r));  
    rdg = rd;
    vec3 p = ro;
    float dO = 0.;

    vec3 d= vec3(0);
    for(float i = 0.; i<STEPS; i++){
        p = ro+rd*dO;
        d = map(p);
        dO+=d.x;
        if(abs(d.x)<0.005){
            break;
        }
        if(dO>MDIST){
            dO = MDIST;
            break;
        }
    }

    {
        vec3 ld = normalize(vec3(0,45,0)-p);
      
        //sss from nusan
        float sss=0.01;
        for(float i=1.; i<20.; ++i){
            float dist = i*0.09;
            sss += smoothstep(0.,1.,map(p+ld*dist).z/dist)*0.023;
        }
        vec3 al = vec3(0.204,0.267,0.373);
        vec3 n = norm(p);
        vec3 r = reflect(rd,n);
        float diff = max(0.,dot(n,ld));
        float amb = dot(n,ld)*0.45+0.55;
        float spec = pow(max(0.,dot(r,ld)),40.0);
        float fres = pow(abs(.7+dot(rd,n)),3.0);     
        //ao from blackle 
        #define AO(a,n,p) smoothstep(-a,a,map(p+n*a).z)
        float ao = AO(.3,n,p)*AO(.5,n,p)*AO(.9,n,p);

        col = al*
        mix(vec3(0.169,0.000,0.169),vec3(0.984,0.996,0.804),mix(amb,diff,0.75))
        +spec*0.3+fres*mix(al,vec3(1),0.7)*0.4;
        col+=sss*hsv(vec3(fract(d.y*0.5+d.y*0.1+0.001)*0.45+0.5,0.9,1.35));
        col*=mix(ao,1.,0.85);
        col = pow(col,vec3(0.75));
    }
    col = clamp(col,0.,1.);
    //col = smoothstep(0.,1.,col);
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}