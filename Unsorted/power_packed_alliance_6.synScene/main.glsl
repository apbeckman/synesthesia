

			//******** Common Code Begins ********


// math
//#define PI 3.141592653589
#define saturate(x) clamp(x,0.,1.)
float hash( vec3 x );
float hash( vec2 p );
float hash( float p );
float hash2Interleaved( vec2 x );
//float noise( vec3 x );
vec4 noised( vec3 x );
mat2 rotate( float t );

// mapping
vec3 randomSphereDir( vec2 rnd );
vec3 randomHemisphereDir( vec3 dir, float i );
vec4 tex3D( sampler2D tex, vec3 p, vec3 n );
vec3 bumpMapping( sampler2D tex, vec3 p, vec3 n, float bf );

// tone mapping
vec3 acesToneMapping( vec3 col );
vec3 filmicToneMapping( vec3 col );





// ---------------------------------------------
// Math
// ---------------------------------------------
float hash( vec3 p )
{
    return fract(sin(dot(p,vec3(127.1,311.7, 74.7)))*43758.5453123);
}

float hash( vec2 p )
{
    return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453123);
}

float hash( float p ) 
{
    return fract(sin(p)*43758.5453123);
}

float hash2Interleaved( vec2 x )
{
    // between random & dithered pattern
    // good for jittering and blur as well as blue noise :)
    // http://www.iryoku.com/next-generation-post-processing-in-call-of-duty-advanced-warfare
    vec3 magic = vec3( 0.06711056, 0.00583715, 52.9829189 );
    return fract( magic.z * fract( dot( x, magic.xy ) ) );
}

vec4 noised( vec3 x )
{
	// https://iquilezles.org/articles/gradientnoise
    vec3 p = floor(x);
    vec3 w = fract(x);
    
    vec3 u = w*w*w*(w*(w*6.0-15.0)+10.0);
    vec3 du = 30.0*w*w*(w*(w-2.0)+1.0);

    float a = hash( p+vec3(0,0,0) );
    float b = hash( p+vec3(1,0,0) );
    float c = hash( p+vec3(0,1,0) );
    float d = hash( p+vec3(1,1,0) );
    float e = hash( p+vec3(0,0,1) );
    float f = hash( p+vec3(1,0,1) );
    float g = hash( p+vec3(0,1,1) );
    float h = hash( p+vec3(1,1,1) );

    float k0 =   a;
    float k1 =   b - a;
    float k2 =   c - a;
    float k3 =   e - a;
    float k4 =   a - b - c + d;
    float k5 =   a - c - e + g;
    float k6 =   a - b - e + f;
    float k7 = - a + b + c - d + e - f - g + h;

    return vec4( -1.0+2.0*(k0 + k1*u.x + k2*u.y + k3*u.z + k4*u.x*u.y + k5*u.y*u.z + k6*u.z*u.x + k7*u.x*u.y*u.z), 
                      2.0* du * vec3( k1 + k4*u.y + k6*u.z + k7*u.y*u.z,
                                      k2 + k5*u.z + k4*u.x + k7*u.z*u.x,
                                      k3 + k6*u.x + k5*u.y + k7*u.x*u.y ) ).yzwx;
}
float noise( vec3 x )
{
	// https://iquilezles.org/articles/gradientnoise
    vec3 p = floor(x);
    vec3 w = fract(x);
    
    vec3 u = w*w*w*(w*(w*6.0-15.0)+10.0);
    vec3 du = 30.0*w*w*(w*(w-2.0)+1.0);

    float a = hash( p+vec3(0,0,0) );
    float b = hash( p+vec3(1,0,0) );
    float c = hash( p+vec3(0,1,0) );
    float d = hash( p+vec3(1,1,0) );
    float e = hash( p+vec3(0,0,1) );
    float f = hash( p+vec3(1,0,1) );
    float g = hash( p+vec3(0,1,1) );
    float h = hash( p+vec3(1,1,1) );

    float k0 =   a;
    float k1 =   b - a;
    float k2 =   c - a;
    float k3 =   e - a;
    float k4 =   a - b - c + d;
    float k5 =   a - c - e + g;
    float k6 =   a - b - e + f;
    float k7 = - a + b + c - d + e - f - g + h;
    return -1.0+2.0*(k0 + k1*u.x + k2*u.y + k3*u.z + k4*u.x*u.y + k5*u.y*u.z + k6*u.z*u.x + k7*u.x*u.y*u.z);
}

mat2 rotate( float t ) {
    float a = cos(t);
    float b = sin(t);
    
    return mat2( a, b, -b, a );
}


// ---------------------------------------------
// Mapping
// ---------------------------------------------
vec3 randomSphereDir( vec2 rnd )
{
    float s = rnd.x*PI*2.;
    float t = rnd.y*2.-1.;
    return vec3(sin(s), cos(s), t) / sqrt(1.0 + t * t);
}

vec3 randomHemisphereDir( vec3 dir, float i )
{
    vec3 v = randomSphereDir( vec2(hash(i+1.), hash(i+2.)) );
    return v * sign(dot(v, dir));
}

vec4 tex3D( sampler2D tex, vec3 p, vec3 n )
{
    n = abs(n);
    
    vec4 c = texture(tex, p.yz) * n.x;
    c += texture(tex, p.xz) * n.y;
    c += texture(tex, p.xy) * n.z;
    
    return c / 3.;
}

vec3 bumpMapping( sampler2D tex, vec3 p, vec3 n, float bf )
{
    // clever code taken from Shane
    // https://www.shadertoy.com/view/MscSDB
    const vec2 e = vec2(0.001, 0);
    
    mat3 m = mat3( tex3D(tex, p - e.xyy, n).rgb,
                   tex3D(tex, p - e.yxy, n).rgb,
                   tex3D(tex, p - e.yyx, n).rgb);
    
    vec3 g = vec3(0.299, 0.587, 0.114) * m;
    g = (g - dot( tex3D(tex,  p , n).rgb, vec3(0.299, 0.587, 0.114)) )/e.x;
    g -= n * dot(n, g);
                      
    return normalize( n + g*bf );
    
}





// ---------------------------------------------
// Tone mapping
// ---------------------------------------------
vec3 acesToneMapping( vec3 col )
{
    // https://www.shadertoy.com/view/XlKSDR
    // Narkowicz 2015, "ACES Filmic Tone Mapping Curve"
    const float a = 2.51;
    const float b = 0.03;
    const float c = 2.43;
    const float d = 0.59;
    const float e = 0.14;
    return (col * (a * col + b)) / (col * (c * col + d) + e);
}

vec3 filmicToneMapping( vec3 col )
{
    // Good reference
    // https://www.shadertoy.com/view/lslGzl
    col = max(vec3(0.), col - vec3(0.004));
    col = (col * (6.2 * col + .5)) / (col * (6.2 * col + 1.7) + 0.06);
    return col;
}

/*

	Tunnel6 a.k.a "Polish up"
    for cooperation of polish demoscene musicians, called:
    
    
    "Power Packed Alliance".
    -----------------------------------

	https://www.youtube.com/watch?v=_lSReW7eRI4
    http://www.pouet.net/prod.php?which=70247

	
    
    also check my chrome extension for Shadertoy:
    https://chrome.google.com/webstore/detail/shadertoy-unofficial-plug/ohicbclhdmkhoabobgppffepcopomhgl?hl=pl

*/

#define getNormal getNormalHex

#define INFINITY 1e32
#define FAR 30.
#define t TIME
#define mt TIME * 1.2 
#define FOV 90.0
#define FOG .7

#define PI 3.14159265
#define TAU (2*PI)
#define PHI (1.618033988749895)

vec3 os;

vec3 pal( in float ta, in vec3 a, in vec3 b, in vec3 c, in vec3 d ){return a + b*cos( 6.28318*(c*ta+d) );}

// 	3D noise function (IQ)
float noise3d(vec3 p)
{
	vec3 ip=floor(p);
    p-=ip; 
    vec3 s=vec3(7,157,113);
    vec4 h=vec4(0.,s.yz,s.y+s.z)+dot(ip,s);
    p=p*p*(3.-2.*p); 
    h=mix(fract(sin(h)*43758.5),fract(sin(h+s.x)*43758.5),p.x);
    h.xy=mix(h.xz,h.yw,p.y);
    return mix(h.x,h.y,p.z); 
}

float yC(float x) {
 	return cos(x * -.134) * 1. * sin(x * .13) * 15.+ noise3d(vec3(x * .01, 0., 0.) * 55.4);
}

float vol = 0.;

vec3 light = vec3(0.0);
vec3 opRep( vec3 p, vec3 c )
{
    return mod(p,c)-0.5*c;
}

void pR(inout vec2 p, float a) {
	p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

// Repeat in three dimensions
vec3 pMod3(inout vec3 p, vec3 size) {
	vec3 c = floor((p + size*0.5)/size);
	p = mod(p + size*0.5, size) - size*0.5;
	return c;
}

float opU2( float d1, float d2 ) {
    if (d1 < d2) return d1;
    return d2;
}

vec3 opU2( vec3 d1, vec3 d2 ) {
    if (d1.x < d2.x) return d1;
    return d2;
}

vec3 opS2( vec3 d1, vec3 d2 )
{	
    if (-d2.x > d1.x) return -d2;
    return d1;
}

float vmax(vec3 v) {
	return max(max(v.x, v.y), v.z);
}

float fBox(vec3 p, vec3 b) {
	vec3 d = abs(p) - b;
	return length(max(d, vec3(0))) + vmax(min(d, vec3(0)));
}

float fCross(vec3 p, vec3 size) {
    float obj = fBox(p, size);
    obj = opU2(obj, fBox(p, size.zxy));
    obj = opU2(obj, fBox(p, size.yzx));               
    return obj;
}

vec3 map(vec3 p) {
    p.y -= yC(p.x);
    vec3 
        obj = vec3(FAR, 0.0, 0.0),
        obj2 = obj,
        obj3 = obj;
    
    vec3 orgP = p;
    vec3 orgP3 = orgP;

	orgP3 = p;
    
    vec3 pp = pMod3(orgP3, vec3(2.4));
    p = orgP3;
    
    obj = vec3(
        fBox(p, vec3(1.05)), 
        1.0, 
        1.0
    );
    
    vec3 orgP2 = orgP;
    
    pR(orgP.zy, orgP.x / 12.);
	
    vec3 size = vec3(0.725 , 1.5, 1.275);
    
    p = opRep(orgP, vec3(0.35, 0.1, .4) + size.x + size.y + size.z);
    
    obj = opS2(
        obj, 
        vec3(                
            fCross(p, size), 
            0.0, 
            1.0
        )
    );
	
    size *= 1.2;
    p = opRep(orgP, vec3(0.35, 0.5, 0.1) + size.x + size.x + size.z);
    
    obj = opS2(
        obj, 
        vec3(                
            fCross(p, size), 
            0.0, 
            1.0
        )
    );
    p = orgP2;
    float n = noise3d(p);
  	pR(p.yz, p.x * .8 + n * 7.);
    p.y += .6;
   
    p = orgP2;

    obj = opS2(obj, vec3(fCross(p, vec3(1e32, .6, .6) ), 1., 1.)); 

	obj3.x = mix(-length(p.zy) + 1., obj3.x,  .6  *  n)- .1;
    obj3 = opU2(obj, vec3(fBox(p, vec3(1.1))));
	os = p;
    
    return obj3;
}

vec3 trace(vec3 prp, vec3 scp) {
    vec3 
        tr = vec3(0., -1., 0.),
        d;
    
    for (int i = 0; i < 164; i++) {	
        d = map(prp + scp * tr.x);
        tr.x += d.x * .4;

        if ((abs(d.x) < .0001) || (tr.x > FAR)) break;
    }
    
    tr.yz = d.yz;
	return tr;
    
}
vec3 traceRef(vec3 ro, vec3 rd) {
    vec3 
        tr = vec3(0., -1., 0.),
        d;
    
    for (int i = 0; i < 50; i++) {
        d = map(ro + rd * tr.x);
        tr.x += d.x;
        
        if (abs(d.x) < 0.0055 || tr.x> FAR) break;
    }
    
    tr.yz = d.yz;
    return tr;
}

float softShadow(vec3 ro, vec3 lp, float k) {
    const int maxIterationsShad = 12;
    vec3 rd = (lp - ro); 

    float shade = .1;
    float dist = 2.2;
    float end = max(length(rd), 0.001);
    float stepDist = end / float(maxIterationsShad);

    rd /= end;
    for (int i = 0; i < maxIterationsShad; i++) {
        float h = map(ro + rd * dist).x;
        shade = min(shade, smoothstep(0.0, 1.0, k * h / dist)); 
        dist += min(h, stepDist * 2.); 
        if (h < 0.001 || dist > end) break;
    }
    return min(max(shade, 0.), 1.0);
}


#define EPSILON .1
vec3 getNormalHex(vec3 pos)
{
	float d=map(pos).x;
	return normalize(
        vec3(
            map(
                pos+vec3(EPSILON,0,0)).x-d,
                map(pos+vec3(0,EPSILON,0)).x-d,
                map(pos+vec3(0,0,EPSILON)).x-d 
        	)
    	);
}

float getAO(in vec3 hitp, in vec3 normal)
{
    float dist = 0.2;
    vec3 spos = hitp + normal * dist;
    float sdist = map(spos).x;
    return clamp(sdist / dist, 0.0, 1.0);
}

vec3 getObjectColor(vec3 p, vec3 n, inout vec2 mat) {
    vec3 col = vec3(.4, 0.6, 1.);    
	
    col = pal( p.x * 0.01, vec3(0.5,0.5,0.5),vec3(0.5,0.5,0.5),vec3(1.0,1.0,1.0),vec3(0.0,0.10,0.20) );
   	col *= 1.+ pow(noise3d( os * 5.) * noise3d(p * 3.2), 1.);
    
    return col;
}

vec3 doColor( in vec3 sp, in vec3 rd, in vec3 sn, in vec3 lp, inout vec2 mat) {
	vec3 sceneCol = vec3(0.0);
    
    vec3 ld = lp - sp; 
    float lDist = max(length(ld), 0.001); 
    ld /= lDist;

    float atten = max(0.1, 2.0 / (1.0 + lDist * .525 + lDist * lDist * 1.05));
    float diff = max(dot(sn, ld), .1);
    float spec = pow(max(dot(reflect(-ld, sn), -rd), .5), 2.0);
   
    vec3 objCol = getObjectColor(sp, sn, mat);
    sceneCol += (objCol * (diff + 0.15) + vec3(.3, .4, .6) * spec * 1.) * atten;

    return sceneCol;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    vec2 uv = fragCoord.xy / RENDERSIZE.xy - .5;
	uv.x /= RENDERSIZE.y / RENDERSIZE.x;
    
    uv *= tan(radians (FOV) / 2.0);
    
    //vol = texture(iChannel0, vec2(.0, .25)).r  * 1.; 
    
    float 
        sk = sin(mt * .2) * 2.0,
        ck = cos(mt * .3) * 2.0,
        
        mat = 0.;
    
    light = vec3(0., 0., 11.);        
    
    vec3 sceneColor = vec3(0.);
    float camx = mt;
    
    vec3 
        vuv = vec3(0., 1., 0.),
    	ro = vec3(camx, yC(camx), 0.),
    	vrp =  vec3(camx + 1., yC(camx + 1.5), 0.),
    	vpn = normalize(vrp - ro),
    	u = normalize(cross(vuv, vpn)),
    	v = cross(vpn, u),
    	vcv = (ro + vpn),
    	scrCoord = (vcv + uv.x * u * RENDERSIZE.x/RENDERSIZE.y + uv.y * v),
    	rd = normalize(scrCoord - ro);        
	
    light = ro;
    
    vec3 lp = vrp;
    
	vec3 orgRO = ro,
         orgRD = rd;
	
    vec3 tr = trace(ro, rd), otr;
    
    float fog = smoothstep(FAR * FOG, 0., tr.x * 1.);
    
    ro += rd * tr.x;
    otr = ro;
    
    vec3 sn = getNormal(ro);	
    float ao = getAO(ro, sn);
   	
    sceneColor += doColor(ro, rd, sn, lp, tr.yz);
    
    float dist = tr.x;
    
    rd = reflect(rd, sn);
    tr = traceRef(ro + rd * .03, rd);
    ro += rd * tr.x;
    sn = getNormal(ro);
    sceneColor += doColor(ro, rd, sn, lp, tr.yz) * .3;

    fragColor = vec4(clamp(sceneColor * 1.3, 0.0, 1.0), tr.x / FAR);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}