

			//******** Common Code Begins ********


// math
#define PI 3.141592653589
#define saturate(x) clamp(x,0.,1.)
float hash( vec3 x );
float hash( vec2 p );
float hash( float p );
float hash2Interleaved( vec2 x );
float noise( vec3 x );
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

// Created by Stephane Cuillerdier - @Aiekick/2018
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

mat3 getRotXMat(float a){return mat3(1.,0.,0.,0.,cos(a),-sin(a),0.,sin(a),cos(a));}
mat3 getRotYMat(float a){return mat3(cos(a),0.,sin(a),0.,1.,0.,-sin(a),0.,cos(a));}
mat3 getRotZMat(float a){return mat3(cos(a),-sin(a),0.,sin(a),cos(a),0.,0.,0.,1.);}

mat3 m1;
mat3 m2;

vec2 path(float t)
{
	return vec2(cos(t*0.08 + cos(t*0.1)*2.), sin(t*0.12 + sin(t*0.05)*2.5)) * 4.;
}

float pattern(vec3 p)
{
	p = abs(fract(p*.3) - 0.5);
	return length(max(abs(p.x), abs(p.y)) - p.z);
}

float map(vec3 p)
{
    vec2 pa = path(p.z); 
	float a = pa.x * pa.y * 0.1;
	
    p.xy -= pa;
    p.xy *= mat2(cos(a),-sin(a),sin(a),cos(a));
    
	vec3 qm1 = p*m1, qm2 = p*m2;
	
    float d0 = min(pattern(qm1), pattern(qm2));
    float d1 = min(pattern(qm1*3.), pattern(qm2*2.));
    
   	float dist0 = (1.-clamp(d0,0.,1.));
	float dist1 = (1.-clamp(d1,0.,1.))*d0;

    float di0 = length(p.xy); // tunnel
    float di1 = abs(cos(p.x*0.3)*p.y); // planes
    float di2 = 1.7 - 0.4 * cos(p.z*0.1) - abs(cos(p.x*0.5)*p.x + sin(p.y*0.5)*p.y); // planes
    
	p.z *= 0.5;
	
	vec2 cs = vec2(cos(p.z), sin(p.z)) * 0.5 + 0.5;
	
    return mix(di2, 1.47-mix(di0, di1, cs.y), cs.x) - dist0*0.75 - dist1*2.25;
}

vec3 nor( vec3 pos, float k )
{
	vec3 eps = vec3( k, 0., 0. );
	vec3 nor = vec3(
	    map(pos+eps.xyy) - map(pos-eps.xyy),
	    map(pos+eps.yxy) - map(pos-eps.yxy),
	    map(pos+eps.yyx) - map(pos-eps.yyx) );
	return normalize(nor);
}

// return color from temperature 
//http://www.physics.sfasu.edu/astro/color/blackbody.html
//http://www.vendian.org/mncharity/dir3/blackbody/
//http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color.html
vec3 blackbody(float Temp)
{
	vec3 col = vec3(255.);
    col.x = 56100000. * pow(Temp,(-3. / 2.)) + 148.;
   	col.y = 100.04 * log(Temp) - 623.6;
   	if (Temp > 6500.) col.y = 35200000. * pow(Temp,(-3. / 2.)) + 184.;
   	col.z = 194.18 * log(Temp) - 1448.6;
   	col = clamp(col, 0., 255.)/255.;
    if (Temp < 1000.) col *= Temp/1000.;
   	return col;
}

// get density of the df at surfPoint // ratio between constant step and df value
float SubDensity(vec3 surfPoint, float prec, float ms) 
{
	vec3 n;
	float s = 0.;
    for (int i=0;i<8;i++)
	{
		n = nor(surfPoint,prec); 
		surfPoint = surfPoint - n * ms; 
		s += map(surfPoint);
	}
	
	return 1.-s/(ms*8.); // s < 0. => inside df
}

float SubDensity(vec3 p, float s) 
{
	vec3 n = nor(p,s);
	return map(p - n * s);
}

// from shane
vec3 tex3D( sampler2D tex, in vec3 p, in vec3 n )
{
    n = max((abs(n) - .2)*7., .001);
    n /= (n.x + n.y + n.z );  
    p = (texture(tex, p.yz)*n.x + texture(tex, p.zx)*n.y + texture(tex, p.xy)*n.z).xyz;
    return p*p;
}

// from shane
vec3 doBumpMap( sampler2D tx, in vec3 p, in vec3 n, float bf)
{
    const vec2 e = vec2(0.001, 0);
    mat3 m = mat3( tex3D(tx, p - e.xyy, n), tex3D(tx, p - e.yxy, n), tex3D(tx, p - e.yyx, n));
    vec3 g = vec3(0.299, 0.587, 0.114)*m; // Converting to greyscale.
    g = (g - dot(tex3D(tx,  p , n), vec3(0.299, 0.587, 0.114)) )/e.x; g -= n*dot(n, g);
    return normalize( n + g*bf ); // Bumped normal. "bf" - bump factor.
}

vec4 shade(vec3 ro, vec3 rd, float d, vec3 lp)
{
	vec3 p = ro + rd * d;											// surface point
	float sb = SubDensity(p, 0.01, 0.1);							// deep subdensity (10 iterations)
	vec3 bb = blackbody(100.*sb+100.);								// bb
	vec3 ld = normalize(lp-p); 										// light dir
	vec3 n = nor(p, .01);											// normal at surface point
	n = doBumpMap(image4, -p*0.5, n, 0.015);
	vec3 refl = reflect(rd,n);										// reflected ray dir at surf point 
	float amb = 0.08; 												// ambiance factor
	float diff = clamp( dot( n, ld ), 0.0, 1.0 ); 					// diffuse
	float fre = pow( clamp( 1. + dot(n,rd),0.0,1.0), 16. ); 			// fresnel
	float spe = pow(clamp( dot( refl, ld ), 0.0, 1.0 ),25.);		// specular
	float sss = 1. - SubDensity(p*0.1, 0.1) * 0.5; 							// one step sub density of df
	return vec4(
        (diff + fre + bb.x * sss) * amb + diff * 0.5, 
        (diff + fre + bb * sb + sss * 0.3) * amb + spe * 0.6 - diff * sss * 0.05	
    );
}

vec4 render(vec3 ro, vec3 rd, float time)
{
    vec3 lp = vec3(path(time + 3.),time + 3.);
    
    float s = 1.;
    float d = 0.;
	
	for(int i=0;i<80;i++)
    {      
        if (d*d/s>1e5) break;
        d += s = map(ro+rd*d) * .6;
	}
	
    vec3 p = ro+rd*d;
    vec3 n = nor(p, 0.01);
        
    vec4 f = shade(ro, rd, d, lp);
	
	f = f.zyww + f.x*0.2;

    float fog = 1.0-exp( -0.01*d*d );
    f = mix( f, vec4(0.8), fog);
       
	f = mix(f, f.grba, sin(fog*5.));
   	return sqrt(f*f*f*1.5);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec4 f = vec4(0);
    
    vec2 g = fragCoord.xy;
    vec2 si = RENDERSIZE.xy;
    
    mat3 mx = getRotXMat(-7.);
	mat3 my = getRotYMat(-5.);
	mat3 mz = getRotZMat(-3.);
	
    m1 = mx * my * mz;
    m2 = m1 * m1;
	
    float time = TIME * 5.;
    
    vec3 cu = vec3(0,1,0);
  	vec3 cv = vec3(path(time + .1),time + .1);
	
  	vec2 uv = (g+g-si)/si.y;
    
    vec3 ro = vec3(path(time),time);
    
	vec2 fov = vec2(0.75,0.9);
  	vec3 z = normalize(cv-ro);
    vec3 x = normalize(cross(cu,z));
  	vec3 y = cross(z,x);
  	vec3 rd = normalize(z + uv.x*x*fov.x + uv.y*y*fov.y);
    
    fragColor = render(ro, rd, time);
}

void mainVR(out vec4 fragColor, vec2 fc, vec3 ro, vec3 rd) 
{
    mat3 mx = getRotXMat(-7.);
	mat3 my = getRotYMat(-5.);
	mat3 mz = getRotZMat(-3.);
	
    m1 = mx * my * mz;
    m2 = m1 * m1;
	
    float time = TIME * 5.;
    
    ro += vec3(path(time), time);

    fragColor = render(ro, rd, time);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}