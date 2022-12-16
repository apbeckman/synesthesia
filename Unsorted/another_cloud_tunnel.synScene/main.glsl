

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

			//******** BuffA Code Begins ********

const int nbSample = 100;
const float invNbSample = 1. / float(nbSample);
const float pi = 3.141592653589;
float time;

float fbm( vec3 p ) {
    float d = 0.;
    float amp = 0.5;
    
    for(int i=0; i<4; i++) {
        vec4 rnd = noised(p) * amp;
        d += rnd.w;
        
        p += rnd.xyz * amp;
        amp *= .5;
        p *= 2.;
    }
    return d;
}


float clouds( vec3 p) {
    p.z += time;
    float d = fbm(p) + dot(cos(p.xyz), sin(p.zxy))*.75;//length(p.xy)-2.;
	
    
    
	return d;
}

vec3 gradientCol(float t) {
    vec3 col = mix( vec3(1.,1.,1.), vec3(1.,.4,.1), 1.-exp(-t*3.) );
    
    col = mix(vec3(0.,.5,1.), col, 1.-exp(-abs(t-.02)*20.));
    
    return col;
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/RENDERSIZE.xy;
    vec2 v = uv*2.-1.;
	v.x *= RENDERSIZE.x/RENDERSIZE.y;
    
    
    // jittering
    float jitt = hash(uv);
    float jittRay = jitt * invNbSample;
    time = TIME+jitt/30.;
    
    // camera ray
    vec3 ro = vec3(-1.5,0.,0.);
    vec3 rd = normalize( vec3(v, 1.) );
    rd.xy = rotate(time*.1) * rd.xy;
    rd.zy = rotate(time*.05) * rd.zy;
    
    
    
    // sky
    vec3 col = vec3(.5,.8,1.)*10.;
    col += vec3(1.) * exp(-uv.y*5.);
    
    
    vec3 ld = normalize(vec3(0.,0.,1.));
        
    for(float i=0.; i<1.; i += invNbSample) {
        vec3 p = ro + rd * (1.-i+jittRay)*20.;
        
        float d = clouds(p);
        vec3 c = gradientCol(d*1.) * exp(-d*5.)*2.;
        
//        c += vec3(1.) * saturate( d - map(p+ld*.3) );
        
        col = mix(col, c, saturate(d*.6));
    }
    
    
    // Output to screen
    fragColor = vec4(col,1.0);
    //fragColor = vec4(gradientCol(uv.x),1.0);
	return fragColor; 
 } 





vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy;
    vec3 col = texture(BuffA, uv).rgb;
    
    
    // vignetting
    col *= saturate(pow( uv.x * uv.y * (1.-uv.x) * (1.-uv.y)*100., .2 ));
    
    // tonemapping
    col = acesToneMapping(col);
    
    
    // Output to screen
    fragColor = vec4(col,1.0);
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