

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

/*
* License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
* Created by bal-khan
*/

// blackbody by aiekick : https://www.shadertoy.com/view/lttXDn

// -------------blackbody----------------- //

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

// -------------blackbody----------------- //

#define I_MAX	150
#define E		0.00001

float	sdTorus( vec3 p, vec2 t );
vec2	march(vec3 pos, vec3 dir);
vec3	camera(vec2 uv);
vec2	rot(vec2 p, vec2 ang);
float	mylength(vec3 p);
float	mylength(vec2 p);
void	rotate(inout vec2 v, float angle);

float	t;
vec3	h;
float	mind;
float	lit;
float	test;
vec3	_color_one, _color_two;

vec4 renderMainImage() {
	vec4 o = vec4(0.0);
	vec2 f = _xy;

    h = vec3(0.);
    t = TIME;
    vec2 R = RENDERSIZE.xy,
          uv  = vec2(f-R/2.) / R.y;
    
    _color_two = vec3(
		sin(.25*TIME*.25+2.08)
        ,
        sin(.25*TIME*.5+2.57)
        ,
        sin(.25*TIME*.75+3.00)
    );
    
	vec3	dir = camera(uv);
    vec4	col = vec4(0.0);
    vec3 pos = vec3(.01, .01, 5.0);

    vec2	inter = (march(pos, dir));
    col.xyz = 1.*blackbody((inter.y*.06125)*100.);
   	col.xyz += h;
    col.xyz += vec3(.2, .2, .2) * inter.x * .0125 * max(.01, 70./max(inter.y, .001) );
    col.xyz += (1.-_color_one) * min(inter.y, 70.) * .003;
    o.xyz = col.xyz;
    return o;
}

float	scene(vec3 p)
{
    mind = 1e5;
    lit  = 1e5;
    p.z-= -30.;
    p.z -= smoothTime*2.;
	vec3	ap = p;
    
    
    p.xy += vec2(cos(p.z*.25), sin(p.z*.25 ) );
    p.xy += vec2(cos(-p.z*.125), sin(-p.z*.125 ) );
    p.xy += vec2(cos(p.z*.0625), sin(p.z*.0625 ) );
    p.xy += vec2(cos(-p.z*.03), sin(-p.z*.03 ) );
    
    p.xy += 
        .175
        *
        vec2(
            cos(TIME*5.+50.*(atan(ap.x, ap.y)*1.-TIME*-.05+p.z*.0251*1.0) )
            ,
            sin(TIME*5.+50.*(atan(ap.x, ap.y)*1.+TIME*-.05+p.z*.0251*1.0) )
        );
    
    float ata = atan(p.x,p.y);
    _color_one = vec3(
    	sin(cos(3.*ata+(p.z*2.)/4.)+7.7+0.00)
        ,
        sin(cos(3.*ata+(p.z*2.)/8.)+7.7-.57)
        ,
        sin(cos(3.*ata+(p.z*2.)/16.)+7.7-1.04)
    );
    
    mind = min(mind
               ,
               max( (length(p.xy)-10.), -(length(p.xy)-9.) )
               );
    
    lit = min(lit, length( (p.xy))-10.);
    lit = max(lit, -(length(p.xy)-9.125) );
    
    lit = max(lit, (length(-1.+1.*cos(p.zz*5.+smoothTimeB*.5) )-.2) );
    h += 1.*_color_two*1./max(lit*lit*10000.1+.0+test*.1251, .1);

    mind = max(mind, -(length(-1.+1.*cos(p.zz*5.+smoothTimeC*.5) )-.2) );
	float tmp = mind;
    
    mind = min(mind, lit);
    
    h += _color_one * 1./max(tmp*tmp*100. + 50. + test, .01);

    return(mind);
}

vec2	march(vec3 pos, vec3 dir)
{
    vec2	dist = vec2(0.0);
    vec3	p = vec3(0.0);
    vec2	s = vec2(0.0);
	vec3	dirr;

    for (int i = -1; i < I_MAX; ++i)
    {
        dirr = dir;
        rotate(dirr.xz, +dist.y*.000+.25*sin(t*.25+dist.y*.005) );
    	p = pos + dirr * dist.y;
        test = dist.y;
        dist.x = scene(p);
        dist.y += dist.x;
        if (dist.x < E)
        {
           break;
	    }
        s.x++;
    }
    s.y = dist.y;
    return (s);
}

void rotate(inout vec2 v, float angle)
{
	v = vec2(cos(angle)*v.x+sin(angle)*v.y,-sin(angle)*v.x+cos(angle)*v.y);
}

float	mylength(vec2 p)
{
	return max(abs(p.x), abs(p.y));
}

float	mylength(vec3 p)
{
	return max(max(abs(p.x), abs(p.y)), abs(p.z));
}

float sdTorus( vec3 p, vec2 t )
{
	vec2 q = vec2(length(p.xy)-t.x,p.z);

    return length(q)-t.y;
}

vec3	camera(vec2 uv)
{
    vec3 dir;
    float   fov = 1.;
	vec3    forw  = vec3(0.0, 0.0, -1.0);
	vec3    right = vec3(1.0, 0.0, 0.0);
	vec3    up    = vec3(0.0, 1.0, 0.0);

    dir = (normalize((uv.x) * right + (uv.y) * up + fov * forw));
	return dir; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}