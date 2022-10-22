//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** Common Code Begins ********

//#define PI 3.141592653589
#define E  2.718281828459
#define EPSILON 0.0000001

#define KEY_LEFT  37
#define KEY_UP    38
#define KEY_RIGHT 39
#define KEY_DOWN  40
#define KEY_SPACE 32
#define KEY_1     49
#define KEY_9     57

const vec3 v0 = vec3(0,0,0);
const vec3 vx = vec3(1,0,0);
const vec3 vy = vec3(0,1,0);
const vec3 vz = vec3(0,0,1);

const vec3 red   = vec3(0.8,0.0,0.0);
const vec3 green = vec3(0.0,0.5,0.0);
const vec3 blue  = vec3(0.2,0.2,1.0);
const vec3 white = vec3(1.0,1.0,1.0);

//  0000000   000       0000000   0000000     0000000   000      
// 000        000      000   000  000   000  000   000  000      
// 000  0000  000      000   000  0000000    000000000  000      
// 000   000  000      000   000  000   000  000   000  000      
//  0000000   0000000   0000000   0000000    000   000  0000000  

struct _text {
    ivec2 size;
    ivec2 adv;
} text;

struct _gl {
    vec2  uv;
    vec2  frag;
    vec2  mouse;
    vec2  mp;
    ivec2 ifrag;
    float aspect;
    vec4  color;
    int   option;
    float time;
} gl;

void initGlobal(vec2 fragCoord, vec3 resolution, vec4 mouse, float time)
{
    //text.size = ivec2(8,16);
    text.size = ivec2(16,32);
    text.adv  = ivec2(text.size.x,0);
    
    mouse.xy = min(_mouse.xy,resolution.xy);
    if (_mouse.z < 1.0)
    {
        if (_mouse.z > -1.0)
        	gl.mouse = resolution.xy*0.5;
       	else
            gl.mouse = _mouse.xy;
       
        
        //gl.mouse.x += sin(time*2.0)*resolution.y/32.0;
    }
    else gl.mouse = _mouse.xy;
    
	gl.mp = (2.0*abs(_mouse.xy)-vec2(resolution.xy))/resolution.y;    

    gl.aspect = resolution.x / resolution.y;
    gl.frag   = fragCoord;
    gl.ifrag  = ivec2(fragCoord);
    gl.uv     = (fragCoord+fragCoord-resolution.xy)/resolution.y;
}

float powi(int a, int b) { return pow(float(a), float(b)); }
float log10(float a) { return log(a)/log(10.0); }
float clamp01(float v) { return clamp(v, 0.0, 1.0); }

// 000   000   0000000    0000000  000   000  
// 000   000  000   000  000       000   000  
// 000000000  000000000  0000000   000000000  
// 000   000  000   000       000  000   000  
// 000   000  000   000  0000000   000   000  

float hash11(float p)
{
    p = fract(p * 0.1031);
    p *= p + 33.33;
    p *= p + p;
    return fract(p);
}

vec3 hash33(vec3 p3)
{
    p3 = fract(p3 * vec3(0.1031, 0.1030, 0.0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy + p3.yxx)*p3.zyx);
}

vec3 hash31(float p)
{
   vec3 p3 = fract(vec3(p) * vec3(.1031, .1030, .0973));
   p3 += dot(p3, p3.yzx+33.33);
   return fract((p3.xxy+p3.yzz)*p3.zyx); 
}

float hash12(vec2 p)
{
    vec3 p3  = fract(vec3(p.xyx) * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

float gradientNoise(vec2 v)
{
    return fract(52.9829189 * fract(dot(v, vec2(0.06711056, 0.00583715))));
}

// 000   000   0000000  000      
// 000   000  000       000      
// 000000000  0000000   000      
// 000   000       000  000      
// 000   000  0000000   0000000  

vec3 hsl2rgb( in vec3 c )
{
    vec3 rgb = clamp( abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );
    return c.z + c.y * (rgb-0.5)*(1.0-abs(2.0*c.z-1.0));
}

vec3 hsl(float h, float s, float l) { return hsl2rgb(vec3(h,s,l)); }

vec3 rgb2hsl( vec3 col )
{
    float minc = min( col.r, min(col.g, col.b) );
    float maxc = max( col.r, max(col.g, col.b) );
    vec3  mask = step(col.grr,col.rgb) * step(col.bbg,col.rgb);
    vec3 h = mask * (vec3(0.0,2.0,4.0) + (col.gbr-col.brg)/(maxc-minc + EPSILON)) / 6.0;
    return vec3( fract( 1.0 + h.x + h.y + h.z ),              
                 (maxc-minc)/(1.0-abs(minc+maxc-1.0) + EPSILON),  
                 (minc+maxc)*0.5 );                           
}

// 00     00   0000000   000000000  00000000   000  000   000  
// 000   000  000   000     000     000   000  000   000 000   
// 000000000  000000000     000     0000000    000    00000    
// 000 0 000  000   000     000     000   000  000   000 000   
// 000   000  000   000     000     000   000  000  000   000  

mat3 alignMatrix(vec3 dir) 
{
    vec3 f = normalize(dir);
    vec3 s = normalize(cross(f, vec3(0.48, 0.6, 0.64)));
    vec3 u = cross(s, f);
    return mat3(u, s, f);
}

// 00000000    0000000   000000000  
// 000   000  000   000     000     
// 0000000    000   000     000     
// 000   000  000   000     000     
// 000   000   0000000      000     

float rad2deg(float r) { return 180.0 * r / PI; }
float deg2rad(float d) { return PI * d / 180.0; }

mat3 rotMat(vec3 u, float angle)
{
    float s = sin(deg2rad(angle));
    float c = cos(deg2rad(angle));
    float i = 1.0-c;
    
    return mat3(
        c+u.x*u.x*i, u.x*u.y*i-u.z*s, u.x*u.z*i+u.y*s,
        u.y*u.x*i+u.z*s, c+u.y*u.y*i, u.y*u.z*i-u.x*s,
        u.z*u.x*i-u.y*s, u.z*u.y*i+u.x*s, c+u.z*u.z*i
        );
}

vec3 rotAxisAngle(vec3 position, vec3 axis, float angle)
{
    mat3 m = rotMat(axis, angle);
    return m * position;
}

//  0000000   000   000   0000000   000000000  
// 000   000  000   000  000   000     000     
// 000 00 00  000   000  000000000     000     
// 000 0000   000   000  000   000     000     
//  00000 00   0000000   000   000     000     

vec4 quatAxisAngle(vec3 axis, float angle)
{ 
    float half_angle = deg2rad(angle*0.5);
    return vec4(axis*sin(half_angle), cos(half_angle));
}

vec4 quatConj(vec4 q)
{ 
    return vec4(-q.x, -q.y, -q.z, q.w); 
}
  
vec4 quatMul(vec4 q1, vec4 q2)
{ 
    vec4 qr;
    qr.x = (q1.w * q2.x) + (q1.x * q2.w) + (q1.y * q2.z) - (q1.z * q2.y);
    qr.y = (q1.w * q2.y) - (q1.x * q2.z) + (q1.y * q2.w) + (q1.z * q2.x);
    qr.z = (q1.w * q2.z) + (q1.x * q2.y) - (q1.y * q2.x) + (q1.z * q2.w);
    qr.w = (q1.w * q2.w) - (q1.x * q2.x) - (q1.y * q2.y) - (q1.z * q2.z);
    return qr;
}

vec3 rotate(vec4 q, vec3 p)
{
    vec4 conj = quatConj(q);
    vec4 q_tmp = quatMul(q, vec4(p, 0));
    return quatMul(q_tmp, conj).xyz;
}

vec3 rotate(vec4 q, vec3 o, vec3 p)
{
    vec4 conj = quatConj(q);
    vec4 q_tmp = quatMul(q, vec4(p-o, 0));
    return o + quatMul(q_tmp, conj).xyz;
}

vec3 rotAxisAngleQuat(vec3 p, vec3 axis, float angle)
{ 
    vec4 qr = quatAxisAngle(axis, angle);
    return quatMul(quatMul(qr, vec4(p, 0)), quatConj(qr)).xyz;
}

vec3 rotRayAngle(vec3 p, vec3 ro, vec3 rd, float angle)
{ 
    return rotAxisAngle(p-ro, rd-ro, angle)+ro;
}

vec3 rotY(vec3 v, float d)
{
    float r = deg2rad(d);
    float c = cos(r);
    float s = sin(r);
    return vec3(v.x*c+v.z*s, v.y, v.z*c+v.x*s);
}

vec3 rotX(vec3 v, float d)
{
    float r = deg2rad(d);
    float c = cos(r);
    float s = sin(r);
    return vec3(v.x, v.y*c+v.z*s, v.z*c+v.y*s);
}

vec3 rotZ(vec3 v, float d)
{
    float r = deg2rad(d);
    float c = cos(r);
    float s = sin(r);
    return vec3(v.x*c+v.y*s, v.y*c+v.x*s, v.z);
}

//  0000000   00000000   0000000   00     00    
// 000        000       000   000  000   000    
// 000  0000  0000000   000   000  000000000    
// 000   000  000       000   000  000 0 000    
//  0000000   00000000   0000000   000   000    


vec3 posOnPlane(vec3 p, vec3 a, vec3 n)
{
    return p-dot(p-a,n)*n;
}

vec3 posOnRay(vec3 ro, vec3 rd, vec3 p)
{
    return ro + max(0.0, dot(p - ro, rd) / dot(rd, rd)) * rd;
}

bool rayIntersectsSphere(vec3 ro, vec3 rd, vec3 ctr, float r)
{
    return length(posOnRay(ro, rd, ctr) - ctr) < r;
}

//  0000000   00000000   
// 000   000  000   000  
// 000   000  00000000   
// 000   000  000        
//  0000000   000        

float opUnion(float d1, float d2)
{
    float k = 0.15;
    float h = clamp(0.5 + 0.5*(d2-d1)/k, 0.0, 1.0);
    return mix(d2, d1, h) - k*h*(1.0-h);
}

float opUnion(float d1, float d2, float k) 
{
    float h = clamp(0.5 + 0.5*(d2-d1)/k, 0.0, 1.0);
    return mix(d2, d1, h) - k*h*(1.0-h);
}

float opDiff(float d1, float d2) 
{
    float k = 0.05;
    float h = clamp(0.5 - 0.5*(d2+d1)/k, 0.0, 1.0);
    return mix(d1, -d2, h) + k*h*(1.0-h); 
}

float opDiff(float d1, float d2, float k) 
{
    float h = clamp(0.5 - 0.5*(d2+d1)/k, 0.0, 1.0);
    return mix(d1, -d2, h) + k*h*(1.0-h); 
}

//  0000000  0000000    
// 000       000   000  
// 0000000   000   000  
//      000  000   000  
// 0000000   0000000    

float sdSphere(vec3 p, vec3 a, float r)
{
    return length(p-a)-r;
}

float sdPlane(vec3 p, vec3 a, vec3 n)
{   
    return dot(n, p-a);
}

float sdPlane(vec3 p, vec3 n)
{   
    return dot(n, p);
}

float sdEllipsoid(vec3 p, vec3 a, vec3 r)
{
    p = p-a;
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float sdCone(vec3 p, vec3 a, vec3 b, float r1, float r2)
{
    vec3 ab = b-a;
    vec3 ap = p-a;
    float t = dot(ab,ap) / dot(ab,ab);
    t = clamp(t, 0.0, 1.0);
    vec3 c = a + t*ab;
    return length(p-c)-(t*r2+(1.0-t)*r1);      
}

float sdCylinder(vec3 p, vec3 a, vec3 b, float r)
{
  vec3  ba = b - a;
  vec3  pa = p - a;
  float baba = dot(ba,ba);
  float paba = dot(pa,ba);
  float x = length(pa*baba-ba*paba) - r*baba;
  float y = abs(paba-baba*0.5)-baba*0.5;
  float x2 = x*x;
  float y2 = y*y*baba;
  float d = (max(x,y)<0.0)?-min(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
  return sign(d)*sqrt(abs(d))/baba;
}

vec3 posOnPlane(vec3 p, vec3 n)
{
    return p-dot(p,n)*n;
}

float sdTorus(vec3 p, vec3 a, vec3 n, float rl, float rs)
{
    vec3 q = p-a;
    return length(vec2(length(posOnPlane(q, n))-rl,abs(dot(n, q))))-rs;
}


#define MOD3 vec3(.1031,.11369,.13787)

vec3 hash353(vec3 p3)
{
	p3 = fract(p3 * MOD3);
    p3 += dot(p3, p3.yxz+19.19);
    return -1.0 + 2.0 * fract(vec3((p3.x + p3.y)*p3.z, (p3.x+p3.z)*p3.y, (p3.y+p3.z)*p3.x));
}

float simplex_noise(vec3 p)
{
    const float K1 = 0.333333333;
    const float K2 = 0.166666667;
    
    vec3 i = floor(p + (p.x + p.y + p.z) * K1);
    vec3 d0 = p - (i - (i.x + i.y + i.z) * K2);
        
    vec3 e = step(vec3(0.0), d0 - d0.yzx);
	vec3 i1 = e * (1.0 - e.zxy);
	vec3 i2 = 1.0 - e.zxy * (1.0 - e);
    
    vec3 d1 = d0 - (i1 - 1.0 * K2);
    vec3 d2 = d0 - (i2 - 2.0 * K2);
    vec3 d3 = d0 - (1.0 - 3.0 * K2);
    
    vec4 h = max(0.6 - vec4(dot(d0, d0), dot(d1, d1), dot(d2, d2), dot(d3, d3)), 0.0);
    vec4 n = h * h * h * h * vec4(dot(d0, hash353(i)), dot(d1, hash353(i + i1)), dot(d2, hash353(i + i2)), dot(d3, hash353(i + 1.0)));
    
    return dot(vec4(31.316), n);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    
    vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;       
    
    float cc = 0.;
    for(float i = 1.; i<1.1; i+=.01){
        vec2 iuv = uv * i*1.1;
        float s = simplex_noise(vec3(iuv*10.*i, script_time*.25+script_bass_time*0.125 + i*5.));        
        s += script_time*.0125 + i*100.;  
        cc += fract(floor(s*5.)/5.) * .1;        
    }    
    
    fragColor = vec4(fract(cc));    
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}