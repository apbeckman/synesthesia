vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** Common Code Begins ********

//#define PI 3.14159265

float sqr(float x)
{
    return x*x;
}

float sdSph( vec3 p, float s )
{
  return length(p)-s;
}

float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

void mengerFold(inout vec4 z) 
{
	float a = min(z.x - z.y, 0.0);
	z.x -= a; z.y += a;
	a = min(z.x - z.z, 0.0);
	z.x -= a; z.z += a;
	a = min(z.y - z.z, 0.0);
	z.y -= a; z.z += a;
}

vec2 minid(vec2 a, vec2 b)
{
    return (a.x<b.x)?a:b;
}

struct mat
{
    vec3 albedo;
    vec3 emiss;
    float rough;
    float metal;
};
    
mat materials[3] = mat[3](mat(vec3(1.,0.02,0), vec3(0), 0.4, 0.2),
    			          mat(vec3(0.9,0.9,1.), vec3(0), 0.6, 0.2),
                          mat(vec3(0.5,0.9,0.2), vec3(0), 0.1, 0.4));

                  
#define MAXD 64.
#define MAXI 128
#define FOV 1.5
float LOD;


mat3 getCamera(vec2 angles)
{
   mat3 theta_rot = mat3(1,   0,              0,
                          0,   cos(angles.y),  sin(angles.y),
                          0,  -sin(angles.y),  cos(angles.y)); 
        
   mat3 phi_rot = mat3(cos(angles.x),   sin(angles.x), 0.,
        		       -sin(angles.x),   cos(angles.x), 0.,
        		        0.,              0.,            1.); 
        
   return theta_rot*phi_rot;
}

vec3 getRay(vec2 angles, vec2 pos)
{
    mat3 camera = getCamera(angles);
    return normalize(transpose(camera)*vec3(FOV*pos.x, 1., FOV*pos.y));
}

float NGGX(vec3 n, vec3 h, float a)
{
    float a2 = sqr(a);
    return a2/(PI*sqr( sqr( max(dot(n,h),0.) )*(a2-1.) + 1.));
}

float GGX(vec3 n, vec3 o, float a)
{
    float ndoto = max(dot(n,o),0.);
    return ndoto/mix(1., ndoto, sqr(a+1.)*0.125);
}

float GS(vec3 n, vec3 i, vec3 o, float a)
{
    return GGX(n,i,a)*GGX(n,o,a);
}

vec3 IR(float D, float k0, vec3 k1)
{
    //interference effect here ->
    return (0.25+ k0*( 1. - cos(2.*PI*pow(vec3(D), -k1)) ))/D ;
}

vec3 BRDF(vec3 i, vec3 o, vec3 n, mat m)
{
    vec3 h = normalize(i + o);
    vec3 F0 = mix(vec3(0.04), m.albedo, m.metal);
    vec3 FS = F0 + (1.0 - F0) * pow(1.0 - max(dot(h, i), 0.0), 5.0);
    vec3 DFG = NGGX(n,h,m.rough)*GS(n,i,o,m.rough)*FS;
    float denom = max(dot(n, i), 0.001) * max(dot(n, o), 0.001);
    return (m.albedo*(1.-FS)/PI +
            DFG*IR(denom, 1., vec3(1.,1.1,1.2)))*max(0., dot(n,o));
}

mat4 rotationMatrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}


vec3 sky(vec3 r)
{
    vec3 c = vec3(.009, .288, .828);
	c = mix(vec3(1.), c, .9);
	c *= .5;
    float atmo = tanh(10.*(r.z-0.05))*0.4 + 0.5 + 0.1*r.z;
    
    vec3 g = vec3(atmo);  
    vec3 A0 = pow(c, g);
    vec3 B0 = 1.-pow(vec3(1.)-c, 1.-g);
    
    vec3 A = A0*(1.-A0);
    vec3 B = B0*(1.-B0);
    
    return mix(A, B, g);
}

//const float PI = 3.14159265359;

float sdBox1( vec3 p, vec3 b )
{
  vec3 d = abs(p) - b;
  return length(max(d,0.0)); 
}


float dBar(vec2 p, float width) {
    vec2 d = abs(p) - width;
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.)) + 0.01 * width;
}

float dCrossBar(vec3 p, float x) {
    float bar_x = dBar(p.yz, x);
    float bar_y = dBar(p.zx, x);
    float bar_z = dBar(p.xy, x);
    return min(bar_z, min(bar_x, bar_y));
}

mat2 Rot(float a) {
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

float sdOctahedron( vec3 p, float s)
{
  p = abs(p);
  return (p.x+p.y+p.z-s)*0.57735027;
}


float dMengerSponge(vec3 p) 
{
 float d = sdBox1(p, vec3(0.6));
 float itt = 5.;
 float one_third = 2. / itt;
 for (float i = 0.0; i < itt; i++) {
  float k = pow(one_third, i);
  float kh = k * 1.;
  d = max(d, -dCrossBar(mod(p + kh, k * 2.) - kh, k * one_third));
 }
 return d;
}
float dMengerSponge2(vec3 p) 
{
 float d = sdBox1(p, vec3(0.5));
 float itt = 4.;
 float one_third = 1. / itt;
 for (float i = 0.0; i < itt; i++) {
  float k = pow(one_third, i);
  float kh = k * 1.;
  d = max(d, -dCrossBar(mod(p + kh, k * 2.) - kh, k * one_third));
 }
 return d;
}

vec2 condmin(in vec2 d1, in vec2 d2) {
return vec2(min(d1.x, d2.x), mix(d1.y, d2.y, step(d2.x, d1.x)));
}



float g1;
vec2 GetDist(vec3 p) {	
    
 float gap = 1.;
 p.xyz = mod(p.xyz + gap,2.0 * gap) - gap;
 vec2 d;
 d=vec2(1.0);  
 vec3 p5 = p;
 vec2 dm1= vec2(dMengerSponge(p5),2);
 p = abs(p-.4);   
 vec3 p1 = p-vec3(-.0,1.,1.0);
 vec3 p2 = p-vec3(1.,0.0,1.0);
 vec3 p4 = p-vec3(1.,1.0,.0);   
 vec2  dm2=vec2(dMengerSponge2(p1),1);
 vec2  dm3=vec2(dMengerSponge2(p2),1);
 vec2  dm4=vec2(dMengerSponge2(p4),1);
 d = condmin( d,dm1);
 d = condmin( d,dm2);
 d = condmin( d,dm3);
 d = condmin( d,dm4);
 g1 +=1./(.001+pow(abs(dm1.x),10.));
 return d;
}

vec2 RayMarch (vec3 ro, vec3 rd) 
 {
 vec2 h, t=vec2( 0.);
 for (int i=0; i<15; i++) 
  {
 h = GetDist(ro + t.x * rd);
 if(h.x<0.001||t.x>4. ) break;
  t.x+=h.x;t.y=h.y;
 }
 if(t.x>4.) t.x=0.;
 return t;
}

vec3 GetNormal(vec3 p){
vec2 e = vec2(.00035, -.00035); 
return normalize(
 e.xyy * GetDist(p + e.xyy).x + 
 e.yyx * GetDist(p + e.yyx).x + 
 e.yxy * GetDist(p + e.yxy).x + 
 e.xxx * GetDist(p + e.xxx).x);
}




vec3 cameraPath(float t) {
 t *= PI *.5 ;
 float t2 =  cos(t)+0.;
 float c = cos(t*2.);
 float x = 0.;
 float y = 0.;
 float z= 0.;
 if (t2<0.){
  x = 1. /1. + 2. +c;
 };       
 if (t2>0.){
  y = 1. /1. + 2. +c;
 };
 vec3 xyz =vec3(x,y,z);
 return xyz;
}




vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

 vec2 uv = (fragCoord-.5*RENDERSIZE.xy)/RENDERSIZE.y;
 vec2 m = iMouse.xy/RENDERSIZE.xy;
 vec3 col = vec3(0);
 float s = mod(smoothTime * 0.125, 1.0);
 float t = 2. * (2.0 * s - s * s);  
 vec3 cameraPos = cameraPath(t);
 vec3 ro = vec3(-cameraPos);
 vec3 rd = normalize(vec3(uv.x, uv.y, 0.4));
 float the = TIME *0.3;
 rd.xz *= mat2(cos(the), -sin(the), sin(the), cos(the));
 rd.yz *= mat2(cos(the), -sin(the), sin(the), cos(the));
 vec2 d = RayMarch(ro, rd);
 float t2;
 t2=d.x;   
 if(t2>10.)
 {
  vec3 p = ro + rd * t2;
  vec3 baseColor = vec3(0.,0.0,0.);
  if(d.y==1.) baseColor=vec3(abs(sin(TIME)),0.,0);
  if(d.y==2.) baseColor=vec3(.0,.9,.0);
  

 }

 float fog = 1. / (1. + d.x * d.x * 1.);
 col *= vec3(fog);  
 vec3 sky = vec3(1., 1., 1.);
 col = mix(sky, col, 1./(d.x*d.x/1./.1*.1+1.)); 
 fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}