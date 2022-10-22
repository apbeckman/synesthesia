

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

                  
#define MAXD 32.
#define MAXI 96
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

#define fracaxis normalize(vec3(2.+sin(smoothTime*0.125),2,3.+cos(smoothTime*0.125)))
#define fracangle 1.
#define fracscale 1.2
#define fracshift vec4(0.05*sin(0.2*smoothTime), 0.1,0.+0.05*cos(0.2*smoothTime),0.)

//fractal matrix
mat4 fmat;

float sdFract(vec3 p)
{
    vec4 cp = vec4(p,1.);
    for(int i = 0; i < 16; i++)
    {
        cp = fmat*cp + fracshift;
        mengerFold(cp);
    }
    return sdBox(cp.xyz + fracshift.xyz, vec3(0.5))/cp.w;
}

vec3 sdFractCol(vec3 p)
{
    p.xy = mod(p.xy + vec2(2.), vec2(4.)) - vec2(2);
    vec4 cp = vec4(p,1.);
    vec3 c = vec3(1,2,1);
    for(int i = 0; i < 6; i++)
    {
        cp = fmat*cp + fracshift;
        mengerFold(cp);
        c = min(c,cp.xyz); 
    }
    return tanh(abs(c));
}

vec2 map(vec3 p)
{
    p.xy = mod(p.xy + vec2(2.), vec2(4.)) - vec2(2);
    vec2 d = vec2(sdFract(p), 0);
    d = minid( vec2(p.z + 1., 1), d);
    return d;
}

vec4 normal(vec3 p, float dx) {
	const vec3 k = vec3(1,-1,0);
	vec4 r = k.xyyx*map(p + k.xyy*dx).x +
			 k.yyxx*map(p + k.yyx*dx).x +
			 k.yxyx*map(p + k.yxy*dx).x +
			 k.xxxx*map(p + k.xxx*dx).x;
    //the normal and the averaged distance
    return vec4(normalize(r.xyz), r.w*0.25);
}

vec4 march(vec3 p, vec3 r)
{
    float td = 0.;
    vec2 d;
    for(int i = 0; i < MAXI; i++)
    {
        d = map(p + td*r);
        if(d.x <= LOD*td) break;
        if(td > MAXD) 
        {
            d.y = -1.;
            break;
        }
        td += d.x;
    }
    td += d.x - LOD*td; //better surface
    return vec4(p + td*r, d.y); //position and ID
}

#define AOscale 30.
float AO(vec3 p, vec3 n, float td)
{
    vec2 de = map(p+AOscale*n*td*LOD);
    return max(0.,tanh(de.x/(AOscale*td*LOD)));
}


vec3 ACESFilm(vec3 x)
{
    float a = 2.51;
    float b = 0.03;
    float c = 2.43;
    float d = 0.59;
    float e = 0.14;
    return (x*(a*x+b))/(x*(c*x+d)+e);
}

vec3 HDRmapping(vec3 color, float exposure)
{
	// Exposure tone mapping
    vec3 mapped = ACESFilm(color * exposure);
    // Gamma correction 
    return pow(mapped, vec3(1.0 / 2.2));
}

vec4 renderMainImage() {
	vec4 col = vec4(0.0);
	vec2 pos = _xy;

    LOD = 1.84/max(RENDERSIZE.x,RENDERSIZE.y);
    fmat = fracscale*rotationMatrix(fracaxis, fracangle);
    
    vec2 uv = (pos - 0.5*RENDERSIZE.xy)/RENDERSIZE.x;
    vec2 angle = (_mouse.z>0.)?(vec2(2.*PI, PI)*_mouse.xy/RENDERSIZE.xy - vec2(0,PI*0.5)):vec2(PI+0.2*TIME, 0.);
	vec3 r = getRay(angle, uv), cr = getRay(angle, vec2(0));
    vec3 p = -0.6*cr;
    
    vec4 res = march(p, r);
    float td = distance(res.xyz,p);
    p = res.xyz;
    
    vec3 L = vec3(0,1,0);
   
    col.xyz = vec3(0.);
    if(res.w < 0.)
    {
        col.xyz = sky(r);
    }
    else
    {
        vec4 n = normal(p, td*LOD);
    	mat a = materials[int(res.w)];
        if(res.w == 0.)
        	a.albedo = sdFractCol(p - n.xyz*n.w);
        float ao = AO(p, n.xyz, td);
        //a few light sources
       	col.xyz += BRDF(-r, vec3(cos(0.5),0,sin(0.5)), n.xyz, a);
        col.xyz += BRDF(-r, vec3(0,cos(1.),sin(1.)), n.xyz, a);
        col.xyz += BRDF(-r, vec3(0,sin(1.),cos(1.)), n.xyz, a)*3.;
        col.xyz += BRDF(-r, vec3(0,-sin(1.),-cos(1.)), n.xyz, a);
        col.xyz*=ao*sky(vec3(1.,0.,0.));
    }
    col = vec4(HDRmapping(col.xyz, 6.),1.);
	return col; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}