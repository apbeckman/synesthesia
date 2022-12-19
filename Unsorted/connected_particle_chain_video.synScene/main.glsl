

			//******** Common Code Begins ********

#define size RENDERSIZE.xy
#define pixel(a, p) texture(a, p/vec2(textureSize(a,0)))
#define texel(a, p) texelFetch(a, ivec2(p), 0)
#define ch0 BuffA
#define ch1 BuffB
#define ch2 BuffC
#define ch3 syn_UserImage
#define PI 3.14159265

#define N ivec2(50,50)
#define dt 1.


const int tot_n = N.x*N.y;

float hash11(float p)
{
    p = fract(p * 15.1031);
    p *= p + 33.33;
    p *= p + p;
    return fract(p) - 0.5;
}

float hash12(vec2 p)
{
	vec3 p3  = fract(vec3(p.xyx) * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}


vec2 hash21(float p)
{
	vec3 p3 = fract(vec3(p) * vec3(.1031, .1030, .0973));
	p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.xx+p3.yz)*p3.zy);

}

vec2 hash22(vec2 p)
{
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+33.33);
    return fract((p3.xx+p3.yz)*p3.zy);
}

const int k = 1664525;  

ivec4 hash( ivec4 x )
{
    x = ((x>>8)^x.wxyz)*k;
    x = ((x>>8)^x.wxyz)*k;
    x = ((x>>8)^x.wxyz)*k;
    x = ((x>>8)^x.wxyz)*k;
    return ivec4(x);
}

ivec2 i2xy(int id)
{
    return ivec2(id%N.x, id/N.x);
}

int xy2i(ivec2 p)
{
    return p.x + p.y*N.x;
}

ivec2 cross_distribution(int i)
{
    return (1<<(i/4)) * ivec2( ((i&2)/2)^1, (i&2)/2 ) * ( 2*(i%2) - 1 );
}

float sdSegment( in vec2 p, in vec2 a, in vec2 b )
{
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}

			//******** BuffA Code Begins ********

//particle buffer

int cid;

ivec4 getA(int id)
{
    return floatBitsToInt(texel(ch1, i2xy(id)));
}

vec4 getparticleA(int id)
{
    return texel(ch0, i2xy(id));
}

float F(float d)
{
    return tanh(0.1*TIME)*(1.*exp(-0.1*d) - 5.*exp(-0.2*d));
}

float imageV(vec2 p)
{
    return 1.-2.*texture(ch3, vec2(1., 1.)*p/size).x;
}

vec2 imageF(vec2 p)
{
    vec3 d = vec3(-1,0,1);
    return vec2(imageV(p+d.zy) - imageV(p+d.xy), imageV(p+d.yz) - imageV(p+d.yx));
}

vec2 Fv(vec4 p0, int pid)
{
    if(pid < 0 || pid >= tot_n || pid == cid) return vec2(0); 
   	vec4 p1 = getparticleA(pid);
    vec2 flocking_force = 0.3*(p1.zw - p0.zw)/(1.+1.*distance(p1.xy,p0.xy));
    vec2 interaction = normalize(p1.xy-p0.xy)*F(distance(p1.xy,p0.xy));
    return flocking_force + interaction;
}

float irad;

vec2 Fspring(vec4 p0, int pid)
{
    if(pid < 0 || pid >= tot_n || pid == cid) return vec2(0); 
   	vec4 p1 = getparticleA(pid);
    vec2 interaction = normalize(p1.xy-p0.xy)*(distance(p1.xy,p0.xy)- 2.*PI*irad/float(tot_n) - 4.*tanh(0.1*TIME));
    return interaction;
}

vec4 renderPassA() {
	vec4 U = vec4(0.0);
	vec2 pos = _xy;

    ivec2 p = ivec2(pos);
    if(p.x < N.x && p.y < N.y)
    {
        irad = 0.3*size.y;
        pos = floor(pos);
        //this pixel value
        U = texel(ch0, pos);
        int id = xy2i(p);
        cid = id;
        
        //this pixel value
        if(FRAMECOUNT<=10)
        {
            float t = 2.*PI*float(id)/float(tot_n);
            U.xy = size*vec2(0.3,0.5) + irad*(1.+0.4*sin(20.*t))*vec2(cos(t),sin(t));
			U.zw = 1.*(hash22(3.14159*pos) - 0.5);
      		return U;
        }
        
        //neighbors
   		ivec4 cp = getA(id);
   	  
        vec2 F = Fv(U, cp.x) +
            	 Fv(U, cp.y) +
            	 Fv(U, cp.z) +
                 Fv(U, cp.w) +
                 Fspring(U, (id + 1)%tot_n) +
                 Fspring(U, (id - 1)%tot_n) +
            	 -4.*imageF(U.xy);
        
        if(_mouse.z > 0.) 
        {
            float d = distance(_mouse.xy, U.xy);
            F += 2.*normalize(_mouse.xy - U.xy)/(sqrt(d)+2.);
        }
        
        U.zw += (-0.002*U.zw + 0.5*F)*dt;
        U.xy += U.zw*dt;
        
        //border conditions
        if(size.x - U.x < 2.) U.z = -abs(U.z);
        if(U.x < 2.) U.z = abs(U.z);
        if(size.y - U.y < 2.) U.w = -abs(U.w);
        if(U.y < 2.) U.w = abs(U.w);
 
        
    }
    else discard;
	return U; 
 } 


			//******** BuffB Code Begins ********

//sorting closest 4 particles in axis directions that make a bounding box
//only in particle space, texture buffer not needed

ivec4 getB(int id)
{
    return floatBitsToInt(texel(ch1, i2xy(id)));
}

vec4 save(ivec4 v)
{
    return intBitsToFloat(v);
}

ivec4 u; //ids
vec4 d; //distances
vec2 pos; //this particle position
int tid;

vec4 getparticleB(int id)
{
    return texel(ch0, i2xy(id));
}

//insertion sort
void sort(int utemp)
{
    if(utemp == tid || utemp < 0) return;
       
   	vec4 part = getparticleB(utemp);
    vec2 dx = part.xy - pos;
    float dtemp = length(dx);
    //sorting
    if(dx.x > abs(dx.y))
    {
        if(d.x > dtemp) 
        {
            d.x = dtemp;
        	u.x = utemp;
        }
    }
    else if(dx.x < -abs(dx.y))
    {
        if(d.y > dtemp) 
        {
            d.y = dtemp;
        	u.y = utemp;
        }
    }
    else if(dx.y > abs(dx.x))
    {
        if(d.z > dtemp) 
        {
            d.z = dtemp;
        	u.z = utemp;
        }
    }
    else if(d.w > dtemp) 
    {
        d.w = dtemp;
        u.w = utemp;
    }
}

void sortneighbor(int id)
{
    ivec4 nb = getB(id);
    for(int j = 0; j < 4; j++)
    {
        sort(nb[j]);
    }
}

vec4 renderPassB() {
	vec4 U = vec4(0.0);
	vec2 fragCoord = _xy;
  
    ivec2 p = ivec2(fragCoord);
    
    if(p.x > N.x || p.y > N.y) discard;
    
    int id = xy2i(p);
     
    u = ivec4(-1); d = vec4(1e10); 
   
    tid = id;
    pos = getparticleB(id).xy;
    
    sortneighbor(id); 
    
    for(int i = 0; i < 8; i++)
    {
        sort(hash(ivec4(p, FRAMECOUNT, i)).x%tot_n); //random sort    
    }
    
    ivec4 nb = getB(id);
    for(int i = 0; i < 4; i++)
    {
        sortneighbor(nb[i]); 
    }
    
    if( any(lessThan(u, ivec4(-1))) || any(greaterThan(u, ivec4(tot_n))))
    {
        u = ivec4(0);
    }
    
    
    U = save(u);
	return U; 
 } 


			//******** BuffC Code Begins ********

//4th order voronoi particle tracking 

ivec4 getC(ivec2 p)
{
    return floatBitsToInt(texel(ch2, p));
}
/*
vec4 save(ivec4 v)
{
    return intBitsToFloat(v);
}
*/
//ivec4 u; //ids
//vec4 d; //distances
//vec2 pos; //pixel position

vec4 getparticleC(int id)
{
    return texel(ch0, i2xy(id));
}

float particleDistance(int id, vec2 p)
{
    return sdSegment(p, getparticleC(id).xy, getparticleC((id+1)%tot_n).xy);
}
/*
//insertion sort
void sort(int utemp)
{
    if(utemp <0) return; 
   	float dtemp = particleDistance(utemp, pos);
    //sorting
    if(d.x > dtemp)
    {
        d = vec4(dtemp, d.xyz);
        u = ivec4(utemp, u.xyz);
    }
    else if(d.y > dtemp && dtemp > d.x)
    {
        d.yzw = vec3(dtemp, d.yz);
        u.yzw = ivec3(utemp, u.yz);
    }
    else if(d.z > dtemp && dtemp > d.y)
    {
        d.zw = vec2(dtemp, d.z);
        u.zw = ivec2(utemp, u.z);
    }
    else if(d.w > dtemp && dtemp > d.z)
    {
        d.w = dtemp;
        u.w = utemp;
    }
}
*/
void sortpos(ivec2 p)
{
    ivec4 nb = getC(p);
    for(int j = 0; j < 4; j++)
    {
        sort(nb[j]);
    }
}

vec4 renderPassC() {
	vec4 U = vec4(0.0);
	vec2 fragCoord = _xy;

    pos = fragCoord;
    ivec2 p = ivec2(pos);
     
    u = ivec4(-1); d = vec4(1e10); 
   
    //jump flood sorting 
    sortpos(p); //resort this position
    for(int i = 0; i < 8; i++)
    {
        sortpos(p+cross_distribution(i)); 
    }
    
    for(int i = 0; i < 4; i++)
    {
        sort(hash(ivec4(p, FRAMECOUNT, i)).x%tot_n); //random sort    
    }
    
    if( any(lessThan(u, ivec4(-1))) || any(greaterThan(u, ivec4(tot_n))) )
    {
        u = ivec4(0);
    }
    
    U = save(u);
	return U; 
 } 


// Fork of "Connected particle chain image" by michael0884. https://shadertoy.com/view/3dXfDN
// 2020-04-30 15:21:35

// Fork of "Large scale flocking" by michael0884. https://shadertoy.com/view/tsScRG
// 2020-04-30 07:24:31

ivec4 getM(ivec2 p)
{
    return floatBitsToInt(texel(ch2, p));
}

ivec4 getb(int id)
{
    return floatBitsToInt(texel(ch1, i2xy(id)));
}

vec4 getparticleM(int id)
{
    return texel(ch0, i2xy(id));
}

vec3 imageC(vec2 p)
{
    return texture(ch3, vec2(1., 1.)*p/size).xyz;
}

float particleDistanceM(int id, vec2 p)
{
    return sdSegment(p, getparticleM(id).xy, getparticleM((id+1)%tot_n).xy);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 pos = _xy;

    ivec4 nb = getM(ivec2(pos));
 	vec4 p0 = getparticleM(nb.x);
   
    
    float d0 = 0.;
    
    d0 += 0.24/(1.+0.3*particleDistance(nb[0], pos));
    d0 += 0.24/(1.+0.3*particleDistance(nb[1], pos));
    
    vec3 v = imageC(pos);
    
    fragColor.xyz = 2.5*v*sin(vec3(1,1,3)*d0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderPassC();
	}
	if(PASSINDEX == 3){
		return renderMainImage();
	}
}