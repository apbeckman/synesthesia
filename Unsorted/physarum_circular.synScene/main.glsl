

			//******** Common Code Begins ********

//simulation variables
float growthFactor = pow((syn_BassLevel*0.5)+(syn_MidLevel*0.35)+(syn_Level*0.15), 2.0);

#define dt 0.25
#define prad 1.4+growthFactor*0.125
#define decay 0.125

//cell speed
float pspeed = 7.;

//sensor distance 
float sdist = 50.*(0.99+growthFactor*0.0125);

//sensor strenght
float sst = 10.*(0.99+growthFactor*0.0125);

//sensor angle
float sangl = 0.2*(1-sin(smoothTime*0.5)*0.5+0.5); //radians

float pdens =  4.*(0.99+growthFactor*0.0125);
//definitions
#define size RENDERSIZE.xy
#define pixel(a, p) texture(a, p/vec2(textureSize(a,0)))
#define texel(a, p) texelFetch(a, ivec2(p-0.5), 0)
#define ch0 BuffA
#define ch1 BuffB
#define ch2 BuffC
#define ch3 BuffD
//#define PI 3.14159265453

//hash functions
//https://www.shadertoy.com/view/4djSRW
float hash11(float p)
{
    p = fract(p * .1031);
    p *= p + 33.33;
    p *= p + p;
    return fract(p);
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


//functions
float gauss(vec2 x, float r)
{
    return exp(-pow(length(x)/r,2.));
}
   

//a rainbow colormap from Matlab
float interpolate(float val, float y0, float x0, float y1, float x1) 
{
    return (val-x0)*(y1-y0)/(x1-x0) + y0;
}

float base(float val) 
{
    if ( val <= -0.75 ) return 0.0;
    else if ( val <= -0.25 ) return interpolate( val, 0.0, -0.75, 1.0, -0.25 );
    else if ( val <= 0.25 ) return 1.0;
    else if ( val <= 0.75 ) return interpolate( val, 1.0, 0.25, 0.0, 0.75 );
    else return 0.0;
}

vec3 jet_colormap(float v)
{
    return vec3(base(v - 0.5),base(v),base(v + 0.5));
}

vec3 jet_range(float v, float a, float b)
{
    return jet_colormap(2.*clamp((v-a)/(b-a),0.,1.) - 1.);
}

//Laplacian operator
vec4 Laplace(sampler2D ch, vec2 p)
{
    vec3 dx = vec3(-1,0.,1);
    return texel(ch, p+dx.xy)+texel(ch, p+dx.yx)+texel(ch, p+dx.zy)+texel(ch, p+dx.yz)-4.*texel(ch, p);
}

//Gradient
vec2 Grad(sampler2D ch, vec2 p)
{
    vec3 dx = vec3(-1,0.,1);
    return vec2(length(texel(ch, p+dx.zy)),length(texel(ch, p+dx.yz))) - length(texel(ch, p));
}


			//******** BuffA Code Begins ********

//voronoi particle tracking 
//simulating the cells

//loop the vector
vec2 loop_d(vec2 pos)
{
	return mod(pos + size*0.5, size) - size*0.5;
}

//loop the space
vec2 loop(vec2 pos)
{
	return mod(pos, size);
}


void Check(inout vec4 U, vec2 pos, vec2 dx)
{
    vec4 Unb = texel(BuffA, loop(pos+dx));
    //check if the stored neighbouring particle is closer to this position 
    if(length(loop_d(Unb.xy - pos)) < length(loop_d(U.xy - pos)))
    {
        U = Unb; //copy the particle info
    }
}

void CheckRadius(inout vec4 U, vec2 pos, float r)
{
    Check(U, pos, vec2(-r,0));
    Check(U, pos, vec2(r,0));
    Check(U, pos, vec2(0,-r));
    Check(U, pos, vec2(0,r));
}

vec4 renderPassA() {
	vec4 U = vec4(0.0);
	vec2 pos = _xy;
    vec2 muv = _mouse.xy/size;
    
    if(length(muv.xy) >0.)
    {
    	sdist *= muv.x;
  		sst *= muv.y; 
    }
    else
    {
        sdist *= 0.8;
  		sst *= 0.05; 
    }
   
    //this pixel value
    U = texel(BuffA, pos);
    
    //check neighbours 
    CheckRadius(U, pos, 1.);
    CheckRadius(U, pos, 2.);
    CheckRadius(U, pos, 3.);
    CheckRadius(U, pos, 4.);
    CheckRadius(U, pos, 5.);
   
    U.xy = loop(U.xy);
    
    //cell cloning 
    if(length(U.xy - pos) > 10.)//change pos for thiccness
    	U.xy += 1.*(hash22(pos)-0.5);

    //sensors
    vec2 sleft = U.xy + sdist*vec2(cos(U.z+sangl), sin(U.z+sangl));
    vec2 sright = U.xy + sdist*vec2(cos(U.z-sangl), sin(U.z-sangl));
    
    float dangl = (pixel(BuffB, sleft).x - pixel(BuffB, sright).x);
    U.z += dt*sst*tanh(3.*dangl)*(1+growthFactor*0.25);
   
    vec2 pvel = (0.9+growthFactor*0.125)*pspeed*vec2(cos(U.z), sin(U.z));
    
    //update the particle
    U.xy += dt*pvel;
    
    U.xy = loop(U.xy*(1+growthFactor*0.00025));
    
    if(length(size*0.5 - U.xy) > 0.5*min(size.x, size.y))
    {
        U.xy = normalize(U.xy- 0.5*size)* 0.5*min(size.x, size.y)+0.5*size;
        U.z += PI;
    }
    
    
    if(FRAMECOUNT <= 1)
    {
        U.xy = vec2(pdens*round(pos.x/pdens),pdens*round(pos.y/pdens));
        U.zw = hash22(U.xy) - 0.5;
    }
	return U; 
 } 


			//******** BuffB Code Begins ********

//depositing and diffusing the pheromone trails 

vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 p = _xy;

    Q = texel(BuffB, p);
   
    //diffusion equation
    Q += dt*Laplace(BuffB, p)*(0.99+growthFactor*0.0125);
    
    vec4 particle = texel(BuffA, p);
    float distr = gauss(p - particle.xy, prad)*(0.9+growthFactor*0.25);
    
    //pheromone depositing
    Q += (dt*distr/1.);
        
    //pheromone decay
    Q += -dt*decay*Q*(1-growthFactor*0.0125);
    
    if(FRAMECOUNT <= 1) Q = vec4(0);
	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 p = _xy;

    Q = texel(BuffB, p);
    
    Q = 0.85*Q + 0.15*texel(BuffA, p); 
    if(FRAMECOUNT <= 1) Q =vec4(0);
	return Q; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 pos = _xy;

	vec4 particle = texel(BuffA, pos);
    float distr = gauss(pos - particle.xy, prad*(1.0+basshits*0.125));
    vec4 pheromone = 2.*texel(BuffB, pos);
    fragColor = vec4(sin(pheromone.xyz*vec3(1,1.2,1.5)), 1.);
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