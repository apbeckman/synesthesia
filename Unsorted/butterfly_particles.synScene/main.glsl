

			//******** Common Code Begins ********
float growthFactor = (pow((syn_BassLevel*0.35)+(syn_MidLevel*0.35)+(syn_Intensity*0.3), 2.0));

//simulation variables
#define dt 0.25
#define prad .24+growthFactor //1.4
#define decay 0.02

//particle speed
float pspeed = 10.;


#define pdens 2.
//definitions
#define size RENDERSIZE.xy
#define pixel(a, p) texture(a, p/vec2(textureSize(a,0)))
#define texel(a, p) texelFetch(a, ivec2(p-0.5), 0)
#define A BuffA
#define B BuffB
#define ch2 iChannel2
#define ch3 iChannel3
//#define PI 3.14159265

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
    return vec3(base(v - 0.5),base(v),base(v + 10.25));
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


float angle_between(vec2 a,vec2 b)
{
    return atan(a.x*b.y-a.y*b.x,dot(a,b));
}

			//******** BuffA Code Begins ********

//voronoi particle tracking 

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
    vec4 Unb = texel(A, loop(pos+dx));
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
   
    //this pixel value
    U = texel(A, pos);
    
    //check neighbours 
    CheckRadius(U, pos, 1.);
    CheckRadius(U, pos, 2.);
    CheckRadius(U, pos, 3.);
    CheckRadius(U, pos, 4.);
    CheckRadius(U, pos, 5.);
   
    U.xy = loop(U.xy);
    
    //cell cloning 
   // if(length(U.xy - pos) > 7.)
    //	U.xy += 1.*(hash22(pos)-0.5);

    //syncronizing the particles with the flow
    
    vec2 vel0 = pspeed*vec2(cos(U.z), sin(U.z)) + 0.1*(hash22(U.xy+smoothTime)-0.5);
    vec4 F =  pixel(B, U.xy);
    vec2 vel1 = F.yz;
    float dangl = (0.8+2.*growthFactor)*angle_between(vel0,vel1);
    U.z += dt*(dangl+0.1*(hash12(U.xy+smoothTime)-0.5)+1.*sin(U.w));
    U.w += dt*0.4;
   
    vec2 pvel = pspeed*vec2(cos(U.z), sin(U.z)) + 1.*(hash22(U.xy+smoothTime)-0.5);;
    
    //update the particle
    U.xy += dt*pvel;
    
    U.xy = loop(U.xy);
    
    
    if(FRAMECOUNT <= 1 || (_mouse.z > 0. && length(_mouse.xy - pos) < 30.) || (length(size.xy*0.5 - pos) < 5.))
    {
        U.xy = vec2(pdens*round(pos.x/pdens),pdens*round(pos.y/pdens));
        U.zw = 2.*PI*(hash22(U.xy) - 0.5);
    }
    
    if(FRAMECOUNT <= 1)
    {
        U = vec4(0);
    }
	return U; 
 } 


			//******** BuffB Code Begins ********

//saving and diffusing the velocity trails 

vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 p = _xy;

    Q = texel(B, p);
   
    //diffusion equation
    Q += 0.05*dt*Laplace(B, p);
    
    vec4 particle = texel(A, p);
    float distr = gauss(p - particle.xy, prad);
    
    vec2 pvel = pspeed*vec2(cos(particle.z), sin(particle.z));
    
    //pheromone depositing
    Q += dt*vec4(1., pvel.x, pvel.y, 1.)*distr;
        
    //pheromone decay
    Q += -dt*decay*Q;
    
    if(FRAMECOUNT <= 1) Q = vec4(0);
	return Q; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 pos = _xy;

	vec4 particle = texel(A, pos);
    float distr = gauss(pos - particle.xy, prad);
    vec4 flow = 0.5*texel(B, pos);
    fragColor = vec4(sin(flow.xxx*vec3(1.5,1.2,1.)), 1.);
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
		return renderMainImage();
	}
}