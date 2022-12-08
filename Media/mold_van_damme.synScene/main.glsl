

			//******** Common Code Begins ********

//simulation variables
float growthFactor = pow((syn_BassLevel*0.35)+(syn_MidLevel*0.35)+(syn_Level*0.2), 2.0)*syn_Intensity;

float dt = 0.25*(0.5+growthFactor/2);
#define prad 1.4 +growthFactor
#define decay 0.3

//cell speed
float pspeed = 6.;

//sensor distance 
float sdist = 10.*(1.0+growthFactor);

//sensor strenght
float sst = 10.*(0.75+syn_Level);

//sensor angle
float sangl = 0.3; //radians

#define pdens 1.
//definitions
#define size RENDERSIZE.xy
#define pixel(a, p) texture(a, p/vec2(textureSize(a,0)))
#define texel(a, p) texelFetch(a, ivec2(p-0.5), 0)
#define ch0 BuffA
#define ch1 BuffB
#define ch2 BuffC
//#define ch3 iChannel3
#define PI 3.14159265

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
    return jet_colormap(12.*clamp((v-a)/(b-a),0.,1.) - 1.);
}

//Laplacian operator
vec4 Laplace(sampler2D ch, vec2 p)
{
    vec3 dx = vec3(-1,0.,1);
    return texel(ch, p+dx.xy)+texel(ch, p+dx.yx)+texel(ch, p+dx.zy)+texel(ch, p+dx.yz)-4.*texel(ch, p);
}

vec4 LaplaceP(sampler2D ch, vec2 p)
{
    vec3 dx = vec3(-5,0.,5);
    return pixel(ch, p+dx.xy)+pixel(ch, p+dx.yx)+pixel(ch, p+dx.zy)+pixel(ch, p+dx.yz)-4.*pixel(ch, p);
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
    vec4 Unb = texel(ch0, loop(pos+dx));
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
    vec2 posCent = _uvc-brushPosition*0.5*vec2(RENDERSIZE.x/RENDERSIZE.y,1.0);
    float mask = smoothstep(0.5,0.5+0.5,length(posCent));
    distFunc = 1.0-mask;



    /*
    if(length(muv.xy) >0.)
    {
    	sdist *= muv.x;
  		sst *= muv.y; 
    }
    else
    {
    }
   */
    //this pixel value


    
    sdist *= 0.8;
    sst *= 0.05; 

    U = texel(ch0, pos);
    //check neighbours 
    CheckRadius(U, pos, 1.);
    //CheckRadius(U, pos, 2.);
    CheckRadius(U, pos, 3.);
   // CheckRadius(U, pos, 4.);
   // CheckRadius(U, pos, 5.);
   
    U.xy = loop(U.xy);
    
    //cell cloning 
    if(length(U.xy - pos) > 15.)
    	U.xy += 1.*(hash22(pos)-0.25);

    //sensors
    vec2 sleft = U.xy + sdist*vec2(cos(U.z+sangl), sin(U.z+sangl));
    vec2 sright = U.xy + sdist*vec2(cos(U.z-sangl), sin(U.z-sangl));
    
    float dangl = (pixel(ch1, sleft).x - pixel(ch1, sright).x);
    U.z += dt*sst*tanh(3.*dangl);
   
    vec2 pvel = pspeed*vec2(cos(U.z), sin(U.z))*(0.8+basshits/2.) + 0.1*(hash22(U.xy+smoothTime*0.3)-0.5);;
    
    //update the particle
    U.xy += dt*pvel;
    
    U.xy = loop(U.xy);
    
    
    if(FRAMECOUNT <= 1)
    {
        U.xy = vec2(pdens*round(pos.x/pdens),pdens*round(pos.y/pdens));
        U.zw = hash22(U.xy) - (0.5);
    }
	return U; 
 } 


			//******** BuffB Code Begins ********

//depositing and diffusing the pheromone trails 

vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 p = _xy;

    Q = texel(ch1, p);
   
    //diffusion equation
    Q += dt*Laplace(ch1, p);
    
    vec4 particle = texel(ch0, p);
    float distr = gauss(p - particle.xy, prad);
    //distr *= (1.0+basshits);
    vec2 ss = vec2(textureSize(ch2,0))/size;
    float video = length(Laplace(syn_UserImage, p*ss)*(1+growthFactor));
    
    //pheromone depositing
       

    Q += dt*(distr+2. * video * 1.0+(highhits*.25*flash)*syn_Intensity);
    //Q += dt*distr/(1.);
    //Q+= basshits/1.;
    //pheromone decay
    Q += -dt*decay*Q;
    
    if(FRAMECOUNT <= 1) Q = vec4(0);
	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 p = _xy;

    Q = texel(ch1, p);
    
    Q = 0.9*Q + 0.1*texel(ch0, p); 
    if(FRAMECOUNT <= 1) Q =vec4(0);
	return Q; 
 } 


// Fork of "Physarum Polycephalum Simulation" by michael0884. https://shadertoy.com/view/tlKGDh
// 2020-04-25 20:06:12

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 pos = _xy;

	vec4 particle = texel(ch0, pos);
    float distr = gauss(pos - particle.xy, prad);
    vec4 pheromone = 1.35*texel(ch1, pos);
    fragColor = vec4(sin(pheromone.xyz*vec3(1,1.2,1.5)), 1.1);
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