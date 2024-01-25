

			//******** Common Code Begins ********
bool fc = FRAMECOUNT <= 1 || Reset != 0.0;
float growthFactor = pow((syn_BassLevel*0.35)+(syn_MidLevel*0.35)+syn_Intensity*0.3, 2.0)*(0.5+0.35*syn_Intensity);
// vec4 mixedEdgeCol = mix( _loadMedia(), mediaEdges, .1+0.5*edge_mix)*media_impact*0.25;
float mediahits_mix =  mix(media_impact, media_impact*(.125+0.875*syn_BassLevel), media_hits);

vec4 mediaEdges = texture(media_pass_fx, _uv);
vec4 mixedEdgeCol = mix(mediaEdges, _loadMedia(), .90-0.6*edge_mix)*(1.0+0.2*syn_Level);
bool mediaOn = _exists(syn_UserImage);
vec4 edgemixed = _brightness( abs(texture(media_pass_fx, _uv)*media_impact), .9 );
float media_lum = sin(_luminance(mixedEdgeCol*mediahits_mix))*media_impact;
float ddist(vec2 p0, vec2 pf) {
    return sqrt((pf.x - p0.x) * (pf.x - p0.x) + (pf.y - p0.y) * (pf.y - p0.y));
}
float d2 = ddist(RENDERSIZE.xy * 0.5, _xy.xy) * 0.004;

//simulation variables
#define dt 0.25*(0.75+0.25*growthFactor)
#define prad 1.4*(0.875+growthFactor*0.125) 
#define decay Decay

//cell speed
float pspeed = 6.+ media_lum*40.*12*(1.0+0.5*syn_BassLevel);

//sensor distance 
float sdist = 10.*(0.875+growthFactor*0.125)+ sin(media_lum)*2.*(1.0+0.2*syn_BassLevel);

//sensor strenght
float sst = 10.+ sin(media_lum)*PI*15.*4*(1.0+0.5*syn_BassLevel);

//sensor angle
float sangl = 0.3; //radians

#define pdens 2.
//definitions
#define size RENDERSIZE.xy
#define pixel(a, p) texture(a, p/vec2(textureSize(a,0)))
#define texel(a, p) texelFetch(a, ivec2(p-0.5), 0)
#define A BuffA
#define B BuffB
#define C BuffC
// #define ch3 iChannel3
#define PI 3.1415926535

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
vec2 rotateCenter(vec2 uvIn, float amount){
  uvIn.y += (RENDERSIZE.x-RENDERSIZE.y)/RENDERSIZE.x;
  uvIn*=vec2(1.0, RENDERSIZE.y/RENDERSIZE.x);
  _uv2uvc(uvIn);
  uvIn = _rotate(uvIn, amount);
  _uvc2uv(uvIn);
  uvIn/=vec2(1.0, RENDERSIZE.y/RENDERSIZE.x);
  uvIn.y -= (RENDERSIZE.x-RENDERSIZE.y)/RENDERSIZE.x;
  return uvIn;
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

vec4 LaplaceP(sampler2D ch, vec2 p)
{
    // vec3 dx = vec3(-5,0.,5);
    vec3 dx = vec3(-6,0.,6);
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
    // pos += _uvc * PI * Fisheye * (1.0 + 0.5 * low) * d2;

    vec2 muv = _mouse.xy/size;
    if ((length(media_lum)) > 0.){
        sst*=(2.25-SST*1.5);
        sdist*=2.5;
    }

    if(invert_paint == 1.0){
        if (_mouse.z> 0. && (length(pos-_mouse.xy))+length(_uvc) > 60.*paint_size){
            sst+=(0.25);
            sdist+= .5;
            sst*=(2.25-SST*1.5);
            sdist*=2.5;
            pspeed *= 1.5;
            sdist-=2.5*media_lum;
            pspeed -= 1.5*media_lum;

        }
    }
    if(invert_paint != 1.0){
    if (_mouse.z> 0. && (length(pos-_mouse.xy))+length(_uvc) < 60.*paint_size){
            sst+=(0.25);
            sdist+= .5;         
            sst*=(2.25-SST*1.5);
            sdist*=2.5;
            sdist-=2.5*media_lum;
            pspeed -= 1.5*media_lum;

        }

    }
    sdist *= dist;
  		// sst *= 0.05*(1.+ SST); 
    sst *= SST; 
    if(invert_paint == 1.0){
        if (_mouse.z> 0. && (length(pos-_mouse.xy))+length(_uvc) > 60.*paint_size){
            sst*=(1.6);
            sdist*=1.5;
        }
    }
        if ((length(media_lum)) > 1){
        sst*=(2.25-SST*1.5);
        sdist*=2.5;
    }

    if(invert_paint != 1.0){
    if (_mouse.z> 0. && (length(pos-_mouse.xy))+length(_uvc) < 60.*paint_size){
            sst*=(1.6);
            sdist*=1.5;
        }

    }
    //this pixel value
    U = texel(A, pos);
    
    //check neighbours 
    CheckRadius(U, pos, 1.);
    //CheckRadius(U, pos, 2.);
    CheckRadius(U, pos, 3.);
   // CheckRadius(U, pos, 4.);
   // CheckRadius(U, pos, 5.);
   
    U.xy = loop(U.xy);
    
    //cell cloning 
    if(length(U.xy - pos) > 10.)
    	U.xy += 1.*(hash22(pos)-0.5);

    //sensors
    vec2 sleft = U.xy + sdist*vec2(cos(U.z+sangl), sin(U.z+sangl));
    vec2 sright = U.xy + sdist*vec2(cos(U.z-sangl), sin(U.z-sangl));
    // sleft = mix(sleft, mixedEdgeCol.xy, Mix*0.75);

    float dangl = (pixel(B, sleft).x - pixel(B, sright).x);
    
    // if (_mouse.z> 0. && (length(pos-_mouse.xy))+length(_uvc) < 60.*paint_size){
    //     // sst*=(2.-SST*1.5);
    //     dangl*=2.5;
    // // pspeed*=.9;
    // // sangl += 0.5;
    // // Q += dt*0.5;
    // }

    U.z += dt*sst*tanh(3.*dangl);
   
    vec2 pvel = pspeed*vec2(cos(U.z), sin(U.z)) + 0.1*(hash22(U.xy+TIME)-0.5);
    
    pvel += mixedEdgeCol.xy;
    //update the particle
    U.xy += dt*pvel;
    
    U.xy = loop(U.xy);
    

    if(fc)
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
    Q = texel(B, p);
    //diffusion equation
    Q += dt*Laplace(B, p);
    p = mix(p, p+PI, Drift);
    p = mix(p, p+PI*_uvc*2, Stretch);
    p = mix(p, p-_uvc, zoom);
    // p -= _uvc*PI*zoom*(1.+0.5*growthFactor); //zoom
    
    vec4 particle = texel(A, p);

    float distr = gauss(p - particle.xy, prad);
    vec2 ss = vec2(textureSize(media_pass_fx,0))/size;


    float video = length(Laplace(media_pass_fx, p*ss));
    if (_mouse.z> 0. && length(p-_mouse.xy) < 40.*paint_size){
    //    distr*=.9;
    // Q += dt*0.5;
    }

    // Q += mix(vec4(0.), mixedEdgeCol, 0.125*int(mediaOn)*media_impact);
    
    //pheromone depositing
    Q += dt*(distr+2.*video*overlay);
    // Q += dt*(distr);

    // Q += mix(vec4(0.), mixedEdgeCol, 0.125*int(mediaOn)*media_impact);
        
    //pheromone decay
    Q += -dt*decay*Q;
    Q = abs(Q);
    if(fc) Q = vec4(0);

	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 p = _xy;

    Q = texel(C, p);

    Q = 0.9*Q + 0.1*texel(B, p); 
    if(fc) Q =vec4(0);
	
    return Q; 
    // Q = mix(Q, mixedEdgeCol, 0.25*int(mediaOn)*media_impact);

 } 


// Fork of "Physarum Polycephalum Simulation" by michael0884. https://shadertoy.com/view/tlKGDh
// 2020-04-25 20:06:12

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 pos = _xy;

	vec4 particle = texel(A, pos);
    float distr = gauss(pos - particle.xy, prad);
    vec4 pheromone = 2.5*texel(C, pos);
    fragColor = vec4(sin(pheromone.xyz*vec3(1,1.2,1.5)), 1.);
    // fragColor += _edgeDetectSobel(BuffC, _uv)*0.1;
    // fragColor+= mix(vec4(0.), mixedEdgeCol, int(mediaOn)*media_impact*2.)*media_impact;

	return fragColor; 
 } 
vec4 mediaPass() {
    vec4 media = vec4(0.);
    vec2 U = _xy;
    media = _brightness(_loadMedia(), 1.0)*media_impact*2.;
    // media = _loadMedia()*media_impact*2.;
    return media;
}
vec4 mediaPassFX() {
    vec4 media = vec4(0.);
    vec2 U = _xy;
    vec4 media_edge_mix = mix(_edgeDetectSobel(media_pass), _loadMedia()*media_impact, 1.0-edge_mix*0.5-0.25);
    return media_edge_mix;
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
        // edgemixed = mix(edgemixed, edgemixed*(0.125+syn_Level*0.875), 1.0-media_hits);

		// return mix(Q, mixedEdgeCol, 0.35*int(mediaOn)*media_impact) ;
		return mediaPass();

	}
	if(PASSINDEX == 4){
        // edgemixed = mix(edgemixed, edgemixed*(0.125+syn_Level*0.875), 1.0-media_hits);

		// return mix(Q, mixedEdgeCol, 0.35*int(mediaOn)*media_impact) ;
		return mediaPassFX();

	}
	if(PASSINDEX == 5){
		return renderMainImage();
	}
}