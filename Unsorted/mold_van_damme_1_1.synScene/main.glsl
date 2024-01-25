

			//******** Common Code Begins ********
bool fc = FRAMECOUNT <= 1 || Reset != 0.0;
float growthFactor = pow((syn_BassLevel*0.35)+(syn_MidLevel*0.35)+syn_Intensity*0.3, 2.0)*(0.5+0.35*syn_Intensity);
// vec4 mixedEdgeCol = mix( _loadMedia(), mediaEdges, .1+0.5*edge_mix)*media_impact*0.25;
float mediahits_mix =  mix(media_impact, media_impact*(.125+0.875*syn_BassLevel), media_hits);

vec4 mediaEdges = texture(media_pass_fx, _uv);
bool mediaOn = _exists(syn_UserImage);
float media_lum = sin(PI*0.75*length(mediaEdges))*mediahits_mix;
// float paintsize = _smin(70.*paint_size, 60.*paint_size, 1.0);
float paintsize = 60.*paint_size;
//simulation variables
#define dt 0.25
#define prad 1.4
#define decay 0.3

//cell speed
float pspeed = 6.;

//sensor distance 
float sdist = 10.;

//sensor strenght
float sst = 10.;

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

    if(mediaOn){(sst *=1.0+2.5*media_lum)*(1.0+syn_MidLevel);}
    if(invert_paint == 1.0){
        if (_mouse.z> 0. && (length(pos-_mouse.xy)) > paintsize){
            sst*=(1.6);
            sdist*=1.5;
            sst-=.25*media_lum;
            sdist-=.25*media_lum;
            pspeed += .25*media_lum;

        }
    }
    if(invert_paint != 1.0){
        if (_mouse.z> 0. && (length(pos-_mouse.xy)) < paintsize){
                sst*=(1.6);
                sdist*=1.5;
            sst-=.25*media_lum;
            sdist-=.25*media_lum;
            pspeed += .25*media_lum;

            }

    }

    sdist *= dist;

    sst *= SST; 
// sdist *= (1.0+0.125*media_lum*(1.0+syn_Level));
// sst *= (1.0+.125*media_lum*(1.0+syn_Level));
    if(invert_paint == 1.0){
        if (_mouse.z> 0. && (length(pos-_mouse.xy)) > paintsize){
            sst*=(1.125-SST*1.05);
            sdist*=1.5;
            sst-=.25*media_lum;
            sdist-=.25*media_lum;
            pspeed += .25*media_lum;
        }
    }
    if(invert_paint != 1.0){
    if (_mouse.z> 0. && (length(pos-_mouse.xy)) < paintsize){
            sst*=(1.25-SST*1.125);
            sdist*=1.5;
            sst-=.25*media_lum;
            sdist-=.25*media_lum;
            pspeed += .25*media_lum;
        }

    }
    pspeed *= 1.0-0.2* growthFactor;
    // if(length(muv.xy) >0.)
    // {
    // 	sdist *= muv.x;
  	// 	sst *= muv.y; 
    // }
    // else
    // {
    //     sdist *= 0.8;
  	// 	sst *= 0.05; 
    // }
   
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
    // sangl*=(1.0-0.025*media_lum);

    //sensors
    vec2 sleft = U.xy + sdist*vec2(cos(U.z+sangl), sin(U.z+sangl));
    vec2 sright = U.xy + sdist*vec2(cos(U.z-sangl), sin(U.z-sangl));
    
    float dangl = (pixel(B, sleft).x - pixel(B, sright).x);
    // dangl*=(1.0+0.025*media_lum);
    U.z += dt*sst*tanh(3.*dangl);
   
    vec2 pvel = pspeed*vec2(cos(U.z), sin(U.z)) + 0.1*(hash22(U.xy+TIME)-0.5);;
    pvel += mediaEdges.xy*0.12;
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
    float dtt = dt;
    if(invert_paint == 1.0){
        if (_mouse.z> 0. && (length(p-_mouse.xy)) > paintsize){
            dtt = dt* .9;

        }
    }
    if(invert_paint != 1.0){
        if (_mouse.z> 0. && (length(p-_mouse.xy)) < paintsize){
            dtt = dt* .9;
            
            }

    }
    dtt -= abs(media_lum) *0.001;
    Q = texel(B, p);
   
    //diffusion equation
    Q += dtt*Laplace(B, p);
    
    vec4 particle = texel(A, p);
    float distr = gauss(p - particle.xy, prad);
    
    vec2 ss = vec2(textureSize(media_pass_fx,0))/size;
    float video = length(Laplace(media_pass_fx, p*ss)*0.125);
    
    //pheromone depositing
    Q += dtt*(distr+video*show_media);
        
    //pheromone decay
    Q += -dtt*decay*Q;
    
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
 } 
vec4 mediaPass() {
    vec4 media = vec4(0.);
    vec2 U = _xy;
    media = _loadMedia();
    // media = _loadMedia()*media_impact*2.;
    return media;
}
vec4 mediaPassFX() {
    vec4 media = vec4(0.);
    vec2 U = _xy;
    vec4 media_edge_mix = mix(_edgeDetectSobel(media_pass), texture(media_pass, _uv), 1.0-edge_mix*0.975-0.0125)*media_impact;
    return media_edge_mix;
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
		return mediaPass();
	}
	if(PASSINDEX == 4){
		return mediaPassFX();
	}
	if(PASSINDEX == 5){
		return renderMainImage();
	}
}