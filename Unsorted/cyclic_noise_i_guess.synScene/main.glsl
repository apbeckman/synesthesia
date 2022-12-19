

			//******** Common Code Begins ********

//simulation variables
#define dt 0.25
#define prad 1.4
#define decay 0.04

//particle speed
float pspeed = 6.;


#define pdens 2.
//definitions
#define size RENDERSIZE.xy
#define pixel(a, p) texture(a, p/vec2(textureSize(a,0)))
#define texel(a, p) texelFetch(a, ivec2(p-0.5), 0)
#define ch0 iChannel0
#define ch1 iChannel1
#define ch2 iChannel2
#define ch3 iChannel3
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

#define pi 3.14159265
#define iterations 13

vec3 function(vec3 x){
    vec3 ca = vec3(0.420,0.827,1.000)*cos(x.xxx)+vec3(1.000,0.902,0.859)*cos(x.yyy);
    vec3 cb = vec3(0.898,1.000,0.922)*cos(x.zzz*1.5);
    return ca*cb-cb;
}

mat2 ROT(float ang)
{
    return mat2(cos(ang), sin(ang), -sin(ang), cos(ang));
}

#define SCALE 1.5
#define SHIFT vec3(2.3, -5.2, 1.0)

vec3 fractal(vec3 x){
    x *= pi;
    vec3 v = vec3(0.0);
    float a = 0.5;
    mat2 rmZ = SCALE*ROT(0.35);
    mat2 rmX = SCALE*ROT(0.46);
    for (int i = 0; i < iterations; i++){
        vec3 F =function(x); 
        v += a*F;
        x.xy = rmZ*x.xy;
        x += 0.3*F;
        x.yz = rmX*x.yz;
        x += SHIFT;
        a /= 1.1*SCALE;
    }
    return v;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = 1.5*(fragCoord-0.5*RENDERSIZE.xy)/max(RENDERSIZE.x, RENDERSIZE.y);
    vec3 color = 0.8*fractal(vec3(uv, 0.05*TIME))+0.4;
    fragColor = vec4(tanh(pow(color, vec3(1.4))), 1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}