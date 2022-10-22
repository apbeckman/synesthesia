

			//******** Common Code Begins ********


#define pi acos(-1.)


#define sint(a) (asin(sin(a))*2. - 1.)

#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))

#define pmod(p,d) mod(p - (d)*0.5, (d)) - 0.5*(d)

float r11(float i){ return fract(sin(i*12.126)*12.6);}

//#define xor(a,b,c) min(max((a),-(b)), max((b),-(a) - c)) 

float ss( float c, float power, float bias){
    c = clamp(c,-0.,1.);
    //c = smoothstep(0.,1.,c);
    
    c = pow(c,1. + bias);
    
    float a = pow( abs(c), power);
    float b = 1.-pow( abs(c - 1.), power);
    
    return mix(a,b,c);
}
float valueNoise(float i, float p){ return mix(r11(floor(i)),r11(floor(i) + 1.), ss(fract(i), p,0.6));}

float valueNoiseStepped(float i, float p, float steps){ return mix(  floor(r11(floor(i))*steps)/steps, floor(r11(floor(i) + 1.)*steps)/steps, ss(fract(i), p,0.6));}


// See: https://www.shadertoy.com/view/ls2Bz1
// Spectral Colour Schemes
// By Alan Zucconi
// Website: www.alanzucconi.com
// Twitter: @AlanZucconi

// Example of different spectral colour schemes
// to convert visible wavelengths of light (400-700 nm) to RGB colours.

// The function "spectral_zucconi6" provides the best approximation
// without including any branching.
// Its faster version, "spectral_zucconi", is advised for mobile applications.


// Read "Improving the Rainbow" for more information
// http://www.alanzucconi.com/?p=6703



float saturate (float x)
{
    return min(1.0, max(0.0,x));
}
vec3 saturate (vec3 x)
{
    return min(vec3(1.,1.,1.), max(vec3(0.,0.,0.),x));
}

// --- Spectral Zucconi --------------------------------------------
// By Alan Zucconi
// Based on GPU Gems: https://developer.nvidia.com/sites/all/modules/custom/gpugems/books/GPUGems/gpugems_ch08.html
// But with values optimised to match as close as possible the visible spectrum
// Fits this: https://commons.wikimedia.org/wiki/File:Linear_visible_spectrum.svg
// With weighter MSE (RGB weights: 0.3, 0.59, 0.11)
vec3 bump3y (vec3 x, vec3 yoffset)
{
	vec3 y = vec3(1.,1.,1.) - x * x;
	y = saturate(y-yoffset);
	return y;
}

// --- Spectral Zucconi 6 --------------------------------------------

// Based on GPU Gems
// Optimised by Alan Zucconi
vec3 spectral_zucconi6 (float x)
{
	// w: [400, 700]
	// x: [0,   1]

	const vec3 c1 = vec3(3.54585104, 2.93225262, 2.41593945);
	const vec3 x1 = vec3(0.69549072, 0.49228336, 0.27699880);
	const vec3 y1 = vec3(0.02312639, 0.15225084, 0.52607955);

	const vec3 c2 = vec3(3.90307140, 3.21182957, 3.96587128);
	const vec3 x2 = vec3(0.11748627, 0.86755042, 0.66077860);
	const vec3 y2 = vec3(0.84897130, 0.88445281, 0.73949448);

	return
		bump3y(c1 * (x - x1), y1) +
		bump3y(c2 * (x - x2), y2) ;
}


			//******** BuffA Code Begins ********


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    vec2 p = fragCoord.xy;
    
    vec2 uv = (p.xy) / RENDERSIZE.xy;
    
    //line thickness
    float th = 0.002;
    
    //gradient detection for r,g and b separatly
    //not sure why I did this actually
    
    float gxr = texture(BuffB, uv.xy+vec2(th, 0.)).r
    - texture(BuffB, uv.xy-vec2(th, 0.)).r ;
    
    float gyr = texture(BuffB, uv.xy+vec2(0., th)).r
    - texture(BuffB, uv.xy-vec2(0., th)).r;
    
    float gxg = texture(BuffB, uv.xy+vec2(th, 0.)).g
    - texture(BuffB, uv.xy-vec2(th, 0.)).g ;
    
    float gyg = texture(BuffB, uv.xy+vec2(0., th)).g
    - texture(BuffB, uv.xy-vec2(0., th)).g;
    
    float gxb = texture(BuffB, uv.xy+vec2(th, 0.)).b
    - texture(BuffB, uv.xy-vec2(th, 0.)).b;
    
    float gyb = texture(BuffB, uv.xy+vec2(0., th)).b
    - texture(BuffB, uv.xy-vec2(0., th)).b;
    
    //hack to concea noise from: https://www.shadertoy.com/view/Mdf3zr
    float gr = gxr*gxr + gyr*gyr;
    float gg = gxg*gxg + gyg*gyg;
    float gb = gxb*gxb + gyb*gyb;
    
    //more noise control
    float g = pow((gr+gg+gb)/1.,1.9);
    
    
    vec3 col = texture(BuffB, p / RENDERSIZE.xy).xyz;
    col = mix(col,  0.5 + 0.5*cos(smoothTimeC*2.+uv.y*6.+vec3(0,2,4)),g*100.);
    
    fragColor = vec4(col,1.);
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

                                                                                                                                                                                                                                                                                        // See Image tab for details, also visit:
//
// https://xemantic.github.io/shader-web-background/
//
// In the original shader-web-background these values are provided as uniforms
// feel free to play with them and if you will find something prettier than
// the equilibrium I established, please send it back to me :)
const vec2  iFeedbackZoomCenter       = vec2(0., 0.);
const float iFeedbackZoomRate         = .001;
const float iFeedbackFadeRate         = .993;
const float iFeedbackColorShiftZoom   = 0.1;
const float iFeedbackColorShiftImpact = -0.00233;
const vec2  iDrawCenter               = vec2(0., 0.);
const float iDrawIntensity            = .975;
const float iBlobEdgeSmoothing        = .09;
const float iBlob1Radius              = .33;
const float iBlob1PowFactor           = 20.;
const float iBlob1ColorPulseSpeed     = .032;
const float iBlob2Radius              = .66;
const float iBlob2PowFactor           = 30.;
const float iBlob2ColorPulseSpeed     = .01234;
const float iBlob2ColorPulseShift     = 2.5;
const float iColorShiftOfRadius       = 2.;
const float iFeedbackMouseShiftFactor = .003;

/*
  Normally it would be provided by texture parameters, but on Mac/iOS the texture REPEAT
  seems to work only for power-of-2 texture sizes.
 */
vec4 repeatedTexture(in sampler2D channel, in vec2 uv) {
    return texture(channel, mod(uv, 1.));
}

float drawBlob(
    in vec2 st,
    in vec2 center,
    in float radius,
    in float edgeSmoothing
) {
    float dist = length((st - center) / radius);
    return dist * smoothstep(1., 1. - iBlobEdgeSmoothing, dist);
}


vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // in shader-web-background provided as uniforms: start
    float iMinDimension = min(RENDERSIZE.x, RENDERSIZE.y);
    vec2 iScreenRatioHalf =
        (RENDERSIZE.x >= RENDERSIZE.y)
            ? vec2(RENDERSIZE.y / RENDERSIZE.x * .5, .5)
            : vec2(.5, RENDERSIZE.x / RENDERSIZE.y);
    vec3 iBlob1Color = spectral_zucconi6(
        mod(smoothTimeB*0.1 * iBlob1ColorPulseSpeed, 1.)
    );
    
    vec3 iBlob2Color = spectral_zucconi6(
        mod(smoothTimeB*0.1 * iBlob2ColorPulseSpeed + iBlob2ColorPulseShift, 1.)
    );
    vec2 iFeedbackShiftVector =
        (_mouse.x > 0. && _mouse.y > 0.)
            ? (_mouse.xy * 2. - RENDERSIZE.xy) / iMinDimension * iFeedbackMouseShiftFactor
            : vec2(0);
    // in shader-web-background provided as uniforms: end
            
    
    vec2 uv = fragCoord / RENDERSIZE.xy;
    vec2 st = (fragCoord * 2. - RENDERSIZE.xy) / iMinDimension;

    vec2  drawDelta = st - iDrawCenter;
    float drawAngle = atan(drawDelta.x, drawDelta.y);
    float drawDist = length(drawDelta);

vec3 feedbk = repeatedTexture(BuffD, uv - st).rgb;
    vec3 colorShift = repeatedTexture(
        BuffC,
        uv - st * iFeedbackColorShiftZoom * iScreenRatioHalf
    ).rgb;

    vec2 stShift = vec2(0);
    stShift += iFeedbackZoomRate * (st - iFeedbackZoomCenter);
    stShift += (feedbk.bg/colorShift.gr - .5) * iFeedbackColorShiftImpact;
    stShift += iFeedbackShiftVector;
    stShift *= iScreenRatioHalf;

    vec3 prevColor = repeatedTexture(BuffB, uv - stShift).rgb;
    prevColor *= iFeedbackFadeRate;

    vec3 drawColor = vec3(0);
   

    float radius =
        1.
        + (colorShift.r + colorShift.g + colorShift.b) * iColorShiftOfRadius;
    drawColor +=
        pow(
          drawBlob(st, iDrawCenter, radius * iBlob1Radius, iBlobEdgeSmoothing),
          iBlob1PowFactor
        ) * iBlob1Color;
    drawColor +=
        pow(
          drawBlob(st, iDrawCenter, radius * iBlob2Radius, iBlobEdgeSmoothing),
          iBlob2PowFactor
        ) * iBlob2Color;

    vec3 color = vec3(0);
    drawColor *= iDrawIntensity;
    prevColor *= iFeedbackFadeRate;
    color += prevColor;
    color += drawColor;

    color = clamp(color, 0., 1.);
    fragColor = vec4(color, 1.);
	return fragColor; 
 } 



			//******** BuffC Code Begins ********

// hueshifter but better https://www.youtube.com/watch?v=CXXmZiVXCxY

#define SWAY 0.01
#define BLURSTRENGTH 1.0

const int samples = 70,
          LOD = 2,         // gaussian done on MIPmap at scale LOD
          sLOD = 1 << LOD; // tile size = 2^LOD
const float sigma = float(samples) * .2;

float gaussian(vec2 i) {
    return exp( -.5* dot(i/=sigma,i) ) / ( 6.28 * sigma*sigma );
}

vec4 blur(sampler2D sp, vec2 U, vec2 scale) {
    vec4 O = vec4(0);  
    int s = samples/sLOD;
    
    for ( int i = 0; i < s*s; i++ ) {
        vec2 d = vec2(i%s, i/s)*float(sLOD) - float(samples)/2.;
        O += gaussian(d) * textureLod( sp, U + scale * d , float(LOD) );
    }
    
    return O / O.a;
}


vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/RENDERSIZE.xy;
    vec2 swayuv = uv * 1.0 + sin(smoothTimeB * 0.5) * SWAY;
    vec3 col = texture(BuffA,swayuv).rgb;
    vec3 blur = blur(BuffA,uv, 1./RENDERSIZE.xy ).rgb;
    blur *= BLURSTRENGTH;
    // col += blur * 0.5;
    bool chroma_shift = false;
    if (chroma_shift) {
        swayuv.x += 0.01;
        col.r *= texture(BuffA,swayuv).r;
        swayuv.x -= 0.02;
        col.b *= texture(BuffA,swayuv).b;
    }
    col = max(col,blur);
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


			//******** BuffD Code Begins ********

// Fork of "spirit aura" by xenn. https://shadertoy.com/view/7lSGzR
// 2021-06-21 18:20:40

/*
	Full Scene Radial Blur
	----------------------

	Radial blur - as a postprocessing effect - is one of the first things I considered doing 
	when the multipass system came out. I've always loved this effect. Reminds me of the early 
	demos from Aardbei et al. 

	Anyway, Shadertoy user, Passion, did a really cool radial blur on a field of spheres that
	inspired me to do my own. Radial blurs are pretty straight forward, but it was still
    helpful to have Passion's version as a guide. 

    As for the radial blur process, there's not much to it. Start off at the pixel position, 
    then radiate outwards gathering up pixels with decreased weighting. The result is a
	blurring of the image in a radial fashion, strangely enough. :)

	Inspired by:

	Blue Dream - Passion
	https://www.shadertoy.com/view/MdG3RD

	Radial Blur - IQ
	https://www.shadertoy.com/view/4sfGRn

	Rays of Blinding Light - mu6k
	https://www.shadertoy.com/view/lsf3Dn

*/

// The radial blur section. Shadertoy user, Passion, did a good enough job, so I've used a
// slightly trimmed down version of that. By the way, there are accumulative weighting 
// methods that do a slightly better job, but this method is good enough for this example.


// Radial blur samples. More is always better, but there's frame rate to consider.
const float SAMPLES = 128.; 


// 2x1 hash. Used to jitter the samples.
float hash( vec2 p ){ return fract(sin(dot(p, vec2(41, 289)))*45758.5453); }


// Light offset.
//
// I realized, after a while, that determining the correct light position doesn't help, since 
// radial blur doesn't really look right unless its focus point is within the screen boundaries, 
// whereas the light is often out of frame. Therefore, I decided to go for something that at 
// least gives the feel of following the light. In this case, I normalized the light position 
// and rotated it in unison with the camera rotation. Hacky, for sure, but who's checking? :)
vec3 lOff(){    
    
    vec2 u = sin(vec2(1.57, 0) - smoothTimeC/2.);
    mat2 a = mat2(u, -u.y, u.x);
    
    vec3 l = normalize(vec3(1.5, 1., -0.5));
    l.xz = a * l.xz;
    l.xy = a * l.xy;
    
    return l;
    
}


vec4 renderPassD() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    // Screen coordinates.
    vec2 uv = fragCoord.xy / RENDERSIZE.xy;

    // Radial blur factors.
    //
    // Falloff, as we radiate outwards.
    float decay = 0.597; 
    // Controls the sample density, which in turn, controls the sample spread.
    float density = 0.725; 
    // Sample weight. Decays as we radiate outwards.
    float weight = 0.1; 
    
    // Light offset. Kind of fake. See above.
    vec3 l = lOff();
    
    // Offset texture position (uv - .5), offset again by the fake light movement.
    // It's used to set the blur direction (a direction vector of sorts), and is used 
    // later to center the spotlight.
    //
    // The range is centered on zero, which allows the accumulation to spread out in
    // all directions. Ie; It's radial.
    vec2 tuv =  uv - .5 - l.xy*.45;
    
    // Dividing the direction vector above by the sample number and a density factor
    // which controls how far the blur spreads out. Higher density means a greater 
    // blur radius.
    vec2 dTuv = tuv*density/SAMPLES;
    
    // Grabbing a portion of the initial texture sample. Higher numbers will make the
    // scene a little clearer, but I'm going for a bit of abstraction.
    vec4 col = texture(BuffC, uv.xy)*0.5;
    
    // Jittering, to get rid of banding. Vitally important when accumulating discontinuous 
    // samples, especially when only a few layers are being used.
    uv += dTuv*(hash(uv.xy + fract(smoothTime))*2. - 1.);
    
    // The radial blur loop. Take a texture sample, move a little in the direction of
    // the radial direction vector (dTuv) then take another, slightly less weighted,
    // sample, add it to the total, then repeat the process until done.
    for(float i=0.; i < SAMPLES; i++){
    
        uv -= dTuv;
        col += texture(BuffA, uv) * weight * col, uv;
        weight *= decay;
        
    }
    
    // Multiplying the final color with a spotlight centered on the focal point of the radial
    // blur. It's a nice finishing touch... that Passion came up with. If it's a good idea,
    // it didn't come from me. :)
    col *= (1. - dot(tuv, tuv)*.75);
    
    // Smoothstepping the final color, just to bring it out a bit, then applying some 
    // loose gamma correction.
    fragColor = sqrt(smoothstep(0., 1., col));
    
    // Bypassing the radial blur to show the raymarched scene on its own.
    //fragColor = sqrt(texture(BuffC, fragCoord.xy / RENDERSIZE.xy));
	return fragColor; 
 } 




// Forked up shit of "Pleasant Yonic" by xenn. https://shadertoy.com/view/ftsSWX
// 2021-07-19 11:02:27


#define AB_SCALE 0.95

vec2 displace(vec2 uv, vec2 offset)
{   
    float d;
    if (media_on == 1.0) {
    d = smoothstep(0.2,2.0,texture(syn_UserImage, (uv*1.0 - vec2(smoothTime /8.0,0)) + offset).r) * 0.25;
}
{
    float d = smoothstep(0.2,2.0,texture(image8, (uv*1.0 - vec2(smoothTime /8.0,0)) + offset).r) * 0.25;
}

    return vec2(d);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
    
    vec2 dr = displace(uv, vec2(0   , 0.02) * AB_SCALE),
         dg = displace(uv, vec2(0.01, 0.01) * AB_SCALE),
         db = displace(uv, vec2(0.01, 0   ) * AB_SCALE);
    
    vec3 color = vec3(0);
    color += vec3(1, 0, 0)*texture(BuffD, uv - dr).r;
    color += vec3(0, 1, 0)*texture(BuffD, uv - dg).g;
    color += vec3(0, 0, 1)*texture(BuffD, uv - db).b;
    
    fragColor = vec4(color, 1.0);
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
		return renderPassD();
	}
	if(PASSINDEX == 4){
		return renderMainImage();
	}
}