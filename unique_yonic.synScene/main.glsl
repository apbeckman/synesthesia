

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
    col = mix(col,  0.5 + 0.5*cos(smoothTimeC*0.1+uv.y*6.+vec3(0,2,4)),g*100.);
    
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
const float iFeedbackColorShiftImpact = .0033;
const vec2  iDrawCenter               = vec2(0., 0.);
const float iDrawIntensity            = .35;
const float iBlobEdgeSmoothing        = .0567;
const float iBlob1Radius              = .33;
const float iBlob1PowFactor           = 20.;
const float iBlob1ColorPulseSpeed     = .2;
const float iBlob2Radius              = .54;
const float iBlob2PowFactor           = 20.;
const float iBlob2ColorPulseSpeed     = .0234;
const float iBlob2ColorPulseShift     = 1.75;
const float iColorShiftOfRadius       = 1.5;
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
        mod(smoothTimeB * iBlob1ColorPulseSpeed, 1.)
    );
    
    vec3 iBlob2Color = spectral_zucconi6(
        mod(smoothTimeB * iBlob2ColorPulseSpeed + iBlob2ColorPulseShift, 1.)
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

vec3 feedbk = repeatedTexture(BuffC, uv - st).rgb;
    vec3 colorShift = repeatedTexture(
        BuffA,
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

// This buffer is the feedback loop

vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

   	vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    
    // Convert the uv's to polar coordinates to scale up  
    vec2 polarUv = (uv * 2.0 - 1.0);
    float angle = atan(polarUv.y, polarUv.x);
    
    // Scale up the length of the vector by a noise function feeded by the angle and length of the vector
    float ll = length(polarUv)*0.4915+(0.1*sin(smoothTime/3.4));
    
    vec3 base = texture(BuffB, uv).rgb;
    
    // Convert the scaled coordinates back to cartesian
    vec2 offs = vec2(cos(angle)*ll + 0.5, sin(angle)*ll + 0.5);
    
    // sample the last texture with uv's slightly scaled up
    vec3 overlay = texture(BuffA,offs).rgb;

        // Additively blend the colors together
    vec4 col = vec4(base + overlay*(0.1+(0.1*cos(smoothTime/2.))), 01.0);
    
    fragColor = col;
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
const float SAMPLES = 32.; 


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
    float decay = 0.97; 
    // Controls the sample density, which in turn, controls the sample spread.
    float density = 0.5; 
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
    vec4 col = texture(BuffA, uv.xy)*0.25;
    
    // Jittering, to get rid of banding. Vitally important when accumulating discontinuous 
    // samples, especially when only a few layers are being used.
    uv += dTuv*(hash(uv.xy + fract(smoothTime))*2. - 1.);
    
    // The radial blur loop. Take a texture sample, move a little in the direction of
    // the radial direction vector (dTuv) then take another, slightly less weighted,
    // sample, add it to the total, then repeat the process until done.
    for(float i=0.; i < SAMPLES; i++){
    
        uv -= dTuv;
        col += texture(BuffB, uv) * weight;
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
    //fragColor = sqrt(texture(BuffA, fragCoord.xy / RENDERSIZE.xy));
	return fragColor; 
 } 




// Fork of "Smooth Soft Colours" by xenn. https://shadertoy.com/view/7tXSzs
// 2021-07-17 10:10:52

// Fork of "soft Colour" by None. https://shadertoy.com/view/-1
// 2021-07-14 22:40:15

// Fork of "soft Colour" by xenn. https://shadertoy.com/view/NtfXD7
// 2021-07-14 22:39:31

// Postprocess copied with some small modifications from Mattias: https://www.shadertoy.com/view/Ms23DR

vec2 curve(vec2 uv)
{
	uv = (uv - 0.5) * 2.0;
	uv *= 1.1;	
	uv.x *= 1.0 + pow((abs(uv.y) / 5.0), 3.0);
	uv.y *= 1.0 + pow((abs(uv.x) / 4.0), 3.0);
	uv  = (uv / 2.0) + 0.5;
	uv =  uv *0.92 + 0.04;
	return uv;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 q = fragCoord.xy / RENDERSIZE.xy;
    vec2 uv = curve( q );
    vec3 oricol = texture( BuffD, vec2(q.x,q.y) ).xyz;
    vec3 col;
    vec3 sha = texture( BuffC, vec2(q.x,q.y) ).xyz;
   
	float x =  sin(0.3*smoothTime+uv.y*21.0)*sin(0.7*smoothTime+uv.y*29.0)*sin(0.3+0.33*smoothTime+uv.y*31.0)*0.0005 ;

    col.r = texture(BuffD,vec2(x+uv.x+0.001,uv.y+0.001)).x+0.05;
    col.g = texture(BuffD,vec2(x+uv.x+0.000,uv.y-0.002)).y+0.05;
    col.b = texture(BuffD,vec2(x+uv.x-0.002,uv.y+0.000)).z+0.05;
    col.r += 0.08*texture(BuffD,0.75*vec2(x+0.025, -0.027)+vec2(uv.x+0.001,uv.y+0.001)).x;
    col.g += 0.05*texture(BuffD,0.75*vec2(x+-0.022, -0.02)+vec2(uv.x+0.000,uv.y-0.002)).y;
    col.b += 0.08*texture(BuffD,0.75*vec2(x+-0.02, -0.018)+vec2(uv.x-0.002,uv.y+0.000)).z;
    
    float vig = (0.0 + 1.0*16.0*uv.x*uv.y*(1.0-uv.x)*(1.0-uv.y));
	col *= vec3(pow(vig,0.15));

    col *= vec3(0.9,1.1,0.9);
    
	float scans = clamp( 0.35+0.05*sin(3.5*smoothTimeC+uv.y*RENDERSIZE.y*1.5), 0.0, 1.0);
	
	float s = pow(scans,0.75);
	col = col *vec3(0.4+0.7*s) + (sha * col);

    col *= 1.0+0.035*sin(110.0*(smoothTime / 4.0));

	col*=1.0-0.75*vec3(clamp((mod(fragCoord.x, 2.0)-1.0)*1.0,0.0,1.0));
	
    float comp = smoothstep( 0.1, 0.9, sin(smoothTime) );

    fragColor = vec4(col,1.0);
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