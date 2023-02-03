/*

    KIFS Flythrough
    ---------------

    After looking at Zackpudil's recent fractal shaders, I thought I'd put together something 
	fractal in nature. It's nothing exciting, just the standard IFS stuff you see here and there. 
	Like many examples, this particular one is based on Syntopia and Knighty's work.

    The construction is pretty similar to that of an infinite sponge, but it has a bit of rotating, 
	folding, stretching, etc, thrown into the mix.

    The blueish environmental lighting is experimental, and based on XT95s environment mapping in 
	his UI example. The idea is very simple: Instead of passing a reflective ray into a cubic 
	texture in cartesian form, convert it to its polar angles, then index into a 2D texture. The 
	results are abstract, and no substitute for the real thing, but not too bad, all things 
	considered.

	The comments are a little rushed, but I'll tidy them up later.

	Examples and references:

    Menger Journey - Syntopia
    https://www.shadertoy.com/view/Mdf3z7

    // Explains the process in more detail.
    Kaleidoscopic (escape time) IFS - Knighty
    http://www.fractalforums.com/ifs-iterated-function-systems/kaleidoscopic-(escape-time-ifs)/

    Ancient Generators - Zackpudil
    https://www.shadertoy.com/view/4sGXzV

*/

#define FAR 40.0
#define TAU 2*PI

float hash( float n ){ return fract(cos(n)*45758.5453); }

// Tri-Planar blending function. Based on an old Nvidia writeup:
// GPU Gems 3 - Ryan Geiss: http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
vec3 tex3D( sampler2D tex, in vec3 p, in vec3 n ){    
    
    // Ryan Geiss effectively multiplies the first line by 7. It took me a while to realize
    // that it's redundant, due to the normalization that follows. I'd never noticed on account
    // of the fact that I'm not in the habit of questioning stuff written by Ryan Geiss. :)
    n = max(abs(n) - 0.2, 0.001); // n = max(abs(n), 0.001), etc.
    n /= (n.x + n.y + n.z ); 
	p = (texture(tex, p.yz)*n.x + texture(tex, p.zx)*n.y + texture(tex, p.xy)*n.z).xyz;
    
    // Rought sRGB to linear RGB conversion in preperation for eventual gamma correction.
    return p*p; 
}


// IFS, or KIFS in particular. The following explains the process in more detail.
//
// Kaleidoscopic (escape time) IFS - Knighty
// http://www.fractalforums.com/ifs-iterated-function-systems/kaleidoscopic-(escape-time-ifs)/
//
// Here's a quick, rushed expanation. Take an input point, and repeat it. After that, rotate
// it, fold it, stretch and translate it about an offset point. Form an object with it and compare
// it to the object formed in the repeat layer. Repeat ad infinitum...
//
// OK, that was a little vague, but it really is a pretty simple process. Playing around with the 
// code will give it more meaning. Change the rotation angles, iteration number, comment some
// things out, etc. I'd also recommend reading Syntopia's blog. He explains things really well...
// and there's so many pretty pictures. :)
//
// Syntopia - http://blog.hvidtfeldts.net/
//
float map(vec3 p){
    
    // I'm never sure whether I should take constant stuff like the following outside the function, 
    // or not. My 1990s CPU brain tells me outside, but it doesn't seem to make a difference to frame 
    // rate in this environment one way or the other, so I'll keep it where it looks tidy. If a GPU
    // architecture\compiler expert is out there, feel free to let me know.
    
    vec3 offs = vec3(1.0*xy_param.x, 0.75*xy_param.y, 0.5*z_param); // Offset point.
    // const vec3 offs = vec3(1.0, 0.75, 0.5);
    const vec2 a = sin(vec2(0, 1.57079632) + 1.57/2.);
    const mat2 m = mat2(a.y, -a.x, a);
    vec2 a2 = sin(vec2(0, 1.57079632) + 1.57/4. + p.z*0.3+0.75*mid_time);
    //a2 +=Warp*_rotate(_uvc, PI*(smoothTimeC*0.01))/PI;
    mat2 m2 = mat2(a2.y, -a2.x, a2);
    
    const float s = 5.; // Scale factor.
    
    float d = 1e5; // Distance.
    
    
    p  = abs(fract(p*.5)*2. - 1.); // Standard spacial repetition.
     
    
    float amp = 1./s; // Analogous to layer amplitude.
    
   
    // With only two iterations, you could unroll this for more speed,
    // but I'm leaving it this way for anyone who wants to try more
    // iterations.
    for(int i=0; i<2; i++){
        
        // Rotating.
        p.xy = m*p.xy;
        p.yz = m2*p.yz;
        
        p = abs(p);

  		// Folding about tetrahedral planes of symmetry... I think, or is it octahedral? 
        // I should know this stuff, but topology was many years ago for me. In fact, 
        // everything was years ago. :)
		// Branchless equivalent to: if (p.x<p.y) p.xy = p.yx;
        p.xy += step(p.x, p.y)*(p.yx - p.xy);
        p.xz += step(p.x, p.z)*(p.zx - p.xz);
        p.yz += step(p.y, p.z)*(p.zy - p.yz);
 
        // Stretching about an offset.
		p = p*s + offs*(1. - s);
        
		// Branchless equivalent to:
        // if( p.z < offs.z*(1. - s)*.5)  p.z -= offs.z*(1. - s);
        p.z -= step(p.z, offs.z*(2. - s)*.5)*offs.z*(1. - s);
        
        // Loosely speaking, construct an object, and combine it with
        // the object from the previous iteration. The object and
        // comparison are a cube and minimum, but all kinds of 
        // combinations are possible.
        p=abs(p);
//standard cubes
        d = min(d, max(max(p.x, p.y), p.z)*amp);

//cool pattern        
        // d = mix(d, min(d, mix(max(max(p.x, p.y), p.z), length(p), sin(p.y*2.0))*amp*0.5), bubbly);

//pointy        
        // d = mix(d, min(d, mix(max(max(p.x, p.y), p.z), length(p), sin(p.z*2.0+TIME))*amp*0.5), pointy);

        // d = min(d, max(max(p.x, p.y), p.z)*amp*1.5);

        amp /= s; // Decrease the amplitude by the scaling factor.
    }
 
 	return d - .035 - size; // .35 is analous to the object size.
}


// Texture bump mapping. Four tri-planar lookups, or 12 texture lookups in total.
// vec3 db( sampler2D tx, in vec3 p, in vec3 n, float bf){
   
//     const vec2 e = vec2(0.001, 0);
//     // Gradient vector, constructed with offset greyscale texture values.
//     //vec3 g = vec3( gr(tpl(tx, p - e.xyy, n)), gr(tpl(tx, p - e.yxy, n)), gr(tpl(tx, p - e.yyx, n)));
    
//     mat3 m = mat3( tex3D(tx, p - e.xyy, n), tex3D(tx, p - e.yxy, n), tex3D(tx, p - e.yyx, n));
    
//     vec3 g = vec3(0.299, 0.587, 0.114)*m;
//     g = (g - dot(tex3D(tx,  p , n), vec3(0.299, 0.587, 0.114)) )/e.x; g -= n*dot(n, g);
                      
//     return normalize( n + g*bf ); // Bumped normal. "bf" - bump factor.
	
// }


// Very basic raymarching equation.
float trace(vec3 ro, vec3 rd){

    
    float t = 0.;//hash(dot(rd, vec3(7, 157, 113)))*0.01;
    for(int i=0; i< 96-liquify*liquify*92; i++){
        
        float d = map(ro + rd*t);
        if (d < 0.00125*(1. + t*.125) || t>FAR*0.5*far_wall*far_wall) break;
        t += d;
    } 
    return t;
}



// Tetrahedral normal, to save a couple of "map" calls. Courtesy of IQ.
// Apart from being faster, it can produce a subtley different aesthetic to the 6 tap version, which I sometimes prefer.
vec3 normal(in vec3 p){

    // Note the slightly increased sampling distance, to alleviate artifacts due to hit point inaccuracies.
    vec2 e = vec2(0.005, -0.005); 
    return normalize(e.xyy * map(p + e.xyy) + e.yyx * map(p + e.yyx) + e.yxy * map(p + e.yxy) + e.xxx * map(p + e.xxx));
}

/*
// Standard normal function.
vec3 normal(in vec3 p) {
	const vec2 e = vec2(0.005, 0);
	return normalize(vec3(map(p + e.xyy) - map(p - e.xyy), map(p + e.yxy) - map(p - e.yxy),	map(p + e.yyx) - map(p - e.yyx)));
}
*/




// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)
// Anyway, I like this one. I'm assuming it's based on IQ's original.
float calculateAO(in vec3 pos, in vec3 nor)
{
	float sca = 1.5, occ = 0.0;
    for( int i=0; i<5; i++ ){
    
        float hr = 0.01 + float(i)*0.5/4.0;        
        float dd = map(nor * hr + pos);
        occ += (hr - dd)*sca;
        sca *= 0.7;
    }
    return clamp( 1.0 - occ, 0.0, 1.0 );    
}


// Cheap shadows are hard. In fact, I'd almost say, shadowing repeat objects - in a setting like this - with limited 
// iterations is impossible... However, I'd be very grateful if someone could prove me wrong. :)
float softShadow(vec3 ro, vec3 lp, float k){

    // More would be nicer. More is always nicer, but not really affordable... Not on my slow test machine, anyway.
    const int maxIterationsShad = 16; 
    
    vec3 rd = (lp-ro); // Unnormalized direction ray.

    float shade = 1.0;
    float dist = 0.05;    
    float end = max(length(rd), 0.001);
    float stepDist = end/float(maxIterationsShad);
    
    rd /= end;

    // Max shadow iterations - More iterations make nicer shadows, but slow things down. Obviously, the lowest 
    // number to give a decent shadow is the best one to choose. 
    for (int i=0; i<maxIterationsShad; i++){

        float h = map(ro + rd*dist);
        //shade = min(shade, k*h/dist);
        shade = min(shade, smoothstep(0.0, 1.0, k*h/dist)); // Subtle difference. Thanks to IQ for this tidbit.
        //dist += min( h, stepDist ); // So many options here: dist += clamp( h, 0.0005, 0.2 ), etc.
        dist += clamp(h, 0.02, 0.25);
        
        // Early exits from accumulative distance function calls tend to be a good thing.
        if (h<0.001 || dist > end) break; 
    }

    // I've added 0.5 to the final shade value, which lightens the shadow a bit. It's a preference thing.
    return min(max(shade, 0.) + 0.25, 1.0); 
}


//////
// Very basic pseudo environment mapping... and by that, I mean it's fake. :) However, it 
// does give the impression that the surface is reflecting the surrounds in some way.
//
// Anyway, the idea is very simple. Obtain the reflected ray at the surface hit point, 
// convert the 3D ray coordinates (rd) to polar coordinates (the angular components) and 
// index into a repeat texture. It can be pretty convincing (in an abstract way) and
// facilitates environment mapping without the need for a cube map, or a reflective pass.
//
// More sophisticated environment mapping:
// UI easy to integrate - XT95    
// https://www.shadertoy.com/view/ldKSDm
vec3 envMap(vec3 rd, vec3 n){
    
    //vec3 col2 = tex3D(image10, rd/4., n).zyx;//*(1.-lod*.8)
    //return smoothstep(.0, 1., col2*2.);
    

    // I got myself a little turned around here, but I think texture size
    // manipulation has to be performed at this level, if you want the angular
    // polar coordinates to wrap... Not sure though... It'll do. :)
    rd /= 8.; 
    
    vec2 uv = vec2(atan(rd.y, rd.x)/6.283, acos(rd.z)/3.14159);
    uv = fract(uv);
   
    vec3 col = texture(image10, uv).zyx;//*(1.-lod*.8)
    return smoothstep(.01, 1., col*col*.85);
    
}



// Pseudo environment mapping. Simlilar to above, but using tri-planar texturing for a more 
// even spread.
// vec3 envMap(vec3 rd, vec3 n){
    
//     vec3 col = tex3D(image10, rd/4., n).zyx;
//     return smoothstep(.1, 1., col*2.);
// }

//////


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    
    // Unit direction ray vector: Note the absence of a divide term. I came across
    // this via a comment Shadertoy user "coyote" made. I'm pretty happy with this.
    vec3 rd = (vec3(2.*fragCoord - RENDERSIZE.xy, RENDERSIZE.y));
    // Barrel distortion;
    // rd = normalize(vec3(rd.xy, sqrt(rd.z*rd.z - dot(rd.xy, rd.xy)*disto_knob)));
     rd = mix(rd, normalize(vec3(_uvc+_uvc*PI*FOV, 1.0)), 1.);
    //rd = normalize(vec3(_uvc+FOV*(_uvc*PI)_uv, 1.0));
    //rd.xy += PI*FOV*_uvc;
    
    
    // Rotating the ray with Fabrice's cost cuttting matrix. I'm still pretty happy with this also. :)
    // vec2 m = sin(vec2(1.57079632, 0) + blur_out*2*PI);
    // rd.xy = rd.xy*mat2(m.xy, -m.y, m.x);
    // rd.xz = rd.xz*mat2(m.xy, -m.y, m.x);
    rd.yz = _rotate(rd.yz, lookXY.y*PI+_uvc.y*Flip*PI*FOV);
    rd.xz = _rotate(rd.xz, lookXY.x*PI);
    rd.xy = _rotate(rd.xy, Spin*PI+spin_time);
    // Ray origin, set off in the YZ direction. Note the "0.5." It's an old lattice trick.
    vec3 ro = vec3(0.0, 0.0, mix(smoothTime*0.35, bass_time*1.0, 1.0));
    
    vec3 lp = ro + vec3(0.0, .25, .25); // Light, near the ray origin.

    // Set the scene color to black.
    vec3 col = vec3(0);

    
    float t = trace(ro, rd); // Raymarch.
    
    // Surface hit, so light it up.
    if(t<FAR){
    //rd.xy += PI*FOV*_uvc;
        
        vec3 sp = ro + rd*t; // Surface position.
        
        vec3 sn = normal(sp); // Surface normal.


        const float sz = 1.; // Texture size.

        // sn = db(image3, sp*sz, sn, .002/(1. + t/FAR)); // Texture bump.

        vec3 ref = reflect(rd, sn); // Reflected ray.


        vec3 oCol = tex3D(image3, sp*sz, sn); // Texture color at the surface point.
        oCol = smoothstep(0., 1., oCol)*standard_lighting;

        float sh = softShadow(sp, lp, 24.); // Soft shadows.
        float ao = calculateAO(sp, sn)*.5 + .5; // Self shadows. Not too much.

        vec3 ld = lp - sp; // Light direction.
        float lDist = max(length(ld), 0.0001); // Light to surface distance.
        ld /= lDist; // Normalizing the light direction vector.

        float diff = max(dot(ld, sn), 0.); // Diffuse component.
        float spec = pow(max(dot(reflect(-ld, sn), -rd), 0.), 16.); // Specular.

        float atten = 1.0 / (1.0 + lDist*0.25 + lDist*lDist*.075); // Attenuation.

        // Combining the elements above to light and color the scene.
        col = oCol*(diff + vec3(1.4, .25, .2)*(1.0+syn_HighHits)) + vec3(.6, .6, .2)*spec*(1+syn_HighLevel)*(standard_lighting);
        
        // Faux environmental mapping.
        col += envMap(reflect(rd, sn), sn)*standard_lighting;
        vec3 texo = tex3D(geom, sp*sz*0.5, vec3(1.0,1.0,0.0));
        if (syn_MediaType>0.5){
            texo = 1.0-tex3D(syn_UserImage, sp*sz, vec3(11.0,1.0,0.0));
        }

        vec3 colPal = _palette(fract(sp.z), vec3(0.720, 0.860, -0.2040), vec3(0.500, 1.000, 0.500), vec3(0.500, 0.345, 0.750), vec3(0.500, 0.440, 0.220));
        col += (1.0-texo)*mix(colPal, vec3(0.1260,0.4,01.8)+syn_HighHits, sp.y/RENDERSIZE.y)*10.0*_pulse(length(ro-sp)*10.0, scan_alt_light*50.0, 1.0*(1.0+alternate_lighting*(50.0)));
        // Environment mapping with a cubic texture, for comparison.


        // Shading the scene color, clamping, and we're done.
        col = min(col*atten*pow(ao,10.0)*(1.0+syn_HighLevel), 1.);
        
        

        //col = clamp(col + hash(dot(rd, vec3(7, 157, 113)))*0.1 - 0.05, 0., 1.);
        
    }
    
    
    // Blend the scene and the background with some very basic, 8-layered smokey haze.
    // float mist = getMist(ro, rd, lp, t);
    // float mist = 0.5;
    // vec3 sky = vec3(.35, .6, 1)* mix(1., .75, mist);//*(rd.y*.25 + 1.);
    
    // // Mix the smokey haze with the object.
    // col = mix(sky, col, 1./(t*t/FAR/FAR*128. + 1.));

      vec3 bg = mix(vec3(1.750, .250, 0.1250), vec3(0.120, 0.0, 1.0), -rd.y*.35 + .35)*0.5*(1.0+syn_HighLevel);


      col = mix(col, bg, smoothstep(0., FAR-15., t));//min(bg.zyx*vec3(1.3, .6, .2)*1.5, 1.)

      float vig = pow(16.*_uv.x*_uv.y*(1.-_uv.x)*(1.-_uv.y), 0.125);
      col *= vec3(1.2, 1.1, .85)*vig;

    // Statistically unlikely 2.0 gamma correction, but it'll do. :)
	fragColor = vec4(sqrt(clamp(col, 0., 1.)), 1);
    float collisionDist = 0.01;
    // if (t < collisionDist){
    //     fragColor.rgb = vec3(t*100);
    // }
    // if (t>FAR*0.05){
    //     fragColor *= 0.5;
    // }
	return fragColor; 
 } 

 vec3 filter_() {
  vec2 delta = 1. / RENDERSIZE;
  
  vec3 val = texture(fractalPass, _uv).xyz;

  vec3 l = texture(fractalPass, _uv + vec2(0., delta.y)).xyz;
  vec3 r = texture(fractalPass, _uv - vec2(0., delta.y)).xyz;
  vec3 u = texture(fractalPass, _uv + vec2(delta.x, 0.)).xyz;
  vec3 d = texture(fractalPass, _uv - vec2(delta.x, 0.)).xyz;

  vec3 n = vec3(_rand(_uvc+fract(TIME))) - 0.5;
  
  vec3 bloom = max(val, max(max(l, r), max( u, d)));
  bloom = bloom  + l + r + u + d;
  bloom /= 5.; // orlando;
  return bloom;

}


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
    else if(PASSINDEX == 1){
        vec4 fractal = vec4(filter_(),1.0);
        return fractal;
    }
}
