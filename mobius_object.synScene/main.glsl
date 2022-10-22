vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


/*

    Mobius Object
    -------------

	I love looking at the Mobius-related renderings that artistically and scientifically inclined people
	put up on places like DeviantArt, etc. There are so many interesting variations and rendering styles 
	out there. There's also a few interesting examples on this site. Anyway, here's yet another one.

	Aesthetically speaking, I started with a clean look that emphasized the Mobius object, then lost focus 
    and got carried away with the surroundings, which for some reason, led to things resembling a grungey, 
	brooding early-2000s demo scene. :) I saved the cleaner version, so I might put that up later.

	Topology was a long time ago for me, so I'm not sure what the object is technically called, but it's 
	easy to see that it's based on the Mobius strip. In descriptive terms, I guess it's a toroidal shape... 
	twisted in the poloidal direction about the toroidal axis? Either way, these things are not much more 
    difficult to code than a torus.

	As for the Mobius object construction, it's pretty straight forward for anyone comfortable with the 
	concepts behind a toroidal distance field. As always, it was helpful to have examples on this site to 
	refer to. One of Dr2's versions helped me correct an annoying quantization error I was making. I'm 
	human, and therefore prone to errors. Thankfully, Dr2 is not. :D

	Other examples:
    
	// Interlinked Mobius strips. Really cool geometry.
	Linked Rings - Dr2    
	https://www.shadertoy.com/view/XsGXR1

	// An Escher recreation: A lot of work went into it.
	Moebius Strip 2 - Dr2
	https://www.shadertoy.com/view/MscXWX

    // I love this one. More variation.
	Twisted Jewelry - vgs
	https://www.shadertoy.com/view/MdjXzG
	Based on:
	Möbius Balls - xTr1m
	https://www.shadertoy.com/view/ldl3zr

*/


#define FAR 40.

// Scene object ID to separate the Mobius object from the terrain and a Mobius object ID to separate 
// the rails from the slats. It was easier to do it this way, but I'll amalgamate them later.
float objID;
float mObjID;

// 2x2 matrix rotation. Note the absence of "cos." It's there, but in disguise, and comes courtesy
// of Fabrice Neyret's "ouside the box" thinking. :)
mat2 r2(float th){ vec2 a = sin(vec2(1.5707963, 0) + th); return mat2(a, -a.y, a.x); }


// IQ's smooth minium function. 
float smin(float a, float b , float s){
    
    float h = clamp( 0.5 + 0.5*(b-a)/s, 0. , 1.);
    return mix(b, a, h) - h*(1.0-h)*s;
}


// Smooth maximum, based on the function above.
float smax(float a, float b, float s){
    
    float h = clamp( 0.5 + 0.5*(a-b)/s, 0., 1.);
    return mix(b, a, h) + h*(1.0-h)*s;
}

// vec2 to vec2 hash.
vec2 hash22(vec2 p) { 

    // Faster, but doesn't disperse things quite as nicely. However, when framerate
    // is an issue, and it often is, this is a good one to use. Basically, it's a tweaked 
    // amalgamation I put together, based on a couple of other random algorithms I've 
    // seen around... so use it with caution, because I make a tonne of mistakes. :)
    float n = sin(dot(p, vec2(41, 289)));
    return fract(vec2(262144, 32768)*n);
    
}

// 2D 2nd-order Voronoi: Obviously, this is just a rehash of IQ's original. I've tidied
// up those if-statements. Since there's less writing, it should go faster. That's how 
// it works, right? :)
//
float Voronoi(in vec2 p){
    
	vec2 g = floor(p), o; p -= g;
	
    // I'm not sure what the largest conceivable closest squared-distance would be, but I think 
    // values as high as 2.5 (center to non-diagonal outer square corner) are possible. Statistically, 
    // it's unlikely, so when scaling back the final value to the zero-to-one range, I divide by 
    // something less than my maximum and cap it to one... It's all a matter of what look you're 
    // trying to achieve.
	vec3 d = vec3(2.5);
    
	for(int y=-1; y<=1; y++){
		for(int x=-1; x<=1; x++){
            
			o = vec2(x, y);
            o += hash22(g + o) - p;
            
			d.z = dot(o, o);
            d.y = max(d.x, min(d.y, d.z));
            d.x = min(d.x, d.z); 
                       
		}
	}

    // Final value, with rough scaling.
    return min((d.y - d.x)*.6, 1.); // Scale: [0, 1].
    
}


// The triangle function that Shadertoy user Nimitz has used in various triangle noise demonstrations.
// See Xyptonjtroz - Very cool. Anyway, it's not really being used to its full potential here.
vec3 tri(in vec3 x){return abs(fract(x)-.5);} // Triangle function.
//vec3 triSmooth(in vec3 x){return cos(x*6.2831853)*0.25+0.25;} // Smooth version. Not used here.


// This is a cheap...ish routine - based on the triangle function - that produces a pronounced jagged 
// looking surface. It's not particularly sophisticated, but it does a surprizingly good job at laying 
// the foundations for a sharp rock face. Obviously, more layers would be more convincing. In fact, 
// I'm disappointed that there weren't enough cycles for one more layer. Unfortunately, this is a 
// GPU-draining distance function. The really fine details have been bump mapped.
float terrain(in vec3 p){
    
    
    // This is just one variation on a common technique: Take a cheap function, then
    // layer it by applying mutations, rotations, frequency and amplitude changes,
    // etc. Feeding the function into itself, folding it, and so forth can also 
    // produce interesting surfaces, patterns, etc.
    //
    // Good examples of the technique include IQ's spiral noise and Nimitz's triangle
    // noise, each of which can be found on Shadertoy. 
    //
    float n = dot(tri(p*0.3 + tri(p.yzx*0.15)), vec3(0.44));
    p = p*1.57-n;//1.5773;// - n; // The "n" mixes things up more.
    p.yz = mat2(.866025, .5, -.5, .866025)*p.yz+smoothTime*0.1;
    p.xz = mat2(.866025, .5, -.5, .866025)*p.xz+smoothTime*0.1;
    n += dot(tri(p*0.545 + tri(p.yzx*0.225)), vec3(0.222));
    
    return smoothstep(0.3, .95, n); // Smoothstep for peaks and troughs. Range [0, 1]


}

// The Mobius object: Take an object, sweep it around in a path (circle) at radius R, and twist (roll) 
// with the axial plane as you do it. Essentially, that's all you're doing.
//
// By the way, I've explained the process in a hurry with a "near enough is good enough" attitude, 
// so if any topology experts out there spot any conceptual errors, mislabling, etc, feel free to 
// let me know.
float Mobius(vec3 q){
 
    //// CONSTANTS ////
    const float toroidRadius = 1.25; // The object's disc radius.
    //const float ringWidth = .15; 
    const float polRot = 4./4.; // Poloidal rotations.
    const float ringNum = 32.; // Number of quantized objects embedded between the rings.
    
    
    //// RAIL SECTION ////
    vec3 p = q;
    
    // Angle of the point on the XZ plane.
    float a = atan(p.z, p.x);
    
    // Angle of the point at the center of 32 (ringNum) partitioned cells.
    //
    // Partitioning the circular path into 32 (ringNum) cells - or sections, then obtaining the angle of 
    // the center position of that cell. The reason you want that angle is so that you can render 
    // something at the corresponding position. In this case, it will be a squared-off ring looking object.  	
    float ia = floor(ringNum*a/6.2831853);  
    // The ".5" value for the angle of the cell center. It was something obvious that I'd overlooked.
    // Thankfully, Dr2 did not. :)
  	ia = (ia + .5)/ringNum*6.2831853; 
    
    // Sweeping a point around a central point at a distance (toroidRadius), more or less. Basically, it's
    // the toroidal axis bit. If that's confusing, looking up a toroidal\poloidal image will clear it up.
    p.xz *= r2(a);
    p.x -= toroidRadius+Radius;
    p.xy *= r2(a*polRot+(smoothTimeC*0.25));  // Twisting about the poloidal direction (controlled by "polRot) as we sweep.
    

    // The rail object. Taking the one rail, then ofsetting it along X and Y, resulting in four rails.
    // This is a neat spacial partitioning trick, and worth knowing if you've never encountered it before.
    // Basically, you're taking the rail, and splitting it into two along X and Y... also along Z, but since 
    // the object is contiunous along that axis, the result is four rails.
    p = abs(abs(p) - .25); // Change this to "p = abs(p)," and you'll see what it does.

    float rail = max(max(p.x, p.y) - .07, (max(p.y-p.x, p.y + p.x)*.7071 - .075)); // Makeshift octagon.
    
    
    //// REPEAT RING SECTION ////
    // The repeat square rings. It's similar to the way in which the rails are constructed, but since the object
    // isn't continous, we need to use the quantized angular positions (using "ia").
    p = q;
    // Another toroidal sweep using the quantized (partitioned, etc) angular position.
    p.xz *= r2(ia); // Using the quantized angle to obtain the position of the center of the corresponding cell.
    p.x -= toroidRadius+Radius;
    p.xy *= r2(a*polRot+smoothTimeC*0.25);  // Twisting about the poloidal direction - as we did with the rails.
    
    // Constructing some square rings.
    p = abs(p);
    float ring = max(p.x, p.y); // Square shape.
    // Square rings: A flat cube, with a thinner square pole taken out.
    ring = max(max(ring - .275, p.z - .03), -(ring - .2));
    
    
    //// WHOLE OBJECT ////
    // Object ID for shading purposes.
    mObjID = step(ring, rail); //smoothstep(0., .07, rail - sqr);
    
    // Smoothly combine (just slightly) the square rings with the rails.
    return smin(ring, rail, .07); 

}


// Combining the Mobius object with the terrain.
float map(vec3 p){
    
    // Mobius object, sitting a bit above the terrain.
    float obj = Mobius(p - vec3(0, .4, 0));
    
    float ter = terrain(p); // The terrain.
 
    float fl = p.y  - ter; // Adding it to a flat plane.
    
    // Creating a flat area to sit the Mobius object on.
    fl =  smax(fl, -max(length(p - vec3(0, 2.5, 0)) - 3., - p.y - .0), .95);
 
    // Object ID.
    objID = step(obj, fl);
    
    // Putting the Mobius object on the terrain.
    return min(fl, obj);
 
}


// Standard raymarching routine.
float trace(vec3 ro, vec3 rd){
   
    float t = 0.0;
    
    for (int i=0; i<96; i++){

        float d = map(ro + rd*t);
        
        if(abs(d)<0.001*(t*.125 + 1.) || t>FAR) break;
        
        t += d*.85;  // Using more accuracy, in the first pass.
    }
    
    return min(t, FAR);
}


// The normal function with some edge detection rolled into it. Sometimes, it's possible to get away
// with six taps, but we need a bit of epsilon value variance here, so there's an extra six.
vec3 getNormal(vec3 p, inout float edge, inout float crv, float ef, float t){ 
	
    vec2 e = vec2(ef/RENDERSIZE.y, 0); // Larger epsilon for greater sample spread, thus thicker edges.

    // Take some distance function measurements from either side of the hit point on all three axes.
	float d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	float d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	float d5 = map(p + e.yyx), d6 = map(p - e.yyx);
	float d = map(p)*2.;	// The hit point itself - Doubled to cut down on calculations. See below.
     
    // Edges - Take a geometry measurement from either side of the hit point. Average them, then see how
    // much the value differs from the hit point itself. Do this for X, Y and Z directions. Here, the sum
    // is used for the overall difference, but there are other ways. Note that it's mainly sharp surface 
    // curves that register a discernible difference.
    edge = abs(d1 + d2 - d) + abs(d3 + d4 - d) + abs(d5 + d6 - d);
    //edge = max(max(abs(d1 + d2 - d), abs(d3 + d4 - d)), abs(d5 + d6 - d)); // Etc.
    
    // Once you have an edge value, it needs to normalized, and smoothed if possible. How you 
    // do that is up to you. This is what I came up with for now, but I might tweak it later.
    edge = smoothstep(0., 1., sqrt(edge/e.x*2.));
    
    
    //crv = clamp((d1 + d2 + d3 + d4 + d5 + d6 - d*3.)*32. + .5, 0., 1.);
	
    // Redoing the calculations for the normal with a more precise epsilon value.
    e = vec2(.002, 0);//*min(1. + t, 5.)
	d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	d5 = map(p + e.yyx), d6 = map(p - e.yyx); 
    
    // Return the normal.
    // Standard, normalized gradient mearsurement.
    return normalize(vec3(d1 - d2, d3 - d4, d5 - d6));
}

// Cheap shadows are hard. In fact, I'd almost say, shadowing particular scenes with limited 
// iterations is impossible... However, I'd be very grateful if someone could prove me wrong. :)
float softShadow(vec3 ro, vec3 lp, float k){

    // More would be nicer. More is always nicer, but not really affordable... Not on my slow test machine, anyway.
    const int maxIterationsShad = 20; 
    
    vec3 rd = (lp-ro); // Unnormalized direction ray.

    float shade = 1.0;
    float dist = 0.025;    
    float end = max(length(rd), 0.001);
    //float stepDist = end/float(maxIterationsShad);
    
    rd /= end;

    // Max shadow iterations - More iterations make nicer shadows, but slow things down. Obviously, the lowest 
    // number to give a decent shadow is the best one to choose. 
    for (int i=0; i<maxIterationsShad; i++){

        float h = map(ro + rd*dist);
        //shade = min(shade, k*h/dist);
        shade = min(shade, smoothstep(0.0, 1.0, k*h/dist)); // Subtle difference. Thanks to IQ for this tidbit.
        //dist += min( h, stepDist ); // So many options here: dist += clamp( h, 0.0005, 0.2 ), etc.
        dist += clamp(h, 0.01, 0.5);
        
        // Early exits from accumulative distance function calls tend to be a good thing.
        if (h<0.001 || dist > end) break; 
    }

    // I've added 0.5 to the final shade value, which lightens the shadow a bit. It's a preference thing.
    return min(max(shade, 0.) + 0.2, 1.0); 
}


// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)
// Anyway, I like this one. I'm assuming it's based on IQ's original.
float cAO(in vec3 pos, in vec3 nor)
{
	float sca = 2.0, occ = 0.0;
    for( int i=0; i<5; i++ ){
    
        float hr = 0.01 + float(i)*0.5/4.0;        
        float dd = map(nor * hr + pos);
        occ += (hr - dd)*sca;
        sca *= 0.7;
    }
    return clamp( 1.0 - occ, 0.0, 1.0 );    
}



 

// The bump mapping function.
float bumpFunction(in vec3 p){
    
    // A reproduction of the lattice at higher frequency. Obviously, you could put
    // anything here. Noise, Voronoi, other geometrical formulas, etc.
    return Voronoi(p.xz*6.);
   
   
}

// Standard function-based bump mapping function with some edging thrown into the mix.
vec3 doBumpMap(in vec3 p, in vec3 n, float bumpfactor, inout float edge, inout float crv){
    
    // Resolution independent sample distance... Basically, I want the lines to be about
    // the same pixel width, regardless of resolution... Coding is annoying sometimes. :)
    vec2 e = vec2(1.5/RENDERSIZE.y, 0); 
    
    float f = bumpFunction(p); // Hit point function sample.
    
    float fx = bumpFunction(p - e.xyy); // Nearby sample in the X-direction.
    float fy = bumpFunction(p - e.yxy); // Nearby sample in the Y-direction.
    float fz = bumpFunction(p - e.yyx); // Nearby sample in the Y-direction.
    
    float fx2 = bumpFunction(p + e.xyy); // Sample in the opposite X-direction.
    float fy2 = bumpFunction(p + e.yxy); // Sample in the opposite Y-direction.
    float fz2 = bumpFunction(p+ e.yyx);  // Sample in the opposite Z-direction.
    
     
    // The gradient vector. Making use of the extra samples to obtain a more locally
    // accurate value. It has a bit of a smoothing effect, which is a bonus.
    vec3 grad = vec3(fx - fx2, fy - fy2, fz - fz2)/(e.x*2.);  
    //vec3 grad = (vec3(fx, fy, fz ) - f)/e.x;  // Without the extra samples.


    // Using the above samples to obtain an edge value. In essence, you're taking some
    // surrounding samples and determining how much they differ from the hit point
    // sample. It's really no different in concept to 2D edging.
    edge = abs(fx + fy + fz + fx2 + fy2 + fz2 - 6.*f);
    edge = smoothstep(0., 1., edge/e.x);
    
    
    // We may as well use the six measurements to obtain a rough curvature value while we're at it.
    //crv = clamp((fx + fy + fz + fx2 + fy2 + fz2 - 6.*f)*32. + .6, 0., 1.);
    
    // Some kind of gradient correction. I'm getting so old that I've forgotten why you
    // do this. It's a simple reason, and a necessary one. I remember that much. :D
    grad -= n*dot(n, grad);          
                      
    return normalize(n + grad*bumpfactor); // Bump the normal with the gradient vector.
	
}

// Tri-Planar blending function. Based on an old Nvidia writeup:
// GPU Gems 3 - Ryan Geiss: http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
vec3 tex3D(sampler2D channel, vec3 p, vec3 n){
    
    n = max(abs(n) - .2, 0.001);
    n /= dot(n, vec3(1));
	vec3 tx = texture(channel, p.zy).xyz;
    vec3 ty = texture(channel, p.xz).xyz;
    vec3 tz = texture(channel, p.xy).xyz;
    
    // Textures are stored in sRGB (I think), so you have to convert them to linear space 
    // (squaring is a rough approximation) prior to working with them... or something like that. :)
    // Once the final color value is gamma corrected, you should see correct looking colors.
    return tx*tx*n.x + ty*ty*n.y + tz*tz*n.z;
}

// Texture bump mapping. Four tri-planar lookups, or 12 texture lookups in total. I tried to 
// make it as concise as possible. Whether that translates to speed, or not, I couldn't say.
vec3 texBump( sampler2D tx, in vec3 p, in vec3 n, float bf){
   
    const vec2 e = vec2(0.002, 0);
    
    // Three gradient vectors rolled into a matrix, constructed with offset greyscale texture values.    
    mat3 m = mat3( tex3D(tx, p - e.xyy, n), tex3D(tx, p - e.yxy, n), tex3D(tx, p - e.yyx, n));
    
    vec3 g = vec3(0.299, 0.587, 0.114)*m; // Converting to greyscale.
    g = (g - dot(tex3D(tx,  p , n), vec3(0.299, 0.587, 0.114)) )/e.x; g -= n*dot(n, g);
                      
    return normalize( n + g*bf ); // Bumped normal. "bf" - bump factor.
	
}

// Cool curve function, by Shadertoy user, Nimitz.
//
// I think it's based on a discrete finite difference approximation to the continuous
// Laplace differential operator? Either way, it gives you the curvature of a surface, 
// which is pretty handy. I used it to do a bit of fake shadowing.
//
// Original usage (I think?) - Cheap curvature: https://www.shadertoy.com/view/Xts3WM
// Other usage: Xyptonjtroz: https://www.shadertoy.com/view/4ts3z2
float curve(in vec3 p, in float w){

    vec2 e = vec2(-1., 1.)*w;
    
    float t1 = map(p + e.yxx), t2 = map(p + e.xxy);
    float t3 = map(p + e.xyx), t4 = map(p + e.yyy);
    
    return 0.125/(w*w) *(t1 + t2 + t3 + t4 - 4.*map(p));
}


// Compact, self-contained version of IQ's 3D value noise function. I have a transparent noise
// example that explains it, if you require it.
float n3D(in vec3 p){
    
	const vec3 s = vec3(7, 157, 113);
	vec3 ip = floor(p); p -= ip; 
    vec4 h = vec4(0., s.yz, s.y + s.z) + dot(ip, s);
    p = p*p*(3. - 2.*p); //p *= p*p*(p*(p * 6. - 15.) + 10.);
    h = mix(fract(sin(h)*43758.5453), fract(sin(h + s.x)*43758.5453), p.x);
    h.xy = mix(h.xz, h.yw, p.y);
    return mix(h.x, h.y, p.z); // Range: [0, 1].
}



// Very basic pseudo environment mapping... and by that, I mean it's fake. :) However, it 
// does give the impression that the surface is reflecting the surrounds in some way.
//
// More sophisticated environment mapping:
// UI easy to integrate - XT95    
// https://www.shadertoy.com/view/ldKSDm
vec3 envMap(vec3 p){
    
    p *= 2.;
    p.xz += smoothTimeB*.5;
    
    float n3D2 = n3D(p*2.);
   
    // A bit of fBm.
    float c = n3D(p)*.157 + n3D2*.28 + n3D(p*4.)*.15;
    c = smoothstep(0.85, 1., c); // Putting in some dark space.
    
    p = vec3(c*c*c*c, c*c, c); // Bluish tinge.
    
    return mix(p.zxy, p, n3D2*.34 + .665); // Mixing in a bit of purple.

}
 

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	
	// Screen coordinates.
	vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;
	
	// Camera Setup.
	vec3 lk = vec3(0, .25, 0);  // "Look At" position.
	vec3 ro = lk + vec3(cos(smoothTime/6.)*2.8, cos(smoothTime/4.)*sin(smoothTime/6.)*.25 + .75, sin(smoothTime/6.)*3.3); // Camera position, doubling as the ray origin.
 
    // Light positioning. One is just in front of the camera, and the other is in front of that.
 	vec3 lp = ro + vec3(.985, 1.75, -1);// Put it a bit in front of the camera.
	


    // Using the above to produce the unit ray-direction vector.
    float FOV = 1.25; // FOV - Field of view.
    vec3 fwd = normalize(lk - ro);
    vec3 rgt = normalize(vec3(fwd.z, 0., -fwd.x )); 
    // "right" and "forward" are perpendicular, due to the dot product being zero. Therefore, I'm 
    // assuming no normalization is necessary? The only reason I ask is that lots of people do 
    // normalize, so perhaps I'm overlooking something?
    vec3 up = cross(fwd, rgt); 

    // rd - Ray direction.
    vec3 rd = normalize(fwd + (rgt*uv.x + up*uv.y)*FOV);

    
    /*   
    // Mouse controls.   
	vec2 ms = vec2(0);
    if (iMouse.z > 1.0) ms = (2.*iMouse.xy - RENDERSIZE.xy)/RENDERSIZE.xy;
    vec2 a = sin(vec2(1.5707963, 0) - ms.x); 
    mat2 rM = mat2(a, -a.y, a.x);
    rd.xz = rd.xz*rM; 
    a = sin(vec2(1.5707963, 0) - ms.y); 
    rM = mat2(a, -a.y, a.x);
    rd.yz = rd.yz*rM;
    */
	 
    
    // Raymarch to the scene.
    float t = trace(ro, rd);
    
    float svObjID = objID;
    float svMObjID = mObjID;
	
    // Initiate the scene color to black.
	vec3 sceneCol = vec3(0);
	
	// The ray has effectively hit the surface, so light it up.
	if(t < FAR){
        
        // Edge and edge-factor. The latter was necessary, to even up the line width between
        // the Mobius object and the terrain... but I ultimately didn't use it. :)
        float edge, crv = 1., ef = 4.; // Curvature variable not used, and commented out in the function.
        
        // Texture scale factor.
        float tSize0 = 1.;
        
        // Texture-based bump mapping factor.
        float bf = .005;
        
        if(svObjID<0.5) { // Different setting for the terrain.
            bf = .02;
            tSize0 = 1./1.;
            ef = 1.5;
        }
    	
    	// Surface position and surface normal.
	    vec3 sp = ro + rd*t;
	    vec3 sn = getNormal(sp, edge, crv, ef, t);
        
        // Texture-based bump mapping. Comment it out, if you prefer a cleaner look.
        // I haven't decided yet. :)
        sn = texBump(image4, sp*tSize0, sn, bf);
        
        // Function-based bump mapping. Note the second edge variable used for bump mapped
        // edging - as opposed to distance field edging. There's also a second curvature
        // variable that isn't used.
        float edge2 = 0., crv2 = 1.;    
        if(svObjID<.5) sn = doBumpMap(sp, sn, .1/(1. + t/FAR), edge2, crv2);    
        

	    
        
        // Obtaining the texel color. 
	    vec3 texCol;        
        
        if(svObjID<0.5) { // Terrain texturing.
            
            texCol = tex3D(image2, sp*tSize0, sn);//*vec3(1, .6, .4);
            texCol = smoothstep(0.05, .5, texCol)*vec3(1, .75, .5);//*vec3(1, .7, .6);//
            texCol *= terrain(sp)*.5 + .5;
            //texCol *= crv*.75 + .25;
        }
        else { // Mobius texturing.
            
            texCol = tex3D(image4, sp*tSize0, sn);
            texCol = smoothstep(0.05, .5, texCol);
            
            // Coloring the Mobius rings. Mixing to avoid nesting an "if" call.
            texCol = mix(texCol, texCol*vec3(1, .4, .2)*2., svMObjID);
        }


    	
    	// Light direction vectors.
	    vec3 ld = lp-sp;

        // Distance from respective lights to the surface point.
	    float lDist = max(length(ld), 0.001);
    	
    	// Normalize the light direction vectors.
	    ld /= lDist;

        
        
        // Shadows and ambient self shadowing.
    	float sh = softShadow(sp, lp, 8.);
    	float ao = cAO(sp, sn); // Ambient occlusion.
	    
	    // Light attenuation, based on the distances above.
	    float atten = 1./(1. + lDist*lDist*0.0125);

    	
    	// Diffuse lighting.
	    float diff = max( dot(sn, ld), 0.0);
        diff = pow(diff, 4.)*1.5; // Ramping up the diffuse.
    	
    	// Specular lighting.
	    float spec = pow(max( dot( reflect(-ld, sn), -rd ), 0.0 ), 16.); 
	    
	    // Fresnel term. Good for giving a surface a bit of a reflective glow.
        float fre = pow( clamp(dot(sn, rd) + 1., .0, 1.), 8.);

        
        
        // I got a reminder looking at XT95's "UI" shader that there are cheaper ways
        // to produce a hint of reflectivity than an actual reflective pass. :)        
        vec3 env = envMap(reflect(rd, sn))*2.; //envMap(refract(rd, sn, 1./1.3))*3.;
        if(svObjID<.5) { // Lowering the terrain settings a bit.
            env *= .25;
            //diff *= .5;
            spec *= .5;            
            fre *= .5;
        }
        

        // Combining the above terms to procude the final color.
        sceneCol += (texCol*(diff + 0.15 + vec3(.7, .9, 1)*fre) + env + vec3(.7, .9, 1)*spec*2.);

        
        // Distance function edging for the Mobius object, and bump mapped edging only
        // for the terrain.
        if(svObjID>.5) sceneCol *= (1. - edge*.8);
        else sceneCol *= (1. - edge2*.7); //*crv2;//(1. - edge*.8)*
        

        // Shading.
        sceneCol *= ao*atten*sh;
        
        
	
	}
    
    // Simple dark fog. It's almost black, but I left a speck of blue in there to account for 
    // the blue reflective glow... Although, it still doesn't explain where it's coming from. :)
    vec3 bg = mix(vec3(.6, .5, 1), vec3(.025, .05, .1), clamp(rd.y + .75, 0., 1.));
    sceneCol = mix(sceneCol, bg/48., smoothstep(0., .95, t/FAR));
    
    // Rought gamma correction.
	fragColor = vec4(sqrt(clamp(sceneCol, 0., 1.)), 1.0);
	
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}