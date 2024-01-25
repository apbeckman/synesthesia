

			//******** BuffA Code Begins ********

/*

    Apollonian Construction
    -----------------------
    
    I'm not sure what random scene this is supposed to be... A moon base 
    under construction? Either way, I love apollonian objects. The fractal 
    formula that produces this particular one is surprisingly succinct for 
    such an intricate object.
    
    I won't go into the details regarding apollonian formation, suffice to
    say that you iteratively subdivide and transform space in a manner 
    similar to a Mobius transform... It's been a while, so I'll leave the
    proper explanation to IQ, Knighty, Kali, MLA et al. :)
    
    I've had this example sitting in my account for way too long, so 
    figured it was time to put in some finishing touches and post it. I
    went for a dark scene with glowy lights, which is what you do when you
    want something to look interesting without having to use creativity. :D
    
    I did at least run some cables through the object to give the lights
    some context. The cables are just a distorted 3D Truchet object. There's 
    a trick to making 3D Truchet objects fast that I wasn't able to utilize 
    here on a account of the glowing lights, so it's in need of some 
    tweaking, which I'll do at a later date.
    
    
    
    Related examples:
    
    // Most of Dave Hoskins's examples are really popular,
    // and this one is no different.
    Fractal Explorer - Dave_Hoskins
    https://www.shadertoy.com/view/4s3GW2
    
    // One of the original examples on here and the one
    // most of us refer to.
    Apollonian - IQ
    https://www.shadertoy.com/view/4ds3zn
    
    // Awesome 2D example.
    It's a Question of Time - Rigel
    https://www.shadertoy.com/view/lljfRD
    
    // Much less code than this one. :)
    [SH17A] Apollonian Structure - Shane
    https://www.shadertoy.com/view/4d2BW1

    
*/

// Maximum ray distance.
#define MAXDIST 15.

// Light color - Reddish purple: 0, Greenish blue: 1 
#define COLOR 0

// Standard 2D rotation formula.
mat2 r2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }


// IQ's vec2 to float hash.
float hash31(vec3 p){  
    return fract(sin(dot(p, vec3(113.619, 57.583, 27.897)))*43758.5453); 
} 

// The path is a 2D sinusoid that varies over time, depending upon the frequencies, and amplitudes.
vec2 path(in float z){ 
 
    //return vec2(0); // Straight path.
    
    // Windy weaved path.
    float c = cos(z*3.14159265/4.);
    float s = sin(z*3.14159265/4.);
    //return vec2(1. + c*2.15 - .15, 1. + s*2.15 - .15); 
    return vec2(1. + c*2. - .0, 1. + s*2. - .0); 
    
}

/*
// Commutative smooth minimum function. Provided by Tomkh, and taken 
// from Alex Evans's (aka Statix) talk: 
// http://media.lolrus.mediamolecule.com/AlexEvans_SIGGRAPH-2015.pdf
// Credited to Dave Smith @media molecule.
float smin(float a, float b, float k){

   float f = max(0., 1. - abs(b - a)/k);
   return min(a, b) - k*.25*f*f;
}
*/

// Commutative smooth maximum function. Provided by Tomkh, and taken 
// from Alex Evans's (aka Statix) talk: 
// http://media.lolrus.mediamolecule.com/AlexEvans_SIGGRAPH-2015.pdf
// Credited to Dave Smith @media molecule.
float smax(float a, float b, float k){
    
   float f = max(0., 1. - abs(b - a)/k);
   return max(a, b) + k*.25*f*f;
}


// Tri-Planar blending function. Based on an old Nvidia writeup:
// GPU Gems 3 - Ryan Geiss: http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
vec3 tex3D(sampler2D channel, vec3 p, vec3 n){
    
    // Matching the texture rotation with the object warping. Matching transformed texture 
    // coordinates can be tiring, but it needs to be done if you want things to line up. :)
    //p.xy *= rot(p.z*ZTWIST); 
    //n.xy *= rot(p.z*ZTWIST); 
    
    n = max(abs(n) - .2, 0.001);
    n /= dot(n, vec3(1));
	vec3 tx = texture(channel, p.yz).xyz;
    vec3 ty = texture(channel, p.xz).xyz;
    vec3 tz = texture(channel, p.xy).xyz;
    
    // Textures are stored in sRGB (I think), so you have to convert them to linear space 
    // (squaring is a rough approximation) prior to working with them... or something like that. :)
    // Once the final color value is gamma corrected, you should see correct looking colors.
    return tx*tx*n.x + ty*ty*n.y + tz*tz*n.z;
}

 
    
// Apollonian sphere packing of sorts, based on IQ's cool example here:
// https://www.shadertoy.com/view/4ds3zn
//
float apollonian(vec3 p){

    // Initial radius and scale.
    const float rad = 1.545;
    float scale = 1.;

    // I'd liken the following to iterative hyperbolic transform subdivision... 
    // Although, I'd definitely look it up for a more formal description. :)
    for( int i=0; i<8; i++ ){
        // Repeat space.
        p = mod(p - 1., 2.) - 1.;
        // Sphere inversion. You perform something similar in 
        // hyperbolic geometry.
        float k = rad/dot(p, p);
        // Scale.
        p *= k;
        scale *= k;
    }
    
    // Rendering tiny spheres at each point. After several iterative space inversions
    // above, the thousands of tiny spheres will form a 3D Apollonian gasket, or whatever
    // you wish to call it. Because we're repeating space from the start, there will be
    // infinite overlapping gasket objects sitting side by side.
    return length(p)/scale/4. - .0025;
}

/*
// Toroidal distance. (Not used).
float lengthT(vec2 p){

    // Circular.
    return length(p);
    
    // Superelliptical.
    //float n = 4.;
    //return pow(dot(pow(abs(p), vec2(n)), vec2(1)), 1./n);
    
    // Square, octagon or dodecahedron.
    //p = abs(p);
    //return max(p.x, p.y);
    //return max(max(p.x, p.y), (p.x + p.y)*.7071);
    //return max(max(dot(p, vec2(.8660254, .5)), p.y), max(dot(p.yx, vec2(.8660254, .5)), p.x));
    
}
*/


// Poloidal distance.
float lengthP(vec2 p){

    // Circular.
    return length(p);
    
    // Square, hexagon or octagon.
    //p = abs(p);
    //return max(p.x, p.y);
    //return max(dot(p, vec2(.8660254, .5)), p.y);
    //return max(max(p.x, p.y), (p.x + p.y)*.7071);
    
}


float gID; // Object ID.
vec3 glow; // Glow.

float map(vec3 p) {
    
    
    // The Truchet cabel object.
    
    // Warping the cables just a bit.
    vec3 q = p - sin(p*8. - cos(p.yzx*8.))*.035; 
    
    // Local cell ID and coordinates. 
    vec3 iq = floor(q/2.);
    q -= (iq + .5)*2.;
    
    // Random cell -- and as such, Truchet cable -- rotation.
    float rnd = hash31(iq + .12);
    if(rnd<.33) q = q.yxz;
    else if(rnd<.66) q = q.xzy;
    
    // Repeat light object along the cable, and angle and cell number.
    float aN = 8., a, n;
    
    // Truchet (cable) and light objects. 
    float lat, light;
    

    // First arc edge point.
    vec3 rq = q - vec3(1, 1, 0); 
    // Torus coordinates around the edge point.
    vec2 tor = vec2(length(rq.xy) - 1., rq.z); 
    tor = abs(abs(tor) - .45); // Repeating the torus four more times.
    lat = lengthP(tor); // Torus (cable) object.
    // Repeat lights on this torus.
    a = atan(rq.y, rq.x)/6.2831;
    n = (floor(a*aN) + .5)/aN;
    rq.xy *= r2(-6.2831*n);
    rq.x -= 1.;   
    rq.xz = abs(abs(rq.xz) - .45); // Repeating the lights four more times.
    //
    // Constructing the light object, if applicable.
    float bx = 1e5;
    if(hash31(iq - vec3(1, 1, 0) + n*.043)<.25) bx = length(rq);
    light = bx;     
    
    // Doing the same for the second arc edge point.
    rq = q - vec3(0, -1, -1); 
    a = atan(rq.z, rq.y)/6.2831;
    tor = vec2(length(rq.yz) - 1., rq.x);
    tor = abs(abs(tor) - .45);
    lat = min(lat, lengthP(tor)); 
    n = (floor(a*aN) + .5)/aN;
    rq.yz *= r2(-6.2831*n);
    rq.y -= 1.;
    rq.xy = abs(abs(rq.xy) - .45);
    bx = 1e5;
    if(hash31(iq - vec3(0, -1, -1) + n*.043)<.25) bx = length(rq);
    light = min(light, bx);
   
    // Doing the same for the third arc edge point.
    rq = q - vec3(-1, 0, 1);
    a = atan(rq.z, rq.x)/6.2831;
    tor = vec2(length(rq.xz) - 1., rq.y);
    tor = abs(abs(tor) - .45);
    lat = min(lat, lengthP(tor));
    n = (floor(a*aN) + .5)/aN;
    rq.xz *= r2(-6.2831*n);
    rq.x -= 1.;
    rq.xy = abs(abs(rq.xy) - .45);    
    bx = 1e5;
    if(hash31(iq - vec3(-1, 0, 1) + n*.043)<.25) bx = length(rq);
    light = min(light, bx);
  

    // Alternative lighting.
    //if(hash31(iq + .03)<.35) ball = 1e5;
    //if(mod(iq.x + iq.y + iq.z, 2.)<.5) ball = 1e5;
    
    // Giving the lattice and lights some thickness.
    lat -= .01;
    light -= .05;
    
    // Use the wire distance to create a diode looking object.
    light = max(light, abs(lat - .015) - .01);
  
  
    // Construct the repeat apollonian object to thread the cables through.
    float ga = apollonian(p);
    
     
    // Repeat sphere space to match the repeat apollonian objects.
    q = mod(p, 2.) - 1.;
    //ga = max(ga, (length(q) - 1.22)); // Hollow out.
    // Smooth combine with the main chamber sphere to give a 
    // slightly smoother surface.
    ga = smax(ga, -(length(q) - 1.203), .005);
    
     
    // Object ID.
    gID = lat<ga && lat<light? 2. : ga<light? 1. : 0.;//
    
    // If we've hit the light object, add some glow.
    if(gID==0.){
      
         glow += vec3(1, .2, .1)/(1. + light*light*256.); // Truchet cable lights.
    }
    

    // Minimum scene object.
    return min(lat, min(ga, light)); 
}


float march(vec3 ro, vec3 rd) {

    // Closest and total distance.
    float d, t = hash31(ro + rd)*.15;// Jittering to alleviate glow artifacts.
    
    glow = vec3(0);
    
    vec2 dt = vec2(1e5, 0); // IQ's clever desparkling trick.
    
    int i;
    const int iMax = 128;
    for (i = 0; i<iMax; i++) 
    {        
        d = map(ro + rd*t);
        dt = d<dt.x? vec2(d, dt.x) : dt; // Shuffle things along.
        if (abs(d)<.001*(1. + t*.1) || t>MAXDIST) break;
        
        // Limit the marching step. It's a bit slower but produces better glow.
        t += min(d*.8, .3); 
    
        
    }
    
     
    // If we've run through the entire loop and hit the far boundary, 
    // check to see that we haven't clipped an edge point along the way. 
    // Obvious... to IQ, but it never occurred to me. :)
    if(i>=iMax - 1) t = dt.y;
   
    return min(t, MAXDIST);
}



// Normal function. It's not as fast as the tetrahedral calculation, but more symmetrical.
vec3 nrm(in vec3 p) {
	
    //const vec2 e = vec2(.001, 0);
    //return normalize(vec3(map(p + e.xyy) - map(p - e.xyy), map(p + e.yxy) - map(p - e.yxy),	
    //                      map(p + e.yyx) - map(p - e.yyx)));
     
    // This mess is an attempt to speed up compiler time by contriving a break... It's 
    // based on a suggestion by IQ. I think it works, but I really couldn't say for sure.
    float sgn = 1.;
    vec3 e = vec3(.001, 0, 0), mp = e.zzz; // Spalmer's clever zeroing.
    for(int i = min(int(FRAMECOUNT), 0); i<6; i++){
		mp.x += map(p + sgn*e)*sgn;
        sgn = -sgn;
        if((i&1)==1){ mp = mp.yzx; e = e.zxy; }
    }
    
    return normalize(mp);
}


// Cheap shadows are hard. In fact, I'd almost say, shadowing particular scenes with limited 
// iterations is impossible... However, I'd be very grateful if someone could prove me wrong. :)
float shadow(vec3 ro, vec3 lp, vec3 n, float k){

    // More would be nicer. More is always nicer, but not always affordable. :)
    const int maxIterationsShad = 48; 
    
    ro += n*.0015; // Coincides with the hit condition in the "trace" function.  
    vec3 rd = lp - ro; // Unnormalized direction ray.

    float shade = 1.;
    float t = 0.; 
    float end = max(length(rd), .0001);
    //float stepDist = end/float(maxIterationsShad);
    rd /= end;

    // Max shadow iterations - More iterations make nicer shadows, but slow things down. Obviously, 
    // the lowest number to give a decent shadow is the best one to choose. 
    for (int i = min(int(FRAMECOUNT), 0); i<maxIterationsShad; i++){

        float d = map(ro + rd*t);
        shade = min(shade, k*d/t);
        //shade = min(shade, smoothstep(0., 1., k*d/t)); // Thanks to IQ for this tidbit.
        // So many options here, and none are perfect: dist += clamp(d, .01, stepDist), etc.
        t += clamp(d*.8, .005, .25); 
        
        
        // Early exits from accumulative distance function calls tend to be a good thing.
        if (d<0. || t>end) break; 
    }

    // Sometimes, I'll add a constant to the final shade value, which lightens the shadow a bit --
    // It's a preference thing. Really dark shadows look too brutal to me. Sometimes, I'll add 
    // AO also just for kicks. :)
    return max(shade, 0.); 
}


// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)
// Anyway, I like this one. I'm assuming it's based on IQ's original.
float cAO(in vec3 pos, in vec3 nor)
{
	float sca = 1.5, occ = 0.0;
    for( int i=0; i<5; i++ ){
    
        float hr = 0.01 + float(i)*0.35/4.0;        
        float dd = map(nor * hr + pos);
        occ += (hr - dd)*sca;
        sca *= 0.7;
    }
    return clamp( 1.0 - occ, 0.0, 1.0 );    
}



vec4 render(vec3 ro, vec3 rd, vec3 lp){ 


    // Raymarch.
	float t = march(ro, rd);
    
    // Saving the object ID.
    float svID = gID;
    // Saving the glow.
    vec3 svGlow = glow;
    
    
    // Initialize to the scene color to zero.
    vec3 col = vec3(0);
    
    // If we've hit the surface color it.
    if (t<MAXDIST){ 
    
        // Surface hit point.
        vec3 sp = ro + rd*t;
    
        // Normal.
        vec3 sn = nrm(sp); 
         
        // Light direction vector and attenuation.
        vec3 ld = (lp - sp);
        float atten = max(length(ld), .001);
        ld /= atten;
        // Light attenuation.
        atten = 1./(1. + atten*.125);
        
        // Diffuse and specular.
        float diff = max(dot(ld, sn), 0.);
        //diff = pow(diff, 4.)*2.;
        float spec = pow(max(dot(reflect(rd, sn), ld), 0.), 8.);
        
        // Ambient occlusion and shadows.
        float ao = cAO(sp, sn);
        float sh = shadow(sp, lp, sn, 8.); 

        
        // Manipulate the glow color.
        float sinF = dot(sin(sp - cos(sp.yzx)*1.57), vec3(1./6.)) + .5;
        //sinF = smoothstep(0., 1., sinF);
        #if COLOR == 0
        svGlow = mix(svGlow, svGlow.zyx, sinF);
        #else
        svGlow = mix(svGlow.zxy, svGlow.zyx, sinF);
        #endif
        
        // Texturing the scene objcts.
        vec3 tx = tex3D(image3, sp*2., sn);
        col = tx;
        

        if (svID==1.) col /= 8.; // Darken the apollonian structure.
        else if (svID==2.) col /= 4.; // Darken the wires a little less.
        
        
        // Applying the diffuse and specular.
        col = col*(diff*sh + .2) + spec*sh*.5;
 
        
        // Adding the glow. In a path tracer, emitters are influenced by and
        // influence everything else, but here, we just add in some glow and hope
        // people don't look into it too much. :)
        col += col*(svGlow)*32.;
    
        /// Attenuation and ambient occlusion.
        col *= atten*ao;

    
    }
    
    
    // Fog out the background.
    col = mix(col, vec3(0), smoothstep(.0, .8, t/MAXDIST));
    
    
    // Color and distance.
   	return vec4(col, t);
}

// The following is based on John Hable's Uncharted 2 tone mapping, which
// I feel does a really good job at toning down the high color frequencies
// whilst leaving the essence of the gamma corrected linear image intact.
//
// To arrive at this non-tweakable overly simplified formula, I've plugged
// in the most basic settings that work with scenes like this, then cut things
// right back. Anyway, if you want to read about the extended formula, here
// it is.
//
// http://filmicworlds.com/blog/filmic-tonemapping-with-piecewise-power-curves/
// A nice rounded article to read. 
// https://64.github.io/tonemapping/#uncharted-2
vec4 uTone(vec4 x){
    return ((x*(x*.6 + .1) + .004)/(x*(x*.6 + 1.)  + .06) - .0667)*1.933526;    
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


	
    vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;
    
    //if(abs(uv.y)>.4){ fragColor = vec4(0); return; }
    
	// Camera Setup.
	vec3 ro = vec3(0, 0, bass_time); // Camera position, doubling as the ray origin.
	vec3 lk = ro + vec3(0, -.05, .25);  // "Look At" position.
    vec3 lp = ro + vec3(0, 0, .75); // Light position.
    
   
	// Using the Z-value to perturb the XY-plane.
	// Sending the camera, "look at," and light vector down the path. The "path" function is 
	// synchronized with the distance function.
    ro.xy += path(ro.z);
	lk.xy += path(lk.z);
	lp.xy += path(lp.z);
    
    

    // Using the above to produce the unit ray-direction vector.
    float FOV = 3.14159/1.5; // FOV - Field of view.
    vec3 fwd = normalize(lk - ro);
    vec3 rgt = normalize(vec3(fwd.z, 0., -fwd.x )); 
    vec3 up = cross(fwd, rgt);

    // rd - Ray direction.
    vec3 rd = normalize(uv.x*rgt + uv.y*up + fwd/FOV);  
   
    // Camera swivel - based on path position.
    vec2 sw = path(lk.z);
    rd.xy *= r2(-sw.x/24.);
    rd.yz *= r2(-sw.y/16.);    
    
   
    // Render.
    vec4 col4 = render(ro, rd, lp);
    
    // Applying tone mapping. I use it to tone down really high intensity lights, 
    // so I want something that does that. A lot of tone mappers completely change
    // the look, which is fine, if that's what you're going for... It's not really
    // my area, but there is a link in the "uTone" function for those interested.
    col4.xyz = uTone(col4).xyz; // Don't apply this to the alpha channel.
    
    #if 1
    // Mix the previous frames in with no camera reprojection.
    // It's OK, but full temporal blur will be experienced.
    vec4 preCol = texelFetch(BuffA, ivec2(fragCoord), 0);
    float blend = (FRAMECOUNT < 2) ? 1. : 1./3.; 
    fragColor = mix(preCol, vec4(clamp(col4.xyz, 0., 0.5), col4.w), blend);
    #else
    // Output to screen.
    fragColor = vec4(clamp(col4.xyz, 0., 1.), col4.w);
    #endif
    
    
	
	return fragColor; 
 } 


/*

    Apollonian Construction
    -----------------------

    Rendering the buffer.
 
    See Buffer A for an explanation.
    
*/

// IQ's vec2 to float hash.
float hash21(vec2 p){  return fract(sin(dot(p, vec2(27.619, 57.583)))*43758.5453); }

// This is an amalgamation of old blur and DOF functions with a couple of borrowed 
// lines from Dave Hoskins's much nicer Fibonacci based "Bokeh disc" function, which 
// you can find here: https://www.shadertoy.com/view/4d2Xzw
//
// This function is only really suitable for this example. If you're interested in 
// bokeh, Dave's function above and some of Shadertoy user, Hornet's, are probably
// the one's you should be looking at. Xor has some cool simple ones on here that I'm
// yet to dig into, but they might worth a look also.
vec4 bloom(sampler2D iCh, vec2 uv){

	vec4 tot = vec4(0);
    
    // UV based DOF. Focused on the horizontal line, then blurring further away.
    //float r = smoothstep(0., 1., abs(uv.y - .57)/.57)*2.;
    // Focal point and circle of confusion.
    const float focD = 2.5, coc = 1.;
    // Linear distance from either side of the focal point.
    float l = abs(texture(iCh, uv).w - focD - coc) - coc;
    // Using it to calculate the DOF.
    float r = clamp(l/coc, 0., 1.);
    
    const int n = 4;
    for (int j = -n; j<=n; j++){
        for (int i = -n; i<=n; i++){
           
            // Random offset contained within a disk or radius n.
            vec2 rnd2 = vec2(hash21(vec2(i, j)), hash21(vec2(i, j) + .1)*6.2831);
            vec2 offs = float(n)*rnd2.x*vec2(cos(rnd2.y), sin(rnd2.y));
            
            vec4 c = texture(iCh, uv + offs/vec2(800, 450)*r, r*RENDERSIZE.y/450.*.7); 
            tot += mix(c, pow(c, vec4(1.25))*3.4, rnd2.x*rnd2.x); //ow(c, vec4(1.5))*4.
            
        }
    }
    
	return tot/float((n*2 + 1)*(n*2 + 1));
}



vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    // Rendering the buffer.
    vec2 uv = fragCoord/RENDERSIZE.xy;
    
    // Retrieving the stored color.
    vec4 col = texture(BuffA, uv);
   
    // Custom DOF bloom-like function.
    col = bloom(BuffA, uv);
    
    

    // Rough gamma correction and screen presentation.
    // "col" should already be above zero, but we're capping it anyway.
    fragColor = pow(max(col, 0.), vec4(1./2.2));
    
	return fragColor; 
 } 




vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}