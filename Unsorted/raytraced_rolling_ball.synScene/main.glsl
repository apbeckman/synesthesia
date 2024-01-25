

			//******** BuffA Code Begins ********

/*

    Raytraced Rolling Ball
    ----------------------
    
    Using a pseudo path tracing technique to produce a simple realtime scene 
    lit up with multiple emitters.
    
    Some time ago, Xor wrote a minimal shader featuring a glitch texture that
    I was pretty taken with, and figured it'd look interesting in some kind of 
    globally illuminated setting, so I quickly repurposed an old example and 
    put this together the same day... then got side tracked (probably by some 
    other Shadertoy post) and forgot about it. :)
    
    Anyway, I came across it today, so quickly tidied it up and added commnents.
    Rolling spheres down a corridor are a bit of a raytracing cliche, but 
    they're visually effective. There's nothing in here that hasn't been covered 
    before, but someone might get something out of it.
    
    I wanted the code to be partly readable, so didn't optimize things as much
    as I probably should have. However, I'll tweak it later.
    
    
    
    Based on:
    
    // Beautiful in its simplicity. Less is more, as they say. 
    Gltch [291 Chars] - Xor
    https://www.shadertoy.com/view/mdfGRs
    
    // A lot of the realtime path tracing demos out there
    // are based on elements from this example.
    past racer by jetlag - w23
    https://www.shadertoy.com/view/Wts3W7

    

*/


// Sample number and blend number: The trick is to find a balance between the
// two, or use a faster computer. :)

// Number of samples: My computer can handle more. If yours is struggling, you 
// can lower this. Naturally, sample number is proportional to noise quality.
#define sampNum 12



// Standard 2D rotation formula.
mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }


// Fabrice's fork of "Integer Hash - III" by IQ: https://shadertoy.com/view/4tXyWN
float hash21(vec2 f){

    // Depending on your machine, this should be faster than
    // the block below it.
    return texture(image17, f*vec2(.2483, .1437)).x;
    /*
    uvec2 p = floatBitsToUint(f);
    p = 1664525U*(p>>1U^p.yx);
    return float(1103515245U*(p.x^(p.y>>3U)))/float(0xffffffffU);
    */
}

// IQ's "uint" based uvec3 to float hash.
float hash31(vec3 f){

    uvec3 p = floatBitsToUint(f);
    p = 1103515245U*((p >> 2U)^(p.yzx>>1U)^p.zxy);
    uint h32 = 1103515245U*(((p.x)^(p.y>>3U))^(p.z>>6U));

    uint n = h32^(h32 >> 16);
    return float(n & uint(0x7fffffffU))/float(0x7fffffff);
}

// Random seed.
vec2 seed = vec2(.183, .257);


// A slight variation on a function from Nimitz's hash collection, here: 
// Quality hashes collection WebGL2 - https://www.shadertoy.com/view/Xt3cDn
vec2 hash22(){

    // I should probably use a "uvec2" seed, but I hacked this from an old
    // example. I'll update it later.
    seed = fract(seed + vec2(.7123, .6457));
    uvec2 p = floatBitsToUint(seed);
    
    // Modified from: iq's "Integer Hash - III" (https://www.shadertoy.com/view/4tXyWN)
    // Faster than "full" xxHash and good quality.
    p = 1103515245U*((p>>1U)^(p.yx));
    uint h32 = 1103515245U*((p.x)^(p.y>>3U));
    uint n = h32^(h32>>16);

    uvec2 rz = uvec2(n, n*48271U);
    // Standard uvec2 to vec2 conversion with wrapping and normalizing.
    return vec2((rz>>1)&uvec2(0x7fffffffU))/float(0x7fffffff);
}

/*
// A slight variation on a function from Nimitz's hash collection, here: 
// Quality hashes collection WebGL2 - https://www.shadertoy.com/view/Xt3cDn
vec3 hash23(vec2 f){

	uvec2 p = floatBitsToUint(f);
    p = 1664525U*(p>>1U^p.yx); 
    
    uint h32 = 1103515245U*((p.x)^(p.y>>3U));
    uint n = h32^(h32>>16);
    
    // See: http://random.mat.sbg.ac.at/results/karl/server/node4.html
    uvec3 rz = uvec3(n, n*16807U, n*48271U); 
    return vec3((rz >> 1) & uvec3(0x7fffffffU))/float(0x7fffffff);
}
*/

// IQ's box routine.
float sBox(in vec2 p, in vec2 b, float r){

  vec2 d = abs(p) - b + r;
  return min(max(d.x, d.y), 0.) + length(max(d, 0.)) - r;
}

 
// A nice random hemispherical routine taken out of one of IQ's examples.
// The routine itself was written by Fizzer.
vec3 cosDir( in float seed, in vec3 n){

    vec2 rnd = hash22();
    float u = rnd.x;
    float v = rnd.y;
    
    // Method 1 and 2 first generate a frame of reference to use with an arbitrary
    // distribution, cosine in this case. Method 3 (invented by fizzer) specializes 
    // the whole math to the cosine distribution and simplfies the result to a more 
    // compact version that does not depend on a full frame of reference.

    // Method by fizzer: http://www.amietia.com/lambertnotangent.html
    float a = 6.2831853*v;
    u = 2.*u - 1.;
    return normalize(n + vec3(sqrt(1. - u*u)*vec2(cos(a), sin(a)), u));
    
}


// Sphere normal.
vec3 sphereNorm(vec3 p, float id, vec4 sph){
   
    return (p - sph.xyz)/sph.w; 
}
 
// Hitting a number of walls from the inside: You could simply raytrace four
// planes, but this is a little more concise. I was too lazy to write my own
// routine, so quickly adapted a working one (sadly, not many of those around) 
// from one of PublicIntI's examples. At some stage, I'll get in amongst it and 
// rewrite one, or find one of my older routines. Alternatively, if someone
// knows of a concise reliable function or sees a way to tidy the following up, 
// feel free to let me know. :)
//
// crystal exhibit(pathtraced) - public_int_i 
// https://www.shadertoy.com/view/wljSRz
//
// Ray-box intersection: The function take in the ray origin (offset if needed)
// the unit direction ray and the box dimensions, then returns the distance and 
// normal.
//
vec4 boxIntersect(vec3 ro, vec3 rd, vec3 dim) {

    const float maxT = 1e8; // Maximum distance.
 
    vec3 minD = (ro + dim)/rd, maxD = (ro - dim)/rd;
	minD = -(minD - step(vec3(-1e-6), minD)*(minD + maxT));
	maxD = -(maxD - step(vec3(-1e-6), maxD)*(maxD + maxT));
	minD = min(minD, maxD);
    
    // Result: Distance and normal.
    vec4 res = vec4(maxT, 0, 0, 0);
    
    // Performing some ray-plane intersections, modified to handle
    // two planes at once. I'd imagine you could cleverly combine this
    // into just one test, but I'm not clever, so I'll leave that to 
    // someone else. :D
    
    // We don't need the left and right walls for this example.
    if (minD.x<maxT){
        float pd = abs(ro.y + rd.y*minD.x) - dim.y;
        if(pd<0.) res = vec4(minD.x, -sign(rd.x), 0, 0);
    }
    
    // Top and bottom surfaces, or ceiling and floor, if you prefer.
    if (minD.y<maxT){
        float pd = abs(ro.x + rd.x*minD.y) - dim.x;
        if(pd<0.) res = vec4(minD.y, 0., -sign(rd.y), 0.);
    }

    
    // Return the distance and normal.
    return res;
}
 
 
// Sphere intersection: Pretty standard, and adapted from one
// of IQ's formulae.
vec2 sphereIntersect(in vec3 ro, in vec3 rd, in vec4 sph){

    vec3 oc = ro - sph.xyz;
	float b = dot(oc, rd);
    if(b > 0.) return vec2(1e8, 0.);
	float c = dot(oc, oc) - sph.w*sph.w;
	float h = b*b - c;
	if(h<0.) return vec2(1e8, 0.);
	return vec2(-b - sqrt(h), 1.); 
    
}


// Sphere position and radius.
const float ballRad = .5;
vec4 sph4 = vec4(0, ballRad - 1., 1., ballRad);

// Hacking in a normal for the box equation.
vec3 boxNrm;

// Scene normal logic: Not that exciting for this example. :)
vec3 getNorm(vec3 p, float id){
    
    return (id<.5)? sphereNorm(p, id, sph4) : boxNrm; 
}


// Intersection logic for all objects.
vec3 intersect(vec3 ro, vec3 rd){
    
    // Containers for two objects. Usually there'd be more.
    vec2[2] q;
    
    // The sphere.
    q[0] = sphereIntersect(ro, rd, sph4);//vec2(1e5);//
    //q[0].x = 1e5;

    // The box tube object, or 4 walls at once, if you prefer. :)
    vec4 bx = boxIntersect(ro - vec3(0, 1.5 - 1., -.5*0.), rd, vec3(2, 1.5, 1e8));
    q[1] = vec2(bx.x, 1);
    boxNrm = bx.yzw; 
   
    
    // Returning the object distance, a hit ID (inside surface, etc, and redundant 
    // for this example) and the object ID used for materials and so forth.
    return q[0].x<q[1].x? vec3(q[0], 0) : vec3(q[1], 1);
    
    /*
    // For more objects, you need to do a little more work.
    vec3 d = vec3(1e5);
    
    for(int i = 0; i<2; i++){
       if(q[i].x< d.x) d = vec3(q[i], i);
    }
        
    return d;
    */
    
}


// The wall and floor pattern, which is just something quick and effective
// that I picked up from Xor's example, here:
//
// Gltch [291 Chars] - Xor
// https://www.shadertoy.com/view/mdfGRs
vec3 distField2(vec2 p, float scl, float rndZ, float oID){
    
    // Edge width.
    float ew = .0125; 
     
    // Outer boundaries, prior to partitioning.
    vec2 pp = abs(fract(p) - .5);
    float sq = abs(max(pp.x, pp.y) - .5) - ew*2.;
    
    // Scale.
    vec2 sc = vec2(1, 1)/scl; 
 
    // Cell ID and local coordinates.
    vec2 ip = floor(p/sc);
    p -= (ip + .5)*sc;
    
    
    // Rounded square.
    float d = sBox(p, sc/2. - ew, .1*min(sc.x, sc.y)*0.);
    //float d = length(p) - sc.x/2. + ew;
   
    // Randomly rotated and scaled lines.
    if(hash21(ip + rndZ*.123 + oID + .055)<.65){
        float f = scl*4.;
        //float n = floor((floor(hash21(ip + .043)*64.)*2.) + 1.)*3.14159/4.;
        float n = floor(hash21(ip + rndZ*.401 + oID + .043)*64.)*3.14159/4.;
        vec2 uv = rot2(n)*(p + vec2(f/2., 0));
        float g = (abs(fract(uv.x*f + n*f*0.) - .5) - .175)/f;
        d = max(d, g);
        
    }
    
    if(oID==0.) d = max(d, -sq);
    
    // Putting a hole in it just to break things up.
    //d = max(d, -(length(p) - .2*sc.x));
    
    // Rings.
    //d = abs(d + .1*sc.x) - .1*sc.x;
    
    
    
    // Returning the distance and local cell ID. Note that the 
    // distance has been rescaled by the scaling factor.
    return vec3(d, ip*sc);
}


// A simple rotated stripey pattern, like above.
vec3 distField(vec2 p, float scl, float oID){
    

    //Noise macro
    #define N(u) hash21(floor(u) + oID*.123) //texture(BuffA, (u)/64.).x
   
    
    p *= scl;
    float rnd = N(p);
    
    
    vec2 p2 = p/(rnd + .1);
    
    vec2 pp = fract(p2) - .5;
    float sq = (abs(max(abs(pp.x), abs(pp.y)) - .5) - 1./40.)/2.*(rnd + .1);;
    rnd = N(p2);
    // Random quarter turn rotation.
    float n = floor(rnd*64.)*3.14159/4.;
    vec2 id = rnd + ceil(p2);
    // Random stripe frequency.
    float f = 1./N(id)/3.14159; 
    //float f = 4.;

    // Random quarter rotation stripes of random frequency.
    //float g = ceil(cos((rot2(n)*p).x*f));
    //float g = cos((rot2(n)*p).x*f);
    float g = (abs(fract((rot2(n)*p).x*f) - .5) - .25)/f; 

    g = max(g, -sq); 
    
    return vec3(g, id);

}


// IQ's signed square formula with some roundness thrown in. 
float sBoxS(in vec2 p, in vec2 b, in float rf){
  
  vec2 d = abs(p) - b + rf;
  return min(max(d.x, d.y), 0.) + length(max(d, 0.)) - rf;
    
}


// mat3 rotation... I did this in a hurry, but I think it's right. :)
// I have a much better version of this that I'll have to find.
mat3 rot(vec3 ang){
    
    vec3 c = cos(ang), s = sin(ang);

    return mat3(c.x*c.z - s.x*s.y*s.z, -s.x*c.y, -c.x*s.z - s.x*s.y*c.z,
                c.x*s.y*s.z + s.x*c.z, c.x*c.y, c.x*s.y*c.z - s.x*s.z,
                c.y*s.z, -s.y, c.y*c.z);
    
}


// Cube mapping - Adapted from one of Fizzer's routines. 
vec4 cubeMap(vec3 p){

    // Elegant cubic space stepping trick, as seen in many voxel related examples.
    vec3 f = abs(p); f = step(f.zxy, f)*step(f.yzx, f); 
    
    ivec3 idF = ivec3(p.x<.0? -1 : 1, p.y<.0? -1 : 1, p.z<0.? -1 : 1);
    
    ivec3 faceID = (idF + 1)/2 + ivec3(0, 2, 4);
    
    return f.x>.5? vec4(p.yz/p.x, idF.x, faceID.x) : 
           f.y>.5? vec4(p.xz/p.y, idF.y, faceID.y) : vec4(p.xy/p.z, idF.z, faceID.z); 
}
 
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;



     
    float sf = 1./RENDERSIZE.y;
        
    // Screen pixel coordinates.
    float iRes = RENDERSIZE.y;
    vec2 seed0 = fract(TIME/vec2(111.13, 57.61))*vec2(-.143, .457);
    vec2 uv0 = (fragCoord - RENDERSIZE.xy*.5)/iRes;
    
  
    float FOV = 1.; // FOV - Field of view.
    vec3 ro = vec3(0, 0, bass_time*2.);
    // "Look At", "forward" and "up" vectors.
    vec3 lk = ro + vec3(0, -.01, .25);
    vec3 fwd = normalize(lk - ro);
    vec3 rgt = normalize(vec3(fwd.z, 0., -fwd.x )); 
    vec3 up = cross(fwd, rgt); 

    
    // Camera.
    mat3 mCam = mat3(rgt, up, fwd);
    // There are ways to rotate the camera all at once, but these will do.
    mCam *= rot(vec3(0, .05, 0)); 
    mCam *= rot(vec3(0, 0, _noise(-sin(TIME/2.))*.125)); 
    mCam *= rot(vec3(_noise(-cos(TIME/2.))*.125, 0, 0)); 
    
    // Mouse driven camera movement.
    if(_mouse.z>0.){
        vec2 ms = (_mouse.xy/RENDERSIZE.xy - .5)*vec2(3.14159/2.);
        mCam *= rot(vec3(0, ms.y/2., -ms.x));
    }
    
    // Positioning the rolling ball.
    sph4.x -= cos(TIME/2.)*.25; // Left to right.
    sph4.z = ro.z + 4.; // In front of the camera.
    
    // Accumulative color.
    vec3 aCol = vec3(0);
    
    float gT = 1e8;
    float avgT = 0.;
    
    
    for(int j = min(0, int(FRAMECOUNT)); j<sampNum; j++){
        
        // Seed value and jitter.
        seed = uv0 + seed0 + vec2(j*57, j*27)/1321.;
        vec2 jit = hash22()*2. - 1.;
        
        // Jittered UV coordinate.
        vec2 uv = uv0 - jit/RENDERSIZE.y;

        // Using the above to produce the unit ray-direction vector.
        vec3 rd = mCam*normalize(vec3(uv, 1./FOV));
        rd.xz = _rotate(rd.xz, _noise(sin(TIME*0.1)*0.125));
        rd.xy = _rotate(rd.xy, spin*PI);
        // Camera position. Initially set to the ray origin.
        vec3 cam = ro;
        // Surface postion. Also initially set to the ray origin.
        vec3 sp = ro;

        vec3 col = vec3(0);
        
        // Emissive, throughput and sample colors.
        vec3 emissive = vec3(0);
        vec3 through = vec3(1);
        vec3 sCol = vec3(0);
        
        // Fog.
        float fogD = 1e8;
       
        
        // Just three bounces. More looks better, but the extra randomess
        // requires more samples. For static scenes, that's not a problem,
        // but this is a realtime one.
        for(int i = min(0, int(FRAMECOUNT)); i<3; i++){

            
            vec3 scene = intersect(sp, rd); // Scene intersection.

            float t = scene.x; // Scene distance.
            float retVal = scene.y; // Redundant here, but used when refraction is involved.
            float id = scene.z;// Object ID.
            
            // Set the fog distance on the first pass.
            if(i==0){ 
                fogD = t;
                avgT += t/float(sampNum);
               if(j==0) gT = fogD;
               
            }

            sp += rd*t; // Advance the ray position.

  
            if(t<1e8){

                
                vec3 sn = getNorm(sp, id); // Normal.

                vec3 oCol = vec3(0), emissive = vec3(0); // Object color, and emissivity.

                emissive = vec3(0);
                float rough = 0.;

               
                if(id<.5) { 
                   
                    // Placing an offset subdivided grid pattern on the sphere,
                    // then randomly lighting up random cells.
 
                    // Texture coordinates.
                    vec3 txP = sp - sph4.xyz;
                    // Rotation.
                    txP.xy *= rot2(sph4.x/(sph4.w)/2.);
                    txP.yz *= rot2(-sph4.z/(sph4.w));
                    
                    // An icosahedral mapping would probably look nicer, but I
                    // wanted to do something different for this example.
                    vec4 q3 = cubeMap(txP);
                    float faceID = q3.w;
                    
                    // Distance field pattern:
                    // Returns the distance field and cell ID.
                    //vec3 d3 = distField((q3.xy/2. + .5), 2., q3.z);
                    vec3 d3 = distField2((q3.xy/2. + .5), 6., q3.z, faceID);
                    // Distance field isoline boundary.
                    d3.x = smoothstep(0., sf, d3.x);

                    float rnd2 = hash21(d3.yz + q3.z*.051 + faceID + .024);;
                    float sRnd = rnd2;
                    //rnd2 = smoothstep(.4, .45, sin(6.2831*rnd2 + smooth_hightime/1.));
                    rnd2 = smoothstep(.4, .45, sin(6.2831*rnd2*syn_HighLevel + smooth_hightime/1.));
                  
 
                    // Render the pattern on the walls, ceiling and floor.
                    vec3 wCol = .5 + .5*cos(6.2831*hash21(d3.yz + q3.z*5.51 + 
                                            faceID + .374)/2. + vec3(0, 1, 2)*1.1 - 0.);
                    oCol = mix(vec3(.9, .95, 1)*(hash21(d3.yz + q3.z*2.035 + 
                                                faceID + .144)*.5 + .5), vec3(.1), d3.x);
           
                    wCol = wCol*vec3(4, 2, 1);  
                   
                    // Classier toned down colors.
                    //wCol = mix(vec3(1), vec3(1, .4, .2), hash21(d3.yz + .14));
                    // Flipping colors.
                    //wCol =  mix(wCol, wCol.zyx, hash21(d3.yz + .14));   
                    
             
                    emissive = mix(wCol*(rnd2*.785 + .015)*3.*vec3(1, .97, .92), 
                                    vec3(.005), d3.x);
                     // Roughness.
                    rough = hash21(d3.yz + q3.z + .11);
                    //rough = smoothstep(.2, .8, rough)*.25;
                    rough = rough*rough*.3 + .025;
                    
                    // Individual pixel roughness.
                    rough *= hash31(sp + .51)*.5 + .75;
                    //rough = min(rough + hash31(sp + .31)*.025, 1.);
                    
                   
                    if(hash21(d3.yz + faceID + .063)<.5){
                        oCol = vec3(1)*dot(oCol, vec3(.299, .587, .114));
                        emissive = vec3(1)*dot(emissive, vec3(.299, .587, .114));
                    }
 
                   
       
          
               }
               else {

                   
                    // Producing a wall and floor pattern, coloring it, and using
                    // parts to act as emitters.
                    
                    // Back wall or not.
                    float sgn = (abs(sn.z)>.5)? 1. : -1.;
                    // Wall ID.
                    float wID = sgn < 0.? sp.y<0.? 0. : 2. : sp.z<0.? 1. : 3.;
                    
                    // UV coordinates for the walls and floors.
                    vec2 q = sgn>.5? sp.xy : abs(sn.x)>.5? sp.yz : sp.xz;
                    
                    // Vertical strips and horizontal wall strips.
                    float strip = abs(mod(sp.z, 4.) - 2.) - 3./6.;
                    float yStrip = abs(sp.y - .5) - 1.5 + 1./6.;
                    
                    
                    // Distance field pattern:
                    // Returns the distance field and cell ID.
                    vec3 d3;
                    // 
                    if(strip<0. && yStrip<0.) d3 = distField2(q, 6., wID, id);
                    else if(abs(sn.y)>.5) d3 = distField(q, 4., wID);
                    else d3 = distField(q, 2., wID);
                    
                    //d3 = distField(q, 2., wID);
                    
                    // Distance field isoline boundary.
                    d3.x = smoothstep(0., sf, d3.x);
             
 
                    // Cell and wall based colors.
                    vec3 wCol = .5 + .5*cos(6.2831*hash21(d3.yz + wID*.054 + .274)/2. +
                                            vec3(0, 1, 2)*1.1 - 1.5*1.);
                    vec3 wCol2 = .5 + .5*cos(6.2831*hash21(d3.yz + wID*.054 + .273)/2. + 
                                             vec3(0, 1, 2)*1.1);
   
                    // Vertical colors and greyscale areas.
                    if(strip>0.) wCol = vec3(1)*dot(wCol2, vec3(.299, .587, .114));
                    else wCol = wCol2*vec3(4, 2, 1);

                    // Greyscale the wall entirely.
                    //wCol = vec3(1)*dot(wCol, vec3(.299, .587, .114));
                    // Flipping wall colors.
                    //wCol =  mix(wCol, wCol.zyx, hash21(d3.yz + .14));
                    // Classier toned down colors.
                    //wCol = mix(vec3(1), vec3(1, .4, .2), hash21(d3.yz + .14));
                    
   
                    
                    // The wall color pattern.
                    oCol = mix(vec3(.9, .95, 1)*(hash21(d3.yz + wID*.054 + .174)*.5 + .5), 
                                                 vec3(.1), d3.x);

                   
                    // Emissivity.
                    float rnd2 = hash21(d3.yz + .067); 
                    // Periodically blinking emissive colors.
                    rnd2 = smoothstep(.4, .47, sin(6.2831*rnd2 + smooth_hightime/4.)*.5);
                  
                    // Pattern based emissivity -- It doesn't always have to be object 
                    // based. Only adding emissivity to the left and right wall strips.
                    if(abs(sn.y)<.5 && yStrip<0.)
                       emissive = mix(wCol*rnd2*4.*vec3(1, .97, .92), vec3(.005), d3.x);
                    // More orange light.
                    //if(abs(mod(sp.z, 4.) - 2.)<.5) 
                    //     emissive *= mix(wCol2*6.*vec3(1, .97, .92), vec3(.0), d3.x); 
                  
                   
                    // Dark strips.
                    if(abs(sn.x)>.5 && abs(yStrip - .015) - .015<0.){ oCol *= 0.; }   
      
                     
                    // Roughness.
                    rough = hash21(d3.yz + wID*.021 + .11);
                    rough = rough*rough*.3 + .025;
                    // Individual pixel roughness.
                    rough *= hash31(sp + .41)*.5 + .75;
                     
                }
                
                // Different emissivity variations.
                //emissive = mix(emissive, emissive.zyx, .7).yxz;
                //emissive = mix(emissive.zyx, emissive, uv0.y*1.5 + .5);
                //emissive = vec3(1)*dot(emissive, vec3(.299, .587, .114));
  
                   
  
                // I definitely like the more natural way in which colors are applied
                // when rendering this way. We only add surface color when it's been
                // hit by a ray that has visited a light source at some point.
                sCol += emissive*through;
                // Applying this bounce's color to future bounces. For instance, if we
                // hit a pink emitter then hit another surface later, that surface will
                // incorporate a bit of pink into it.
                through *= oCol;

 
                vec3 ref = reflect(rd, sn); // Purely reflected vector.
                vec3 rrd = cosDir(0., sn); // Random half hemisphere vector.
                //vec3 rrd = normalize(hash23() - .5); // Less evenly distributed.

         
                // Mimicking surface inconsistancies with fuzzy reflections.
                // Rougher surfaces have a greater chance of randomly reflecting at any 
                // direction and smoother surfaces are more likely to purely reflect.
                //float rChance = step(rough*2. + .4, hash21(uv + vec2(i*277, j*113) + 
                //                     fract(TIME*.977 + .137)));
                //rd = (mix(rrd, ref, rChance));
                //rd = normalize(mix(ref, rrd, rough));
             
               
                //rd = normalize(ref + rrd*rough);
                rd = normalize(mix(ref, rrd, rough));
                // Not sure this line matters too much, but with the fake random
                // bounce above, I guess rays could head into the surface, so 
                // it's there, just in case.
                if(dot(rd, sn)<0.) rd = -rd; 


                sp += sn*1e-5;
                //rd = ref; // Pure reflection override.

            } 
            
            
             if(aCol.x>1e5) break; // Attempting to reduce compile time. 
        }
      
        // Applying some fog, if necessary. You don't actually see this, but
        // I want it there for completeness.
        sCol = mix(vec3(0), sCol, 1./(1. + fogD*fogD*.02));

        
        // Accumulate the sample color.
        aCol += sCol;
        
        if(sCol.x>1e5) break; // Attempting to reduce compile time.
        
        
    }
    
    // Average color over all samples.
    aCol /= float(sampNum);
    
    
  
    
    // Mix the previous frames in with no camera reprojection.
    // It's OK, but full temporal blur will be experienced.
    vec4 preCol = texelFetch(BuffA, ivec2(fragCoord), 0);
    float blend = (int(int(FRAMECOUNT)) < 2) ? 1. : 1./2.;//1./(1. + length(uv0)*3.); 
    fragColor = mix(preCol, vec4(max(aCol, 0.), avgT), blend);
    
    // No temporal blur, for comparisson.
    //fragColor = vec4(max(aCol, 0.), 1);
    
	return fragColor; 
 } 


/*

    Raytraced Rolling Ball
    ----------------------

    See "Buffer A" for an explanation.

*/


 
// Just a very basic depth of field routine -- I find a lot of it is
// common sense. Basically, you store the scene distance from the camera 
// in the fourth channel, then use it to determine how blurry you want
// your image to be at that particular distance.
//
// For instance, in this case, I want pixels that are 2.25 units away from 
// the camera to be in focus (not blurred) and for things to get more 
// blurry as you move away from that point -- aptly named the focal point 
// for non camera people. :)
//
// I based this on old code of mine, but adopted things that I found in 
// IQ and Nesvi7's examples, which you can find here:
//
// Ladybug - IQ
// https://www.shadertoy.com/view/4tByz3
//
// Cube surface II - Nesvi7
// https://www.shadertoy.com/view/Mty3DV
//
vec3 DpthFld(sampler2D iCh, vec2 uv){
	
    // Focal point and circle of confusion.
    const float focD = 3.5, coc = 2.;
    // Linear distance from either side of the focal point.
    float l = abs(texture(iCh, uv).w - focD - coc) - coc;
    // Using it to calculate the DOF.
    float dof = clamp(l/coc, 0., 1.);
    //dof *= smoothstep(0., .25, length(uv - .5));

    
    // Combine samples. Samples with a larger DOF value are taken further 
    // away from the original point, and as such, appear blurrier.
    vec3 acc = vec3(0);

    for(int i = 0; i<25; i++){
        // Accumulate samples.
        acc += texture(iCh, uv + (vec2(i/5, i%5) - 2.)/vec2(800, 450)*dof).xyz;
        //acc.x *= dof/2.;
    }

    // Return the new variably blurred value.
    return acc /= 25.;
    // Visual debug representation of DOF value.
    //return vec3(dof);
}

 
/*
// Standard 2D rotation formula.
//mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }


// IQ's float to float hash. I've added an extra sine wrapping modulo to
// cater for my annoying AMD based system, which can't wrap sine with a 
// proper degree of accuracy.
float hash11B(float x){ return fract(sin(mod(x, 6.2831853))*43758.5453); }

// This is an amalgamation of old blur and DOF functions with a heap of borrowed 
// lines from Dave Hoskins's much nicer Fibonacci based "Bokeh disc" function, which 
// you can find here: https://www.shadertoy.com/view/4d2Xzw
//
// If you're interested in bokeh, Dave's function above and some of Shadertoy user, 
// Hornet's, are probably the one's you should be looking at. Xor has some cool simple 
// ones on here too.
//
vec4 bokeh(sampler2D iCh, vec2 uv){

    
    vec4 colOrig = texture(iCh, uv);

    // UV based DOF. Focused on the horizontal line, then blurring further away.
    //float r = smoothstep(0., 1., abs(uv.y - .57)/.57)*2.;
    // Focal point and circle of confusion.
    const float focD = 3.5, coc = 1.25;
    // Linear distance from either side of the focal point.
    float l = abs(colOrig.w - focD - coc) - coc;
 
    // Using it to calculate the DOF.
    float ra = clamp(l/coc, 0., 2.);
    //float ra = smoothstep(.3, 1., abs(uv.y - .5)*2.);
    //float ra = mix(clamp(l/coc, 0., 2.), smoothstep(.3, 1., abs(uv.y - .5)*2.), .5);
    ra *= smoothstep(0., .1, length(uv - .5) - .25);

    // Standard Fibonacci distribution calculations, compliments of Dave Hoskins.
    const int iter = 96;
    float rad = 1.6;//max(2.*ra, .5); // Bokeh radius.
    float r = 1.;
	vec4 tot = vec4(0), div = tot;
    vec2 vangle = vec2(0., rad*.01/sqrt(float(iter)));
    #define GA 2.3999632 // Golden angle.
    const mat2 rot = mat2(cos(GA), sin(GA), -sin(GA), cos(GA));

    // Aspect ratio.
    vec2 aspect = vec2(RENDERSIZE.y/RENDERSIZE.x, 1);
    
    
	for (int i = 0; i<iter; i++){
        
        #if 1
        
        // Dave Hoskin's Fibonacci based scattering. Cheaper and much nicer, so
        // it's set as the default.
        // The approx increase in the scale of sqrt(0, 1, 2, 3...).
        r += 1./r;
	    vangle = rot*vangle;
        vec4 col = texture(iCh, uv - (r - 1.)*vangle*aspect, RENDERSIZE.y/450.*1.5);
        
        
        #else
        
        // A hash based random distribution, for anyone who wants to try it.
        //int ii = i%10; // For square bokeh.
        //int jj = i/10;
    
        // Random offset contained within a disk or radius n.
        float fi = float(i) + fract(TIME);
        //vec2 fi = vec2(ii, jj) - 5. + fract(TIME);
        vec2 rnd2 = vec2(hash11B(fi), hash11B(fi + .1)*6.2831);
        vec2 offs = 6.*sqrt(rnd2.x)*vec2(cos(rnd2.y), sin(rnd2.y));
        ////////
         
        // Polygons, if desired. Comment out the line above and comment in
        // the "rot2" formula above, if using it.
        //const float N = 6.;
        //float ra = rnd2.y;
        //float a = (floor(ra*N) + .5)*6.2831859/N;
        //vec2 offs  = mix(rot2(a)*vec2(0, 1), rot2(a + 6.2831859/N)*vec2(0, 1), 
        //                 fract(ra*N));
        //offs *= 6.*sqrt(rnd2.x);
        ////////
        offs *= rad;
        //offs = rad*(offs + (vec2(hash11B(fi), hash11B(fi + .21)) - .5));
        vec4 col = texture(iCh, uv - offs/RENDERSIZE.xy, RENDERSIZE.y/450.*1.5);  
  
        #endif
         
        // Thanks to Dave for figuring out how to tweak the colors to produce brighter 
        // contrast. It's common sense... once someone figures it out for you. :D 
        vec4 bokeh = pow(col, vec4(2));
		tot += bokeh*col*col;
		div += bokeh;
        
	}
    
    
    // Mixing the original value with the bokeh tweaked value according
    // to the depth of field.
    
    // Not entirely correct, but no one will notice here. :)
	return mix(colOrig, colOrig*.25 + tot/div*4., ra);
}

*/

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;



    // Screen oordinates.
    vec2 uv = fragCoord/RENDERSIZE.xy;
    
     // Retrieving the stored color only.
    //vec4 col = texture(BuffA, uv);

    // Depth of field.
    vec4 col = DpthFld(BuffA, uv).xyzz;
    
    
   
    // Mixing in some bokeh. I thought it was a bit much for
    // this example (The "bokeh" function would need to be 
    // commented back in first).
    //vec4 col2 = bokeh(BuffA, uv);
    //col = mix(col, col2, .35);
    
     
    /*
    // Chromatic aberration. Not for this example, but interesting.
    #if 0
    const float focD = 3.5, coc = 1.25;
    // Linear distance from either side of the focal point.
    float l = abs(texture(BuffA, uv).w - focD - coc) - coc;
    // Using it to calculate the DOF.
    float ra = clamp(l/coc, 0., 1.);
    vec2 e = ra*4./RENDERSIZE.xy;
    #else
    vec2 e = 4./RENDERSIZE.xy;
    #endif
    vec4 colX = texture(BuffA, uv + e.xy);
    vec4 col = texture(BuffA, uv);
    vec4 colZ = texture(BuffA, uv - e.xy);
    col = vec4(colX.x, col.y, colZ.z, 1);
    */

    // Subtle vignette. Designers use them to frame things and guide
    // the viewer's eyes toward the center... or something like that.
    //col *= pow(16.*uv.x*uv.y*(1. - uv.x)*(1. - uv.y) , 1./32.);
    // Colored vignette.
    col = mix(col, col.yzxw, smoothstep(.35, .6, length(uv - .5)));


    // Rough gamma correction and screen presentation.
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