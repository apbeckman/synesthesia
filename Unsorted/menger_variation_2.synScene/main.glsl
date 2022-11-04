

/*

    Menger Sponge Variation
    -----------------------

	This is just a variation on my Menger Sponge example. I put in some rounded corners,
	bracing, better lighting and some post processing effects. It was just an excuse to
	make something artsy, so not that exciting.

	It seems to take a while to compile. The reflection pass would add to things, but it
	feels more like a WebGL 2 thing... Yeah, I know, only bad programmers blame their tools. :D	

	For anyone who's never put a Menger Sponge together, here's a very, very short, overly 
	generalized explanation:
	
	Construct a Void Cube (or repeat Void Cubes, as the case may be), which is analogous 
	to a Rubix Cube with the center mechanism removed. Create Void Cubes from the 20 cubies 
	(the remaining smaller cubes), and continue to iterate ad infinitum.
	
	In code:

	// Repeat Void Cubes - A Void Cube is a Level-1 Menger Sponge.
	float map(vec3 p){
    	p = abs(mod(p, 3.) - 1.5); // Repeat space.
    	return min(max(p.x, p.y), min(max(p.y, p.z), max(p.x, p.z))) - 1.; // Void Cube.
	}

	// More than one level Menger Sponge - Infinitely repeated, in this case.
	float map(vec3 q){
        
		vec3 p; float d = 0.;
        
        // One Void Cube.
    	p = abs(mod(q, 3.) - 1.5);
    	d = max(d, min(max(p.x, p.y), min(max(p.y, p.z), max(p.x, p.z))) - 1.);

        // Subdividing into more Void Cubes.    
    	p = abs(mod(q, 1.) - 0.5); // Dividing above by 3.
    	d = max(d, min(max(p.x, p.y), min(max(p.y, p.z), max(p.x, p.z))) - 1./3.);
        
        // And so on.
    	p = abs(mod(q, 1./3.) - 0.5/3.); // Dividing above by 3.
    	d = max(d, min(max(p.x, p.y), min(max(p.y, p.z), max(p.x, p.z))) - 1./3./3.);
        
		// Continue on in this manner. For more levels, you'll want to loop it. There's
		// a commented out example in the code somewhere. Also, you can experiment with 
		// the code to create more interesting variants.

		return d;
	}
	
	For a more formal explanation, look up "Menger Sponge," "Cantor Sets," "Void Cube," 
	etc., on the web, or better yet, refer to the countless Menger Sponge examples
	on this site.	

	Examples:

	Menger Journey - Syntopia (A favorite of mine, and everyone else.)
	https://www.shadertoy.com/view/Mdf3z7

*/


#define FAR 50. // Maximum ray distance. Analogous to the far plane.

#define BRACING // Chrome bracing or no bracing. Obviously, faster without it.


// Scene object ID. Either the Menger object (1) or the chrome bracing (0).
float objID;
float svObjID; // Global ID to keep a copy of the above from pass to pass.

float hash(float n){ return fract(sin(n)*43758.5453); }

// Fabrice's consice, 2D rotation formula.
//mat2 r2(float th){ vec2 a = sin(vec2(1.5707963, 0) + th); return mat2(a, -a.y, a.x); }
// Standard 2D rotation formula.
mat2 r2(in float a){ float c = cos(a), s = sin(a); return mat2(c, s, -s, c); }

// The path is a 2D sinusoid that varies over time, depending upon the frequencies, and amplitudes.
vec2 path(in float t){ 

    return vec2(0);
    //float s = sin(t/24.)*cos(t/12.);
    //return vec2(s*4., 0.);
    
    float a = sin(t*.11);
    float b = cos(t*.14);
    return vec2(a*4. -b*1.5, b*1.7 + a*1.5);
    
    //return vec2(sin(t*.15)*2.4, cos(t*.25)*1.7*.5); 
}

// Tri-Planar blending function. Based on an old Nvidia tutorial.
vec3 tex3D( sampler2D t, in vec3 p, in vec3 n ){ 
     
    //p.xy *= r2(p.z*zRot);
    //n.xy *= r2(p.z*zRot);    
    
    n = max(abs(n), 0.001);
    n /= dot(n, vec3(1));
	vec3 tx = texture(t, p.yz).xyz;
    vec3 ty = texture(t, p.zx).xyz;
    vec3 tz = texture(t, p.xy).xyz;
    
    // Textures are stored in sRGB (I think), so you have to convert them to linear space 
    // (squaring is a rough approximation) prior to working with them... or something like that. :)
    // Once the final color value is gamma corrected, you should see correct looking colors.
    return (tx*tx*n.x + ty*ty*n.y + tz*tz*n.z);
}

// Compact, self-contained version of IQ's 3D value noise function.
float n3D(vec3 p){
    
	const vec3 s = vec3(7., 157., 113.);
	vec3 ip = floor(p); p -= ip; 
    vec4 h = vec4(0., s.yz, s.y + s.z) + dot(ip, s);
    //p = p*p*(3. - 2.*p);
    p *= p*p*(p*(p * 6. - 15.) + 10.);
    h = mix(fract(sin(h)*43758.5453), fract(sin(h + s.x)*43758.5453), p.x);
    h.xy = mix(h.xz, h.yw, p.y);
    return mix(h.x, h.y, p.z); // Range: [0, 1].
}



// Smooth maximum, based on IQ's smooth minimum.
float smax(float a, float b, float s){
    
    float h = clamp(.5 + .5*(a - b)/s, 0., 1.);
    return mix(b, a, h) + h*(1. - h)*s;
}



// Variation on a Menger Sponge (See the formula above). This one has four layers. The 
// easiest way to understand this is to comment out layers, then add them back in to 
// see what each does.
float map(vec3 q){
    

    q.xy -= path(q.z); // Wrap the object around the path.

    
    // Layer one. The ".05" on the end varies the hole size.
 	vec3 p = abs(fract(q/3.)*3. - 1.5);
 	float d = min(max(p.x, p.y), min(max(p.y, p.z), max(p.x, p.z))) - 1. + .03;
   // d+= sin(TIME);
    ///////////
    // Chrome bracing section.
    #ifdef BRACING
    float tb;
    
    float rp = 3.;
    p = abs(fract(q/rp + .5)*rp - rp/2.);
    float x1 = sqrt(min(dot(p.xy, p.xy),min(dot(p.yz, p.yz),dot(p.xz, p.xz))))-rp/2.*.475; // EQN 2

    
    // Repeat field entity two, which is just an abstract object repeated every half unit. 
    rp = 3./15.;
    p = abs(fract(q/rp)*rp - rp/2.);
    p = abs(p - rp/2.);
    float x2 = min(p.x, min(p.y,p.z)) - .0125; // EQN 1
    

    // Combining the two entities above.
    tb = smax(abs(x1), abs(x2), .1) - .035; 
    #endif
    /////////////
    
 
    
    // Layer two.
    p =  abs(fract(q) - .5);    
 	d = smax(d, min(max(p.x, p.y), min(max(p.y, p.z), max(p.x, p.z))) - 1./3. + .05, .05);
    
   
    // Layer three. 3D space is divided by two, instead of three, to give some variance.
    p =  abs(fract(q*2.)*.5 - .25);
 	d = max(d, min(smax(p.x, p.y, .125), min(smax(p.y, p.z, .125), smax(p.x, p.z, .125))) - .5/3. - .025); 

    // Layer four. The little holes, for fine detailing.
    p =  abs(fract(q*3./.5)*.5/3. - .5/6.);
 	d = max(d, min(max(p.x, p.y), min(max(p.y, p.z), max(p.x, p.z))) - 1./18. - .015);
    //d = max(d, max(max(p.x, p.y), p.z) - 1./18. - .024);
    //d = max(d, length(p) - 1./18. - .048);
    
    #ifdef BRACING
    objID = step(d, tb);
    return min(d, tb);
    #else
    objID = 0.;
    return d;
    #endif
    
}


// Standard raymarching routine.
float trace(vec3 ro, vec3 rd){
   
    float t = 0., d;
    
    for (int i=0; i<96; i++){

        d = map(ro + rd*t);
        
        if(abs(d)<.001*(t*.125 + 1.) || t>FAR) break;//.001*(t*.125 + 1.)
        
        t += d; // Using slightly more accuracy in the first pass.
    }
    
    return min(t, FAR);
}

// Second pass, which is the first, and only, reflected bounce. 
// Virtually the same as above, but with fewer iterations and less 
// accuracy.
//
// The reason for a second, virtually identical equation is that 
// raymarching is usually a pretty expensive exercise, so since the 
// reflected ray doesn't require as much detail, you can relax things 
// a bit - in the hope of speeding things up a little.
float traceRef(vec3 ro, vec3 rd){
    
    float t = 0., d;
    
    for (int i=0; i<48; i++){

        d = map(ro + rd*t);//*rDir;
        
        if(abs(d)<0.001*(t*.125 + 1.) || t>FAR) break;
        
        t += d;
    }
    
    return min(t, FAR);
}

// Another pass, which is the first, and only, refracted bounce. 
// Virtually the same as above, but uses a slimmed down distance function - due
// to the fact that the water plane doesn't need to be included.
float traceRefr(vec3 ro, vec3 rd){
    
    float t = 0., d;
    
    for (int i=0; i<48; i++){

        d = map(ro + rd*t);
        
        if((d<0. && abs(d)<0.001*(t*.25 + 1.)) || t>FAR) break;
        
        t += -d;
    }
    
    return min(t, FAR);
}


// Cheap shadows are the bain of my raymarching existence, since trying to alleviate artifacts is an excercise in
// futility. In fact, I'd almost say, shadowing - in a setting like this - with limited  iterations is impossible... 
// However, I'd be very grateful if someone could prove me wrong. :)
float softShadow(vec3 ro, vec3 lp, float k, float t){

    // More would be nicer. More is always nicer, but not really affordable... Not on my slow test machine, anyway.
    const int maxIterationsShad = 24; 
    
    vec3 rd = lp-ro; // Unnormalized direction ray.

    float shade = 1.;
    float dist = .0025*(t*.125 + 1.);  // Coincides with the hit condition in the "trace" function.  
    float end = max(length(rd), 0.0001);
    //float stepDist = end/float(maxIterationsShad);
    rd /= end;

    // Max shadow iterations - More iterations make nicer shadows, but slow things down. Obviously, the lowest 
    // number to give a decent shadow is the best one to choose. 
    for (int i=0; i<maxIterationsShad; i++){

        float h = map(ro + rd*dist);
        //shade = min(shade, k*h/dist);
        shade = min(shade, smoothstep(0.0, 1.0, k*h/dist)); // Subtle difference. Thanks to IQ for this tidbit.
        // So many options here, and none are perfect: dist += min(h, .2), dist += clamp(h, .01, stepDist), etc.
        dist += clamp(h, .01, .2); 
        
        // Early exits from accumulative distance function calls tend to be a good thing.
        if (h<0.0 || dist > end) break; 
    }

    // I've added a constant to the final shade value, which lightens the shadow a bit.
    return min(max(shade, 0.) + .2, 1.); 
}

/*
// Standard normal function. It's not as fast as the tetrahedral calculation, but more symmetrical. 
// Due to the intricacies of this particular scene, it's kind of needed to reduce jagged effects.
vec3 getNormal(in vec3 p) {
	const vec2 e = vec2(0.0025, 0);
	return normalize(vec3(map(p + e.xyy) - map(p - e.xyy), map(p + e.yxy) - map(p - e.yxy),	
                          map(p + e.yyx) - map(p - e.yyx)));
}
 

// Normal calculation, with some edging and curvature bundled in.
vec3 getNormal(vec3 p, inout float edge, inout float crv) { 
	
    // Roughly two pixel edge spread, regardless of resolution.
    vec2 e = vec2(3./RENDERSIZE.y, 0);

	float d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	float d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	float d5 = map(p + e.yyx), d6 = map(p - e.yyx);
	float d = map(p)*2.;

    edge = abs(d1 + d2 - d) + abs(d3 + d4 - d) + abs(d5 + d6 - d);
    //edge = abs(d1 + d2 + d3 + d4 + d5 + d6 - d*3.);
    edge = smoothstep(0., 1., sqrt(edge/e.x*2.));
     
    // Wider sample spread for the curvature.
    //e = vec2(12./450., 0);
	//d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	//d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	//d5 = map(p + e.yyx), d6 = map(p - e.yyx);
    //crv = clamp((d1 + d2 + d3 + d4 + d5 + d6 - d*3.)*32. + .5, 0., 1.);
 
    
    e = vec2(.0015, 0); //RENDERSIZE.y - Depending how you want different resolutions to look.
	d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	d5 = map(p + e.yyx), d6 = map(p - e.yyx);
	
    return normalize(vec3(d1 - d2, d3 - d4, d5 - d6));
}
*/

// Standard normal function. It's not as fast as the tetrahedral calculation, but more symmetrical.
float getEdge(in vec3 p, in vec2 e) { 

     
    // This mess is an attempt to speed up compiler time by contriving a break... It's 
    // based on a suggestion by IQ. I think it works, but I really couldn't say for sure.
    float sgn = 1.;
    float mp[6];
    vec3[3] e6 = vec3[3](e.xyy, e.yxy, e.yyx);
    for(int i = min(0, int(FRAMECOUNT)); i<6; i++){
		mp[i] = map(p + sgn*e6[i/2]);
        sgn = -sgn;
        if(sgn>2.) break; // Fake conditional break;
    }
        
    float d = map(p)*2.;

    float edge = abs(mp[0] + mp[1] - d) + abs(mp[2] + mp[3] - d) + abs(mp[4] + mp[5] - d);
    //edge = abs(mp[0] + mp[1] + mp[2] + mp[3] + mp[4] + mp[5] - d*3.);
    edge = smoothstep(0., 1., sqrt(edge/e.x*2.));
    
    return edge;
}

// Standard normal function. It's not as fast as the tetrahedral calculation, but more symmetrical.
vec3 getNrm(in vec3 p, in vec2 e) {
    
    //vec3 n = normalize(vec3(map(p + e.xyy) - map(p - e.xyy),
    //map(p + e.yxy) - map(p - e.yxy),	map(p + e.yyx) - map(p - e.yyx)));
    
    // This mess is an attempt to speed up compiler time by contriving a break... It's 
    // based on a suggestion by IQ. I think it works, but I really couldn't say for sure.
    float sgn = 1.;
    float mp[6];
    vec3[3] e6 = vec3[3](e.xyy, e.yxy, e.yyx);
    for(int i = min(0, int(FRAMECOUNT)); i<6; i++){
		mp[i] = map(p + sgn*e6[i/2]);
        sgn = -sgn;
        if(sgn>2.) break; // Fake conditional break;
    }
    
    return normalize(vec3(mp[0] - mp[1], mp[2] - mp[3], mp[4] - mp[5]));
}


// Normal calculation, with some edging and curvature bundled in.
//
// Addendum: I've rewritten this in a very contrived and ugly form to 
// appease the compiler. It seems to work, but I still don't like it. :)
vec3 getNormal(vec3 p, inout float edge, inout float crv, float t) { 
	
    
    // Roughly two pixel edge spread, regardless of resolution.
    vec2 e = vec2(3./RENDERSIZE.y*(1. + t*.5), 0);
    
    edge = getEdge(p, e);
/*   
	float d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	float d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	float d5 = map(p + e.yyx), d6 = map(p - e.yyx);
	float d = map(p)*2.;
    
    edge = abs(d1 + d2 - d) + abs(d3 + d4 - d) + abs(d5 + d6 - d);
    //edge = abs(d1 + d2 + d3 + d4 + d5 + d6 - d*3.);
    edge = smoothstep(0., 1., sqrt(edge/e.x*2.));
    */
/*    
    // Wider sample spread for the curvature.
    e = vec2(12./450., 0);
	d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	d5 = map(p + e.yyx), d6 = map(p - e.yyx);
    crv = clamp((d1 + d2 + d3 + d4 + d5 + d6 - d*3.)*32. + .5, 0., 1.);
*/
    
    e = vec2(.0015, 0); //RENDERSIZE.y - Depending how you want different resolutions to look.
    /*
    d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	d5 = map(p + e.yyx), d6 = map(p - e.yyx);
	
    return normalize(vec3(d1 - d2, d3 - d4, d5 - d6));
    */
    
    return getNrm(p, e);
}

/*
// Normal calculation, with some edging and curvature bundled in.
vec3 getNormalRefr(vec3 p, inout float edge, inout float crv) { 
	
    // Roughly two pixel edge spread, regardless of resolution.
    vec2 e = vec2(3./RENDERSIZE.y, 0);

	float d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	float d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	float d5 = map(p + e.yyx), d6 = map(p - e.yyx);
	float d = map(p)*2.;

    edge = abs(d1 + d2 - d) + abs(d3 + d4 - d) + abs(d5 + d6 - d);
    //edge = abs(d1 + d2 + d3 + d4 + d5 + d6 - d*3.);
    edge = smoothstep(0., 1., sqrt(edge/e.x*2.));
   
    // Wider sample spread for the curvature.
    //e = vec2(12./450., 0);
	//d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	//d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	//d5 = map(p + e.yyx), d6 = map(p - e.yyx);
    //crv = clamp((d1 + d2 + d3 + d4 + d5 + d6 - d*3.)*32. + .5, 0., 1.);

    
    e = vec2(.0015, 0); //RENDERSIZE.y - Depending how you want different resolutions to look.
	d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	d5 = map(p + e.yyx), d6 = map(p - e.yyx);
	
    return normalize(vec3(d1 - d2, d3 - d4, d5 - d6));
}
*/

// Coloring\texturing the scene objects, according to the object IDs.
vec3 getObjectColor(vec3 p, vec3 n){
    
    // Object texture color.

    // Contorting the texture coordinates to math the contorted scene.
    vec3 txP = p - vec3(path(p.z), 0.);
    
    // Grab the texel for the hit point.
    vec3 tx = tex3D(image3, txP/1., n);
    tx = mix(tx, vec3(1)*dot(tx, vec3(.299, .587, .114)), .25);
    tx = smoothstep(0., 1., tx); // Ramp up the color.

    vec3 lCol = vec3(1.25, 1, .75);
    
    
    // Darken part of the sponge object, and color up the rest.
    vec3 q = abs(mod(txP, 3.) - 1.5);
    if (max(max(q.x, q.y), q.z) > 1.063) tx *= lCol;
    else tx *= vec3(1., .75, .5);
    
    // Shine up the bracing.
    #ifdef BRACING
    if(svObjID<.5) tx *= 5./vec3(1.25, 1, .75);
    #endif
    
    tx *= vec3(1.2, 1, .8); // Golden hue. 

    
    return tx; // Return the object color.
    
}

// Using the hit point, unit direction ray, etc, to color the scene. Diffuse, specular, falloff, etc. 
// It's all pretty standard stuff.
vec3 doColor(in vec3 sp, in vec3 rd, in vec3 sn, in vec3 lp, float edge, float crv, float ao, float t){
    
    // Initiate the scene (for this pass) to zero.
    vec3 sceneCol = vec3(0);
    
    if(t<FAR){ // If we've hit a scene object, light it up.
    
        vec3 ld = lp - sp; // Light direction vector.
        float lDist = max(length(ld), 0.001); // Light to surface distance.
        ld /= lDist; // Normalizing the light vector.

        // Attenuating the light, based on distance.
        float atten = 2./(1. + lDist*0.125 + lDist*lDist*0.025);

        // Standard diffuse term.
        float diff = max(dot(sn, ld), 0.);
        //diff = pow(diff, 2.)*.66 + pow(diff, 4.)*.34;
        // Standard specualr term.
        float spec = pow(max( dot( reflect(-ld, sn), -rd ), 0.0 ), 32.0);
        float fres = clamp(1. + dot(rd, sn), 0., 1.);
        float Schlick = pow( 1. - max(dot(rd, normalize(rd + ld)), 0.), 5.0);
        float fre2 = mix(.5, 1., Schlick);  //F0 = .5.

        // Coloring the object. You could set it to a single color, to
        // make things simpler, if you wanted.        
        vec3 objCol = getObjectColor(sp, sn);

        // Combining the above terms to produce the final scene color.
        sceneCol = objCol*(diff + .6*ao + fres*fres) + vec3(1, .7, .5).zyx*spec*2.;
        
        // Edges and curvature.
        //sceneCol *= clamp(crv, 0., 1.);
        //sceneCol += (sceneCol*.75 + .25)*edge;
        sceneCol *= 1. - edge*.9;
        

        // Attenuation only. To save cycles, the shadows and ambient occlusion
        // from the first pass only are used.
        sceneCol *= atten;
    
    }
    
  
    // Return the color. Done once for each pass.
    return sceneCol;
    
}

// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)
// Anyway, I like this one. I'm assuming it's based on IQ's original.
float calculateAO(in vec3 pos, in vec3 nor)
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


// Simple environment mapping. Pass the reflected vector in and create some
// colored noise with it. The normal is redundant here, but it can be used
// to pass into a 3D texture mapping function to produce some interesting
// environmental reflections.
//
// More sophisticated environment mapping:
// UI easy to integrate - XT95    
// https://www.shadertoy.com/view/ldKSDm
vec3 eMap(vec3 rd, vec3 sn){
    
    return getObjectColor(rd, sn);
    
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    // Screen coordinates.
	vec2 uv = (fragCoord.xy - RENDERSIZE.xy*.5) / RENDERSIZE.y;
    
	
	// Camera Setup.
	vec3 ro = vec3(0, 0, smoothTime); // Camera position, doubling as the ray origin.
	vec3 lk = ro + vec3(0, 0, .25);  // "Look At" position.

   
    // Light position. Set in the vicinity the ray origin.
    vec3 lp = ro + vec3(0, 1, 0);
    
   
	// Using the Z-value to perturb the XY-plane.
	// Sending the camera, "look at," and light vector down the path, which is 
	// synchronized with the distance function.
    /*
    ro.xy += path(ro.z);
	lk.xy += path(lk.z);
	lp.xy += path(lp.z);
 */
    
/* replacing with below code for ease of adding controls AB
    // Using the above to produce the unit ray-direction vector.
    float FOV = 3.14159/2.; // FOV - Field of view.
    vec3 forward = normalize(lk-ro);
    vec3 right = normalize(vec3(forward.z, 0., -forward.x )); 
    vec3 up = cross(forward, right);

    // rd - Ray direction.
    vec3 rd = normalize(forward + FOV*uv.x*right + FOV*uv.y*up);
*/  

    float FOV = 3.14159/(2/FOV);
    vec3 fwd = normalize(ro); // o = ray origin vector
    vec3 rgt = normalize(vec3(fwd.z, 0, -fwd.x ));
    vec3 up = cross(fwd, rgt);
    // Unit direction ray.
    vec3 rd = normalize(fwd + FOV*(uv.x*rgt + uv.y*up));
    rd.yz = _rotate(rd.yz, lookXY.y*PI) ;
    rd.xy = _rotate(rd.xy, lookXY.x*PI);




    //rd.xy *= r2(ro.z*zRot);
    /*
    rd.xy *= r2(TIME/4.);
    rd.xz *= r2(TIME/4.);
    */
    // Edge and curvature variables. Passed to the normal functions... The refraction
    // pass has seperate normal function.
    float edge = 0., crv = 1.;

    
    
    // FIRST PASS.
    //
    float t = trace(ro, rd); // Trace.

    // Save the object IDs after the first pass.
    svObjID = objID;
    float oSvObjID = svObjID;
    
    // Advancing the ray origin, "ro," to the new hit point.
    vec3 sp = ro + rd*t;
    
    // Retrieving the normal at the hit point, plus the edge and curvature values.
    vec3 sn = getNormal(sp, edge, crv, t);

    
    // Fresnel. Handy for all kinds of aesthetic purposes. Not used here.
    //float fr = clamp(1. + dot(rd, sn), 0., 1.);
    
    // Shading. Shadows, ambient occlusion, etc. We're only performing this on the 
    // first pass. Not accurate, but faster, and in most cases, not that noticeable.
    // In fact, the shadows almost didn't make the cut, but it didn't quite feel 
    // right without them.
    float sh = softShadow(sp + sn*.002, lp, 16., t); // Set to "1.," if you can do without them.
    float ao = calculateAO(sp, sn);
    sh = (sh + ao*.3)*ao;
    
    // Retrieving the color at the initial hit point.
    vec3 sceneColor = doColor(sp, rd, sn, lp, edge, crv, ao, t);
    

    // Fog - based off of distance from the camera. This will be applied at the end.
    float fog = smoothstep(0., .95, t/FAR);
    
    
    // Cheap and nasty alternative, if not using the reflection pass.
    //sceneColor += eMap(reflect(rd, sn)/2., sn);
   
   
    // SECOND PASS
    
    // Reflected and refracted rays.
    vec3 refl = reflect(rd, sn); // Standard reflection.
    //vec3 refr = refract(rd, sn, 1./1.33); // Water refraction. Note the inverted index.
    
    // We're branching off from the same spot in two directions, so we'll use this so as
    // not to interfere with the original surface point vector, "sp." It was a style
    // choice on my part, but there are other ways.
    vec3 refSp; 
    
    // REFLECTED PASS
    //
    // Standard reflected ray, which is just a reflection of the unit
    // direction ray off of the intersected surface. You use the normal
    // at the surface point to do that. Hopefully, it's common sense.


    // The ray is edged off the surface, as required, but note that it has to be enough
    // to avoid conflict with the break condition in the "reflected" trace algorithm.
    t = traceRef(sp + refl*0.005*(t*.125 + 1.), refl);

    // Save the object IDs after the second pass.
    svObjID = objID;
    oSvObjID = svObjID; 
    
    
    // Advancing the ray from the new origin, "sp," to the new reflected hit point.
    refSp = sp + refl*t;
    
    // Retrieving the normal at the reflected hit point.
    sn = getNormal(refSp, edge, crv, t);
    //sn = getNormal(refSp);//*rDir;
    //edge = 0.;
 
    // Color at the reflected hit point.
    vec3 reflColor = doColor(refSp, refl, sn, lp, edge, crv, 1., t);
    sceneColor = sceneColor + reflColor*.5 ;
    //sceneColor = sceneColor*.5 + mix(reflColor, sceneColor, fr*fr*.66 + .34);
    

    
    // APPLYING SHADOWS
    //
    // Multiply the shadow from the first pass by the final scene color. Ideally, you'd check to
    // see if the reflected point was in shadow, and incorporate that too, but we're cheating to
    // save cycles and skipping it. It's not really noticeable anyway.
    sceneColor *= sh;
    
    
    // APPLYING FOG
    // Blend in a bit of light fog for atmospheric effect. I really wanted to put a colorful, 
    // gradient blend here, but my mind wasn't buying it, so dull, blueish grey it is. :)
    vec3 fogCol = vec3(.7, .8, 1.)*(rd.y*.5 + .5)*2.5;
    sceneColor = mix(sceneColor, fogCol, fog); // exp(-.002*t*t), etc. fog.zxy //pow(fogCol, vec3(1.33))*1.66
    
    
    // POSTPROCESSING
    // Interesting red to blueish mix.
    sceneColor = mix(sceneColor, pow(min(vec3(1.5, 1, 1)*sceneColor, 1.), vec3(1, 1.5, 8.)), uv.y);
    //sceneColor = pow(max(sceneColor, 0.), vec3(1.33))*1.66; // Adding a bit of contrast.
    //sceneColor *= vec3(1.2, 1, .8);
    
    //vec2 u2 = uv*r2(3.14159/6.);
    //float overlay = 1. + .4*sin(u2.x*3.14159*RENDERSIZE.y/1.5);
    //overlay *= 1. + .4*sin(u2.y*3.14159*RENDERSIZE.y/1.5); 
    //sceneColor *= overlay;
    
    // Subtle vignette.
    uv = fragCoord/RENDERSIZE.xy;
    sceneColor *= pow(16.*uv.x*uv.y*(1. - uv.x)*(1. - uv.y) , .125)*.5 + .5;
    // Colored varation.
    //sceneColor = mix(pow(min(vec3(1.5, 1, 1)*sceneColor, 1.), vec3(1, 2.5, 12.)).zyx, sceneColor, 
                    // pow( 16.0*uv.x*uv.y*(1.0-uv.x)*(1.0-uv.y) , .125)*.5 + .5);
    


    // Clamping the scene color, then presenting to the screen.
	fragColor = vec4(sqrt(clamp(sceneColor, 0.0, 1.0)), 1.0);
	return fragColor; 
 } 




vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}