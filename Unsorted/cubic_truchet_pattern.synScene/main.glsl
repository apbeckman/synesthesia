//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick);


/*

	Cubic Truchet Pattern
	---------------------

	I have a more elaborate 3D Truchet related example to post, but wanted to put up a
	simpler version first. There are a few cubic Truchet examples on here, including one
	of my own, so I'm not bringing anything particularly new to the table. :)

	I like 3D cubic Truchet patterns, because they're geometrically interesting, and are
	reasonably easy to produce. A standard 3D Truchet tile consists of three toroids
	centered on three edges of a cube, positioned to enter and exit six cube faces... Look
	one up on the internet, and that diatribe should make more sense. :) The idea is to
	interconnect the tiles in a 3D grid - by randomly rotating each around one of the axes -
	to produce an interesting spaghetti looking pattern.

	Constructing the individual tiles is as simple as breaking space into a cubic grid then
	positioning three torii in each cell. If you can position, rotate and render a torus,
	then it should be rudimentary.

	Detailed scenes look cool, but sometimes it can be difficult to separate the main point
	of the code from the details themselves. The Truchet code requires but a few lines.
	Adding a few decorations complicates things slightly, and trying to do it in a more
	efficient way can confuse things even further. I've tried my best to mitigate this.
    However, you'll see some esoteric object ID code scattered throughout the example, which
	consists of last minute logic that I put together to get the job done... You can safely
	ignore it. :)

	I code most things on a fast computer these days, which doesn't exactly help when it
	comes to optimization. As a rough guage, I'll try to get it running as smoothly as
	possible in fullscreen. This one seems to run OK in fullscreen, but it could do with
	a few extra tweaks. Aesthetically speaking, I had the 800 by 450 canvas in mind.

	I'm not really sure what kind of look I was going for, but I wanted something simple
	and clean looking that clearly showed the Truchet pattern. The lighting is very basic,
	but reasonably effective.

	Anyway, I'll put up my more interesting example next... provided I don't get
	sidetracked. :)


	Other cubic Truchet examples:

	// The first cubic Truchet example I saw on here. Really cool, and shiny. :)
	Truchet Tentacles - WAHa_06x36
	https://www.shadertoy.com/view/ldfGWn

	// Dr2's interpretation.
	Truchet Flythrough 2 - Dr2
	https://www.shadertoy.com/view/4dsBWf

	// Another one I did a while back. Less detailed, so probably easier to understand.
	Twisted Tubes - Shane
	https://www.shadertoy.com/view/lsc3DH

*/

// Maximum ray distance.
#define FAR Fog * Fog

// Scene object ID: Main tube, colored inner tube or colored band.
float objID; // Global ID to keep a copy of the above from pass to pass.

// Storage vectors to help determine the above and produce the blinking light effect. They're
// used inside the distance function. I'm not fond of using globals inside distance field
// functions, but felt it was necessary in this case.
vec3 vObjID;
float gID;


// Standard 2D rotation formula.
mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }

// Tube: Cylindrical tube, square tube, etc. In this case, it's a squarish tube with some
// beveled sides.
float tube(vec2 p){

    // Normally needed, but in this example "p" is passed to the function in its absoluted form.
    //p = abs(p);

    return max(max(p.x, p.y), (p.x + p.y)*.5773+Slim); // .7071 for an octagon, etc.
    //return max(p.x, p.y); // Square tube, etc.
}



// The toroidal tube objects. Each consist of a white squarish outer tube, a similar colored
// inner one (only visible through the holes) and some colored bands.
vec4 torTube(vec3 p){


    // Tube width.
    float rad2 = .065 + (Inner_Width * Inner_Width); // MS: Width of inner tube


    // Main tube. If it were not for the additional tube decorations, the following
    // would be all that'd be required.
    //
    // Note that we're converting one of the coordinates to its circular form. That way,
    // we're rendering a circular tube, instead of a straight one. It's an oversimplification,
    // but that's basically all a torus is. By the way, it doesn't have to be circular,
    // converting "p.xy" to an octagonal form, etc, would work also.
    float tb = tube(abs(vec2(length(p.xy) - .5, p.z))) - rad2;



    // Adding some details to the tube.
    //

    // Inner tube for colored lights.
    float innerTb = tb + .015 + (Shed_Outer * Shed_Outer * -1); // MS: Might want to mess w/ this at some point for a totally different scene


    // Tube segments - for the bands and holes.
    //
    // Number of tube segments. Breaking a circle into 8 lots of 3. Other combinations can
    // work though.
    float aNum = Texture_Detail + 24.;

    // To place things evenly around the tube, you need to obtain the angle subtended to the center,
    // partition it into the required number of cells (aNum), then obtain the angle at the center.
    float a = atan(p.y, p.x);
    float ia = floor(a/6.283*aNum) + .5; // .5 to move to the cell center. MS: This controls the windows and bands on tubes

    // Converting to polar coordinates - In effect: Radial position, "p.x," and angular position, "p.y."
    p.xy = rot2(ia*6.283/aNum)*p.xy;
    // The radial coordinate effective starts off at the center, so to spread the objects out, you have
    // to advance them  along the radial coordinate in the radial direction. In this case, we want the
    // objects to have the same radius as the torus we're attaching them to, which is ".5."
    p.x -= .5;

    // Drawing the objects within each of the partitioned cells. In this case, we're rendering some
    // colored sleeves (or bands), and boring out some holes.

    p = abs(p);

    float band = 1e5;

    // Group the 24 cell partitions into groups of 3 - in order to cover every third cell with the
    // colored band and bore some holes out in the others... I figured it'd break up the monotony. :)
    // On a side note, I try to avoid "if" statements inside distance functions when I can, but I
    // figured this would be the best way for readability. Although, I might rework it later.
    if(mod(ia + 1., 3.)>2.){

        band = max(tube(p.xz) - rad2 - .01, p.y - .04);
    		band = max(band, min(band + .005, -p.y + .015));
    }
    else {

        // Cute trick to break the cell into four - in order to bore out four holes in each cell.
        // Comment it out to produce just one hole.
        if (Just_One_Hole == 0) { p = abs(p - .02); }

        // Cut out two cross sections from the main tube.
        tb = max(tb, -min(tube(p.xy) - rad2 + .055, tube(p.yz) - rad2 + .055));

    }


    // Return the tube, bands, and inner tube objects.
		vec4 returnVec = vec4(tb, band, innerTb, ia);
		if (Remove_Texture == 1) { returnVec = vec4(1., band, innerTb, ia); }
		return returnVec;
}


/*

	The Truchet pattern:

	A standard 3D Truchet tile consists of three toroids centered on three edges of a cube,
    positioned to enter and exit six cube faces... Look one up on the internet, and that
	diatribe will make more sense. :) The idea is to connect the tiles in a 3D grid, then
	randomly rotate each around one of the axes to produce an interesting spaghetti looking
	pattern.

	Constructing the individual tiles is as simple as breaking space into a cubic grid then
	positioning three torii in each cell. If you can position, rotate and render a torus,
	then it should be rudimentary.

	On a side note, if you're one of those people who have trouble with the torus concept,
    you're basically rendering a straight tube that has had one of its coordinates warped
    into a circle first:

	float torus(vec3 p, vec3 center){

	    // Position the torus.
        p -= center;
		// Warp - Comment out below, and you're left with a straight tube.
    	vec2 q = vec2(length(p.xy) - .5, p.z); // vec2(length(p.xz) - .5, p.y), etc.
    	// Render a circular tube.
		float dist = length(q) - rad;
	}

*/


// I can thank Mattz for reminding me of this. You don't need to call all three decorated tubes,
// then determine the minimum. You can determine the minimum main tube, then call the function
// for the tube containing the more elaborate detailing that corresponds to it. And by that I
// mean return the unique oriented point that corresponds to the nearest tube segment distance.
//
vec4 torTubeTest(vec3 p){

    vec2 v = vec2(length(p.xy) - .5, p.z);

    // Main tube distance squared. Note: If a + c < b + c, then a*a<b*b.
    // Ie: we don't need to test length(v) - r, just dot(v, v);
    return vec4(p, dot(v, v));
}


float map(vec3 p)
{

    // Random ID for each grid cube.
    float rnd = fract(sin(dot(floor(p + vec3(111, 73, 27)), vec3(7.63, 157.31, 113.97)))*43758.5453);

    // Partition space into a grid of unit cubes - centered at the origin and ranging from
    // vec3(-.5, -.5, -.5) to vec3(-.5, -.5, -.5).
    p = fract(p) - .5;

    // Use each cube's random ID to rotate it in such a way that another one of its faces is
    // facing forward. In case you're not aware, the swizzling below is a cheap trick used to
    // achieve this. By the way, there may be a faster way to write the conditionals - using
    // ternary operators, or something to that effect, but I'm leaving it this way for now...
    // However, if a GPU expert feels that it's unnecessarily slow, then feel free to let me
    // know, and I'll change it.
    if(rnd>.833) p = p.xzy;
    else if(rnd>.666) p = p.yxz;
    else if(rnd>.5) p = p.yzx;
    else if(rnd>.333) p = p.zxy;
    else if(rnd>.166) p = p.zyx;

    // I can thank Mattz for reminding me of this step. Each Truchet tile contains three decorated
    // tubes. However, you only need to find the closest tube, "not" the closest decorated tube, which
    // requires a lot more GPU power. Each of these return the closest point and the distance...
    // Actually, the squared distance, which for comparisson purposes, is the same thing.
    vec4 tb1 = torTubeTest(vec3(p.xy + .5, p.z));
    vec4 tb2 = torTubeTest(vec3(p.yz - .5, p.x));
    vec4 tb3 = torTubeTest(vec3(p.xz - vec2(.5, -.5), p.y));

    // Sort the distances, then return the closest point.
    p = tb1.w<tb2.w && tb1.w<tb3.w ? tb1.xyz : tb2.w<tb3.w ? tb2.xyz : tb3.xyz;

    // Render the randomly aligned Truchet block. Ie, the three torii - plus bells and whistles.
    // Each quarter torus consists of three separate objects: A white tube with some holes in it,
    // some bracing (the colored sleeve looking things) and a colored inner tube. That's nine
    // objects returned in all. If it were not for the need to sort objects and attain a segment
    // identifier (tb.w), only a float would be necessary.
    vec4 tb = torTube(p);


    /// A unique angular segment identifier - used to produce the blinking lights.
    gID = tb.w;


    // Each torus segment contains three individual objects. Here, we're finding the minimum in
    // each category. We're keeping a global copy here that will be sorted for object identification
    // outside the raymarching loop. The reason this step is necessary is because the line below
    // finds the closest object, but doesn't tell us which object that is. That requires sorting,
    // which is best done outside the loop, for speed reasons.
    vObjID = tb.xyz - Tube_Width - syn_Hits*0.005*Tube_Gain;

    // Finding the minimum of the above to determine the overall minimum object in the scene.
    return min(min(vObjID.x, vObjID.y), vObjID.z);


}



// Recreating part of the distance function to obtain the segment IDs, which in turn is used
// to create the blink effect.
// MS: Removed the randomizer, instead the lights are just mapped to bass hits
float lightBlink(vec3 p, float gID){

    // Unique identifier for the cubic grid cell.
    // float rnd = fract(sin(dot(floor(p + vec3(111, 73, 27)), vec3(7.63, 157.31, 113.97)))*43758.5453);

    // Reusing "rnd" to produce a new random number, then using that
    // random number to create lights that blink at random intervals.
    //rnd = fract(rnd + gID*43758.54571);

    // Blink at random.
    // return smoothstep(0.33, .66, sin(rnd*6.283 + TIME*3.)*.5 + .5);
		//return smoothstep(0.33, .66, sin(syn_OnBeat)*.5 + .5);
		return syn_BassHits;


}

//float glow;
// Standard raymarching algorithm.
float trace(vec3 o, vec3 r){

    //glow = 0.;

    // Total ray distance travelled, and nearest distance at the current ray position.
    float t = 0., d;

		float numSteps = Globular * Globular;

    for (int i = 0; i<numSteps; i++) {

        // Surface distance.
        d = map(o + r*t) * 0.5 * Flip;

        //if(abs(d)<.05) glow += (.05 - abs(d))/(1. + d*d);
        // If the ray position is within the surface threshold ("abs" means either side of the
        // surface), or if we've traversed beyond the maximum, exit the loop.
        if(abs(d)<.001*(t*.125 + 1.) || t>FAR) break;


        // Standard jump.
        t += d + Remove_Surface;

        // Shortening the ray jump right near the camera to alleviated near-camera artifacts.
        //t += t<.125 ? d*.7 : d;
    }

    // Clamp the total distance to "FAR." It can sometimes get rid of far surface artifacts.
    return min(t, FAR);
}

// Cheap shadows are the bain of my raymarching existence, since trying to alleviate artifacts is an excercise in
// futility. In fact, I'd almost say, shadowing - in a setting like this - with limited  iterations is impossible...
// However, I'd be very grateful if someone could prove me wrong. :)
float shadow(vec3 ro, vec3 lp, float k, float t){

    // More would be nicer. More is always nicer, but not really affordable... Not on my slow test machine, anyway.
    const int maxIterationsShad = 32;

    vec3 rd = lp-ro + syn_MidHighLevel; // Unnormalized direction ray. MS: Adding some audioreactivity

    float shade = 1.;
    float dist = .001*(t*.125 + 1.);  // Coincides with the hit condition in the "trace" function.
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

    // I sometimes add a constant to the final shade value, which lightens the shadow a bit. It's a preference
    // thing. Really dark shadows look too brutal to me. Sometimes, I'll also add AO, just for kicks. :)
    return min(max(shade, 0.) + .0 + Brightness, 1.);
}

// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)
// Anyway, I like this one. I'm assuming it's based on IQ's original.
float cAO(in vec3 pos, in vec3 nor)
{
	float sca = 1., occ = 0.0;
    for( int i=0; i<5; i++ ){

        float hr = 0.01 + float(i)*0.35/4.0;
        float dd = map(nor * hr + pos);
        occ += (hr - dd)*sca;
        sca *= 0.7;
    }
    return clamp(1.0 - occ, 0.0, 1.0 );
}

// Normal calculation, with some edging and curvature bundled in.
// MS: I think this is all just lighting
vec3 nrm(vec3 p, inout float edge, inout float crv, float t) {

    // It's worth looking into using a fixed epsilon versus using an epsilon value that
    // varies with resolution. Each affects the look in different ways. Here, I'm using
    // a mixture. I want the lines to be thicker at larger resolutions, but not too thick.
    // As for accounting for PPI; There's not a lot I can do about that.
    vec2 e = vec2(1./mix(400., RENDERSIZE.y, .5)*(1. + t*.5), 0);

	float d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	float d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	float d5 = map(p + e.yyx), d6 = map(p - e.yyx);
	float d = map(p)*2.;

    edge = abs(d1 + d2 - d) + abs(d3 + d4 - d) + abs(d5 + d6 - d);
    //edge = abs(d1 + d2 + d3 + d4 + d5 + d6 - d*3.);
    edge = smoothstep(0., 1., sqrt(edge/e.x*2.));
/*
    // Wider sample spread for the curvature.
    e = vec2(12./450., 0);
	d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	d5 = map(p + e.yyx), d6 = map(p - e.yyx);
    crv = clamp((d1 + d2 + d3 + d4 + d5 + d6 - d*3.)*32. + .5, 0., 1.);
*/

    e = vec2(.002, 0); //RENDERSIZE.y - Depending how you want different resolutions to look.
	d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	d5 = map(p + e.yyx), d6 = map(p - e.yyx);

    return normalize(vec3(d1 - d2, d3 - d4, d5 - d6));
}


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // Aspect correct screen coordinates.
	vec2 uv = (fragCoord.xy - RENDERSIZE.xy*.5)/RENDERSIZE.y;
    
    // Ray origin, or camera - Moving along the Z-axis.
    vec3 o = vec3(0., 0., fly_time);
    // Light. Situated near the camera whilst moving along with it.
	vec3 lp = vec3(-1, 3, -.25) + o;

    // Unit ray vector.
		// MS: Adding toggle for bulbous scene warp
		float isWarp = Fish_Eye_Lens;
		vec3 r = normalize(vec3(uv, 1));
		if (isWarp == 1) { r = (vec3(uv, 1.1)); }
    //r = normalize(vec3(r.xy, r.z - length(r.xy)*.1));



    // Rotating "r" back and forth along various axes for some cheap camera movement.
		// MS: Removing this, since we want to have control over the camera with XY controller rather than leaving it up to a rotation based on time
    //r.xz *= rot2(sin(TIME/2.) * 0.4);
    //r.xy *= rot2(cos(TIME/2.) * 0.2);
        r.xz = _rotate(r.xz, LookXY.x*PI);
        r.yz = _rotate(r.yz, LookXY.y*PI);
		
        // MS: Adding camera movement here
		//r.xy += LookXY.xy;
		//o.xy += LookXY.xy*Cam_Toggle; // MS: Just testing to see what happens if we map o as well

    // Trace out the scene.
    float t = trace(o, r);

    // Determining the object ID. Sorting the three different objects outside the loop
    // is a little less readable, but usually faster. See the distance function.
    objID = (vObjID.x<vObjID.y && vObjID.x<vObjID.z) ? 0. : (vObjID.y<vObjID.z) ? 1. : 2.;

    // Segment ID: Sorting the segments to determine the unique ID. This ID is fed
    // into a function to give the blinking light effect.
    float svGID = gID;


	// Initiate the scene color to zero.
    vec3 sc = vec3(0);

		// MS: Check for Synesthesia media
		if (_exists(syn_UserImage)) { sc = _loadUserImage().rgb; }

    // An object in the scene has been hit, so light it.
    if(t<FAR){

        // Hit position.
        vec3 sp = Cartoonize + o + r*t;

        // Normal, plus edges and curvature. The latter isn't used.
        float edge = 0., crv = 1.;
        vec3 sn = nrm(sp, edge, crv, t);

        // Producing a gradient color based on position. Made up on the spot.
        vec3 oCol = mix(vec3(1, .1, .3), vec3(1, .5, .1), dot(sin(sp*8. - cos(sp.yzx*4.)), vec3(.166)) + .5);
        oCol = mix(oCol, oCol.yzx, smoothstep(.3, 1., dot(sin(sp*4. + cos(sp.zxy*4. + TIME)), vec3(.166*.6)) + .3));

        // Color the individual objects, based on object ID.
        if(objID<.5)oCol = mix(oCol, vec3(Red, Green, Blue), 0.97); // The whitish tube.
        else if(objID<1.5) oCol = mix(oCol, vec3(1), .05); // The colorful bands.
        else {

            oCol = mix(oCol, vec3(1), .05); // Inner tube color. Same as above, but you could change it.
          	//oCol = mix(oCol, oCol.zyx, dot(cos(sp*32. + sin(sp.yzx*16.)), vec3(.166)) + .5);

            // The blinking light effect. In effect, the number varies color intensity is periodically
            // ramped right up. The individual segment ID is responsible for the randomness.
            oCol *= lightBlink(sp, svGID)*7.5 + .5;
        }


        // Ambient occlusion and shadows.
        float ao = cAO(sp, sn);
        float sh = shadow(sp + sn*.002, lp, 16., t);


        // Point light direction vector.
        vec3 ld = lp - sp;
        float dist = max(length(ld), 0.001); // Distance.
        ld /= dist; // Using the distance to nomalize the point light direction vector.


        // Attenuation - based on light to surface distance.
        float atten = 3.5/(1. + dist*0.05 + dist*dist*0.05);

        // Diffuse light.
        float diff = max(dot(ld, sn), 0.);



        // Combining the above terms to produce the final color.
    	sc = oCol*(diff + ao*.35);


        // Fake caustic lighting... It didn't sit right with the scene, so it didn't make the cut. :)
        //sc += .02/max(abs(.05 - map(sp*1.5)), .01)*oCol*vec3(1, .7, .5);
        //sc += oCol*abs(tan(t*1.5 + TIME/2.))*vec3(.1, .2, 1)*.05;

        // Applying the dark edges, attenuation, shadows and ambient occlusion.
        if (Remove_Shadows == 0) { sc *= (1. - edge*.7)*atten*(sh + ao*.3)*ao; }


    }


    // Applying some basic camera distance fog. Not to be confused with the light
    // to surface attenuation.
    float fog = 1./(1. + t*.125 + t*t*.05);
    sc = mix(vec3(0), sc, fog);
    //sc = mix(sc, vec3(0), smoothstep(0.0, .2, t/FAR));


    // Subtle vignette.
    uv = fragCoord/RENDERSIZE.xy;
    //sc *= pow(16.*uv.x*uv.y*(1. - uv.x)*(1. - uv.y) , .125);
    // Colored variation.
    sc = mix(pow(min(vec3(1.5, 1, 1).zyx*sc, 1.), vec3(1, 3, 16).zyx), sc,
             pow(16.*uv.x*uv.y*(1. - uv.x)*(1. - uv.y) , .125)*.75 + .25);


	fragColor = vec4(sqrt(max(sc, 0.)), 1);
	return fragColor;
 }


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}
