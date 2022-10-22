/*

	Dry Rocky Gorge
	---------------

	This is a very simple terrain example, and not much different to the many others on here.
	Overall, it's not that inspiring - After all, there are only so many ways you can render a
	terrain scene; Terrain, sky, fake normal-based ground coloring. Not much effort went into
	it. In fact, I probably spent more time choosing the ground texture and sky color than
	creating the terrain. :)

    I've been playing around with low-poly terrain lately, but believe it or not, using a modern
	machine to emulate the way things looked on old machines isn't always as easy as you'd think. :)
	Therefore, I took a break and coded up a very cliche medium-level terrain fly-though.

	The lighting is fudged in a lot of places - especially where the clouds are concerned, so I
	wouldn't pay too much attention to it. In fact, you could ignore most of the code and just
	look at the camera setup and distance function.

	The geometry is very basic. Render a plane, carve out a squarish tube, wrap it around the
	camera path, then add some noise layers to the result. Not much to it. For the terrain itself,
	I'd hoped to make use of IQ's gradient noise derivatives code, but speed was an issue, so I let
	it go. Maybe next time. :)

	There wouldn't be much code here if it were not for the cloud layering routine that I dropped
    in. It's only used for four layers, which meant the aesthetic returns were rather diminished,
	so it was hardly worth the effort.

	Anyway, I have a lot of more interesting examples than this that I hope to add at some point.


*/

#define FAR 80. // Maximum ray distance. Analogous to the far plane.


// Fabrice's concise, 2D rotation formula.
//mat2 r2(float th){ vec2 a = sin(vec2(1.5707963, 0) + th); return mat2(a, -a.y, a.x); }
// Standard 2D rotation formula - See Nimitz's comment.
mat2 r2(in float a){ float c = cos(a), s = sin(a); return mat2(c, s, -s, c); }

// vec3 to float hash.
float hash31( vec3 p ){ return fract(sin(dot(p, vec3(157, 113, 7)))*45758.5453); }

// vec3 to float hash.
float hash21( vec2 p ){ return fract(sin(dot(p, vec2(41, 289)))*45758.5453); }

// Non-standard vec3-to-vec3 hash function.
vec3 hash33(vec3 p){

    float n = sin(dot(p, vec3(7, 157, 113)));
    return fract(vec3(2097152, 262144, 32768)*n);
}

// Compact, self-contained version of IQ's 3D value noise function. I put this together, so be
// careful how much you trust it. :D
float n3D(vec3 p){

	const vec3 s = vec3(7, 157, 113);
	vec3 ip = floor(p); p -= ip;
    vec4 h = vec4(0., s.yz, s.y + s.z) + dot(ip, s);
    p = p*p*(3. - 2.*p); //p *= p*p*(p*(p * 6. - 15.) + 10.);
    h = mix(fract(sin(h)*43758.5453), fract(sin(h + s.x)*43758.5453), p.x);
    h.xy = mix(h.xz, h.yw, p.y);
    return mix(h.x, h.y, p.z); // Range: [0, 1].
}


// Cheap and nasty 2D smooth noise function, based on IQ's original. Very trimmed down. In fact,
// I probably went a little overboard. I think it might also degrade with large time values. I'll
// swap it for something more robust later.
float n2D(vec2 p) {

	vec2 i = floor(p); p -= i; p *= p*(3. - p*2.); //p *= p*p*(p*(p*6. - 15.) + 10.);

	return dot(mat2(fract(sin(vec4(0, 41, 289, 330) + dot(i, vec2(41, 289)))*43758.5453))*
               vec2(1. - p.y, p.y), vec2(1. - p.x, p.x));

}



// Tri-Planar blending function. Based on an old Nvidia tutorial.
vec3 tex3D( sampler2D t, in vec3 p, in vec3 n ){

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

// The path is a 2D sinusoid that varies over time, depending upon the frequencies, and amplitudes.
vec2 path(in float z){

    vec2 straightPathVec = vec2(0);
  //  if (Straight_Path < 0.5) {
      vec2 curvedPathVec = vec2(sin(z*.075)*8., cos(z*.1)*.75); // Windy path.
      float morphTime = Morph_Amount*sin(syn_BassTime*Morph_Speed);
      curvedPathVec += morphTime;
    //}

    return mix(curvedPathVec, straightPathVec, Straight_Path);
    //return pathVec;

}


// The triangle function that Shadertoy user Nimitz has used in various triangle noise demonstrations.
// See Xyptonjtroz - Very cool. Anyway, it's not really being used to its full potential here.
//vec2 tri(in vec2 x){return abs(x-floor(x)-.5);} // Triangle function.

// Smooth maximum, based on IQ's smooth minimum.
float smax(float a, float b, float s){

    float h = clamp(.5 + .5*(a - b)/s, 0., 1.);
    return mix(b, a, h) + h*(1. - h)*s;
}


/*
// IQ's smooth minium function.
float smin(float a, float b , float s){

    float h = clamp( 0.5 + 0.5*(b-a)/s, 0. , 1.);
    return mix(b, a, h) - h*(1.0-h)*s;
}
*/

// Basic terrain - in the form of noise layering. How it's approached depends on the desired outcome, but
// in general you add a layer of noise, vary the amplitude and frequency, then add another layer, and so
// on. Other additions include skewing between layers (which is done here) and attenuating the amplitudes
// with the noise derivatives (For cost reasons, I left that out). Of course, you can do this with functions
// other than noise, like Voronoi, sinusoidal variations, triangle noise, etc.
float terrain(vec2 p){

    p /= 8.; // Choosing a suitable starting frequency.

    p *= Frequency;

    // Edging the terrain surfacing into a position I liked more. Not really necessary though.
    p += .5;

    // Amplitude, amplitude total, and result variables.
    float a = 1., sum = 0. + Amplitude + clamp(Amplitude_Time*sin(syn_Time*.25),0,3), res = 0.;

    // Only five layers. More layers would be nicer, but cycles need to be taken into
    // consideration. A simple way to give the impression that more layers are being added
    // is to increase the frequency by a larger amount from layer to layer.
    for (int i=0; i<5; i++){

        res += n2D(p)*a; // Add the noise value for this layer - multiplied by the amplitude.
        //res += abs(n2D3(p) - .5)*a; // Interesting variation.
        //res += n2D3(p)*abs(a)*.8; // Another one.

        // Scaling the position and doing some skewing at the same time. The skewing isn't
        // mandatory, but it tends to give more varied - and therefore - interesting results.
        // IQ uses this combination a bit, so I'll assume he came up with the figures. I've
        // tried other figures, but I tend to like these ones as well.
        p = mat2(1, -.75, .75, 1)*p*2.72;
        //p *= 3.2; // No skewing. Cheaper, but less interesting.

        sum += a; // I reasoned that the sum will always be positive.

        // Tempering the amplitude. Note the negative sign - a less common variation - which
        // was thrown in just to mix things up.
        a *= -.5/1.7;
    }


    return res/sum; // Return the noisy terrain value.

}



// The gorge terrain setup up: It's just a flat plane with a channel cut out of it, which is
// wrapped around the camera path. Then, a few layers of simple 2D noise is added to it.
float map(vec3 p){


    // The noise layers.
    float trSf = terrain(p.xz);

    p.xy -= path(p.z); // Wrap the gorge around the path.
    p.y += sin(p.z + syn_Time*Wave_Speed)*Wave_Height; // Waaavy
    // The canyon - or gorge, which consists of a mixed circle and square shape, extruded along
    // the path. It's been stretched, lowered, then subtracted from the flat plane, before adding
    /// the noise layers.
    vec2 ca = abs(p.xy*vec2(1, .7) + vec2(0, -2.75)); // Stretch and lower space.

    // Smoothly carve out the gorge from the plane, then add the noise to the result.
    float n = smax(6. - mix(length(ca), max(ca.x, ca.y), .25), p.y - 1.75, 2.) + (.5 - trSf)*4.;


    return n*.7; // Return the minimum hit point.

}




// Texture bump mapping. Four tri-planar lookups, or 12 texture lookups in total. I tried to
// make it as concise as possible. Whether that translates to speed, or not, I couldn't say.
vec3 texBump( sampler2D tx, in vec3 p, in vec3 n, float bf){

    const vec2 e = vec2(.001, 0);

    // Three gradient vectors rolled into a matrix, constructed with offset greyscale texture values.
    mat3 m = mat3( tex3D(tx, p - e.xyy, n), tex3D(tx, p - e.yxy, n), tex3D(tx, p - e.yyx, n));

    vec3 g = vec3(.299, .587, .114)*m; // Converting to greyscale.
    g = (g - dot(tex3D(tx,  p , n), vec3(.299, .587, .114)) )/e.x;

    // Adjusting the tangent vector so that it's perpendicular to the normal -- Thanks to
    // EvilRyu for reminding why we perform this step. It's been a while, but I vaguely recall
    // that it's some kind of orthogonal space fix using the Gram-Schmidt process. However,
    // all you need to know is that it works. :)
    g -= n*dot(n, g);

    return normalize( n + g*bf ); // Bumped normal. "bf" - bump factor.

}



// Standard raymarching routine.
float trace(vec3 ro, vec3 rd){

    float t = 0., d;

    int numSteps = int(160 - 145*pow(Ray_Steps,.08));
    for (int i=0; i<numSteps; i++){

        d = map(ro + rd*t);

        if(abs(d)<.001*(t*.125 + 1.) || t>FAR) break;

        t += d;
    }

    return min(t, FAR);
}


// Cheap shadows are the bain of my raymarching existence, since trying to alleviate artifacts is an excercise in
// futility. In fact, I'd almost say, shadowing - in a setting like this - with limited  iterations is impossible...
// However, I'd be very grateful if someone could prove me wrong. :)
float softShadow(vec3 ro, vec3 lp, float k, float t){

    // More would be nicer. More is always nicer, but not really affordable... Not on my slow test machine, anyway.
    const int maxIterationsShad = 48;

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
        dist += clamp(h, .02, .25);

        // Early exits from accumulative distance function calls tend to be a good thing.
        if (h<0. || dist > end) break;
    }

    // I've added a constant to the final shade value, which lightens the shadow a bit. It's a preference thing.
    // Really dark shadows look too brutal to me. Sometimes, I'll add AO also just for kicks. :)
    return min(max(shade, 0.) + .15, 1.);
}

/*
// Standard normal function. It's not as fast as the tetrahedral calculation, but more symmetrical. Due to
// the intricacies of this particular scene, it's kind of needed to reduce jagged effects.
vec3 getNormal(in vec3 p) {
	const vec2 e = vec2(.002, 0);
	return normalize(vec3(map(p + e.xyy) - map(p - e.xyy), map(p + e.yxy) - map(p - e.yxy),	map(p + e.yyx) - map(p - e.yyx)));
}
*/


// Tetrahedral normal, to save a couple of "map" calls. Courtesy of IQ.
vec3 getNormal( in vec3 p ){

    // Note the slightly increased sampling distance, to alleviate
    // artifacts due to hit point inaccuracies.
    vec2 e = vec2(0.002, -0.002);
    return normalize(e.xyy*map(p + e.xyy) + e.yyx*map(p + e.yyx) + e.yxy*map(p + e.yxy) + e.xxx*map(p + e.xxx));
}




// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)
// Anyway, I like this one. I'm assuming it's based on IQ's original.
float calcAO(in vec3 p, in vec3 nor){

	float sca = 1.5, occ = 0.;

    for(float i=0.; i<5.; i++){
        float hr = .01 + i*.5/4.;
        float dd = map(nor*hr + p);
        occ += (hr - dd)*sca;
        sca *= .7;
    }

    return clamp(1. - occ, 0., 1.);
}




// Distance function.
float fmap(vec3 p){

    // Three layers of noise. More would be nicer.
    p *= vec3(1, 4, 1)/400.;

    return n3D(p)*0.57 + n3D(p*4.)*0.28 + n3D(p*8.)*0.15;
}

// Used in one of my volumetric examples. With only four layers, it's kind of going to waste
// here. I might replace it with something more streamlined later.
vec4 cloudLayers(vec3 ro, vec3 rd, vec3 lp, float far){

    // The ray is effectively marching through discontinuous slices of noise, so at certain
    // angles, you can see the separation. A bit of randomization can mask that, to a degree.
    // At the end of the day, it's not a perfect process. Note, the ray is deliberately left
    // unnormalized... if that's a word.
    //
    // Randomizing the direction.
    rd = (rd + (hash33(rd.zyx)*0.004-0.002));
    // Randomizing the length also.
    rd *= (1. + fract(sin(dot(vec3(7, 157, 113), rd.zyx))*43758.5453)*0.04-0.02);

    // Some more randomization, to be used for color based jittering inside the loop.
    //vec3 rnd = hash33(rd+311.);

    // Local density, total density, and weighting factor.
    float ld=0., td=0., w=0.;

    // Closest surface distance, and total ray distance travelled.
    float d=1., t=0.;


    // Distance threshold. Higher numbers give thicker clouds, but fill up the screen too much.
    const float h = .5;


    // Initializing the scene color to black, and declaring the surface position vector.
    vec3 col = vec3(0), sp;

    vec4 d4 = vec4(1, 0, 0, 0);



    // Particle surface normal.
    //
    // Here's my hacky reasoning. I'd imagine you're going to hit the particle front on, so the normal
    // would just be the opposite of the unit direction ray. However particles are particles, so there'd
    // be some randomness attached... Yeah, I'm not buying it either. :)
    vec3 sn = normalize(hash33(rd.yxz)*.03-rd);

    // Raymarching loop.
    for (int i=0; i<4; i++) {

        // Loop break conditions. Seems to work, but let me
        // know if I've overlooked something.
        if(td>1. || t>far)break;


        sp = ro + rd*t; // Current ray position.
        d = fmap(sp); // Closest distance to the surface... particle.
        //d4 = fmap(sp); // Closest distance to the surface... particle.

        //d = d4.x;
        //sn = normalize(d4.yzw);

        // If we get within a certain distance, "h," of the surface, accumulate some surface values.
        // The "step" function is a branchless way to do an if statement, in case you're wondering.
        //
        // Values further away have less influence on the total. When you accumulate layers, you'll
        // usually need some kind of weighting algorithm based on some identifying factor - in this
        // case, it's distance. This is one of many ways to do it. In fact, you'll see variations on
        // the following lines all over the place.
        //
        ld = (h - d) * step(d, h);
        w = (1. - td) * ld;

        // Use the weighting factor to accumulate density. How you do this is up to you.
        //td += w*w*8. + 1./60.; //w*w*5. + 1./50.;
        td += w*.5 + 1./65.; // Looks cleaner, but a little washed out.


        // Point light calculations.
        vec3 ld = lp-sp; // Direction vector from the surface to the light position.
        float lDist = max(length(ld), 0.001); // Distance from the surface to the light.
        ld/=lDist; // Normalizing the directional light vector.

        // Using the light distance to perform some falloff.
        float atten = 100./(1. + lDist*0.005 + lDist*lDist*0.00005);

        // Ok, these don't entirely correlate with tracing through transparent particles,
        // but they add a little anglular based highlighting in order to fake proper lighting...
        // if that makes any sense. I wouldn't be surprised if the specular term isn't needed,
        // or could be taken outside the loop.
        float diff = max(dot( sn, ld ), 0.);
        float spec = pow(max( dot( reflect(-ld, sn), -rd ), 0. ), 4.);

        // Accumulating the color. Note that I'm only adding a scalar value, in this case,
        // but you can add color combinations.
        //col += w*(1. + diff*.5 + spec*.5)*atten;

        // Try this instead, to see what it looks like without the fake contrasting. Obviously,
        // much faster.
        col += w*(diff + vec3(1, .75, .5)*spec + .5)*atten;//*1.25;

        // Optional extra: Color-based jittering. Roughens up the grey clouds that hit the camera lens.
        //col += (fract(rnd*289. + t*41.)-.5)*0.02;;



        // Enforce minimum stepsize. This is probably the most important part of the procedure.
        // It reminds me a little of of the soft shadows routine.
        t += max(d4.x*.5, 0.25)*100.; //* 0.75
        // t += 0.2; // t += d*0.5;// These also work, but don't seem as efficient.

    }

    //t = min(t, FAR); //24.

    return vec4(col, t);

}

// Pretty standard way to make a sky.
vec3 getSky(in vec3 ro, in vec3 rd, vec3 lp, float t){


	float sun = max(dot(rd, normalize(lp - ro)), 0.0); // Sun strength.
	float horiz = pow(1.0-max(rd.y, 0.0), 3.)*.25; // Horizon strength.

	// The blueish sky color. Tinging the sky redish around the sun.
	vec3 col = mix(vec3(.25, .5, 1)*.8, vec3(.8, .75, .7), sun*.5);//.zyx;
    // Mixing in the sun color near the horizon.
	col = mix(col, vec3(1, .5, .25), horiz);

    //vec3 col = mix(vec3(1, .7, .55), vec3(.6, .5, .55), rd.y*.5 + .5);

    // Sun. I can thank IQ for this tidbit. Producing the sun with three
    // layers, rather than just the one. Much better.
	col += 0.25*vec3(1, .7, .4)*pow(sun, 5.0);
	col += 0.25*vec3(1, .8, .6)*pow(sun, 64.0);
	col += 0.15*vec3(1, .9, .7)*max(pow(sun, 512.0), .25);

    // Add a touch of speckle. For better or worse, I find it breaks the smooth gradient up a little.
    col = clamp(col + hash31(rd)*0.04 - 0.02, 0., 1.);

    //return col; // Clear sky day. Much easier. :)

	// Clouds. Render some 3D clouds far off in the distance. I've made them sparse and wispy,
    // since we're in the desert, and all that.

    // Mapping some 2D clouds to a plane to save some calculations. Raytrace to a plane above, which
    // is pretty simple, but it's good to have Dave's, IQ's, etc, code to refer to as backup.

    // Give the direction ray a bit of concavity for some fake global curvature - My own dodgy addition. :)
    //rd = normalize(vec3(rd.xy, sqrt(rd.z*rd.z + dot(rd.xy, rd.xy)*.1) ));

    // If we haven't hit anything and are above the horizon point (there for completeness), render the sky.

    // Raytrace to a plane above the scene.
    float tt = (1000. - ro.y)/(rd.y + .2);

    if(t>=FAR && tt>0.){

        // Trace out a very small number of layers. In fact, there are so few layer that the following
        // is almost pointless, but I've left it in.
        vec4 cl = cloudLayers(ro + rd*tt, rd, lp, FAR*3.);
        vec3 clouds = cl.xyz;

        // Mix in the clouds.
        col = mix( col, vec3(1), clouds); // *clamp(rd.y*4. + .0, 0., 1.)
    }

    return col;

}



// Coloring\texturing the scene objects, according to the object IDs.
vec3 getObjectColor(vec3 p, vec3 n){

    //p.xy -= path(p.z);

    // Object texture color.
    vec3 tx = tex3D(image45, p/8., n ); // Texture value. Pinkish limestone.

    // Hinting that there's some dry vegetation below. The flatter the surface (based on n.y), the greater
    // the chance that something's growing on it. Physical trees would be much nicer, and I'm working on that,
    // but for now, cheap trickery will have to suffice. :) By the way, take a look at IQ's "Rainforest"
    // example for an amazing looking compromise.
    vec3 gr = mix(vec3(1), vec3(.8, 1.3, .2), smoothstep(.5, 1., n.y));
    return mix(tx, tx*gr, smoothstep(.7, 1., (n.y)));

}

// Using the hit point, unit direction ray, etc, to color the scene. Diffuse, specular, falloff, etc.
// It's all pretty standard stuff.
vec3 doColor(in vec3 ro, in vec3 rd, in vec3 lp, float t){

    // Initiate the scene (for this pass) to zero.
    vec3 sceneCol = vec3(0);

    if(t<FAR){ // If we've hit a scene object, light it up.

            // Advancing the ray origin, "ro," to the new hit point.
        vec3 sp = ro + rd*t;

        // Retrieving the normal at the hit point.
        vec3 sn = getNormal(sp);

        vec3 tx = sp;
        //tx.xy -= path(tx.z);
        sn = texBump(image45, tx/2., sn, .15);

        // Shading. Shadows, ambient occlusion, etc.
        float sh = softShadow(sp + sn*.002, lp, 16., t); // Set to "1.," if you can do without them.
        float ao = calcAO(sp, sn);
        sh = (sh + ao*.25)*ao;


        vec3 ld = lp - sp; // Light direction vector.
        float lDist = max(length(ld), 0.001); // Light to surface distance.
        ld /= lDist; // Normalizing the light vector.

        // Attenuating the light, based on distance.
        float atten = 3./(1. + lDist*0.005 + lDist*lDist*0.00005);

        // Standard diffuse term.
        float diff = max(dot(sn, ld), 0.);
        //diff = pow(diff, 2.)*.66 + pow(diff, 4.)*.34;
        // Standard specualr term.
        float spec = pow(max( dot( reflect(-ld, sn), -rd ), 0.0 ), 64.0);
        //float fres = clamp(1. + dot(rd, sn), 0., 1.);
        //float Schlick = pow( 1. - max(dot(rd, normalize(rd + ld)), 0.), 5.0);
        //float fre2 = mix(.5, 1., Schlick);  //F0 = .5.

        // Coloring the object. You could set it to a single color, to
        // make things simpler, if you wanted.
        vec3 objCol = getObjectColor(sp, sn);//mix(sn, oSn, .75)

        // Combining the above terms to produce the final scene color.
        sceneCol = objCol*(diff + ao*.5 + vec3(1, .7, .5)*spec);

        // Apply the attenuation and shadows.
        sceneCol *= atten*sh;

    }


    // Return the color. Done once for each pass.
    return sceneCol;

}


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    // Screen coordinates.
	vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;

	// Camera Setup.
	vec3 ro = vec3(0, Cam_Height, fly_time); // Camera position, doubling as the ray origin.
	vec3 lk = ro + vec3(0, -.04 - Cam_Height*.0025, .25);  // "Look At" position.
  lk.xy += vec2(LookXY.x, LookXY.y*.25);

    // Light position. Set reasonably far away in the background somewhere. A sun is usually so far
    // away that direct light is called for, put I like to give it just a bit of a point light feel.
    vec3 lp = ro + vec3(8, FAR*.26, FAR*.52)*3.;
    //vec3 lp = ro + vec3(0., 0, 4);


	// Using the Z-value to perturb the XY-plane.
	// Sending the camera, "look at," and light vector down the path. The "path" function is
	// synchronized with the distance function.
    ro.xy += path(ro.z);
	lk.xy += path(lk.z);
	lp.xy += path(lp.z);


    // Using the above to produce the unit ray-direction vector.
    float FOV = 3.14159/FOV; // 3 FOV - Field of view.
    vec3 forward = normalize(lk-ro);
    vec3 right = normalize(vec3(forward.z, 0., -forward.x ));
    vec3 up = cross(forward, right);

    // rd - Ray direction.
    vec3 rd = normalize(forward + FOV*uv.x*right + FOV*uv.y*up);
    //rd = normalize(vec3(rd.xy, sqrt(max(rd.z*rd.z - dot(rd.xy, rd.xy)*.15, 0.)) ));

    // Camera swivel - based on path position.
    vec2 sw = path(lk.z);
    rd.xy *= r2(-sw.x/24.);
    rd.yz *= r2(-sw.y/16.);

    // Trace the scene.
    float t = trace(ro, rd);


    // Retrieve the background color.
    vec3 sky = getSky(ro, rd, lp, t);




    // Retrieving the color at the initial hit point.
    vec3 sceneColor = doColor(ro, rd, lp, t);


    // APPLYING FOG
    // Fog - based off of distance from the camera.
    float fog = smoothstep(0., .95, t/FAR); // t/FAR;

    // Blend in the sky. :)
    vec3 fogCol = sky;
    //mix(vec3(.6, .9, 1).zyx, vec3(.62, .68, 1).zyx, rd.y*.5 + .5);
    sceneColor = mix(sceneColor, fogCol, fog); // exp(-.002*t*t), etc. fog.zxy


    // POSTPROCESSING


    // Subtle vignette.
    uv = fragCoord/RENDERSIZE.xy;
    sceneColor *= pow(16.*uv.x*uv.y*(1. - uv.x)*(1. - uv.y) , .125)*.75 + .25;
    // Colored varation.
    //sceneColor = mix(pow(min(vec3(1.5, 1, 1)*sceneColor, 1.), vec3(1, 3, 16)), sceneColor,
                     //pow( 16.0*uv.x*uv.y*(1.0-uv.x)*(1.0-uv.y) , .125)*.5 + .5);

    // A very simple overlay. Two linear waves - rotated at 60 degree angles - to give a dot-matrix vibe.
    //uv = sin(uv*r2(3.14159/6.)*3.14159*RENDERSIZE.y/1.5)*.1 + 1.;
    //sceneColor *= uv.x*uv.y;


    // Clamping the scene color, then presenting to the screen.
	fragColor = vec4(sqrt(clamp(sceneColor, 0.0, 1.0)), 1.0);
	return fragColor;
 }


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}
