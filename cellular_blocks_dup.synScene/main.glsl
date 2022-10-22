vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick);





/*



	Fast Cellular Blocks

	--------------------



	After looking at IQ and Aiekick's 2D Voronoi triangle-metric examples (the ones that

	look like blocks), I thought it'd be interesting to raymarch the pattern onto a back

	plane to see just how 3D it looked. The result was an interesting voxelized look with

	more depth than the regular, or bump mapped, examples you see around. After that, I

	figured I'd add some variance by rounding things off, combining layers, etc, but ran

	into speed problems.



	At that stage, I wondered whether the much faster cellular tiled approach could emulate

	the pattern, and thankfully, it could. I came up with the layered repeat-tile approach

    a while back, and have explained it in other examples. The result, here, is an

	artifact-free, animated triangle metric at a small fraction of the cost.



    Unlike my other cellular tiled examples, this one is animated. It was a proof of concept,

	which was put in as an afterthought, so probably needs some tuning. Nevertheless, it

	works well enough for the purpose of this demonstration.



	At best, a normal animated Voronoi algorithm utilizing the triangle metric would have

	nine iterations (many use more) consisting of several operations. The algorithm used

	here uses just four, each with a fraction of the operations.



	I haven't utilized the speed to its full potential. However, there are still two

	raymarched layers - each containing a more expensive distance metric - and my slowest

	computer can handle it easily. I'll do something more interesting with it later.



	By the way, I made a very basic example with very little code for anyone who'd like to

	look at the algorithm without having to sift through a bunch of window dressing. The

	link is below.



	// The main inspiration for this, and practically all the other examples.

	Blocks -IQ

    https://www.shadertoy.com/view/lsSGRc





	// Just the cellular block algorithm.

    Fast, Minimal Animated Blocks - Shane

	https://www.shadertoy.com/view/MlVXzd



*/





// Far plane. Redundant here, but included out of habit.

#define FAR 10.



float objID = 0.; // Object ID. Used to identify the large block and small block layers.





/*



// Just the function itself. Everything else can be ignored, if you just want to produce repeat

// block tile patterns.



// Distance metric. A slightly rounded triangle is being used, which looks a little more organic.

float dm(vec2 p){



    p = fract(p) - (syn_Level*5*Gain);



    //return max(abs(p.x)*.866025 + (p.y)*.5, -(p.y)); // Regular triangle metric.



    return (dot(p, p)*6.*syn_BassLevel*.25 + .75)*max(abs(p.x)*.866025 + p.y*syn_Level*5*Gain, -p.y);

    //return (1.-(dot(p, p)*4.*.5 + .5)*.5)*max(abs(p.x)*.866025 + p.y*.5, -p.y);

    //return (1.-dot(p, p)*4.*.166)*max(abs(p.x)*.866025 + p.y*.5, -p.y);



    //return max(-max(abs(p.x)*.866025 - p.y*.5, p.y) + .25, max(abs(p.x)*.866025 + p.y*.5, -p.y));

    //return (length(p) + .5)*max(abs(p.x)*.866025 + p.y*.5, -p.y);



}





// Very cheap wrappable cellular tiles. This one produces a block pattern on account of the

// metric used, but other metrics will produce the usual patterns.

//

// Construction is pretty simple: Plot two points in a wrappble cell and record their distance.

// Rotate by a third of a circle then repeat ad infinitum. Unbelievably, just one rotation

// is needed for a random looking pattern. Amazing... to me anyway. :)

//

// Note that there are no random points at all, no loops, and virtually no setup, yet the

// pattern appears random anyway.

float cell(vec2 p){





    // Matrix to rotate the layer by TAU/3. radians - or 120 degrees -

    // which makes sense when dealing with the equilateral triangle metric.

    const mat2 m = mat2(-.5, .866025, -.866025, -.5);



    // Abstract varitation, just for fun.

    //const mat2 m = mat2(-1, 1, -1, -1)*.7071;

    //const mat2 m = mat2(.5, .866025, -.866025, .5);



    // Two rotating point plus offset

    const float offs = .666 - .166;

    vec2 a = sin(vec2(1.93, 0) + TIME)*.166;

    float d0 = dm(p + vec2(a.x, 0));

    float d1 = dm(p + vec2(0, offs + a.y));



    // Rotate the layer, and plot another two points.

    p = m*(p + .5);

    float d2 = dm(p + vec2(a.x, 0));

    float d3 = dm(p + vec2(0, offs + a.y));



    // Find the distance to the nearest point.

    // It works with just one rotation and four points very well.

    return min(min(d0, d1), min(d2, d3))*2.;





    // Add another two points and a rotation, and the pattern looks even more random.

    //p = m*(p + .5);

    //float d4 = dm(p +  vec2(a.x, 0));

    //float d5 = dm(p +  vec2(0, offs + a.y));



    //return min(min(min(d0, d1), min(d2, d3)), min(d4, d5))*2.;





}



*/



// Distance metric. A slightly rounded triangle is being used, which looks a little more organic.

float dm(vec2 p){



    p = fract(p) - .5;



    //return max(abs(p.x)*.866025 + (p.y)*.5, -(p.y)); // Regular triangle metric.



    return (dot(p, p)*4.*.25 + .75)*max(abs(p.x)*.866025 + p.y*.5, -p.y);

    //return (1.-(dot(p, p)*4.*.5 + .5)*.5)*max(abs(p.x)*.866025 + p.y*.5, -p.y);

    //return (1.-dot(p, p)*4.*.166)*max(abs(p.x)*.866025 + p.y*.5, -p.y);



    //return max(-max(abs(p.x)*.866 - p.y*.5, p.y) + .25, max(abs(p.x)*.866025 + p.y*.5, -p.y));

    //return (length(p) + .5)*max(abs(p.x)*.866025 + p.y*.5, -p.y);



}



// Distance metric for the second pattern. It's just a reverse triangle metric.

float dm2(vec2 p){



    p = fract(p) - .5;



    //return max(abs(p.x)*.866025 + -p.y*.5, p.y);

    return (dot(p, p)*4.*.25 + .75)*max(abs(p.x)*.866025 - p.y*.5, p.y);



    //return max(-max(abs(p.x)*.866025 + p.y*.5, -p.y) + .2, max(abs(p.x)*.866025 - p.y*.5, p.y));

    //return (length(p)*1. + .5)*max(abs(p.x)*.866025 - p.y*.5, p.y);



}



// Very cheap wrappable cellular tiles. This one produces a block pattern on account of the

// metric used, but other metrics will produce the usual patterns.

//

// Anyway, plot two points in a wrappble cell and record the minimum distance, rotate by a third

// of a circle whilst storing the overall minimum, then repeat ad infinitum. In this case just

// one rotation is needed for a random looking pattern. Amazing.

//

// Note that there are no random points at all, no loops, and virtually no setup, yet the

// pattern appears random anyway.

//

// By the way, this particular function combines two patterns for the large and small blocks,

// but the original is commented out above, for anyone interested.

float cell(vec2 q){





    // SETUP.

    //

    // Matrix to rotate the layer by TAU/3. radians - or 120 degrees -

    // which makes sense when dealing with the equilateral triangle metric.

    const mat2 m = mat2(-.5, .866025, -.866025, -.5);



    // Abstract varitation, just for fun.

    //const mat2 m = mat2(-1, 1, -1, -1)*.7071;

    //const mat2 m = mat2(.5, .866025, -.866025, .5);



    // FIRST PATTERN.

    // Two rotating points plus offset

    vec2 p = q;

    const float offs = .666 - .166;

    vec2 a = sin(vec2(1.93, 0) + smoothTime*.5)*.166*(syn_BassLevel*BassPower*.5);

    float d0 = dm(p + vec2(a.x, syn_Level*.03*(Gain-1)));

    float d1 = dm(p + vec2(0, offs + a.y));



    // Rotate the layer, and plot another two points.

    p = m*(p + .5);

    float d2 = dm(p + vec2(a.x, 0));

    float d3 = dm(p + vec2(0, offs + a.y));



    // Find the distance to the nearest point.

    // It works with just one rotation and four points very well.

    float l1 = min(min(d0, d1), min(d2, d3))*ShapeSlide;





    // SECOND PATTERN... The small blocks, just to complicate things. :)

    p = q;

    d0 = dm2(p + vec2(a.x, 0));

    d1 = dm2(p + vec2(0, offs + a.y));



    // Rotate the layer, and plot another two points.

    p = m*(p + .5);

    d2 = dm2(p + vec2(a.x, 0));

    d3 = dm2(p + vec2(0, offs + a.y));



    // Find the distance to the nearest point.

    // It works with just one rotation and four points very well.

    float l2 = min(min(d0, d1), min(d2, d3))*(2.+pow(syn_Hits*.1+syn_MidHits*0.2, 2.0));





    // COMBINING PATTERNS.

    objID = step(l1, -(l2 - .4)); // Object, or pattern, ID.



    // Combine layers.

    return max(l1, -(l2 - .4 + Combiner));







}



// The heightmap. We're combining two patterns, each with their own distance metric, so the cell

// function is more complicated than it normally would be. However, there's another that has been

// commented out that people can refer to, if they're interested.

float heightMap(vec3 p){



    return cell(p.xy*2.); // Just one layer.



}



// The distance function. Just a heightmap function applied to a plane. Pretty standard stuff.

float map(vec3 p){



    float tx = heightMap(p);



    return Zoom - p.z + (.5 - tx)*.125;



}



// Normal calculation, with some edging and curvature bundled in.

vec3 nr(vec3 p, inout float edge, inout float crv) {



    // Roughly two pixel edge spread, regardless of resolution.

    vec2 e = vec2(2./RENDERSIZE.y, 0);



	float d1 = map(p + e.xyy), d2 = map(p - e.xyy);

	float d3 = map(p + e.yxy), d4 = map(p - e.yxy);

	float d5 = map(p + e.yyx), d6 = map(p - e.yyx);

    float highHits = pow(syn_Level*0.25 + syn_HighLevel*0.75, 1.25);

	float d = map(p)*(2.+((highHits+Shimmer)/10));


    edge = abs(d1 + d2 - d) + abs(d3 + d4 - d) + abs(d5 + d6 - d);

    edge = smoothstep(0., 1., sqrt(edge/e.x*2.+((highHits+Shimmer)/5)));



    // Wider sample spread for the curvature.

    e = vec2(12./450., 0);

	d1 = map(p + e.xyy), d2 = map(p - e.xyy);

	d3 = map(p + e.yxy), d4 = map(p - e.yxy);

	d5 = map(p + e.yyx), d6 = map(p - e.yyx);

    crv = clamp((d1 + d2 + d3 + d4 + d5 + d6 - d*(3. + smoothTimeC))*32. + .5, 0., 1.);





    e = vec2(2./450., 0); //RENDERSIZE.y - Depending how you want different resolutions to look.

	d1 = map(p + e.xyy), d2 = map(p - e.xyy);

	d3 = map(p + e.yxy), d4 = map(p - e.yxy);

	d5 = map(p + e.yyx), d6 = map(p - e.yyx);



    return normalize(vec3(d1 - d2, d3 - d4, d5 - d6));

}



// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)

// Anyway, I like this one. I'm assuming it's based on IQ's original.

float cao(in vec3 p, in vec3 n){



    float sca = 2., occ = 0.;

    for(float i=0.; i<6.; i++){



        float hr = .01 + i*.75/5.;

        float dd = map(n * hr + p);

        occ += (hr - dd)*sca;

        sca *= 0.8;

    }



    return clamp(1.0 - occ, 0., 1.);





}



// Compact, self-contained version of IQ's 3D value noise function.

float n3D(vec3 p){



	const vec3 s = vec3(7, 157, 113);

	vec3 ip = floor(p); p -= ip;

    vec4 h = vec4(0., s.yz, s.y + s.z) + dot(ip, s);

    p = p*p*(3. - 2.*p); //p *= p*p*(p*(p * 6. - 15.) + 10.);

    h = mix(fract(sin(h)*43758.5453), fract(sin(h + s.x)*43758.5453), p.x);

    h.xy = mix(h.xz, h.yw, p.y);

    return mix(h.x, h.y, p.z); // Range: [0, 1].

}



// Simple environment mapping. Pass the reflected vector in and create some

// colored noise with it. The normal is redundant here, but it can be used

// to pass into a 3D texture mapping function to produce some interesting

// environmental reflections.

vec3 eMap(vec3 rd, vec3 sn){



    vec3 sRd = rd; // Save rd, just for some mixing at the end.



    // Add a time component, scale, then pass into the noise function.

    //rd.xy -= TIME*.25;
    rd.xy -= smoothTime*.1;

    rd *= 3.;



    float c = n3D(rd)*.57 + n3D(rd*2.)*.28 + n3D(rd*4.)*.15; // Noise value.

    c = smoothstep(0.5, 1., c); // Darken and add contast for more of a spotlight look.



    //vec3 col = vec3(c, c*c, c*c*c*c).zyx; // Simple, warm coloring.

    vec3 col = vec3(min(c*1.5, 1.), pow(c, 2.5), pow(c, 12.)).zyx; // More color.



    // Mix in some more red to tone it down and return.

    return mix(col, col.yzx, sRd*.25+.25);



}



vec4 renderMainImage() {

	vec4 fragColor = vec4(0.0);

	vec2 u = _xy;





    // Unit direction ray, ray origin (camera position), and light.

    vec3 rd = normalize(vec3(u - RENDERSIZE.xy*.5, RENDERSIZE.y)),

         ro = vec3(-TIME*.125, -TIME*.05, 0), l = ro + vec3(.5, -1.5, -1.);



    /*

    // Mild perspective and lens effects.

    r = normalize(vec3(r.xy, (r.z - length(r.xy)*.2)*1.2));

    vec2 a = sin(vec2(1.57, 0) - 3.14159/32.);

    r.yz = mat2(a, -a.y, a.x) * r.yz;

    r.xz = r.xz*mat2(a, -a.y, a.x);

    */





    // Raymarching against a back plane usually doesn't require many iterations -

    // nor does it require a far-plane break - buy I've given it a few anyway.

    float d, t = 0.;



    for(int i=0; i<64;i++){



        d = map(ro + rd*t); // Distance the nearest surface point.

        if(abs(d)<0.001 || t>FAR) break; // The far-plane break is redundant here.

        t += d*.86; // The accuracy probably isn't needed, but just in case.

    }



    //t = min(t, FAR); // Capping "t" to the far plane. Not need here.



    float svObjID = objID; // Store the object ID just after raymarching.



    vec3 sCol = vec3(0); // Scene color.



    // Edge and curvature variables. Passed into the normal function.

    float edge = 0., crv = 1.;



    if(t<FAR){



        vec3 p = ro + rd*t, n = nr(p, edge, crv);//normalize(fract(p) - .5);



        l -= p; // Light to surface vector. Ie: Light direction vector.

        d = max(length(l), 0.001); // Light to surface distance.

        l /= d; // Normalizing the light direction vector.





        // Attenuation and extra shading.

        float atten = 1./(1. + d*d*.05);

        float shade = heightMap(p);





        // Texturing. Because this is a psuedo 3D effect that relies on the isometry of the

        // block pattern, we're texturing isometrically... groan. :) Actually, it's not that

        // bad. Rotate, skew, repeat. You could use tri-planar texturing, but it's doesn't

        // look quite as convincing in this instance.

        //

        // By the way, the blocks aren't perfectly square, but the texturing doesn't seem to

        // be affected.

        vec2 tuv = vec2(0);

        vec3 q = p;

        const mat2 mr3 = mat2(.866025, .5, -.5, .866025); // 60 degrees rotation matrix.

        q.xy *= mr3; // Rotate by 60 degrees to the starting alignment.

        if((n.x)>.002) tuv = vec2((q.x)*.866 - q.y*.5, q.y); // 30, 60, 90 triangle skewing... kind of.

        q.xy *= mr3*mr3; // Rotate twice for 120 degrees... It works, but I'll improve the logic at some stage. :)

        if (n.x<-.002) tuv = vec2((q.x)*.866 - q.y*.5, q.y);

        q.xy *= mr3*mr3; // Rotate twice.

        if (n.y>.002) tuv = vec2((q.x)*.866 - q.y*.5, q.y);



        // Pass in the isometric texture coordinate, roughly convert to linear space (tx*tx), and

        // make the colors more vibrant with the "smoothstep" function.

        //vec3 tx = texture(image3, tuv*2.).xyz;

        vec3 tx = vec3(.85, .2, .5);



        if(svObjID>.5) tx *= vec3(2, .9, .3); // Add a splash of color to the little blocks.





        float ao = cao(p, n); // Ambient occlusion. Tweaked for the this example.





        float diff = max(dot(l, n), 0.); // Diffuse.

        float spec = pow(max(dot(reflect(l, n), rd), 0.), 6.); // Specular.

        //diff = pow(diff, 4.)*0.66 + pow(diff, 8.)*0.34; // Ramping up the diffuse.





        // Cheap way to add an extra color into the mix. Only applied to the small blocks.

        if(svObjID>.5) {

        	float rg = dot(sin(p*6. + cos(p.yzx*4. + 1.57/3.)), vec3(.333))*.5 + .5;

        	tx = mix(tx, tx.zxy, smoothstep(0.6, 1., rg));

        }





        // Applying the lighting.

        sCol = tx*(diff + .5) + vec3(1, .6, .2)*spec*3.;





        // Alternative, mild strip overlay.

        //sCol *= clamp(sin(shade*6.283*24.)*3. + 1., 0., 1.)*.35 + .65;





        // Adding some cheap environment mapping to help aid the illusion a little more.

        sCol += (sCol*.75 + .25)*eMap(reflect(rd, n), n)*3.; // Fake environment mapping.



        //sCol = pow(sCol, vec3(1.25))*1.25; More contrast, if you were going for that look.



        // Using the 2D block value to provide some extra shading. It's fake, but gives it a

        // more shadowy look.

        sCol *= (smoothstep(0., .5, shade)*.75 + .25);



        // Applying curvature, edging, ambient occlusion and attenuation. You could apply this

        // in one line, but I thought I'd seperate them for anyone who wants to comment them

        // out to see what effect they have.

        sCol *= min(crv, 1.)*.7 + .3;

        sCol *= 1. - edge*.85;

        sCol *= ao*atten;





    }





    // Rough gamma correction.

    fragColor = vec4(sqrt(clamp(sCol, 0., 1.)), 1.);





	return fragColor;

 }





vec4 renderMain(){

	if(PASSINDEX == 0){

		return renderMainImage();

	}

}
