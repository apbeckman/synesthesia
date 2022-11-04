/*

	Rounded Voronoi Blocks
	----------------------

    Applying some bump mapping and lighting to a Voronoi variation to produce an 
	animated, pseudo-lit, rounded block effect.

	I had originally intended to do a standard square block version using a first order
    triangular distance metric, but noticed that IQ, Aeikick and others had already done 
	it, so I switched to a second order distance setup and experimented with a few different 
	metrics. The end result was the oddly distributed round-looking blocks with defined 
	edges you see.

	To add a little more to the illusion, I bump mapped it, put in some fake occlusion,
	and some subtle edging. As you can see, there's not a lot of code.


	Related examples:

	// Uses the more standard, triangle metric. So stylish.
	IQ - Blocks 
	https://www.shadertoy.com/view/lsSGRc

	// Aiekick's tunnelized take on the above. 
	Voro Tri Tunnel - aiekick
	https://www.shadertoy.com/view/XtGGWy

*/





// vec2 to vec2 hash.
vec2 hash22(vec2 p) { 

    // Faster, but doesn't disperse things quite as nicely. However, when framerate
    // is an issue, and it often is, this is a good one to use. Basically, it's a tweaked 
    // amalgamation I put together, based on a couple of other random algorithms I've 
    // seen around... so use it with caution, because I make a tonne of mistakes. :)
    float n = sin(dot(p, vec2(41, 289)));
    //return fract(vec2(262144, 32768)*n); 
    
    // Animated.
    p = fract(vec2(262144, 32768)*n); 
    // Note the ".25," insted of ".5".
    return sin( p*6.2831853 + TIME )*.25 + .5; 
    
}


//float tri(float x){ return abs(fract(x) - .5)*2.; }

// 2D 2nd-order Voronoi: Obviously, this is just a rehash of IQ's original. I've tidied
// up those if-statements. Since there's less writing, it should go faster. That's how 
// it works, right? :)
//
float Voronoi(in vec2 p){
    
	vec2 g = floor(p), o; p -= g; // Cell ID, offset variable, and relative cell postion.
	
	vec3 d = vec3(1); // 1.4, etc. "d.z" holds the distance comparison value.
    
	for(int y=-1; y<=1; y++){
		for(int x=-1; x<=1; x++){
            
			o = vec2(x, y); // Grid cell ID offset.
            o += hash22(g + o) - p; // Random offset.
			
            // Regular squared Euclidean distance.
            d.z = dot(o, o); 
            // Adding some radial variance as we sweep around the circle. It's an old
            // trick to draw flowers and so forth. Three petals is reminiscent of a
            // triangle, which translates roughly to a blocky appearance.
            d.z *= cos(atan(o.y, o.x)*3. - 3.14159/2.)*.333 + .667;
            //d.z *= (1. -  tri(atan(o.y, o.x)*3./6.283 + .25)*.5); // More linear looking.
            
            d.y = max(d.x, min(d.y, d.z)); // Second order distance.
            d.x = min(d.x, d.z); // First order distance.
                      
		}
	}

    // A bit of science and experimentation.
    return d.y*.5 + (d.y-d.x)*.5; // Range: [0, 1]... Although, I'd check. :)
    
    //return d.y; // d.x, d.y - d.x, etc.
    
    
}


// Bump mapping function. Put whatever you want here. If you wish to do some
// fake shadowing, it's usually helpful to keep the range between zero and one.
float bumpFunc(vec2 p){ 
    
    return Voronoi(p*4.); // Range: [0, 1] 
    
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    // Screen coordinates.
	vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;
    
    
    // PLANE ROTATION
    //
    // Rotating the canvas back and forth. I don't feel it adds value, in this case,
    // but feel free to uncomment it.
    //vec2 a = sin(vec2(1.57, 0) + sin(TIME*0.1)*sin(TIME*0.12)*2.);
    //uv *= mat2(a, -a.y, a.x);
  
    // Position offset.
    vec3 offs = vec3(-TIME/16., -TIME/8., 0);

    // VECTOR SETUP - surface postion, ray origin, unit direction vector, and light postion.
    //
    // Setup: I find 2D bump mapping more intuitive to pretend I'm raytracing, then lighting a bump mapped plane 
    // situated at the origin. Others may disagree. :)  
    vec3 sp = vec3(uv, 0) + offs; // Surface posion. Hit point, if you prefer. Essentially, a screen at the origin.
    vec3 rd = normalize(vec3(uv, 1.)); // Unit direction vector. From the origin to the screen plane.
    vec3 lp = vec3(cos(TIME)*0.5, sin(TIME)*0.2, -2.) + offs; // Light position - Back from the screen.
	vec3 sn = vec3(0., 0., -1); // Plane normal. Z pointing toward the viewer.
 
     

    
    // BUMP MAPPING - PERTURBING THE NORMAL
    //
    // Setting up the bump mapping variables and edge calcultion. Normally, you'd amalgamate a lot 
    // of the following, and roll it into a single function, but I wanted to show the workings.
    //
    // f - Function value
    // fx - Change in "f" in in the X-direction.
    // fy - Change in "f" in in the Y-direction.
    vec2 eps = vec2(3./RENDERSIZE.y, 0.);
    
    float f = bumpFunc(sp.xy); // Function value.
    float fx = bumpFunc(sp.xy+eps.xy); // Nearby sample in the X-direction.
    float fy = bumpFunc(sp.xy+eps.yx); // Nearby sample in the Y-direction.
  
    float fx2 = bumpFunc(sp.xy-eps.xy); // Sample on the other side in the X-direction.
    float fy2 = bumpFunc(sp.xy-eps.yx); // Same on the other side in the Y-direction.
    
    // Using the samples to provide an edge measurement. How you do it depends on the
    // look you're going for.
    //float edge = abs(fx + fy + fx2 + fy2 - 4.*f); //abs(fx - f) + abs(fy - f);
    float edge = abs(fx + fy + fx2 + fy2) - 4.*f; //abs(fx - f) + abs(fy - f);
    //float edge = abs(fx + fx2) + abs(fy + fy2) - 4.*f; //abs(fx - f) + abs(fy - f);

    edge = smoothstep(0., 8., edge/eps.x*4.);
   
 	// Controls how much the bump is accentuated.
	const float bumpFactor = 0.35;
    
    // Using the above to determine the dx and dy function gradients.
    fx = (fx-fx2)/eps.x/2.; // Change in X
    fy = (fy-fy2)/eps.x/2.; // Change in Y.
    // Using the gradient vector, "vec3(fx, fy, 0)," to perturb the XY plane normal ",vec3(0, 0, -1)."
    // By the way, there's a redundant step I'm skipping in this particular case, on account of the 
    // normal only having a Z-component. Normally, though, you'd need the commented stuff below.
    //vec3 grad = vec3(fx, fy, 0);
    //grad -= sn*dot(sn, grad);
    //sn = normalize( sn + grad*bumpFactor ); 
    sn = normalize( sn + vec3(fx, fy, 0)*bumpFactor );           
     
    
    // LIGHTING
    //
	// Determine the light direction vector, calculate its distance, then normalize it.
	vec3 ld = lp - sp;
	float lDist = max(length(ld), 0.001);
	ld /= lDist;

    // Light attenuation.    
    float atten = 1./(1.0 + lDist*lDist*0.15);
    
    // Using the bump function, "f," to darken the crevices. Completely optional, but I
    // find it gives extra shadowy depth.
    //atten *= smoothstep(0.1, 1., f)*.9 + .1; // Or... f*f*.7 + .3; //  pow(f, .75); // etc.
    atten *= ((1.-f)*.9 + .1); // Or... f*f*.7 + .3; //  pow(f, .75); // etc.
    

	
    // TEXTURE COLOR
    //
	// Some fake tri-plannar mapping. I wouldn't take it too seroiusly.    
    const float ts = 2.;
    vec3 nsn = max(abs(sn)-.2, .001);
    nsn /= dot(nsn, vec3(1));
    vec3 texCol = vec3(0);
    sp.z += -f*1.; // Pretending that the Z-value is not sitting flat on the plane.
    // Tri-planar.
    texCol += texture(image3, sp.xy*ts).xyz*nsn.z;
    texCol += texture(image3, sp.xz*ts).xyz*nsn.y;
    texCol += texture(image3, sp.yz*ts).xyz*nsn.x;
    // sRGB to linar with processing. Basically, minipulating the color a bit.
    texCol = smoothstep(0.075, .5, texCol*texCol)*1.;
    
    
    // Diffuse value.
	float diff = max(dot(sn, ld), 0.);  
    // Enhancing the diffuse value a bit. Made up.
    diff = pow(diff, 4.)*0.66 + pow(diff, 8.)*0.34; 
    // Specular highlighting.
    float spec = pow(max(dot( reflect(-ld, sn), -rd), 0.), 8.); 
    // Fresnel reflection. It's one line, so why not.
    float fres = pow(clamp(dot(rd, sn) + 1., 0., 1.), 16.);
    
    
    // Apply the edging here. You could do it later, but I wanted the specular lighting to 
    // supercede it. Not for a scientific reason. I just thought it looked nicer this way. :)
    texCol *= 1.-edge*.5;
    
    // I did this by accident, but found it added to the depth, so I kept it.
    texCol *= smoothstep(0.1, .6, atten);
	
    
    // FINAL COLOR
    // Using the values above to produce the final color.   
    vec3 col = (texCol*(diff*2. + 0.2 + vec3(1, 1.4, 2)*fres) + vec3(1., .7, .3)*spec*2.)*atten;

    // Doing the edging here is fine, but it overrides the specular, which wasn't the look I was going for.
    //col *= 1. - edge*.5;
    
    // Postprocesing - A subtle vignette, contrast and coloring.
    uv = fragCoord/RENDERSIZE.xy;
    col *= pow(16.*uv.x*uv.y*(1.-uv.x)*(1.-uv.y), 0.125); // Vignette.
  
    col = vec3(1.4, 1.2, .9)*pow(max(col, 0.), vec3(1, 1.2, 1.5))*1.5; // Contrast, coloring.
    

    // Perform some statistically unlikely (but close enough) 2.0 gamma correction. :) 
	fragColor = vec4(sqrt(clamp(col, 0., 1.)), 1.);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}