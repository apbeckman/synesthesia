

/*
	Apollonian Structure
	--------------------

	Overall, there's nothing particularly exciting about this shader, but I've 
	always liked this structure, and thought it'd lend itself well to the two 
	tweet environment.

	I couldn't afford shadows, AO, etc, so applied a bit of fakery to at least
	convey that feel.


*/

// Apollonian based fractal: I couldn't find the original source, but it's
// been around for a while. IQ has a really cool variation here: 
// Apollonian: https://www.shadertoy.com/view/4ds3zn
//
// I'm guessing the original was posted by someone on a fractal forum somewhere 
// at some stage.
//
float m(vec3 p){
    
    // Moving the scene itself forward, as opposed to the camera.
    // IQ does it in one of his small examples.
    p.z += TIME;
    // Loop counter and variables.
    float i = 0., s = 1., k;

	// Repeat Apollonian distance field. It's just a few fractal related 
    // operations. Break up space, distort it, repeat, etc. More iterations
    // would be nicer, but this function is called a hundred times, so I've
    // used the minimum to give just enough intricate detail.
    while(i++ < 9.) p *= k = 1.5/dot(p = mod(p - 1., 2.) - 1., p), s *= k;
		
	// Render numerous little spheres, spread out to fill in the 
    // repeat Apollonian lattice-like structure you see.
    //
    // Note the ".01" at the end. Most people make do without it, but
    // I like the tiny spheres to have a touch more volume, especially
    // when using low iterations.
    return length(p)/s - .0025; 
    
}


vec4 renderMainImage() {
	vec4 c = vec4(0.0);
	vec2 u = _xy;

    // Direction ray and origin. By the way, you could use "o = d/d" (Thanks, Fabrice),
    // then do some shuffling around in the lighting calculation, but I didn't quite 
    // like the end result, so I'll leave it alone, for now anyway.
    vec3 d = vec3(u/RENDERSIZE.y - .5, 1)/4., o = vec3(1, 1, 0);
    d.xy = _rotate(d.xy, Rotate*PI);
    //d.xy+= -.16;

    // Initialize to zero.
    c -= c;
    
    // Raymarching loop -- sans break, which always makes me cringe. :)
    while(c.w++<1e2) o += m(o)*d;

    
    // Lame lighting - loosely based on directial derivative lighting and the 
    // way occlusion is performed, but mostly made up. It'd be nice to get rid 
    // of that "1.1," but it's kind of necessary.  Note that "o.z" serves as  
    // a rough distance estimate, and to give a slight volumetric light feel. 
    //c += (m(o - .01)*m(o - d)*4e1 + o.z*1.1 - 2.)/o.z;
    // I stared at the line above for ages and got nothing. Fabrice looked at it
    // instantly, and saw the obvious. :)
    c += (m(o - .01)*m(o - d)*4e1 - 2.)/o.z + 1.1;
    
	return c; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}