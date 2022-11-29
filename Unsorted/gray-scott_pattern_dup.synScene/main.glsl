

			//******** BuffA Code Begins ********

/*
	Gray-Scott Pattern
	------------------

    Diffusion reaction: The process is simple enough to understand, but whenever I'm 
    working through it, I often wonder which brainiac nerd poured a solution from one 
    beaker into another, then thought to themselves, "That'd be a great way to make 
    procedural giraffe textures." :) ... Actually, I think the initial idea came from
	the famous Engima code breaking guy, Alan Turing, but don't quote me on that. :)
       
    At its very core, diffusion textures are just a visual representation over time of 
    a particular kind of partial differential equation, nothing more, so in essence, you 
    could ignore all explanations and simply apply it. However, if you're like me and 
    you require a better understanding, refer to the articles below. If you're after a 
	decent summary, refer to the following article: 
     
    http://www.karlsims.com/rd.html

	Anyway A lot of the code here is just some standard noise routines for anyone who
	wants to experiment with different intial conditions. The actual code to produce the 
	diffusion texture is contained in the "mainImage" function and there isn't much of 
	that at all.


	// Articles that may be helpful:

	Reaction Diffusion: The Gray-Scott Algorithm
 	http://www.algosome.com/articles/reaction-diffusion-gray-scott.html
	
	Gray Scott Model of Reaction Diffusion
	https://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/
    
    Karl Sims - http://www.karlsims.com/rd.html
	
	Robert Munafo - http://mrob.com/pub/comp/xmorphia/

	Reaction–Diffusion System
	https://en.wikipedia.org/wiki/Reaction%E2%80%93diffusion_system

*/
float growthFactor = normalize(pow((syn_BassLevel*0.5)+(syn_MidLevel*0.35)+(syn_Level*0.15), 2.0));

float hash21(vec2 p){
    p = mod(p, 64.);
    return fract(sin(dot(p, vec2(1.27, 113.93)))*43758.5453);
}

// 2x2 hash algorithm.
vec2 hash22(vec2 p) { 

    p = mod(p, 128.);///vec2(RENDERSIZE.y/RENDERSIZE.x, 1)
    // More concise, but wouldn't disperse things as nicely as the block below.
    float n = sin(dot(p,vec2(113, 1))); 
    return fract(vec2(2097152, 262144)*n);
    
    // Animation.
    p = fract(vec2(2097152, 262144)*n);
    return sin(p*6.2831853)*0.5 + 0.5;
}

/*
vec3 hash33(in vec2 p){ 
    float n = sin(dot(p, vec2(41, 289)));    
    return fract(vec3(2097152, 262144, 32768)*n); 
}
*/


vec3 sTexture(sampler2D smp, vec2 uv) {
 
    vec2 textureResolution = RENDERSIZE.yy;
	uv = uv*textureResolution + 0.5;
	vec2 iuv = floor( uv );
	uv -= iuv;
	uv = iuv + smoothstep(0., 1., uv); 
    //uv = iuv +  uv*uv*uv*(uv*(uv*6. - 15.) + 10.);
	uv = (uv - .5)/textureResolution;
    return texture(smp, uv).xyz;
    
}


/*
// IQ's value noise formula.
float noise( in vec2 p ){
   
    vec2 i = floor(p); p -= i; 
    p *= p*p*(p*(p*6. - 15.) + 10.);
    //p *= p*(3. - p*2.);  

    return mix( mix( hash21(i + vec2(0, 0)), 
                     hash21(i + vec2(1, 0)), p.x),
                mix( hash21(i + vec2(0, 1)), 
                     hash21(i + vec2(1, 1)), p.x), p.y);
}
*/

// 2D 2nd-order Voronoi: Obviously, this is just a rehash of IQ's original. I've tidied
// up those if-statements. Since there's less writing, it should go faster. That's how 
// it works, right? :)
//
float Voronoi(vec2 p){
    
	vec2 g = floor(p), o;
	p -= g;// p = fract(p);
	
	vec2 d = vec2(2); // 1.4, etc.
    
	for(int y = -1; y <= 1; y++){
		for(int x = -1; x <= 1; x++){
            
			o = vec2(x, y);
            o += hash22(g + o) - p;
            
			float h = dot(o, o*(1.0+cos(TIME*0.1)*0.001));
            d.y = max(d.x, min(d.y, h)); 
            d.x = min(d.x, h);            
		}
	}
	
	return min(sqrt(d.x), 1.);
    //return min(d.y - d.x, 1.); // etc.
}




// Shorthand, so that the texture lines read a little better.
vec4 tx(vec2 p){ return texture(BuffA, p); }


// 25 (or 9) tap Laplacian -- Gaussian Laplacian, to be more precise. I think of it as taking
// the sum of the partial second derivatives of a blurry 2D height map... in each channel...
// I think I'm making things more confusing, but it works anyway. :D Seriously, just look
// up the Laplacian operator of a 2D function.
vec4 Laplacian(vec2 p) {
    

	// Kernel matrix dimension, and a half dimension calculation.
    const int mDim = 5, halfDim = (mDim - 1)/2;

    // You can experiment with different Laplacian related setups here. Obviously, when 
    // using the 3x3, you have to change "mDim" above to "3." There are widely varying 
    // numerical variances as well, so that has to be considered also.
    float kernel[mDim*mDim] = float[mDim*mDim](
    
 		//The Laplacian-Gaussian... I can't remember where I found this,
        // but I've seen it in more than one place, and it seems about right.
		0.,  0., .25,  0.,  0.,
        0., .25,  .5, .25,  0.,
        .25, .5,  -4., .5, .25,
        0., .25,  .5, .25,  0.,
        0.,  0., .25,  0.,  0.);
/*       
        // Another variation -- Might need scaling first.
        1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1.,
        1., 1.,-24.,1., 1.,
        1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1.);
*/

/*      

		// 3x3 variation.
        1., 2., 1.,
        2., -12., 2.,
        1., 2., 1.);
*/
    
    
/*      
		// 3x3 variation. Slightly different to above, and scaled differently,
		// so that has to be taken into account.
        .05, .2, .05,
        .2,  -1., .2,
        .05, .2, .05);
*/
    
    // Initiate the color. In this example, we only want the XY values, but just
    // in case you'd like to apply this elsewhere, I've included all four texture
    // channnels.
    vec4 col = vec4(0);
    
    // We're indexing neighboring pixels, so make this a pixel width. In fact, these
    // texture Laplacian calculations are annoyingly touchy, so it has to be one
    // pixel width. Not two, not a half... Computers are tools. :D
    float px = 1./RENDERSIZE.y; 

    
    // There's a million boring ways to apply a kernal matrix to a pixel, and this 
    // is one of them. :)
    for (int j=0; j<mDim; j++){
        for (int i=0; i<mDim; i++){ 
            
            col += kernel[j*mDim + i]*tx(p + vec2(i - halfDim, j - halfDim)*px);
        }
    }


    return col;
}



// Nine tap Laplacian.
vec4 Laplacian9(vec2 p) {
    
    vec2 px = 1.*_uvc/RENDERSIZE.yy;
    // Four spots aren't required in this case, but are when the above isn't aspect correct.
    vec4 e = vec4(px, 0, -px.x);
 
    // Laplacian with some Gaussian smoothing applied... I'm guessing, so I probably wouldn't
    // quote me on it. :)
    return (tx(p - e.xy)*.5 +  tx(p - e.zy ) + tx(p - e.wy)*.5 + // First row.
			tx(p - e.xz) - tx(p)*6. + tx(p + e.xz) + 		     // Seond row.
			tx(p + e.wy)*.5 + tx(p + e.zy) +  tx(p + e.xy)*.5);  // Third row
  
    
    // Laplacian with no Gaussian component... I think.
    //return (tx(p - e.xy) +  tx(p - e.zy ) + tx(p - e.wy) +   // First row.
	//		tx(p - e.xz) - tx(p)*8. + tx(p + e.xz) + 		 // Seond row.
	//		tx(p + e.wy) + tx(p + e.zy) +  tx(p + e.xy))/2.; // Third row
 
}



// Five tap Laplacian -- The simplest Laplacian filter... unless there's a more minimalistic one.
vec4 Laplacian5( vec2 p) {
    
    vec3 e = vec3(1./RENDERSIZE.yy, 0);

	return tx(p - e.zy) + tx(p - e.xz) - tx(p)*4. + tx(p + e.xz) + tx( p +  e.zy);
}




// Analytic Laplacian of sorts.
vec4 grad(in vec2 uv) {
    
    vec2 e = vec2(1, 0)/RENDERSIZE.y/2.;
	vec2 dfdx = tx(uv + e).xy - tx(uv - e).xy;
    vec2 dfdy = tx(uv + e.yx).xy - tx(uv - e.yx).xy;
    return vec4(dfdx, dfdy)*RENDERSIZE.y;
}

// Divergence of gradient.
vec2 LaplacianA(in vec2 uv) {
 
    vec2 e = vec2(1, 0)/RENDERSIZE.y/2.;
    vec2 dgdx = grad(uv + e).xy - grad(uv - e).xy;
    vec2 dgdy = grad(uv + e.yx).zw - grad(uv - e.yx).zw;
    vec2 lap = (dgdx + dgdy)*RENDERSIZE.y;
    
    return lap/RENDERSIZE.y/RENDERSIZE.y;
}





vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    // Screen coordinates. I can't remember why I didn't want aspect correctness here...
    // I'm sure there was a reason. Fixed size square buffers would make life a lot
    // easier -- I know that much. :)
    vec2 p = fragCoord.xy/RENDERSIZE.xy;
    

    // Grab the reaction diffusion values from the previous frame. These values are
    // are representative of concentrations of hypothetical liquid, gas, etc, solutions,
    // which are denoted by A and B. A is stored in the X channel, and B is stored in 
    // the Y channel.
	vec4 rdVal = texture(BuffA, p);
    
    // The concentrations of elements A and B will tend to spread out in a smoother -- and
    // sometimes, faster -- way if we smooth out the underlying values. Making that happen 
    // is as simple as blurring the A and B concoction (the X and Y texture values) every 
    // pass. The wider the blur, the better. You could achieve that by by adding 
    // an extra buffer and performing the blur in conjunction with the reaction 
    // diffusion calculations performed here, or you could combine a blur matrix with 
    // the Lapacian matrix all in one step, which is what we're doing here. In other 
    // words we're combining a 5x5 Gaussian matrix with a Laplacian matrix and 
    // using the resultant matrix in the Laplacian step.
    //
    // And yes, more buffers could be added to decrease equilibrium time, but I 
    // wanted to keep things simple.
      
    
    // In regard to the form of the reaction-diffusion equation we're using, there's some 
	// kind of physical flow involved, which tends to require second derivative measurements.
	// To get those from an underlying pixelated map, pixel samples from a spread of neighbors 
	// will be required, so it's not a stretch to imagine that some kind of pixel matrix will 
	// be involved. Hence, the Laplacian step below:
    //
    // 25 tap Laplacian to really smoothen things out.

    vec2 lap;// = Laplacian(p).xy;
    // Other Laplacian functions to experiment with. The different pixels arrangements
    // produce subtle differences. Each function would need to be uncommented above.
    //vec2 lap = Laplacian9(p).xy;
    //vec2 lap = LaplacianA(p);
    //vec2 lap = Laplacian5(p).xy;
    int LPM = int(Laplacian_Mode); //lap mode

    if(LPM==0.){
    lap = Laplacian(p).xy; //original
    }

    if(LPM==1.){
    lap = Laplacian9(p).xy; //
    }
    if(LPM==2.){
    lap = LaplacianA(p).xy; //
    }
    if(LPM==3.){
    lap = Laplacian5(p).xy; //
    }

    // AB: refactored to allow each of the alternate sets of feed/kill values to be used via the Feed_Kill_Mode slider
    // Original code below
    float feed, kill;
    int FKM = int(Feed_Kill_Mode); //feed kill mode integer

    if(FKM==0.){
    feed = 0.0545, kill = 0.062; //original constants
    }

    if(FKM==1.){
        feed = 0.058, kill = 0.065; // Lines and dots.
    }

    if(FKM==2.){
        feed = 0.046, kill = 0.063; // Lines and dots2.
    }

    if(FKM==3.){
        feed = 0.098, kill = 0.057; // Solid cells.
    }

    if(FKM==4.){
        feed = 0.037, kill = 0.06; // Random joined lines.
    }

    if(FKM==5.){
        feed = 0.058, kill = 0.065; // Lines and dots.
    }

    if(FKM==6.){
        feed = 0.03, kill = 0.0565; // Maze.
    }
    /*
    if(FKM==6.){
        feed = 0.03, kill = 0.063; // Self replicating spots.
    }
    */
    if(FKM==7.){
        feed = 0.082, kill = 0.06; // Self replicating spots.
    }

    //AB: fine adjust
    feed += Feed;
    kill += Kill;

    // Feed and kill rates. These are close to standard values that you see everywhere.
    // You can change them, but they're sensitive.
    //float feed = 0.0545, kill = 0.062; //original constants
    // More constants to try. A lot of them are presets that I found using a handy pattern
    // generator at this link: http://mrob.com/pub/comp/xmorphia/ogl/index.html

    //const float feed = 0.058, kill = 0.065; // Lines and dots.
    //const float feed = 0.046, kill = 0.063; // Lines and dots.
    //const float feed = 0.098, kill = 0.057; // Solid cells.
    //const float feed = 0.037, kill = 0.06; // Random joined lines.
    //float feed = 0.0353, kill = 0.06; // Random joined lines.
    //float feed = 0.03, kill = 0.063; // Self replicating spots.

    //feed *= (0.99975+basshits*0.00025);
    //const float feed = 0.03, kill = 0.0565; // Maze.
    
    // Diffusion rates for concentrations A and B. The first component needs to be higher
    // than the second. The amplitude depends on the size of the Laplacian. I've noticed 
    // that if "dAB.x" is exactly twice "dAB.y," the pattern will form quicker, but won't
    // vary greatly. Tiny changes in these values will widely change the pattern, but can
    // also bring about a blank screen. I guess that's one of the downsides of diffusion 
    // patterns.
	//vec2 dAB = vec2( ((.2)*(0.9975+growthFactor*0.0025))+DiffA, .106+DiffB);
    vec2 dAB = vec2(.2, .1);
    dAB += vec2 (DiffA, DiffB);
    // Time component. Kind of redundant here, but it can be used to control reaction rate.
    // Unfortunately, like all the figures used here, just a tiny change can ruin the 
    // reaction, which results in no pattern at all. Either way, you should try to set 
    // this number as high as you can without trashing the pattern.
    float t = 1.;//*(0.99875+growthFactor*0.00125); 
    //float t = growthFactor;
    // The diffusion term: Just the constant diffusion rates multiplied by the Laplacian
    // for each concentration.
    vec2 diffusion = (dAB)*lap;
    diffusion *= (1.+growthFactor*0.0025);
    // The reaction term: The "rdVal.x*rdVal.y*rdVal.y" value is representative of the
    // chance that a particle from concentration A will react with two particles from
    // concentration B, which from probability theory is: A*B*B. This results in a decrease
    // of concentration A and an increase in concentration B by that particular amount.
    // Hence, the negative and positive vector terms on the end.
    rdVal.y*=(0.9975+growthFactor*0.0035);
    //rdVal.y *=(0.999875+growthFactor*0.00025);
    vec2 reaction = vec2(rdVal.x*rdVal.y*rdVal.y)*vec2(-1, 1);
    reaction.x *=1.-(Reaction*0.025);
    reaction.y *=1.+(Reaction*0.005);
    reaction.xy *= vec2(1.0+growthFactor*0.0001, 1.-growthFactor*0.0001);
    // Feed and kill rates. Substance A is added at a feed rate, and substance B is taken
    // away at a kill rate. Hence, the positive and negative vector terms on the end.
    // It took me a while to figure out why the "1. - rdVal.x" term was necessary. It's 
    // necessary to reduce the amount that is fed into the system as the concentration of
    // A "rdVal.x" builds up, otherwise, we'd never reach equilibrium.
    vec2 feedKill = vec2(feed*(1. - rdVal.x), (kill + feed)*rdVal.y)*vec2(1, -1);
    // Try the following with just the "rdVal.x" and the initial condition:
    // fragColor.xy = vec2(1, 0) + (hash22(p*64.).xy - .5)*.75;
    //
    // Interesting, but equilibrium is never attained.
    //vec2 feedKill = vec2(feed*rdVal.x, (kill + feed)*rdVal.y)*vec2(1, -1);

    // Calculating the change in concentration of A and B. This calculation the crux of the
    // whole thing. It's just an applied partial differential equation... which is a little 
    // difficult to write in ASCII form, but easy to apply. To see the actual equations in
    // mathematical or physical form, you can read about it in the sources I've provided above.
    
    // New u value: du\dt = rU*LapU - u*v*v + f*(1. - u);
    // New v value: dv\dt = rV*LapV + u*v*v + (f + k)*v;  
	vec2 delta = diffusion + reaction + feedKill;
    
    // Updating the old A and B concentrations by adding the change in concentration
    // over time.
    fragColor.xy = clamp(rdVal.xy + delta*t, 0., 1.);
    
    // Cute trick to allow reinitialization whenever the user switches between window
    // sizes. However, if all four channels were needed, you'd have to load in another
    // buffer and store it there.
    fragColor.zw = RENDERSIZE.xy; 
    
    
    //AB: trying to get fancy, transplanting some code from glassier
    /*float distFunc = 0.0;

    vec3 logoCol = vec3(0.0);
    if (_exists(syn_UserImage)){
        logoCol = _loadUserImageAsMask().rgb;
    }

    if (image_block < 0.5){

        distFunc = distFunc-length(logoCol)*10.5;
    }


    if (image_block>=0.5){

        if (logoCol.r>0.5){
            fragColor.r = 1.0*syn_BassLevel;
        }
    }

    else{
        distFunc = 0.;
    }*/
    // Initializing the substances A and B: It requires a little finesse.
    //
    // You could literally spend weeks playing around with concentration initialization and 
    // never quite understand what works and what doesn't... OK, maybe that's just me, but
    // I've never been able to determine a general routine that enables a pattern to form.
    // Either way, here are a bunch of different initial conditions that do produce patterns.
    //
    // Sometimes, the application won't recognize the first frame -- or something, so it's 
    // necessary to initialize for the first few frames to give it a chance to catch on.
    // Also, check for changes in canvas size.
    //AB: added reinitialization trigger for 'Reset' control
    if(FRAMECOUNT<10 || abs(rdVal.w - RENDERSIZE.y)>.001|| Reset > 0.)  {

        //AB: adding some simple variance to the initial state so resetting multiple times looks better

        //fragColor.xy =  vec2(1, 0) +  (vec2(Voronoi(p*64.), Voronoi(p*64. + vec2(2.93*(1.0+cos(smoothTime)), 19.37+sin(smoothTime)))) - .5); 

        
        //AB: adding option for user media
        if (_exists(syn_UserImage)){
        //fragColor.xy =  vec2(1, 0) +  texture(syn_UserImage, fragCoord.xy/RENDERSIZE.y).xy *.5;   
        fragColor.xy =  vec2(1, 0) +  (vec2(Voronoi(p*64.), Voronoi(p*64. + vec2(2.93, 19.37))) - .5)*0.1+vec2(-1, 1)*(sTexture(syn_UserImage, p).xy*(1.+basshits) - .5);
        }
        else{
        // Required multi-channel noise texture in "iChannel1."
        //fragColor.xy =  vec2(1, 0) +  texture(image30, fragCoord.xy/RENDERSIZE.y).xy -.5;   
        fragColor.xy =  vec2(1, 0) +  (vec2(Voronoi(p*64.), Voronoi(p*64. + vec2(2.93, 19.37))) - .5)+basshits*  (vec2(Voronoi(p*64.), Voronoi(p*64. + vec2(2.93, 19.37))) - .5);
        }
        //fragColor.xy =  vec2(1, 0) +  vec2(-1, 1)*(Voronoi(p*64.)*1.3 - .5);
        // Some of these require functions above which need to be uncommented first.
        //fragColor.xy =  vec2(1, 0) +  vec2(-1, 1)*(sTexture(image30, p).xy - .5);
        //fragColor.xy =  vec2(1, 0) +  vec2(-1, 1)*(noise(p*64.) - .5);
        //fragColor.xy =  vec2(1, 0) +  vec2(-1, 1)*(abs(noise(p*64.) - .5)*2. - .5);
        //fragColor.xy =  vec2(1, 0) + (hash22(p*64.).xy - .5); //(hash22(p*64.).xy - .5)*1.1, etc.
        //fragColor.xy =  vec2(.5, -.5)*.75 + vec2(Voronoi(p*64.), Voronoi(p*64. + vec2(3, 7))); 
        //fragColor.xy =  vec2(.45, -.53) + dot(sin(p*9.*6.283 - cos(p.yx*7.*6.283 + 1.15)*3.14159), vec2(.25)) + .5;
        
    }
	return fragColor; 
 } 




/*
	Gray-Scott Pattern
	------------------

	Of all the standard texturing algorithms, I probably like the look of reaction-diffusion textures most. 
	Unfortunately, I find them a little cumbersome to work with inside realtime shader environments, so I tend 
	to stick with the classics, like smooth noise and Voronoi.

	Anyway, this is just a basic example that I put together for my own amusement. For the really interesting 
    possibilities, take a look at some of Flexi and Cornusammonis's shaders. 

	I restricted the pattern to something pretty simple, since I didn't want people having to wait for 20 years
	for something legible to form. This particular pattern is recognizable in under ten seconds and tends to reach
    a state of equilibrium in about 30 seconds -- Not instant gratification, but not too bad for just one extra
	buffer. Obviously, outside of Shadertoy, the idea would be to prerender to a wrappable texture, then load it 
	in.	

	There isn't just one way to perform diffusion, but there are some common ones that graphics programmers like
    to employ. The Gray-Scott model is pretty popular, and that's the one I'm utilizing here. It's been done many 
    times before, but I wanted to post my own version. I went out of my way to make the formula as straight forward 
    as possible and to explain it. Thanks to the people on Shadertoy who posted examples before me, because having
    working shader examples to refer to always makes life easier.

	In regard to the rendering of the pattern, this is a just a bit of bump mapping. The lighting, shadows, etc,
	are all layered. 3D rendering can get a little repetitive, so just for fun, I occasionally enjoy applying 2D 
    techniques to imagery in an attempt to make them look 3D -- We all need our hobbies. :D

	Anyway, I have a proper 3D version that I'll put up at some stage.

	// Other diffusion examples:

    // I have a 3D version of this that I might put up at some stage.
	Viscous Fingering - cornusammonis
	https://www.shadertoy.com/view/Xst3Dj

	// I remember Flexi's original from years ago. :)
	expansive reaction-diffusion - Flexi
    https://www.shadertoy.com/view/4dcGW2

	// A nice all round example.
	ReactDiff Experiment 1 - aiekick
	https://www.shadertoy.com/view/XtlcDl

	// More interactive and fun than my example. :)
	Gray-Scott Explorer 2 - Dr2
	https://www.shadertoy.com/view/MtlfDN

    // Fabrice has a pretty sharp mind. He whipped this up not long after reading a pretty esoteric paper. 
    // I remember perusing the paper briefly, then I put it on the "I'll look at it later" list, and we all 
    // know what that means. :D
	Seashells reaction-diffusion - FabriceNeyret2
	https://www.shadertoy.com/view/MdGGRD

*/

// Global resolution variable: Kind of hacky, but necessary to accommodate a repeat texture created on a 
// variable sized buffer that needs to coincide with canvas size changes... Working with variable sized
// buffers can be a little tiring. :)
float gResY;

// 2x2 matrix rotation. Note the absence of "cos." It's there, but in disguise, and comes courtesy
// of Fabrice Neyret's "ouside the box" thinking. :)
mat2 rot2( float a ){ vec2 v = sin(vec2(PI/2., 0) + a);	return mat2(v, -v.y, v.x); }


// Tri-Planar blending function. Based on an old Nvidia writeup:
// GPU Gems 3 - Ryan Geiss: https://developer.nvidia.com/gpugems/GPUGems3/gpugems3_ch01.html
vec3 tex3D(sampler2D t, in vec3 p, in vec3 n){

    
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


// Cheap and nasty 2D smooth noise function with inbuilt hash function - based on IQ's 
// original. Very trimmed down. In fact, I probably went a little overboard. I think it 
// might also degrade with large time values. I'll swap it for something more robust later.
float n2D(vec2 p) {

	vec2 i = floor(p); p -= i; 
    //p *= p*p*(p*(p*6. - 15.) + 10.);
    p *= p*(3. - p*2.);  
    
	return dot(mat2(fract(sin(vec4(0, 1, 113, 114) + dot(i, vec2(1, 113)))*43758.5453))*
                vec2(1. - p.y, p.y), vec2(1. - p.x, p.x) );

}

// Gradient noise fBm.
float fBm(in vec2 p){
    
    return n2D(p)*.57 + n2D(p*2.)*.28 + n2D(p*4.)*.15;
    
}
/*
// 2x2 hash algorithm.
vec2 hash22(vec2 p) { 

    //p = mod(p, 24.);
    // More concise, but wouldn't disperse things as nicely as the block below.
    float n = sin(dot(p,vec2(1, 113))); 
    //return fract(vec2(2097152, 262144)*n);
    
    // Animation.
    p = fract(vec2(2097152, 262144)*n);
    return sin(p*6.2831853 + TIME)*0.5 + 0.5;
}

// 2D 2nd-order Voronoi: Obviously, this is just a rehash of IQ's original. I've tidied
// up those if-statements. Since there's less writing, it should go faster. That's how 
// it works, right? :)
//
float Voronoi(vec2 p){
    
	vec2 g = floor(p), o;
	p -= g;// p = fract(p);
	
	vec2 d = vec2(1); // 1.4, etc.
    
	for(int y = -1; y <= 1; y++){
		for(int x = -1; x <= 1; x++){
            
			o = vec2(x, y);
            o += hash22(g + o) - p;
            
			float h = dot(o, o);
            d.y = max(d.x, min(d.y, h)); 
            d.x = min(d.x, h);            
		}
	}
    
    
	
	//return sqrt(d.y) - sqrt(d.x);
    return (d.y - d.x); // etc.
}

*/

float ggs;

// Bump mapping function. Put whatever you want here. In this case, 
// we're returning the length of the sinusoidal warp function.
float bumpFunc(vec2 p){ 

    //vec3 tx = texture(iChannel1, p*2.).xyz; tx *= tx;
    //float bump = dot(tx, vec3(.299, .587, .114));
    
 
    p = p/2.25*450./gResY;//*vec2(RENDERSIZE.y/ RENDERSIZE.x, 1)
	float gs = 1. - (texture(BuffA, (p)).x);//*.9 + bump*.1; // Range: [0, 1]
    

    
    ggs = gs;
    
   
    
    // Trick to smoothly flatten out the top a bit.
    gs = mix(gs, smoothstep(0., .1, gs - .275), .75);
    
    //gs += smoothstep(0., .1, ggs - .125)*.05 - .025;
    gs = mix(gs, smoothstep(0., .05, ggs - .525), .015);
     
    
    // Subtly blend in a more detailed pattern over the top.
    float difPat = texture(BuffA, p*6., -100.).x;
    difPat = smoothstep(0., .5, difPat - .65);
    //difPat *= (1.0+basshits);
    
    difPat *= smoothstep(0., .1, ggs - .525);
    
    gs = mix(gs, difPat, .015); //difPat*gs
    
    gs += smoothstep(0., .05, ggs - .525)*.05 - .025;

    
    return clamp(gs, 0., 1.);//*.95 + bump*.05;

}

vec3 bump(vec3 sp, vec3 sn, float bumpFactor, inout float edge){
    
      // BUMP MAPPING - PERTURBING THE NORMAL
    //
    // Setting up the bump mapping variables. Normally, you'd amalgamate a lot of the following,
    // and roll it into a single function, but I wanted to show the workings.
    //
    // f - Function value
    // fx - Change in "f" in in the X-direction.
    // fy - Change in "f" in in the Y-direction.
    vec2 eps = vec2(1./RENDERSIZE.y, 0.);
    
    float f = bumpFunc(sp.xy); // Sample value multiplied by the amplitude.
    float fx = bumpFunc(sp.xy - eps.xy); // Same for the nearby sample in the X-direction.
    float fy = bumpFunc(sp.xy - eps.yx); // Same for the nearby sample in the Y-direction.
    
    
    
    // Using the above to determine the dx and dy function gradients.
    fx = (fx - f)/eps.x; // Change in X
    fy = (fy - f)/eps.x; // Change in Y.
    // Using the gradient vector, "vec3(fx, fy, 0)," to perturb the XY plane normal ",vec3(0, 0, -1)."
    // By the way, there's a redundant step I'm skipping in this particular case, on account of the 
    // normal only having a Z-component. Normally, though, you'd need the commented stuff below.
    //vec3 grad = vec3(fx, fy, 0);
    //grad -= sn*dot(sn, grad);
    //sn = normalize( sn + grad*bumpFactor ); 
      
    sn = normalize( sn + vec3(fx, fy, 0)*bumpFactor ); 
    
    
    eps = vec2(6./RENDERSIZE.y, 0);
    fx = bumpFunc(sp.xy - eps.xy);// - bumpFunc(sp.xy + eps.xy); // Same for the nearby sample in the X-direction.
    fy = bumpFunc(sp.xy - eps.yx);// - bumpFunc(sp.xy + eps.yx);
    edge = (abs(f*1. - fx) + abs(f*1. - fy))*(0.75+highhits); // Edge value.
 	
 
    
    return sn;
    
}


// Standard ray-plane intersection.
float rayPlane(vec3 p, vec3 o, vec3 n, vec3 rd) {
    
    float dn = dot(rd, n), t = 1e8;
    
    if (abs(dn)>.0001){
        t = dot(p - o, n)/dn;
        t += float(t<0.)*1e8;
    }
    
    return t;
}


//vec3 smoothFract(vec3 x){ x = fract(x); return min(x, x*(1.-x)*12.); }



vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    // The following hack is necessary to accommodate a repeat texture created on a variable 
    // sized buffer. The reason it's capped between  350 and 800 pixels is because the 
    // pattern looks too bulky when moving to fullscreen and too busy on smaller canvases. 
    // Of course, on a phone, the  PPX can be double that of a laptop, so this figure would 
    // have to be set to double for the viewer to see the same size image I'm viewing... 
    // Unfortunately, I can't cater for every situation, so I don't. :)
    gResY = clamp(RENDERSIZE.y, 350., 800.);
    
    // Centered, aspect-correct, screen coordinates. Only one line is necessary to
    // accomplish this.
	vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/gResY;

    // BASIC SETUP - surface postion, ray origin, unit direction vector, and light postion.
    //
	// Setup for a raytraced plane. Not all that necessary, since you could effect a tiled
    // plane without any trace, but I didn't see the harm. :)
    //
	vec3 rd = normalize(vec3(_uvc, 1.)); // Unit direction vector.
    rd = normalize(vec3(rd.xy, sqrt(rd.z*rd.z - dot(rd.xy, rd.xy)*_uvc*.25))); // Very subtle lens warping.
    
	vec3 ro = vec3(0., 0., -2/FOV); // Camera position, ray origin, etc.

	// Plane normal -- Titled very slightly.
    vec3 sn = normalize(vec3(.01+cos(smoothTime*0.1)*0.005, .01+sin(smoothTime*0.1)*0.005, -1)); 
	// Hit point on the plane.
	float t = rayPlane(vec3(0), ro, sn, rd);
    vec3 sp = ro + rd*t;
    float FOVlp = FOV*2.0;
    vec3 lp = vec3(.125+cos(smoothTimeB*0.3), .35+sin(smoothTimeB*0.3), -1); // Light position - Back from the screen.
    lp.xy/= vec2(FOVlp);
    sp.xy += Camera_Pos.xy;
    sp.xy -= vec2(-5*sin(smoothTime*0.0025+TIME*0.01), -5.5*cos(smoothTime*0.0025+TIME*0.01));
    lp.xy -= vec2(- 5*sin(smoothTime*0.0025+TIME*0.01), -5.5*cos(smoothTime*0.0025+TIME*0.01));
	
    //vec3 sn = vec3(0., 0., -1); // Plane normal. Z pointing toward the viewer.
    
    // Last term controls how much the bump is accentuated.
    float edge = 0.;
    sn = bump(sp, sn, 0.25, edge);
    
    // The main bump mapped pattern. It has already been calculated in the "bump" function
    // above, so this double handling, but it's a pretty cheap example.
    float mainPat = bumpFunc(sp.xy);      
    
    
    // LIGHTING
    //
	// Determine the light direction vector, calculate its distance, then normalize it.
	vec3 ld = lp - sp;
	float lDist = max(length(ld), 0.001);
	ld /= lDist;

    // Light attenuation.    
    float atten = min(1./(1. + lDist*0.125 + lDist*lDist*0.05), 1.);
	//float atten = min(1./(lDist*lDist*1.), 1.);
    
    // Using the bump function, "f," to darken the crevices. Completely optional, but I
    // find it gives extra depth.
    atten *= smoothstep(0., 1., mainPat*.5 + .5 - .25);//*.95 + .05; // Or... f*f*.5 + .5; //  pow(f, .75); // etc.

	

	// Diffuse value.
	float diff = max(dot(sn, ld), 0.);  
    // Enhancing the diffuse value a bit. Made up.
    diff = pow(diff, 2.)*0.66 + pow(diff, 4.)*0.34; 
    // Specular highlighting.
    float spec = pow(max(dot( reflect(-ld, sn), -rd), 0.), 64.)*(1.+highhits); 

	
    // TEXTURE COLOR
    //
    vec3 texCol = mix(vec3(.5), vec3(.6, .4, .2), dot(sin(sp.xy*4. - cos(sp.yx*4.)), vec2(.25)) + .5);//doColor(vor);
    
    
    // Lightening the ground a little.
    texCol = mix(vec3(1), texCol, mainPat);
    
    // Using the diffusion texture to create a more refined diffusion pattern to decorate the main object.
    float difPat2 = texture(BuffA, sp.xy/2.25*450./gResY*6.).x;
    difPat2 = smoothstep(0., .05, -difPat2 + .5);
    
    // Small pattern edges. Comment it out to see what it does.
    texCol = mix(texCol, vec3(0), (1. - smoothstep(0., .05, max(ggs - .525 - .025, -(ggs - .525 + .025))))*.35);
    
    // Coloring the main object yellow, and the ground brownish with some fine diffusion crevices,
    // or something to that effect.
    texCol = mix(vec3(.6, .48, .42)*(difPat2*.25 + .75), vec3(1, .45, .05)*1.5, mainPat)*texCol;
    
    // Darkening the fine grade diffusion crevices on the yellow object.
    texCol *= mix(vec3(1), vec3(0), difPat2*smoothstep(0., .05, ggs - .525)*.6);//.525
    //texCol *= mix(vec3(1), vec3(0), clamp(difPat2*mainPat*smoothstep(0., .05, ggs - .525)*.75, 0., 1.));
    
    // Applying some orangey patches.
    texCol = mix(texCol*1.5, texCol.xzy/4., dot(sin(sp.xy*12. - cos(sp.yx*8.)), vec2(.2*.5)) + .2);
    
    //texCol = mix(texCol, texCol.yxz*mainPat*.7, fBm(sp.xy*3.));
    // Greenish moldy weathering.
    texCol = mix(texCol, texCol.yxz*mix(.2, .7, mainPat), smoothstep(.35, 1., fBm(sp.xy*7.))*.7);
    
    // Extra combinations.
    //texCol = mix(texCol, texCol.xzy, smoothstep(.6, .8, fBm(sp.xy*3. + .5))*mainPat*.7);
    //texCol = mix(texCol, texCol.zyx, smoothstep(.6, .8, fBm(sp.xy*3. - .5))*mainPat*.9);
    //texCol = mix(vec3(1)*dot(texCol, vec3(.299, .587, .114)), texCol.zyx, smoothstep(.5, .7, fBm(sp.xy*3.))*mainPat*.8);

    
    // FINAL COLOR
    // Using the values above to produce the final color.   
    vec3 col = (texCol*(diff*2. + 0.25) + (texCol*.5 + .5)*vec3(.25, .5, 1)*spec)*atten;
    
    // Lightening the edges a bit -- Very subtle, but it enhances the pseudo 3D look a fraction.
    col += col*vec3(1)*diff*edge;
    
    // Extra dark edging. Not needed here.
    //col *= vec3(1)*(1. - edge*.75);
    
    
    // Fake shadows. Just a bit of sampling and masking trickery. I made it up on the spot,
    // but it works well enough. Obviously, the fake shadows are there to add false depth.
    vec3 lp2 = lp + vec3(.5, 1.4, -.5);
    vec2 dir = normalize(lp2 - sp).xy*1./min(800., RENDERSIZE.y);
    float shad = max(1. - bumpFunc(sp.xy + dir*16.), mainPat);
    float shad2 = max(bumpFunc(sp.xy - dir*16.), 1. - mainPat);
    shad = min(shad, shad2);
    // Applying the shadow layer. Comment this line out to see the difference. The image becomes
    // very flat looking.
    col *= vec3(1)*(shad*.75 + .25);

    
    // If you want to see just the plane diffusion pattern.
    //vec4 reactDiff = texture(BuffA, sp.xy/2.25*450./gResY);
    //float sm = 1. - reactDiff.x;
    //col = vec3(1)*smoothstep(0., .075, sm - .475);
    
     
    // Subtle vignette.
    //uv = fragCoord/RENDERSIZE.xy;
    //col *= pow(16.*uv.x*uv.y*(1. - uv.x)*(1. - uv.y) , .125) + .1;
    // Colored variation.
    //col = mix(col.zyx/2., col, 
              //pow(16.*uv.x*uv.y*(1. - uv.x)*(1. - uv.y) , .125));

    // Rough gamma correction, and we're done.
	fragColor = vec4(sqrt(min(col, 1.)), 1.);
    
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