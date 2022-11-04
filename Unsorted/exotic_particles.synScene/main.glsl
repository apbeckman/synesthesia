vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********



// Start off with a function, warp it, and accumulate color along the way.
// This one is just a more mutated version of a simple sine warp function,
// of which there are plenty of examples on Shadertoy.
vec3 warp(vec2 u, float ph1, float ph2){

    // Initializing the warped UV coordinates. This gives it a bit 
    // of a worm hole quality. There are infinitly other mutations.
    vec2 v = u - log(1./max(length(u), .001))*vec2(-1, 1);
    
    // Scene color.
    vec3 col = vec3(0.);
    
    // Number of iterations.
    const int n = 9;
    
    for (int i = 0; i<n; i++){
    
        // Warp function.
        v = cos(v.y - vec2(0, 1.57))*exp(sin(v.x + ph1) + cos(v.y + ph2));
        v -= u;
        
        // Color via IQ's cosine palatte and shading.
        vec3 d = (.5 + .45*cos(vec3(i)/float(n)*3. + vec3(0, 1, 2)*1.5))/max(length(v), .001);
        // Accumulation.
        col += d*d/32.;
        
        // Adding noise for that fake path traced look. 
        // Also, to hide speckling in amongst noise. :)
        /*col += fract(sin(u.xyy*.2 + u.yxx + dot(u + fract(TIME), 
                    vec2(113.97, 27.13)))*45758.5453)*.01 - .005;*/    }
    
    return col;
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    // Aspect correct UV coordinates.
    vec2 u = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y*2.;

   
    // Angular offsets.
    float ph1 = smoothTime*.1;
    float ph2 = sin(smoothTimeC*0.25)*.125;
    
    // Adding two warp functions phase shifted by a certain amount was
    // Jolle's interesting addition. Just the one would work, but isn't
    // as interesting.
    vec3 col = warp(u, ph1, ph2) + warp(u, ph1, ph2 + 1.57);
    
    // Toning things down slightly.
    col = mix(col, col.zyx, .1);
    
    // Noise, for that fake path traced feel. :)
    //col.xyz += fract(sin(u.xyy*.7 + u.yxx + dot(u + fract(TIME), 
    //                 vec2(113.97, 27.13)))*45758.5453)*.1 - .05;    
    
    // Mix the previous frames in.
    vec4 preCol = texelFetch(BuffA, ivec2(fragCoord), 0);
    float blend = (FRAMECOUNT < 2) ? 1. : .25; 
    col = mix(preCol.xyz, col, blend);
    
    
    // Clamp and add to Buffer A.
    fragColor = vec4(clamp(col, 0., 1.), 1);
	return fragColor; 
 } 


/*

    Exotic Particles
    ----------------
    
    Accumulating color values via transcental function warping to create some 
    pretty moving imagery that resembles colliding exotic particles in a 
    chamber... of paint... I actually have no idea what this looks like. :D
    
    Function warping is nothing new, and this particular example is just a 
    slightly dressed up version of Jolle and Jarble's previous work, which in 
    turn was very loosely based on one of Lomateron's recent examples -- The 
    respective links are below.
    
    Anyway, I've commented the code. However, there's definitely nothing 
    difficult to grasp here. The simple color imagery was produced in 
    "Buffer A", which was blended with previous frames for a bit of temporal 
    blurring. The result (Image tab) was then used to take two 3x3 blurred 
    samples in order to add some highlights.
    
    
    
    
    Uses elements from the following shaders:
    
    Glass bubble lamp - Jarble: https://www.shadertoy.com/view/ttcfD7
    
    Glass bubble lamp fork - Jolle: https://www.shadertoy.com/view/WtdBDM
    
    Mount Mask - lomateron: https://www.shadertoy.com/view/WdsfRf
    
    
*/

// Things look cleaner without highlights, and in some ways I prefer it.
// However, it's less interesting... I think? :)
#define HIGHLIGHTS

// Serves no other purpose than to save having to write this out all the time. I'm using this 
// on a buffer texture, so no sRGB to linear operation needs to be performed. I'm also
// using (and prefer to use) aspect correct pixel coordinates, so it's necessary to stretch 
// out the X values before retrieving them. It's also possible to stretch out the UV coordinates
// first, then use a stretched sample spread, which is faster... Yeah, it's confusing, but it 
// doesn't matter, just so long as you have a method you're happy with. :)
//
vec4 tx(in vec2 p){ 
     p *= vec2(RENDERSIZE.y/RENDERSIZE.x, 1);
     return texture(BuffA, p + .5/RENDERSIZE.y); 
}

// Blur function. Pretty standard.
vec4 bTx(in vec2 p){
    
    // Sample spread -- Measured in pixels.
    float px = 2.;
    
    // Result.
	vec4 c = vec4(0);
    
    // Standard equally weighted 3x3 blur.
    for(int i = 0; i<9; i++) c += tx(p + (vec2(i/3, i%3) - 1.)*px/RENDERSIZE.y);
 
    // Normalizing the return value.
    return c/9.;  
    
    /*
    // NxN blur.
    const int N = 5;
    for(int i = 0; i<N*N; i++) c += tx(p + (vec2(i/N, i%N) - float(N - 1)/2.)*px/RENDERSIZE.y);
    return c/float(N*N); 
    */
}
 
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    // Aspect correct pixel coordinates.
    vec2 uv = fragCoord/RENDERSIZE.y;
   
    // A 3x3 blurred texture sample. The generated warped imagery contains a few
    // high frequency speckles, so blurred samples mitigate that somewhat. Denoising
    // would be better, but this will do.
    vec4 col = bTx(uv);
    //vec4 col = tx(uv); // Standard single sample.
     
    #ifdef HIGHLIGHTS
    // Bump mapping via cheap, directional derivative-based highlighting.
    vec2 px = 4./RENDERSIZE.yy; // Sample spread.
    vec4 col2 = bTx(uv - px); // Seperate sample.
    float b = max(dot(col2 - col, vec4(.299, .587, .114, 0)), 0.)/length(px); // Bump.
    col += col2.yzxw*col2.yzxw*b/12.; // Add the colored highlights.
    #endif
    
    // Toning down the lower half slightly.
    col = mix(col, col.zyxw, max(.3 - uv.y, 0.));
    

    // Rough gamma correction.
    fragColor = sqrt(max(col, 0.));
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