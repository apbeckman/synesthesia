

			//******** BuffA Code Begins ********

/*

    Random Asymmetric Quad Grid
    ---------------------------
    
    I coded up an asymmetric quad cell grid some time ago, just for
    the fun of it, then forgot about it. I came across it the other 
    day and decided to add some bump mapping, highlights, etc, and 
    wound up with whatever this is. :) It's a 2D example, so should 
    run at a reasonable pace on most machines.
    
    Out of sheer boredome, I put some letters in the cells, then 
    proceeded to hack away at it until it looked like... a picture of
    a keyboard that's been run through a poorly written AI art 
    generating algorithm or something. :)
    
    I originally used static cells, which is a little more believable,
    since hardened graphite material doesn't flow, but I thought the
    animation made it a little more interesting.
    
    There's nothing ground breaking in here, but the grid code might
    be interesting to some. The code was originally written with 
    extruded grids in mind. I have a few examples, and will post one 
    or two at some stage. There's a custom bloom-like postprocessing 
    routine in the "Image" tab that someone might find useful also.
    

*/


// Animate the blocks or not. The animation reduces the beveled 
// letter illusion a bit, but it's more interesting.
#define ANIMATE

// Different color configurations.
// Orange: 0, Blue: 1, Green 2, Orange and blue: 3.
#define COLOR 0


// IQ's vec2 to float hash.
float hash21(vec2 p){ return fract(sin(dot(p, vec2(117.619, 57.623)))*43758.5453); } 


// Standard 2D rotation formula.
mat2 r2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }


// Commutative smooth maximum function. Provided by Tomkh, and taken 
// from Alex Evans's (aka Statix) talk: 
// http://media.lolrus.mediamolecule.com/AlexEvans_SIGGRAPH-2015.pdf
// Credited to Dave Smith @media molecule.
float smax(float a, float b, float k){

   float f = max(0., 1. - abs(b - a)/k);
   return max(a, b) + k*.25*f*f;
}
 

// Entirely based on IQ's signed distance to a 2D triangle -- Very handy.
// I have a generalized version somewhere that's a little more succinct,
// so I'll track that down and drop it in later.
float quad(in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2, in vec2 p3){

    //////////////////
    // Extra calculations here.
    //////////////////
    
    // Unfortunately, to render rounded quads on a grid, the vertices need to 
    // be moved inward in order to add a rounding factor and have the quad still
    // fit within the cell. Contracting all vertices by a simple shrinking factor 
    // doesn't preserve edge width, so some extra trigonometry is required.
    //
    // By the way, I can't but help think there must be a simpler way to go about
    // this that I'm not seeing, so if someone can figure out a way to produce 
    // the same without the following addition, feel free to let me know. :)


    // Nudge factor to move the vertices inward. 
    const float ndg = .02;
    
    // Normalized edge lines.
	vec2 e0 = normalize(p1 - p0);
	vec2 e1 = normalize(p2 - p1);
	vec2 e2 = normalize(p3 - p2);
	vec2 e3 = normalize(p0 - p3);
    
    // Angle between vectors. Calculating four at a time, which one would
    // assume is less draining on the GPU.
    vec4 ang, sl;
    ang = acos(vec4(dot(e0, -e3), dot(e1, -e0), dot(e2, -e1), dot(e3, -e2)));
    sl = ndg/tan(ang/2.);
    // Using the above to move the vertices in by the nudge factor.
    p0 += sl.x*e0 + ndg*e0.yx*vec2(1, -1);
    p1 += sl.y*e1 + ndg*e1.yx*vec2(1, -1);
    p2 += sl.z*e2 + ndg*e2.yx*vec2(1, -1);
    p3 += sl.w*e3 + ndg*e3.yx*vec2(1, -1);
    
    // Setting the new edge vectors.
    e0 = p1 - p0;
    e1 = p2 - p1;
    e2 = p3 - p2;
    e3 = p0 - p3;
    
    //////////////////
    //////////////////

    // Continuing with IQ's quad calculations... Actually, this is a 
    // slighlty modified version of his triangle function. 
	vec2 v0 = p - p0;
	vec2 v1 = p - p1;
	vec2 v2 = p - p2;
	vec2 v3 = p - p3;
   
	vec2 pq0 = v0 - e0*clamp(dot(v0, e0)/dot(e0, e0), 0., 1.);
	vec2 pq1 = v1 - e1*clamp(dot(v1, e1)/dot(e1, e1), 0., 1.);
	vec2 pq2 = v2 - e2*clamp(dot(v2, e2)/dot(e2, e2), 0., 1.);
	vec2 pq3 = v3 - e3*clamp(dot(v3, e3)/dot(e3, e3), 0., 1.);
    
    float s = sign(e0.x*e3.y - e0.y*e3.x);
    vec2 d = min( min( vec2(dot(pq0, pq0), s*(v0.x*e0.y - v0.y*e0.x)),
                       vec2(dot(pq1, pq1), s*(v1.x*e1.y - v1.y*e1.x))),
                       vec2(dot(pq2, pq2), s*(v2.x*e2.y - v2.y*e2.x)));
    
    d = min(d, vec2(dot(pq3, pq3), s*(v3.x*e3.y - v3.y*e3.x)));
 
    // Signed, rounded quad distance.
	return -sqrt(d.x)*sign(d.y) - ndg;
}


// IQ's vec2 to float hash.
vec2 hash22TA(vec2 p){ 

    p = texture(image30, p/63.2453).xy;
    return sin(p*6.2831853 + TIME);//p;//
    
  
}

// IQ's vec2 to float hash.
vec2 hash22T(vec2 p){ return texture(image30, p/65.6217).xy*2. - 1.; }


vec4 gVal; // Storage for the cell contents.
vec2 gSID; // Static ID.

// Global normal. It's a weird place to put it, but I'm doing the 
// bump mapping inside the "map" function to save calculations.
vec3 gN; 

vec2 scale = vec2(1.5, 1)/10.; // Scale.

// The asymmetric quad grid. 
vec4 map(vec2 p){
    
    // Distance and edge width.
    float d = 1e5;
    float ew = .005;
    
    // Scale.
    vec2 sc = scale;
    
    
    // Centers for all four tiles.
    const mat4x2 cntr = mat4x2(vec2(-.5), vec2(-.5, .5), vec2(.5), vec2(.5, -.5)); 
    
    // Saving the local coordinates and four vertices of the nearest quad. 
    // Only used for bump mapping.
    vec2 svQ;
    mat4x2 svV;
    
    // Because the asymmetric boundaries of the quads overlap neighboring cells, 
    // neighbors need to be considered. In this case four cell renders will do.
    const int n = 2;
    //const float m = floor(float(n)/2. + .001) - .5;
    for(int i = 0; i<n*n; i++){

        // Local coordinates and ID.
        vec2 q = p.xy;
        vec2 iq = floor(q/sc - cntr[i]) + .5; 
        q -= (iq)*sc;
        
        // The four vertices for this cell.
        mat4x2 v = cntr;
        
        // Offset the vertices.
        #ifdef ANIMATE
        v[0] += hash22TA(iq + v[0])*.2;
   		v[1] += hash22TA(iq + v[1])*.2;
        v[2] += hash22TA(iq + v[2])*.2; 
        v[3] += hash22TA(iq + v[3])*.2;
        #else
        v[0] += hash22T(iq + v[0])*.2;
   		v[1] += hash22T(iq + v[1])*.2;
        v[2] += hash22T(iq + v[2])*.2; 
        v[3] += hash22T(iq + v[3])*.2;
        #endif
        // Scale.
        v[0] *= sc; v[1] *= sc; v[2] *= sc; v[3] *= sc;
        

        // Render the quad. This is a modified version of one of IQs
        // polygon functions to render rounded general quads that
        // fit on a grid.
        float d2 = quad(q, v[0], v[1], v[2], v[3]) + ew;
        
        // If this quad distance is nearer, update.
        if(d2<d){
            
            // New distance, static ID and moving ID.
            // The zero field is an unused height value holder.
            d = d2;
            gSID = iq;
            iq += (v[0] + v[1] + v[2] + v[3])/4./sc;
            gVal = vec4(d2, 0., iq);
            
            // Saving the local coordinates and vertices to perform 
            // some bump mapping. These were hacked in as an afterthough.
            // Normally, you wouldn't need to do this.
            svQ = q;
            svV = v;
      
        }
    
    }
    
    // Set the normal to the back plane.
    gN = vec3(0, 0, -1);
    
    // If we've hit a quad, bump map the normal.
    if(d<1e5){
    
        // Saving some calculations by constructing the normal here,
        // instead of calling the entire map function two more times.
        //
        vec2 e = vec2(.005, 0); // Sample spread.

        vec2 q = svQ;
        mat4x2 v = svV;

        // Nearby X and Y samples -- Standard gradient calculations.
        float d2X = quad(q - e, v[0], v[1], v[2], v[3]) + ew;// - .02;
        float d2Y = quad(q - e.yx, v[0], v[1], v[2], v[3]) + ew;// - .02;
        // Modifying the samples to give a slightly rounded beveled look.
        d = mix(d, smax(d, smax(d*2., -.01, .01), .1), .75); 
        d2X = mix(d2X, smax(d2X, smax(d2X*2., -.01, .01), .1), .75);  
        d2Y = mix(d2Y, smax(d2Y, smax(d2Y*2., -.01, .01), .1), .75);  
        /*
        // Alternative flatter tops, which are more button like.
        d = smax(smax(d*6., -.03, .01), smax(d*2., -.02, .01), .03); 
        d2X = smax(smax(d2X*6., -.03, .01), smax(d2X*2., -.02, .01), .03); 
        d2Y = smax(smax(d2Y*6., -.03, .01), smax(d2Y*2., -.02, .01), .03); 
        */
        // Applying the gradient to the flat plane normal, vec3(0, 0, -1).
        float bf = .7; // Bump factor.
        gN = normalize(vec3((d - d2X)/e.x, (d - d2Y)/e.x, -1./bf));
    }
    
    
    // Combining the floor with the extruded object.
    return  gVal;
 
}

// Obtaining an ASCII character from the font texture. This is pretty standard.
// I think I based it on one of Fabrice's early examples. I'll find it and put
// in the link later.
vec4 chTx(vec2 p, int n){

    vec2 i = vec2(n&15, 15 - (n>>4))/16.;
    return (p.x<0. || p.x>1./16. || p.y<0. || p.y>1./16.)? 
            vec4(0, 0, 0, 1e5) : texture(image49, p + i);
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;



    // Coordinates.
    vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;
    
    // Moving coordinates.
    vec2 p = r2(-3.14159/16.)*uv + vec2(0, 2)*TIME/64.;//r2(-3.14159/12.)*
     
    // Unit direction ray and light direction.
    vec3 rd = normalize(vec3(uv, 1)); 
    vec3 ld = normalize(vec3(-.5, .5, -3) - vec3(uv, 0));

    // Obtaining the asymmetric quad cell information.
    vec4 d4 = map(p);
    
    // Set the normal to the one calculated in the "map" function.
    vec3 sn = gN;
    
    // Smoothing factor.
    float sf = 1./RENDERSIZE.y;
    
    vec3 col = vec3(0);
    
    // Random graphite cell color.
    float rnd = hash21(gSID); // Random object value.
    vec3 oCol = vec3(rnd*rnd*.04 + .05);
    
    // Random color for blinking light cells.
    vec3 oColL = .5 + .45*cos(6.2831*rnd/4. + vec3(0, 1, 2)*1.2); // Orange.
    //vec3 oColL = .5 + .45*cos(6.2831*rnd/15. + vec3(0, 1, 2.2)*1.1 + 1.15);
    //vec3 oColL = .5 + .38*cos(6.2831*rnd/7.5 + vec3(0, 1, 2.2)/1.2 + .7);
    
    // Mixing in another color.
    oColL = mix(oColL, oColL.xzy, hash21(gSID + .42)*.25);
    
    #if COLOR == 1
    oColL = oColL.zyx; // Blue.
    #elif COLOR == 2
    oColL = mix(oColL, oColL.yxz*.8, .8); // Green.
    #elif COLOR == 3
    // Random compliment colors.
    oColL = mix(oColL.zyx, oColL, step(.5, hash21(gSID + .44)));
    #endif
    

    
    // Random lighter grey.
    if(hash21(gSID + .17)<.5) {
        oColL = vec3(.88, .9, 1)*(rnd*.12 + .18);//dot(oColL, vec3(.299, .587, .114));
    }
         
    // Random blinking colors. I overuse this trick a bit, but it's effective.    
    float rnd2 = hash21(gSID + .23);
    oCol = mix(oCol, oColL*1.5, smoothstep(.93, .96, sin(6.2831*rnd2 + TIME/4.)*.5 + .5));

    
    ///// FONTS /////
    // Adding text to the cells.
    vec2 sc = scale;
    vec2 tPos = p - d4.zw*sc - vec2(1, 2)*32./1024.*sc + vec2(sc.x/4., sc.y/2.4);
    int nc = int(floor(hash21(gSID)*93.) + 33.);
    float sD = .003;
    float bF = 2.;
    // Texture samples.
    vec4 dC = chTx(tPos, nc);
    vec4 dCX = chTx(tPos - vec2(sD, 0), nc);
    vec4 dCY = chTx(tPos - vec2(0, sD), nc);
    // Gradient indentation -- I was in a hurry, so there are probably neater
    // ways to do this... Probably a question for Fabrice. :)
    dC.x = smoothstep(0., .25, -(dC.w - .53));
    dCX.x = smoothstep(0., .25, -(dCX.w - .53));
    dCY.x = smoothstep(0., .25, -(dCY.w - .53));
    
    // Applying the font indentation bump to the normal. It's not perfect, but it'll do.
    float bf = .5;
    sn = normalize(sn + vec3((dC.x - dCX.x)*bf, (dC.x - dCY.x)*bf, 0)/sD*bF);
   
    // Applying a font overlay color to the font indentation.
    vec3 fCol = oCol/(.5 + dot(col, vec3(.299, .587, .114)));  
    float gr = dot(fCol, vec3(.299, .587, .114));
    fCol = max(vec3(1) - gr*gr*8., .0)*.9 + .05;
    //svC = mix(oColL, svC*, smoothstep(.93, .96, sin(6.2831*rnd2 + TIME/4.)*.5 + .5));
    oCol = mix(oCol, fCol, (1. - smoothstep(0., sf*3., dC.w - .51))*.9);
    //////////////
    
    
    // Diffuse and specular calculations. 
    float dif = max(dot(ld, sn), 0.); 
    float spe = pow(max(dot(reflect(ld, sn), rd), 0.), 8.);
    //float speR = pow(max(dot(normalize(ld - rd), sn), 0.), 16.);

    
    // Overall coloring and lighting for the cell.
    oCol *= (dif*.35 + .3 + vec3(1, .97, .92)*spe*.5);
    
    // Using the cube map for some fake specular reflections.
    //vec3 tx = texture(image47, reflect(rd, sn)).xyz; tx *= tx;
    //oCol += (oCol*.75 + .25)*mix(tx.zyx, tx.yyy, .3)*4.*spe;
 
    // Applying the cell to the background.
    col = mix(col, oCol, 1. - smoothstep(0., sf, d4.x));
 
 
  
    #if 1
    // Mix the previous frames in with no camera reprojection.
    // It's OK, but full temporal blur will be experienced.
    vec4 preCol = texelFetch(BuffA, ivec2(fragCoord), 0);
    float blend = (FRAMECOUNT < 2) ? 1. : 1./3.; 
    fragColor = mix(preCol, vec4(clamp(col, 0., 1.), d4.x), blend);
    #else
    // Output to screen.
    fragColor = vec4((max(col, 0.)), d4.x);
    #endif
	return fragColor; 
 } 


/*

    Random Asymmetric Quad Grid
    ---------------------------

    Rendering the buffer.
 
    See Buffer A for an explanation.
    
*/

// IQ's vec2 to float hash.
//float hash21(vec2 p){  return fract(sin(dot(p, vec2(27.619, 57.583)))*43758.5453); }

// This is an amalgamation of old blur and DOF functions with a couple of borrowed 
// lines from Dave Hoskins's much nicer Fibonacci based "Bokeh disc" function, which 
// you can find here: https://www.shadertoy.com/view/4d2Xzw
//
// This function is only really suitable for this example. If you're interested in 
// bokeh, Dave's function above and some of Shadertoy user, Hornet's, are probably
// the one's you should be looking at. Xor has some cool simple ones on here that I'm
// yet to dig into, but they might worth a look also.
vec4 bloom(sampler2D iCh, vec2 uv, float radius){

	vec4 tot = vec4(0);
    
    // UV based DOF. Focused on the horizontal line, then blurring further away.
    float r = smoothstep(0., 1., abs(uv.y - .57)/.57)*3.;
    
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
    col = bloom(BuffA, uv, 1.);

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