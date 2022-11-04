

			//******** Common Code Begins ********


// Display the individual hexagon grid cells.
//#define SHOW_HEX_GRID
    
// Flat shading. No lighting.
//#define FLAT_SHADING
 
// Apply hatching
#define HATCHING

// Tile scale.
#define GSCALE vec2(1./5.)*vec2(.5, .8660254)


//  Vertices and edge midpoints: Clockwise from the bottom left. -- Basically, the ones 
// above rotated anticlockwise. :)
vec2[6] vID = vec2[6](vec2(-.5, -1./3.)/vec2(.5, 1), vec2(-.5, 1./3.)/vec2(.5, 1), vec2(0, 2./3.)/vec2(.5, 1), 
                      vec2(.5, 1./3.)/vec2(.5, 1), vec2(.5, -1./3.)/vec2(.5, 1), vec2(0, -2./3.)/vec2(.5, 1));
//vec2[6] eID = vec2[6](vec2(-.5, 0)/vec2(.5, 1), vec2(-.25, .5)/vec2(.5, 1), vec2(.25, .5)/vec2(.5, 1), vec2(.5, 0)/vec2(.5, 1), 
                      //vec2(.25, -.5)/vec2(.5, 1), vec2(-.25, -.5)/vec2(.5, 1));

 
// Standard 2D rotation formula.
mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }


// IQ's vec2 to float hash.
float hash21(vec2 p){ 
    
    p = mod(p, 256.);
    // An annoying, but necessary, hack for systems with less sin
    // function accuracy. If anyone knows a way around it, feel 
    // free to let me know.
    //p = floor(p*1048576.)/1048576.;
    return fract(sin(dot(p, vec2(27.649, 57.583)))*43758.5453); 
}

/*
// UE4 random function: I like this because it incorporates a modulo
// 128 wrap, so in theory, things shouldn't blow up with increasing input.
// Also, in theory, you could tweak the figures by hand to get a really
// scrambled output... When I'm feeling less lazy, I might do that.
float hash21(vec2 p) {

    p = fract(p/128.)*128. - vec2(59.340627, 73.465623);
    return fract(dot(p.xyx*p.xyy, vec3(20.390625, 80.703127, 12.4281203)));
    
}
*/


/*
// This is IQ's WebGL 2 hash formula: I've modified it slightly to 
// take in the normal decimal floats that we're used to passing. It 
// works here, I think, but I'd consult the experts before using it.
//
// I remember reading through a detailed explanation of the C++ hash 
// we all used to use many years ago (which the following would be
// similar to), but have long since forgotten why it works. By the 
// way Nimitz and Dave Hoskins have formulae on Shadertoy that's worth
// looking up.
//
// Integer Hash - III - IQ
// https://www.shadertoy.com/view/4tXyWN
float hash21(vec2 p){
    
    uvec2 q = uvec2(ivec2(p*1023.));
	q = 1103515245U*((q>>1U)^q.yx);
    uint n = 1103515245U*(q.x^(q.y>>3U));
    return float(n)/float(0xffffffffU);
}
*/


// IQ's distance to a regular pentagon, without trigonometric functions. 
// Other distances here:
// https://iquilezles.org/articles/distfunctions2d
//
#define NV2 4
//
float sdPoly4(in vec2 p, in vec2[NV2] v){

    int num = v.length();
    float d = dot(p - v[0],p - v[0]);
    float s = 1.0;
    for( int i = 0, j = num - 1; i < num; j = i, i++){
    
        // distance
        vec2 e = v[j] - v[i];
        vec2 w =    p - v[i];
        vec2 b = w - e*clamp(dot(w, e)/dot(e, e), 0., 1. );
        d = min( d, dot(b,b) );

        // winding number from http://geomalgorithms.com/a03-_inclusion.html
        bvec3 cond = bvec3( p.y>=v[i].y, p.y<v[j].y, e.x*w.y>e.y*w.x );
        if( all(cond) || all(not(cond)) ) s*=-1.0;  
    }
    
    return s*sqrt(d);
}


// Determines which side of a line a pixel is on. Zero is the threshold.
float line(vec2 p, vec2 a, vec2 b){
     return ((b.x - a.x)*(p.y - a.y) - (b.y - a.y)*(p.x - a.x));
}

/*
// IQ's unsigned box formula.
float sBoxS(in vec2 p, in vec2 b, in float sf){

  vec2 d = abs(p) - b + sf;
  return min(max(d.x, d.y), 0.) + length(max(d, 0.)) - sf;
}
*/
// IQ's standard box function.
float sBox(in vec2 p, in vec2 b){
   
    vec2 d = abs(p) - b;
    return min(max(d.x, d.y), 0.) + length(max(d, 0.));
}

// This will draw a box (no caps) of width "ew" from point "a "to "b". I hacked
// it together pretty quickly. It seems to work, but I'm pretty sure it could be
// improved on. In fact, if anyone would like to do that, I'd be grateful. :)
float lBox(vec2 p, vec2 a, vec2 b, float ew){
    
    float ang = atan(b.y - a.y, b.x - a.x);
    p = rot2(ang)*(p - mix(a, b, .5));

    vec2 l = vec2(length(b - a), ew);
    return sBox(p, (l + ew)/2.) ;
}

 


/*

	Isometric Weave Pattern
    -----------------------

	Using a subdivided hexagon grid to produce a random weave, rendered
    in an impossible geometry style.

	Since I hadn't seen a randomized isometric weave before, I thought 
    it'd be interesting to programatically produce one. There's a bit of 
    code here, but a lot of it is window dressing. The principle behind
	its constuction is pretty simple:

    Partition space into a hexagonal grid. With each hexagonal cell, 
    partition it into three rhomboids which will each make up a cube face. 
    Choosing the SHOW_HEX_GRID define will make it more clear, since you 
    can see the three rhomboid partitionings in the background.

    With each rhomboid quad, render two straight lines between the mid 
    side edges, or one line over the other to form a cross. Cut each line 
    down the middle and assign a normal, depending upon whether it's 
    facing in the X, Y or Z directions.

    The process is very similar to a hexagonal Truchet weave, but with 
	the lines having a straight edge isometric quality. This is just one
	of many styles, and just one way to produce this particular 
    arrangement -- This is a weave, but you could include cross beams and
    so forth for way more variation. I'd imagine there'd be more elegant 
    solutions... which I'll leave to anyone bored enough to give it a go. :)



    Related examples:

    // Cool, clever and concise.
	isometric textured 3-story-map - FabriceNeyret2 
    https://www.shadertoy.com/view/WdsXW4

    // An oldschool isometric tiling example. One of my favorites on here.
    Isometric City 2.5D - knarkowicz
    https://www.shadertoy.com/view/MljXzz

    // I'm pretty sure BigWIngs was the first to put a hexagonal weave
    // on Shadertoy. He also has other weave examples that are worth
    // looking up.
    Hexagonal Truchet Weaving - BigWIngs
    https://www.shadertoy.com/view/llByzz

    // Awesome impossible geometry example... I probably should have 
    // consulted it before writing mine, as the code is much cleaner. :)
    Impossible Chainmail - BigWIngs
    https://www.shadertoy.com/view/td23zV
 
    // Another impossible geometry example done in a different way.
    Impossible Geometric Lattice - Shane
    https://www.shadertoy.com/view/wd3XRj

*/ 
 


// Hexagon cell routine. It's a shortened version of a 4-tap 3D routine, so it's
// not as consise as it should be. If you're after that, try this version:
//
// Minimal Hexagonal Grid - Shane
// https://www.shadertoy.com/view/Xljczw
vec4 hexCell(vec2 q){

    // Block dimension: Length to height ratio with additional scaling.
    //const vec2 dim = GSCALE;
    // A helper vector, but basically, it's the size of the repeat cell.
	const vec2 s = GSCALE*2.;
    
    // Two initial tile positions for a hexagon with a pointed top.
    //vec2[2] ps4 = vec2[2](vec2(-1, .5), vec2(0, -.5));
    vec2[2] ps4 = vec2[2](vec2(-.5, .5), vec2(.5, -.5));
  
    //Local coordinates of the closest hexagon cell and cell ID
    vec2 lP, id;
    
        // Distance.
    float minD = 1e5;
    
    // You definitely don't need a loop here, but as mentioned above, this
    // has been cut down from a 3D routine.
    for(int i = 0; i<2; i++){

        // Cell center.
        vec2 cntr = ps4[i]/2.; 
         
        vec2 p = q.xy; // Local coordinates.
        vec2 ip = floor(p/s - cntr) + .5; // Local tile ID.
        p -= (ip + cntr)*s; // New local position.
        
           
        // The hexagon cell distance or bound.... Kind of... Technically, it's 
        // "max(abs(p.x)*.8660256 + abs(p.y)*.5, abs(p.y))," but this will work.
        float d = dot(p, p);
        
        // If applicable, update the overall minimum distance value,
        // ID, and object ID. 
        if(d<minD){
            
            minD = d;
            // Setting the local coordinates and ID.
            lP = p;
            // Individual tile ID.
            id = ip + cntr; // Multiply by "s" for scaling.
        }
        
    }
    
    // Return the distance, position-based ID and triangle ID.
    return vec4(lP, id);
}


// A swap without the extra declaration, but involves extra operations -- 
// It works fine on my machine, but if it causes trouble, let me know. :)
#define swap(a, b){ a = a + b; b = a - b; a = a - b; }
/*
// Object swapping.
void swap(inout float a, inout float b){ float tmp = a; a = b; b = tmp; }
void swap(inout vec2 a, inout vec2 b){ vec2 tmp = a; a = b; b = tmp; }
void swap(inout vec3 a, inout vec3 b){ vec3 tmp = a; a = b; b = tmp; }
*/


vec3 colorSurf(vec3 rd, vec3 ld, vec3 sn){
   
    vec3 col = vec3(.6, .25, .1)*(max(dot(sn, ld), 0.) + .25);
    col += vec3(1, .6, .2)*pow(max(dot(reflect(sn, ld), rd), 0.), 4.); 
    //col += abs(sn)/2.;
    return col;
}

vec3 UVToXYZ(vec2 p) {
    return vec3(p.x + p.y, p.x, p.y);
}
 
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // Screen coordinates.
	vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;
    
     // Global scale. Rendundant here.
    const float gSc = 1.;
    // Smoothing factor.
    float sf = gSc/RENDERSIZE.y;
    
    // Scaled than transalated coordinates.
    vec2 p = uv*gSc + vec2(0, .5)*TIME/8.;
  
    // The hexagon information for the current cell. This returns the
    // local coordinates (h4.xy), and the unique central position ID (h4.zw).
    vec4 h4 = hexCell(p);
    vec2 svP = h4.xy; // Local coordinates.
    
    vec2 ap = abs(svP);
    float hexDist = max(ap.y*.8660256 + ap.x*.5, ap.x) - (GSCALE).x;
    
    // The offset vertex information.
    // Hexagon vertices with scaling to enable rendering back in normal space. 
    vec2 dim = GSCALE;
    vec2[6] svV = vec2[6](vID[0]*dim, vID[1]*dim, vID[2]*dim, vID[3]*dim, vID[4]*dim, vID[5]*dim);

    
    // Containers for the midpoint edges and corresponding
    // normals for the string and string shadows.
    vec4 svEP[6];
    
    
    
    // The two line distance fields.
    float ln = 1e5, ln2 = 1e5;
    // Line edge width.
    #ifdef SHOW_HEX_GRID
    // Thinner lines when the hexagonal cell grid is displayed, just to 
    // make it easier to see things.
    const float ew = .02;
    #else
    const float ew = .03;
    #endif
    
    // The quad distace field. There are 3 rhomboidal quads per hexagon cell.
    float qDist = 1e5;
    // The minimum quad ID, based on the minimum quad central position.
    vec2 qID = vec2(0);
    // The box quad ID.
    int bID = 0;
    
    
    // The face surface normals. The first is for upward facing surfaces, and the others are
    // for left and right facing. As an aside, it'd be possible to code hexagon, octagon
    // surface cross sections too, but that'd be pretty painful, so I'm sticking to squares. :)
    vec3 nr1 = vec3(0, 1, 0), nr2 = normalize(vec3(0, 0, -1)), nr3 = normalize(vec3(1, 0, 0));
    //swap(nr1, nr2); // Swapping face normals.
    // The two surface normals -- Initially to the left edge, but probably don't need to be.
    vec3 sn = nr2, sn2 = nr2; 
   
    
    // For each hexagon, partition into three rhomboids then render two lines through 
    // each whilst assigning normals.
    
    for(int j = 0; j<6; j+=2){
        
        // Constructing the edge midpoints and normals at those
        // points for both the string and corresponding shadows.
        vec2[4] v = vec2[4](svV[(j + 1)%6], svV[(j + 2)%6], svV[(j + 3)%6], vec2(0));
 
        // Quad center and local quad ID.
        vec2 ctr = (v[0] + v[1] + v[2] + v[3])/4.;
        vec2 lID = h4.zw + ctr;
        
        // Some random numbers for this particular quad.
        float rndI = hash21(lID + .327);
        float rndI2 = hash21(lID + .493);
       
        // Mid edge points. You can offset the lines from the mid point
        // slightly (<.25) to create warped patterns. You could also 
        // offset on an edge by edge basis, but I wanted to keep things
        // simple.
        float offs = 0.;//sin(TIME)*.05;
        vec2 e0 = mix(v[0], v[1], .5 + offs);
        vec2 e1 = mix(v[1], v[2], .5 - offs);
        vec2 e2 = mix(v[2], v[3], .5 + offs);
        vec2 e3 = mix(v[3], v[0], .5 - offs);
     

        float quad = sdPoly4(svP, v);
        
        // If this particular quad is closer, calculate some lines. By the way, 
        // the "quad<1e-3" line is there, so that we mostly calculate the 
        // following just once. Normally, you'd just find the nearest quad, 
        // save the edge points and perform this outside the loop, but I'm trying
        // to keep things more compact.
        if(quad<qDist && quad<1e-3){
            
            
            qDist = quad; // Minimum distance.
            qID = lID; // Minimum quad ID.
            bID = j/2; // Cube or box ID -- There are 3 in all.
 
            // Swapping end points is a way to switch from two parallel quad lines
            // to ones that cross over.
            if(rndI<.5) { 
                 swap(e2, e3); 
            }
            
            // Line one and two tangents.
            vec2 tn = normalize(e2 - e0);
            vec2 tn2 = normalize(e3 - e1);
            
            // The distance fields for line one and two. Note that I've extended the
            // mid edge points beyond the quad boundaries. This will get rid of
            // border seam lines.
            const float ndg = .5; // Change to "-.02," or something to shorten the lines.
            ln = lBox(svP, e0 - tn*ndg, e2 + tn*ndg, ew);
            ln2 = lBox(svP, e1 - tn2*ndg, e3 + tn2*ndg, ew); 
            
            // Determine which side of the center line we're on for each
            // line segment.
            float line02 = line(svP, e0, e2);
            float line13 = line(svP, e1, e3);
            
            // Setting the face normal: This is dependent on the quad we're
            // in (top or two sides) which side of the central line of each
            // line segment we're in and the normal orientation (up, sloping
            // down, etc). By the way, there's probably a more elegant way
            // to go about setting the face normals, but these calculations
            // are only performed once, so it doesn't matter too much... I'll
            // leave the more elegant solution to Fabrice. :)
            //
            if(j==0){ // Top rhomboidal quad.
            
                 // Straight up and down line segment one.
                if(abs(tn.x)<.1) {
                    // Left side of the line, X-dominant normal.
                    if(line02<0.) sn = nr3;
                    else sn = nr2; // Right side of the line, Z-dominant normal.
                }
                else{ // Left or right line one.
                    if(line02<0.) sn = nr3; // Bottom side of the line, X-dominant normal.
                    else sn = nr1; // Top side of the line, Y-dominant normal.
                }
                
                // Straight up and down line segment two.
                if(abs(tn2.x)<.1) {
                    if(line13<0.) sn2 = nr3; // Left side of the line, X-dominant normal.
                    else  sn2 = nr2; // Right side of the line, Z-dominant normal.
                }
                else{  // Left or right line two.
                    if(line13<0.) sn2 = nr1;  // Top side of the line, Y-dominant normal.
                    else sn2 = nr2;  // Bottom side of the line, Z-dominant normal.
                }
                
            }
             
            if(j==2){  // Bottom right rhomboidal quad.
            
                // Similar logic to above (See "j==0").
                if(abs(tn.x)<.2) {
                    
                    if(line02<0.) sn = nr2;
                    else sn = nr3;
                }
                else{
                    
                    if(line02<0.) sn = nr1;
                    else {                        
                        if(tn.y>.0) sn = nr3;
                        else sn = nr2;
                    }
                }
                
                if(abs(tn2.x)<.2) {
                    if(line13<0.) sn2 = nr2;
                    else sn2 = nr3;
                }
                else{
                    if(line13<0.) sn2 = nr1;
                    else sn2 = nr3;
                }
                
            }  
            
            
             if(j==4){ // Bottom right rhomboidal quad.
            
                // Similar logic to above (See "j==0").
                if(abs(tn.x)<.2){
                    if(line02<0.) sn = nr2;
                    else sn = nr3;
                }
                else{
                    if(line02<0.) sn = nr2;
                    else sn = nr1;
                }
                 
                
                if(abs(tn2.x)<.2) {
                    if(line13<0.) sn2 = nr2;
                    else sn2 = nr3;
                }
                else{
                    if(line13<0.){
                        if(tn2.y>.0) sn2 = nr2;
                        else sn2 = nr3;
                    }
                    else sn2 = nr1;
                }
                
            }           
 
            // Randomly swapping the rendering order of the lines for
            // more variation.
            if(rndI2<.5){
                swap(ln, ln2);
                swap(sn, sn2);
            }
 
            
            
        }
 
    }
    
    
    // Rendering.
 
    // At this point, you have access to the two line distance field values, and
    // their respective normals. You also have the quad distance field values and
    // normals, so at this point, you can render whatever you want in as simple
    // or as complex a manner as you'd like.
    
    
    // Unit direction ray.
    vec3 rd = normalize(vec3(uv, 1));
    // Fake point light... Very fake, but it's only used to create a bit of
    // gradient shine within the scene, so it doesn't matter.
    vec3 ld = normalize(vec3(.5, 3, -3) - UVToXYZ(uv));
    //vec2 xy = rot2(3.14159/3.)*uv;
    //vec3 ld = normalize(vec3(.5, 3, -3) - vec3(xy.x, uv.y, xy.y));
    vec3 ld2 = -ld.zxy; // Another light.

    
        
    // Initializing the scene color to the background color.
    vec3 col = vec3(1);
    
    // The two line distance field colors.
    vec3 lnCol = vec3(1);
    vec3 lnCol2 = vec3(1);
    
    
    
    // Face IDs. One for each pylon, or line.
    int faceID1 = abs(sn.x)>.8? 0 : abs(sn.y)>.8? 1 : 2;
    int faceID2 = abs(sn2.x)>.8? 0 : abs(sn2.y)>.8? 1 : 2;
 
    #ifdef FLAT_SHADING
        // Flat shading option.
    
        // The three face colors. I've arranged for the lightest side to face up.
        const vec3 cc = vec3(1, .8, .55);
        //vec3 aCol[3] = vec3[3](cc*cc*cc*.5, cc, cc*cc*.8);
        //vec3 aCol[3] = vec3[3](vec3(1, .15, .3), vec3(1, .5, .85), vec3(.6, .1, 1));
        vec3 aCol[3] = vec3[3](vec3(1, .1, .01), vec3(1, .6, .1), vec3(1, .2, .0));
        // The two line colors.
        lnCol = aCol[faceID1];
        lnCol2 = aCol[faceID2];
    #else    
        // Using the face normal to lighting the surface.
    
        lnCol = colorSurf(rd, ld, sn); 
        lnCol2 = colorSurf(rd, ld, sn2); 
        //lnCol = vec3(1)*dot(lnCol, vec3(.299, .587, .114));
        //lnCol2 = vec3(1)*dot(lnCol2, vec3(.299, .587, .114));
        //lnCol += pow(colorSurf(rd, ld2, sn).zyx, vec3(2))/16.; 
        //lnCol2 += pow(colorSurf(rd, ld2, sn2).zyx, vec3(2))/16.; 
    #endif

    // The cross hatching pattern on each face. Oriented to match the face.
    vec3 aR = vec3(-3.14159/3., 0, 3.14159/3.) + 3.14159/1.;
    vec2 qUV = rot2(aR[faceID1])*p;
    vec2 qUV2 = rot2(aR[faceID2])*p;
    float freq = 150.;
    float pat = abs(fract(qUV.y*freq) - .5)*2. - .1;
    float pat2 = abs(fract(qUV2.y*freq) - .5)*2. - .1;
    pat = smoothstep(0., sf*freq*2., pat);
    pat2 = smoothstep(0., sf*freq*2., pat2);
    #ifdef HATCHING
    // Applying the hatch pattern to each line.
    lnCol *= pat*.4 + .7;
    lnCol2 *= pat2*.4 + .7;
    #endif
    
    // The background cubes consist of the three face quads. The top quad has
    // an upward facing normal, and the sides have X or Z faceing normals. 
    vec3 cN = bID == 0? nr3 : bID == 1? nr2 : nr1;

    // The cube cross hatch pattern.
    freq = 150.;
    vec2 qUV3 = rot2(aR[(bID + 1)%3])*p;
    float pat3 = abs(fract(qUV3.y*freq) - .5)*2. - .1;
    pat3 = smoothstep(0., sf*freq*2., pat3);
    
    #ifndef FLAT_SHADING
    //vec3 cCol = colorSurf(rd, ld, cN);
    vec3 cCol = vec3(1, .95, .9);
    // Add some blue light to the cube faces.
    cCol += pow(colorSurf(rd, ld2, cN).zyx, vec3(2))/2.; 
    #else 
    vec3 cCol = vec3(.25, .2, .15);
    col = cCol;
    #endif
    
    // Random blinking.
    float rnd = hash21(qID + .13);
    rnd = smoothstep(.9, .95, sin(6.2831*rnd + TIME*2.)*.5 + .5);
    // Shading the background with one of the normals.
    float sh3 = dot(colorSurf(rd, ld, nr1), vec3(.299, .587, .114));
    
    // Hexagonal cell background shade.
    //float sh = max(.5 - qDist/.1, 0.);
    float sh = max(.5 - hexDist/.1, 0.);
    
    #ifdef FLAT_SHADING
    sh = 1., sh3 = 1.; rnd = 0.;
    #endif
    // Border highlight color.
    vec3 hiCol = cCol*vec3(.3, .4, 1)*4.*sh3;
   
    // Applying shade and hatching to the cube background.
    cCol *= sh3*sh;
    #ifdef HATCHING
    cCol *= (pat3*.4 + .7);
    #endif
   
    // Rendering the hexagonal quad segment in the background -- Dark edge and coloring.
    col = mix(col, cCol, (1. - smoothstep(0., sf, qDist + .0035)));
    #ifndef FLAT_SHADING
    col = mix(col, vec3(0), (1. - smoothstep(0., sf*6., abs(qDist) - .0035/2.*6.))*.5);
    col = mix(col, vec3(0), (1. - smoothstep(0., sf, abs(qDist) - .0035/2.*6.))); 
    col = mix(col, mix(vec3(sh3), hiCol, rnd), (1. - smoothstep(0., sf, abs(qDist) - .0035/2.*3.))); 
    col = mix(col, vec3(0), (1. - smoothstep(0., sf, abs(qDist) - .0035/2.)));
    #else
    //col = mix(col, vec3(0), (1. - smoothstep(0., sf, abs(qDist) - .0035/2.*5.)));
    //col = mix(col, vec3(.25, .2, .15)*2., (1. - smoothstep(0., sf, abs(qDist) - .0035/2.*3.)));
    col = mix(col, vec3(0), (1. - smoothstep(0., sf, abs(qDist) - .0035/2.)));
    #endif
   
    // Shading the lines.
    sh = max(.5 - ln/.03, 0.);
    float sh2 = max(.5 - ln2/.03, 0.);
    
    #ifdef FLAT_SHADING
    sh2 = 1.;
    #endif

     
    // Rendering the first line layer -- Drop shadow, edges, color and center line.
    col = mix(col, vec3(0), (1. - smoothstep(0., sf*10., ln))*.5);
    col = mix(col, vec3(0), 1. - smoothstep(0., sf, ln));
    //col = mix(col, vec3(.5)*sh, 1. - smoothstep(0., sf, ln + .005));
    //col = mix(col, vec3(0), 1. - smoothstep(0., sf, ln + .005 + .0035));
    col = mix(col, lnCol*sh, 1. - smoothstep(0., sf, ln + .005)); // + .005 + .0035
    col = mix(col, vec3(0), 1. - smoothstep(0., sf, ln + ew - .0015));
 
    // Doing the same for the second line layer 
    col = mix(col, vec3(0), (1. - smoothstep(0., sf*10., ln2))*.5);
    col = mix(col, vec3(0), 1. - smoothstep(0., sf, ln2)); 
    //col = mix(col, vec3(.5)*sh2, 1. - smoothstep(0., sf, ln2 + .005));
    //col = mix(col, vec3(0), 1. - smoothstep(0., sf, ln2 + .005 + .0035));
    col = mix(col, lnCol2*sh2, 1. - smoothstep(0., sf, ln2 + .005)); // + .005 + .0035
    col = mix(col, vec3(0), 1. - smoothstep(0., sf, ln2 + ew - .0015));
 
 
    #ifdef SHOW_HEX_GRID
    col = mix(col, vec3(0), (1. - smoothstep(0., sf*8., abs(hexDist) - .005/2.*3.))*.5); 
    col = mix(col, vec3(0), (1. - smoothstep(0., sf, abs(hexDist) - .005/2.*3.))); 
    col = mix(col, vec3(1, .95, .9), (1. - smoothstep(0., sf, abs(hexDist) - .005/2.))); 
    #endif
 
    // Mixing in a bit of gradient color.
    col = mix(col, col.xzy, length(uv)*.25);
    
    //col = mix(col.zyx, col.yxz, .8); // Other colors.
    
    // Applying a subtle silhouette, for art's sake.
	uv = fragCoord/RENDERSIZE.xy;
    col *= pow(16.*(1. - uv.x)*(1. - uv.y)*uv.x*uv.y, 1./16.)*1.05;
                          

    // Rough gamma correction, and we're done.
	fragColor = vec4(sqrt(max(col, 0.)), 1);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}