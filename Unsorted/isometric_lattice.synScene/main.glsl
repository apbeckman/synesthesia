

/*

	Isometric Lattice
    -----------------

    I thought it'd be fun to quickly code up a standard impossible geometry 
    lattice with a couple of extra elements -- In theory, it involves
	strategically placing some overlapping quads inside the cells of a 
    triangle grid, which should be pretty easy... for a halfway competent
	coder, which I used to be, but I struggle with everything these days. :D

	The code was slapped together with a view to getting the job done, so 
	you can safely ignore it as I'd imagine there'd be way more elegant
	ways to produce the same. 

	Anyway, this is just one of countless impossible geometry examples out 
    there. It'd be great to see others on Shadertoy. All involve placing 2D 
    elements in various ways on a grid or the whole canvas.




    Other examples:


    // Oldschool isometric maze. I have a version of this somewhere.
    Isometric Maze - fizzer
    https://www.shadertoy.com/view/Md2XRd

	// Flopine has a few interesting isometric examples worth looking at.
    // This particular one was coded in a few minutes... It took me that
	// long just to decide what colors I wanted to use. :)
	Flopine - I'm trapped 
	https://www.shadertoy.com/view/WtsfWr


*/


// Triangle grid cell borders. If you like spoiling illusions like I do,
// this is the command to comment in. :)
//#define SHOW_CELLS

// Display T-joins instead of cube joins.
//#define TJOINS

// A singe join, instead of four. It produces a cleaner design, but a less 
// interesting one, I believe.
//#define ONE_JOIN

// Display the mid edge green block.
#define MID_BLOCK
    


// A swap without the extra declaration, but involves extra operations -- 
// It works fine on my machine, but if it causes trouble, let me know. :)
#define swap(a, b){ a = a + b; b = a - b; a = a - b; }


// Standard 2D rotation formula.
mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }


// IQ's vec2 to float hash.
float hash21(vec2 p){  return fract(sin(dot(p, vec2(27.609, 57.583)))*43758.5453); }


// Compact, self-contained version of IQ's 2D value noise function.
float n2D(vec2 p){
   
    // Setup.
    // Any random integers will work, but this particular
    // combination works well.
    const vec2 s = vec2(1, 113);
    // Unique cell ID and local coordinates.
    vec2 ip = floor(p); p -= ip;
    // Vertex IDs.
    vec4 h = vec4(0., s.x, s.y, s.x + s.y) + dot(ip, s);
   
    // Smoothing.
    p = p*p*(3. - 2.*p);
    //p *= p*p*(p*(p*6. - 15.) + 10.); // Smoother.
   
    // Random values for the square vertices.
    h = fract(sin(h)*43758.5453);
   
    // Interpolation.
    h.xy = mix(h.xy, h.zw, p.y);
    return mix(h.x, h.y, p.x); // Output: Range: [0, 1].
}

// FBM -- 4 accumulated noise layers of modulated amplitudes and frequencies.
float fbm(vec2 p){ return n2D(p)*.533 + n2D(p*2.)*.267 + n2D(p*4.)*.133 + n2D(p*8.)*.067; }


vec3 pencil(vec3 col, vec2 p){
    
    // Rough pencil color overlay... The calculations are rough... Very rough, in fact, 
    // since I'm only using a small overlayed portion of it. Flockaroo does a much, much 
    // better pencil sketch algorithm here:
    //
    // When Voxels Wed Pixels - Flockaroo 
    // https://www.shadertoy.com/view/MsKfRw
    //
    // Anyway, the idea is very simple: Render a layer of noise, stretched out along one 
    // of the directions, then mix similar, but rotated, layers on top. Whilst doing this,
    // compare each layer to it's underlying greyscale value, and take the difference...
    // I probably could have described it better, but hopefully, the code will make it 
    // more clear. :)
    // 
    // Tweaked to suit the brush stroke size.
    vec2 q = p*4.;
    const vec2 sc = vec2(1, 12);
    q += (vec2(n2D(q*4.), n2D(q*4. + 7.3)) - .5)*.03;
    q *= rot2(-3.14159/2.5);
    // I always forget this bit. Without it, the grey scale value will be above one, 
    // resulting in the extra bright spots not having any hatching over the top.
    col = min(col, 1.);
    // Underlying grey scale pixel value -- Tweaked for contrast and brightness.
    float gr = (dot(col, vec3(.299, .587, .114)));
    // Stretched fBm noise layer.
    float ns = (n2D(q*sc)*.66 + n2D(q*2.*sc)*.34);
    //
    // Repeat the process with a couple of extra rotated layers.
    q *= rot2(3.14159/2.); q += 2.5;
    float ns2 = (n2D(q*sc)*.66 + n2D(q*2.*sc)*.34);
    q *= rot2(-3.14159/5.);  q += 4.5;
    float ns3 = (n2D(q*sc)*.66 + n2D(q*2.*sc)*.34);
    //
    // Compare it to the underlying grey scale value.
    //
    // Mix the two layers in some way to suit your needs. Flockaroo applied common sense, 
    // and used a smooth threshold, which works better than the dumb things I was trying. :)
    const float contrast = 1.;
    ns = (.5 + (gr - (max(max(ns, ns2), ns3)))*contrast); // Same, but with contrast.
    //ns = smoothstep(0., 1., .5 + (gr - max(max(ns, ns2), ns3))); // Different contrast.
    // 
    // Return the pencil sketch value.
    return vec3(clamp(ns, 0., 1.));
    
}


// IQ's distance to a regular polygon, without trigonometric functions. 
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



// IQ's 2D box function.
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
   return sBox(p, (l + ew)/2.);
    
}

 
// Skewing coordinates. "s" contains the X and Y skew factors.
vec2 skewXY(vec2 p, vec2 s){ return mat2(1, -s.y, -s.x, 1)*p; }


// Unskewing coordinates. "s" contains the X and Y skew factors.
vec2 unskewXY(vec2 p, vec2 s){ return inverse(mat2(1, -s.y, -s.x, 1))*p; }



const float scale = 1.;
// Rectangle stretch.
const vec2 rect = vec2(1, 1./.8660254)*scale; 
//const vec2 rect = vec2(.85, 1.15)*scale; 
// Skewing half way along X, and not skewing in the Y direction.
const vec2 sk = vec2(0, rect.y*.5); // 12 x .2
// Irregular skewing is possible too, since it's all just math.
//const vec2 sk = vec2(rect.x*.4, rect.y*.2); // 12 x .2


/*
//#define FLAT_TOP

#ifdef FLAT_TOP
//vec2 dim = vec2(1.5, 1)*scale;
//const vec2 rect = vec2(1., 1.5)*scale; // 12 x .2
const vec2 rect = (vec2(1./.8660254, 1))*scale;

// Skewing half way along X, and not skewing in the Y direction.
const vec2 sk = vec2(rect.x*.5, 0); // 12 x .2
#else
const vec2 rect = (vec2(1, 1./.8660254))*scale; 

// Skewing half way along X, and not skewing in the Y direction.
const vec2 sk = vec2(0, rect.y*.5); // 12 x .2

#endif
*/



float gTri;
vec4 getTri(vec2 p){
    
     p = skewXY(p, sk);
    
    // Unique position-based ID for each cell. Technically, to get the central position
    // back, you'd need to multiply this by the "rect" variable, but it's kept this way
    // to keep the calculations easier. It's worth putting some simple numbers into the
    // "rect" variable to convince yourself that the following makes sense.
	vec2 id = floor(p/rect) + .5; 
    // Local grid cell coordinates -- Range: [-rect/2., rect/2.].
	p -= id*rect; 
    
    
    // Equivalent to: 
    //float tri = p.x/rect.x < -p.y/rect.y? 1. : 0.;
    // Base on the bottom (0.) or upside down (1.);
    gTri = dot(p, 1./rect)<0.? 0. : 1.;
   
    p = unskewXY(p, sk);
    
    return vec4(p, id);
    
    
}


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    // SETUP.
    
    // Aspect correct pixel coordinates.
    vec2 uv = fragCoord/RENDERSIZE.y;
    
    // Scale, or zoom variable, if you prefer.
    const float sc = 4.;
    
    // Scaling and translation.
    vec2 p = uv*sc - vec2(1)*TIME/2.;
    
    // Perturbing the coordinates for imperfect hand drawn lines.
    vec2 offs = vec2(fbm(p*16.), fbm(p*16. + .35));
    const float oFct = .01;
    p -= (offs - .5)*oFct;
    
    // Copy of the original coordinates.
    vec2 oP = p;
    
    // Scale based smoothing factor. Good as a general scene formula, but for trickier
    // fields, derivatives, etc, should be taken into account.
    float sf = sc/RENDERSIZE.y;


    
    // Shadow coordinate, ID and triangle orientation id.
    vec2 shOff = normalize(vec2(1, -1.25))*.175;
    vec4 p4Sh = getTri(p - shOff);
    vec2 pSh = p4Sh.xy;
    vec2 idSh = p4Sh.zw;
    float triSh = gTri;
    
    // Cell coordinate, ID and triangle orientation id.
    vec4 p4 = getTri(p);
    p = p4.xy;
    vec2 id = p4.zw;
    float tri = gTri;

    
    
    // Using the unique cell ID to produce a random number.
    float rnd = hash21(id);
     
   
    // DISTANCE FIELD CALCULATION.
    
       
    // Grid vertices, clockwise from the bottom left:
    //
    //vec2[4] vert = vec2[4](vec2(-rect.x, -rect.y)/2., vec2(-rect.x, rect.y)/2.,
                            //vec2(rect.x, rect.y)/2., vec2(rect.x, -rect.y)/2.));
    
    // However, what you'll note is that the above can be written in the following
    // manner.
    //
    // These are just the vertices of a unit grid cell, and when coupled with the
    // cell ID, can be used for all kinds of things. You can also construct
    // midway points between vertices for unique edge identifiers, which enables
    // the construction of weaves, jigsaw patterns, etc.
    //
    // By the way, and this can get confusing, all IDs in this example will work,
    // and will be unique. However, if you want them be position bases, like reading
    // into a heightmap or whatever, all IDs need to be multiplied by the cell
    // dimensions, "rect."
    vec2[4] vertID = vec2[4](vec2(-.5, .5), vec2(.5), vec2(.5, -.5), vec2(-.5));
    
    // The rectangular grid vertices.
    vec2[4] vert = vec2[4](vertID[0]*rect, vertID[1]*rect, vertID[2]*rect, vertID[3]*rect);
    
    for(int i = 0; i<4; i++) vert[i] = unskewXY(vert[i], sk); // Unskew.
    
    
   
 
    
    // RENDERING.
    
     // Using the unique ID to give the cells a unique background color
    // via IQ's versatile palette formula.
    vec3 bg = vec3(.3, .6, 1);
    //bg = .6 + .4*cos(6.283*hash21(id + tri)/4. + vec3(0, 1, 2)); // Random cell color.
    
    /*
    // Background lines.
    vec2 qUV2 = rot2(-3.14159/3.)*p;
    float freq2 = 30.;
    float pat2 = abs(fract(qUV2.y*freq2) - .5)*2.;
    pat2 = smoothstep(0., sf*freq2*2., pat2 - .2);    
    bg *= pat2*.5 + .7;
    */
    
    // Set the scene to the cell background.
    vec3 col = bg;
    
    /*==
    // Color the upside down triangles the opposing color.
    if(tri<.5) {
        col = 1. - col*.5;
    }
    */
    
 
    
    // Cell vertices and vertex ID.
    vec2[3] v, vID;
    //vec2 ctr = vec2(0);

    
    if(tri>=.5){
        vID = vec2[3](vertID[0], vertID[2], vertID[1]);
    	v = vec2[3](vert[0], vert[2], vert[1]);
        //ctr = (v[0] + v[1] + v[2])/3.;
    }
    else {
        vID = vec2[3](vertID[0], vertID[2], vertID[3]);
        v = vec2[3](vert[0], vert[2], vert[3]); 
        //ctr = (v[0] + v[1] + v[2])/3.;
    }
    
      
    // Shadow vertices -- They need to be handled seperately... Sigh! :)
    vec2[3] vSh, vIDSh;
    
    if(triSh>=.5){
        vIDSh = vec2[3](vertID[0], vertID[2], vertID[1]);
    	vSh = vec2[3](vert[0], vert[2], vert[1]);
    }
    else {
        vIDSh = vec2[3](vertID[0], vertID[2], vertID[3]);
        vSh = vec2[3](vert[0], vert[2], vert[3]);
    }    
    
    
     

    // Assigned face normals, colors, and pattern.
    vec3 vN[3] = vec3[3](vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1));
    vec3 vCol[3] = vec3[3](vec3(.5), vec3(1), vec3(.2));
    vec3 aR = vec3(-3.14159/3., 0, 3.14159/3.);
    vec3 vPat = vec3(1);
    
    /*
    // Cross hatching pattern, depending on face orientation -- Not used here.
    for(int i = 0; i<3; i++){
        // The cross hatching pattern on each face. Oriented to match the face.
        vec2 qUV = rot2(aR[i])*p;
        float freq = 30.;
        float dir = tri<.5? -1. : 1.; 
        
        if(i==2) dir = 0.;

        float pat = abs(fract(qUV.y*freq + dir*.25) - .5)*2.;
        pat = min(pat, abs(fract(qUV.x*freq + dir*.5) - .5)*2.);
        //pat = smoothstep(0., sf*freq*2., pat - .2);
        // Applying the hatch pattern to each line.
        vPat[i] = 1.;;pat*.4 + .75;
        
        //if(i==2) vCol[i] *= 4.;
       
    }
    */
    
    /*
    // Swapping face normals.
    swap(vN[0], vN[1]);
    swap(vCol[0], vCol[1]); 
    swap(vPat[0], vPat[1]);
    //swap(vCol[1], vCol[2]);
    
    // Face normals, colors, etc, face in different directions for opposing grid
    // triangles... It can get confusing, but that's triangle grids for you. :)
     if(tri<.5) {
        swap(vN[0], vN[1]);
        swap(vCol[0], vCol[1]);
        swap(vPat[0], vPat[1]);
    } 
    */
    
    // Does the same as above, but only for this particular examples.
    if(tri>.5) {
        swap(vN[0], vN[1]);
        swap(vCol[0], vCol[1]);
        swap(vPat[0], vPat[1]);
    }
     
    // Background triangle cell lines.
    for(int i = 0; i<3; i++){
         float lnBord = lBox(p, v[i], v[(i + 1)%3], .00625);
         col = mix(col, vec3(0), (1. - smoothstep(0., sf, lnBord)));
     
    }   
    
   
    
    float oLn = 1e5;
    float ln = 1e5, ln2 = 1e5;
    
    // Main block width.
    float hexW = scale*.8660254/2.; //(.8660254*3.);//
 
    // Middle block width.
    float mBW = hexW/6.;
    
    // Joiner line width and edge width.
    #ifdef ONE_JOIN
    const float lw = .8660254/8.; 
    #else
    const float lw = .8660254/12.; 
    #endif
    const float ew = .02;
    
    #ifdef ONE_JOIN
    float ndg = hexW/((3. + lw)*.8660254);
    #else
    float ndg = hexW/((5. + lw)*.8660254);
    #endif
    
    
    // Drop shadow distance.
    float gSh = 1e5;
    
    
    // Render the drop shadows.
    for(int i = min(0, int(FRAMECOUNT)); i<3; i++){
        
        // Random edge variable.
        float rndI = hash21(idSh + mix(vIDSh[i], vIDSh[(i + 1)%3], .5));

        // Edge tangent vectors.
        vec2 tgnt0Sh = normalize(vSh[i] - vSh[(i + 1)%3]);
        vec2 tgnt1Sh = normalize(vSh[(i + 1)%3] - vSh[(i + 2)%3]);
        vec2 tgnt2Sh = normalize(vSh[(i + 2)%3] - vSh[(i + 0)%3]);
        
         // Central join line vertices.        
        vec2[4] qvSh = vec2[4](vSh[i] - tgnt0Sh*ndg, vSh[(i + 1)%3] + tgnt0Sh*ndg, 
                             vSh[(i + 1)%3] + tgnt0Sh*ndg + tgnt2Sh*lw, vSh[i] - tgnt0Sh*ndg + tgnt2Sh*lw);
        
        // Edge join line vertices.
        vec2[4] qvBSh = vec2[4](v[i] - tgnt0Sh*ndg, v[(i + 1)%3] + tgnt0Sh*ndg, 
                            v[(i + 1)%3] + tgnt0Sh*ndg + tgnt1Sh*lw, v[i] - tgnt0Sh*ndg + tgnt1Sh*lw);
        
        // Central join line.
        ln = sdPoly4(pSh, qvSh);
        
        // Outside joins line - Top.
        float ln2Offs = hexW/3.;
        float ln2 = sdPoly4(pSh - tgnt2Sh*ln2Offs, qvBSh);///////
        // Outside joins line - Side.
        float ln2B = sdPoly4(pSh - tgnt2Sh*ln2Offs, qvBSh);//////
        
        // If more than one join, join the line shadows.
        #ifndef ONE_JOIN
        ln = min(ln, min(ln2, ln2B));
        #endif
        
        // Shadow main block vertices.
        qvSh = vec2[4](vSh[i], vSh[i] + tgnt2Sh*hexW, 
                    vSh[i] + tgnt2Sh*hexW + tgnt1Sh*hexW, vSh[i] + tgnt1Sh*hexW);
        
        // Main block.
        float bl = sdPoly4(pSh, qvSh);
        
        // Main block T-Joins -- Shadow.
        #ifdef TJOINS
        // T - intersection.
        //if((i==1 && triSh>.5)|| (i==2 && triSh<.5)){// 
            vec2 ap = abs(pSh - vSh[i] + tgnt0Sh*.7);
            float blN = max(ap.y*.8660254 + ap.x*.5, ap.x) - hexW*.8660254;
            bl = max(bl, -blN);
        //}
        #endif  
       
        // Middle block.
        float blockDist = (1. - hexW/length(vSh[i] - vSh[(i + 1)%3]))/2.;// + mBW;
        vec2 eM = mix(vSh[i], vSh[(i + 1)%3], blockDist);
        blockDist *= length(vSh[i] - vSh[(i + 1)%3]);
        float hexW2 = hexW;//*.85;
        
        // Middle block vertices.
        qvSh = vec2[4](vSh[i], vSh[i] + tgnt2Sh*hexW2, 
                     vSh[i] + tgnt2Sh*hexW2 + tgnt1Sh*hexW2, vSh[i] + tgnt1Sh*hexW2);
        
        // Middle block face.
        float bl2B = sdPoly4(pSh + tgnt0Sh*(blockDist + mBW), qvSh);
        
      
        // Middle block edge vertices.
        qvSh = vec2[4](eM - tgnt0Sh*mBW, eM + tgnt0Sh*mBW, 
                     eM + tgnt0Sh*mBW + tgnt2Sh*hexW2, eM - tgnt0Sh*mBW + tgnt2Sh*hexW2);
        
        // Middle block edge distance.
        float bl2 = sdPoly4(pSh, qvSh);
        
        // Don't show various random edge shadows.
        if(rndI<1./2.){
            bl2 = 1e5;
            bl2B = 1e5;
            ln = 1e5; 
        }
        
        // If we're not displaying the middle green block.
        #ifndef MID_BLOCK
        bl2 = bl2B = 1e5;
        #endif
               
        // Accumulate all the object shadows into the generalized shadow variable.
        gSh = min(gSh, min(min(bl, ln), min(bl2, bl2B))); 
        
        // Dummy break to avoid unrolling... which may or may not work.
        if(rndI>1e6) break;
        
        
    }
    
    // Render the drop shadows.
    col = mix(col, vec3(0), (1. - smoothstep(0., sf, gSh - ew/3.))*.5);

    
    
    // Produce the main scene.
    
    // Render various shapes at each of the triangle cells vertices, edges, etc.
    for(int i = min(0, int(FRAMECOUNT)); i<3; i++){
        
         
        // Unique random edge value.
        float rndI = hash21(id + mix(vID[i], vID[(i + 1)%3], .5));
 
        // Edge tangents.
        vec2 tgnt0 = normalize(v[i] - v[(i + 1)%3]);
        vec2 tgnt1 = normalize(v[(i + 1)%3] - v[(i + 2)%3]);
        vec2 tgnt2 = normalize(v[(i + 2)%3] - v[(i + 0)%3]);
        
        // Normal -- Not used.
        //vec2 n = tgnt0.yx*vec2(1, -1);
        //if(tri<.5){  n *= -1.;  }
         
        
        // Central join line vertices.
        vec2[4] qv = vec2[4](v[i] - tgnt0*ndg, v[(i + 1)%3] + tgnt0*ndg, 
                             v[(i + 1)%3] + tgnt0*ndg + tgnt2*lw, v[i] - tgnt0*ndg + tgnt2*lw);
        
        // Edge join line vertices.
        vec2[4] qvB = vec2[4](v[i] - tgnt0*ndg, v[(i + 1)%3] + tgnt0*ndg, 
                            v[(i + 1)%3] + tgnt0*ndg + tgnt1*lw, v[i] - tgnt0*ndg + tgnt1*lw);
        
        // Central join line.
        ln = sdPoly4(p - tgnt0*.0, qv);
        // Central join line - Side -- Redundant here, but included anyway.
        float lnB = sdPoly4(p - tgnt0*.0, qvB);
       
        
        // Outside joins line - Top.
        float ln2Offs = hexW/3.;
        ln2 = sdPoly4((p) - tgnt2*ln2Offs, qv);///////
        // Outside joins line - Side.
        float ln2B = sdPoly4((p) - tgnt2*ln2Offs, qvB);//////
        
        #ifdef ONE_JOIN
        ln2 = ln2B = 1e5;
        #endif
        
        // Save the old lines -- to deal with overlap on the last vertex. 
        if(i==0) oLn = min(min(ln, lnB), min(ln2, ln2B)); 
      
        
         
        // Main block vertices.
        qv = vec2[4](v[i], v[i] + tgnt2*hexW, 
                     v[i] + tgnt2*hexW + tgnt1*hexW, v[i] + tgnt1*hexW);
        // Main block (The pink one).
        float bl = sdPoly4(p, qv);
        
        
        
        // Choppin bits off the main block to create T-Joins.
        #ifdef TJOINS
        // T - intersection.
        //if((i==1 && tri>.5)|| (i==2 && tri<.5)){// 
            vec2 ap2 = abs(p - v[(i + 0)%3] + tgnt0*.7);
            float blN = max(ap2.y*.8660254 + ap2.x*.5, ap2.x) - hexW*.8660254;
            bl = max(bl, -blN);
        //}
        #endif
        
       
        // Middle block (The green one).
        float blockDist = (1. - hexW/length(v[i] - v[(i + 1)%3]))/2.;// + mBW;
        vec2 eM = mix(v[i], v[(i + 1)%3], blockDist);
        blockDist *= length(v[i] - v[(i + 1)%3]);
        
        float hexW2 = hexW;//*.85;
        // Middle block face vertices.
        qv = vec2[4](v[i], v[i] + tgnt2*hexW2, 
                     v[i] + tgnt2*hexW2 + tgnt1*hexW2, v[i] + tgnt1*hexW2);
        // Middle block face distance.
        float bl2B = sdPoly4(p + tgnt0*(blockDist + mBW), qv);
        
        
        float bl2BLine = sdPoly4(p + tgnt0*(blockDist + mBW) + tgnt0*ndg, qv);
        bl2B = max(bl2B, -max(min(ln, lnB), bl2BLine));
        float bl2BLine2 = sdPoly4(p + tgnt0*(blockDist + mBW) + tgnt0*ndg - tgnt2*ln2Offs, qv);
        bl2B = max(bl2B, -max(min(ln2, ln2B), bl2BLine2));
        
     
        
        
        // Middle block side vertices.
        qv = vec2[4](eM - tgnt0*mBW, eM + tgnt0*mBW, 
                     eM + tgnt0*mBW + tgnt2*hexW2, eM - tgnt0*mBW + tgnt2*hexW2);
        // Middle block side distance.
        float bl2 = sdPoly4(p, qv);
        
        // If we're not displaying the middle green block.
        #ifndef MID_BLOCK
        bl2 = bl2B = 1e5;
        #endif
        
        /*
        qvB = vec2[4](eM - tgnt1*mBW, eM + tgnt1*mBW, 
                     eM + tgnt1*mBW + tgnt0*hexW, eM - tgnt1*mBW + tgnt0*hexW);
        */  
        
        
        // Overlap occurs when producing the final vertex objects, so you need
        // to perform CSG with the previous objects... Sometimes, it's easier
        // just doing things with conventional 3D. :)
        if(i == 2) {
            
            // Next vertex block
            vec2 ap = abs(p - v[(i + 1)%3]);
            float blN = max(ap.y*.8660254 + ap.x*.5, ap.x) - hexW*.8660254;
            
            // Combine it with previous line elements.
            blN = min(blN, oLn);

            // Take the above away from the final blocks and lines so that
            // the overlap doesn't appear... Yeah, I find 2D layering 
            // confusing too. :D
            ln = max(ln, -blN);
            ln2 = max(ln2, -blN);
            lnB = max(lnB, -blN);
            ln2B = max(ln2B, -blN);
            bl2 = max(bl2, -blN);
            bl2B = max(bl2B, -blN);
 
        }  
      
        // Don't show random connections -- The open space kind of looks 
        // more interesting.
        if(rndI<.5){
            
            ln = 1e5;
            ln2 = 1e5;
            lnB = 1e5;
            ln2B = 1e5;
            bl2 = 1e5;
            bl2B = 1e5;
            if(i==0) oLn = 1e5;
        }
        
        
        // Line colors.
        vec3 lnCol = vCol[i]; // Top face line.
        vec3 lnColB = vCol[(i + 1)%3]; // Side face line.
        vec3 lnCol2B = vCol[(i + 1)%3]; // Side face line.
        
        // Block colors.
        vec3 blCol = vCol[(i + 2)%3];
        vec3 blColB = vCol[i];
        
        
        lnCol *= vec3(.7, .85, 1)*.75;
        blCol = min(blCol*vec3(1, .2, .4)*2.8, 1.);
        blColB = min(blColB*vec3(1, .2, .4)*2.8, 1.);
        
        lnColB *= vec3(.7, .85, 1)*.7;
        lnCol2B *= vec3(.7, .85, 1)*.7;
    
        
       
        // Main block render for this pass.
        //float sh = max(.75 - bl/.35, 0.);
        col = mix(col, vec3(0), 1. - smoothstep(0., sf, bl - ew/3.));
        col = mix(col, blCol*vPat[(i + 2)%3], 1. - smoothstep(0., sf, bl + ew*2./3.));    


        // Side joiner line.
        col = mix(col, vec3(0), (1. - smoothstep(0., sf*8., ln2 - ew/3.))*.35);
        col = mix(col, vec3(0), 1. - smoothstep(0., sf, ln2 - ew/3.));
        col = mix(col, lnCol*vPat[i], 1. - smoothstep(0., sf, ln2 + ew*2./3.));
        col = mix(col, vec3(0), 1. - smoothstep(0., sf, ln2B - ew/3.));
        col = mix(col, lnCol2B*vPat[(i + 1)%3], 1. - smoothstep(0., sf, ln2B + ew*2./3.));

        // Central joiner line.
        col = mix(col, vec3(0), (1. - smoothstep(0., sf*8., ln - ew/3.))*.35);
        col = mix(col, vec3(0), 1. - smoothstep(0., sf, ln - ew/3.));
        col = mix(col, lnCol, 1. - smoothstep(0., sf, ln + ew*2./3.));
        col = mix(col, vec3(0), 1. - smoothstep(0., sf, lnB - ew/3.));
        col = mix(col, lnColB*vPat[(i + 1)%3], 1. - smoothstep(0., sf, lnB + ew*2./3.));

        // Middle block.
        //sh = max(.75 - bl2B/.35, 0.); 
        col = mix(col, vec3(0), 1. - smoothstep(0., sf, bl2B - ew/3.));
        col = mix(col, blCol.yxz*vPat[(i + 2)%3], 1. - smoothstep(0., sf, bl2B + ew*2./3.)); 
        //sh = max(.75 - bl2/.35, 0.);
        col = mix(col, vec3(0), 1. - smoothstep(0., sf, bl2 - ew/3.));
        col = mix(col, blColB.yxz*vPat[i], 1. - smoothstep(0., sf, bl2 + ew*2./3.)); 

        // Fake AO on the green block.
        col = mix(col, vec3(0), (1. - smoothstep(0., sf*5., max(ln, bl2B + .03) - ew/3.))*.25);
        col = mix(col, vec3(0), (1. - smoothstep(0., sf*5., max(ln2, bl2B + .03) - ew/3.))*.25);  
        
        // Dummy break to avoid unrolling... which may or may not work. If you think appeasing
        // compilers instead of making code faster is stupid, you'll get no argument from me. :)
        if(rndI>1e6) break;

        
    } 
    
      
    // Subtle pencil overlay... It's cheap and definitely not production worthy,
    // but it works well enough for the purpose of the example. The idea is based
    // off of one of Flockaroo's examples.
    vec2 q = oP*2.;
    vec3 colP = pencil(col, q*RENDERSIZE.y/450.);
    #if 0
    // Just the pencil sketch. The last factor ranges from zero to one and 
    // determines the sketchiness of the rendering... Pun intended. :D
    col = mix(dot(col, vec3(.299, .587, .114))*vec3(1), colP, .7);
    #else
    col = mix(col, 1. - exp(-(col*2.)*(colP + .25)), .9); 
    #endif
    //col = mix(col, colP, .5);
    //col = mix(min(col, colP), max(col, colP), .5); 
    
    
    #ifdef SHOW_CELLS
    // Triangle cell borders.
    float gLnBrd = 1e5;     
    for(int i = 0; i<3; i++){
        gLnBrd = min(gLnBrd, lBox(p, v[i], v[(i + 1)%3], .005));
    }
    col = mix(col, vec3(0), (1. - smoothstep(0., sf, gLnBrd - .0175)));
    col = mix(col, vec3(1, .9, .5), (1. - smoothstep(0., sf, gLnBrd)));
    #endif
    

    // Rough gamma correction and output.
    fragColor = vec4(sqrt(max(col, 0.)), 1);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}