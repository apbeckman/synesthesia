//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick);


			//******** Common Code Begins ********

// The cubemap texture resultion.
#define cubemapRes vec2(1024)

// If you use all four channels of one 1024 by 1024 cube face, that would be
// 4096000 storage slots (1024*1024*4), which just so happens be 160 cubed.
// In other words, you can store the isosurface values of a 160 voxel per side
// cube into one cube face of the cubemap.
//
// The voxel cube dimensions -- That's the one you'd change, but I don't really
// see the point, since setting it to the maximum resolution makes the most
// sense. For demonstrative purposes, dropping it to say, vec3(80), will show
// how a decrease in resolution will affect things. Increasing it to above the
// allowable resolution (for one cube face) to say, vec3(200), will display the
// wrapping issues.
//
// On a side note, I'm going to put up an example later that uses four of the
// cubemap faces, which should boost the resolution to 256... and hopefully,
// not add too much to the complexity, and consequent lag that would follow.
const vec3 dimsVox = vec3(160, 160, 160);
const vec3 scale = vec3(4, 1, 1);
const vec3 dims = dimsVox/scale;



// Reading into one of the cube faces, according to the face ID. To save on cycles,
// I'd hardcode the face you're after into all but the least costly of situations.
// This particular function is used just once for an update in the "CubeA" tab.
//
// The four cube sides - Left, back, right, front.
// NEGATIVE_X, POSITIVE_Z, POSITIVE_X, NEGATIVE_Z
// vec3(-.5, uv.yx), vec3(uv, .5), vec3(.5, uv.y, -uv.x), vec3(-uv.x, uv.y, -.5).
//
// Bottom and top.
// NEGATIVE_Y, POSITIVE_Y
// vec3(uv.x, -.5, uv.y), vec3(uv.x, .5, -uv.y).
vec4 tx(samplerCube tx, vec2 p, int id){

    //vec4 rTx;

    vec2 uv = fract(p) - .5;
    // It's important to snap to the pixel centers. The people complaining about
    // seam line problems are probably not doing this.
    //p = (floor(p*cubemapRes) + .5)/cubemapRes;

    vec3[6] fcP = vec3[6](vec3(-.5, uv.yx), vec3(.5, uv.y, -uv.x), vec3(uv.x, -.5, uv.y),
                          vec3(uv.x, .5, -uv.y), vec3(-uv.x, uv.y, -.5), vec3(uv, .5));


    return texture(tx, fcP[id]);
}


vec4 texMapCh(samplerCube tx, vec3 p){

    p *= dims;
    int ch = (int(p.x*4.)&3);
    p = mod(floor(p), dims);
    float offset = dot(p, vec3(1, dims.x, dims.x*dims.y));
    vec2 uv = mod(floor(offset/vec2(1, cubemapRes.x)), cubemapRes);
    // It's important to snap to the pixel centers. The people complaining about
    // seam line problems are probably not doing this.
    uv = fract((uv + .5)/cubemapRes) - .5;
    return vec4(1)*texture(tx, vec3(-.5, uv.yx))[ch];

}

// Used in conjunction with the function below. When doing things eight times over, any
// saving is important. If I could trim this down more, I would, but there's wrapping
// and pixel snapping to consider. Having said that, I might take another look at it,
// at some stage.
vec4 txChSm(samplerCube tx, in vec3 p, in int ch){

    p = mod(floor(p), dims);
    //vec2 uv = mod(floor(dot(p, vec3(1, dims.x, dims.x*dims.y))/vec2(1, cubemapRes.x)), cubemapRes);
    // I think the fract call below already wraps things, so no "mod" call needed.
    vec2 uv = floor(dot(p, vec3(1, dims.x, dims.x*dims.y))/vec2(1, cubemapRes.x));
    // It's important to snap to the pixel centers. The people complaining about
    // seam line problems are probably... definitely not doing this. :)
    uv = fract((uv + .5)/cubemapRes) - .5;
    return vec4(1)*texture(tx, vec3(-.5, uv.yx))[ch];

}

// Smooth texture interpolation that access individual channels: You really need this -- I
// wish you didn't, but you do. I wrote it a while ago, and I'm pretty confident that it works.
// The smoothing factor isn't helpful at all, which surprises me -- I'm guessing it molds things
// to the shape of a cube. Anyway, it's written in the same way that you'd write any cubic
// interpolation: 8 corners, then a linear interpolation using the corners as boundaries.
//
// It's possible to use more sophisticated techniques to achieve better smoothing, but as you
// could imagine, they require more samples, and are more expensive, so you'd have to think about
// it before heading in that direction -- Perhaps for texturing and bump mapping.
vec4 texMapSmoothCh(samplerCube tx, vec3 p){

    // Voxel corner helper vector.
	const vec3 e = vec3(0, 1, 1./4.);

    // Technically, this will center things, but it's relative, and not necessary here.
    //p -= .5/dimsVox.x;

    p *= dimsVox;
    vec3 ip = floor(p);
    p -= ip;


    int ch = (int(ip.x)&3), chNxt = ((ch + 1)&3);  //int(mod(ip.x, 4.))
    ip.x /= 4.;

    vec4 c = mix(mix(mix(txChSm(tx, ip + e.xxx, ch), txChSm(tx, ip + e.zxx, chNxt), p.x),
                     mix(txChSm(tx, ip + e.xyx, ch), txChSm(tx, ip + e.zyx, chNxt), p.x), p.y),
                 mix(mix(txChSm(tx, ip + e.xxy, ch), txChSm(tx, ip + e.zxy, chNxt), p.x),
                     mix(txChSm(tx, ip + e.xyy, ch), txChSm(tx, ip + e.zyy, chNxt), p.x), p.y), p.z);


 	/*
    // For fun, I tried a straight up average. It didn't work. :)
    vec4 c = (txChSm(tx, ip + e.xxx*sc, ch) + txChSm(tx, ip + e.yxx*sc, chNxt) +
             txChSm(tx, ip + e.xyx*sc, ch) + txChSm(tx, ip + e.yyx*sc, chNxt) +
             txChSm(tx, ip + e.xxy*sc, ch) + txChSm(tx, ip + e.yxy*sc, chNxt) +
             txChSm(tx, ip + e.xyy*sc, ch) + txChSm(tx, ip + e.yyy*sc, chNxt) + txChSm(tx, ip + e.yyy*.5, ch))/9.;
 	*/

    return c;

}




// If you want things to wrap, you need a wrapping scale. It's not so important
// here, because we're performing a wrapped blur. Wrapping is not much different
// to regular mapping. You just need to put "p = mod(p, gSc)" in the hash function
// for anything that's procedurally generated with random numbers. If you're using
// a repeat texture, then that'll have to wrap too.
float gSc;


// IQ's vec2 to float hash.
float hash21(vec2 p){
    p = mod(p, gSc);
    return fract(sin(dot(p, vec2(27.609, 157.583)))*43758.5453);
}


// Commutative smooth maximum function. Provided by Tomkh, and taken
// from Alex Evans's (aka Statix) talk:
// http://media.lolrus.mediamolecule.com/AlexEvans_SIGGRAPH-2015.pdf
// Credited to Dave Smith @media molecule.
float smax(float a, float b, float k){

   float f = max(0., 1. - abs(b - a)/k);
   return max(a, b) + k*.25*f*f;
}


// Commutative smooth minimum function. Provided by Tomkh, and taken
// from Alex Evans's (aka Statix) talk:
// http://media.lolrus.mediamolecule.com/AlexEvans_SIGGRAPH-2015.pdf
// Credited to Dave Smith @media molecule.
float smin(float a, float b, float k){

   float f = max(0., 1. - abs(b - a)/k);
   return min(a, b) - k*.25*f*f;
}

/*
// IQ's exponential-based smooth maximum function. Unlike the polynomial-based
// smooth maximum, this one is associative and commutative.
float smaxExp(float a, float b, float k){

    float res = exp(k*a) + exp(k*b);
    return log(res)/k;
}
*/

// IQ's exponential-based smooth minimum function. Unlike the polynomial-based
// smooth minimum, this one is associative and commutative.
float sminExp(float a, float b, float k){

    float res = exp(-k*a) + exp(-k*b);
    return -log(res)/k;
}


// With the spare cycles, I thought I'd splash out and use Dave's more reliable hash function. :)
//
// Dave's hash function. More reliable with large values, but will still eventually break down.
//
// Hash without Sine.
// Creative Commons Attribution-ShareAlike 4.0 International Public License.
// Created by David Hoskins.
// vec3 to vec3.
vec3 hash33G(vec3 p){


    p = mod(p, gSc);
	p = fract(p * vec3(.10313, .10307, .09731));
    p += dot(p, p.yxz + 19.1937);
    p = fract((p.xxy + p.yxx)*p.zyx)*2. - 1.;
    return p;

    /*
    // Note the "mod" call. Slower, but ensures accuracy with large time values.
    mat2  m = rot2(mod(TIME, 6.2831853));
	p.xy = m * p.xy;//rotate gradient vector
    p.yz = m * p.yz;//rotate gradient vector
    //p.zx = m * p.zx;//rotate gradient vector
	return p;
    */

}

// Cheap vec3 to vec3 hash. I wrote this one. It's much faster than others, but I don't trust
// it over large values.
vec3 hash33(vec3 p){


    p = mod(p, gSc);
    //float n = sin(dot(p, vec3(7, 157, 113)));
    //p = fract(vec3(2097152, 262144, 32768)*n)*2. - 1.;

    //mat2  m = rot2(TIME);//in general use 3d rotation
	//p.xy = m * p.xy;//rotate gradient vector
    ////p.yz = m * p.yz;//rotate gradient vector
    ////p.zx = m * p.zx;//rotate gradient vector
	//return p;

    float n = sin(dot(p, vec3(57, 113, 27)));
    return fract(vec3(2097152, 262144, 32768)*n)*2. - 1.;


    //float n = sin(dot(p, vec3(7, 157, 113)));
    //p = fract(vec3(2097152, 262144, 32768)*n);
    //return sin(p*6.2831853 + TIME)*.5;
}

// This is a variation on a regular 2-pass Voronoi traversal that produces a Voronoi
// pattern based on the interior cell point to the nearest cell edge (as opposed
// to the nearest offset point). It's a slight reworking of Tomkh's example, which
// in turn, is based on IQ's original example. The links are below:
//
// On a side note, I have no idea whether a faster solution is possible, but when I
// have time, I'm going to try to find one anyway.
//
// Voronoi distances - iq
// https://www.shadertoy.com/view/ldl3W8
//
// Here's IQ's well written article that describes the process in more detail.
// http://www.iquilezles.org/www/articles/voronoilines/voronoilines.htm
//
// Faster Voronoi Edge Distance - tomkh
// https://www.shadertoy.com/view/llG3zy
//
//
vec3 cellID;
int gIFrame;
//
vec3 Voronoi(in vec3 p, in vec3 rd){

    // One of Tomkh's snippets that includes a wrap to deal with
    // larger numbers, which is pretty cool.


    vec3 n = floor(p);
    p -= n + .5;


    // Storage for all sixteen hash values. The same set of hash values are
    // reused in the second pass, and since they're reasonably expensive to
    // calculate, I figured I'd save them from resuse. However, I could be
    // violating some kind of GPU architecture rule, so I might be making
    // things worse... If anyone knows for sure, feel free to let me know.
    //
    // I've been informed that saving to an array of vectors is worse.
    //vec2 svO[3];

    // Individual Voronoi cell ID. Used for coloring, materials, etc.
    cellID = vec3(0); // Redundant initialization, but I've done it anyway.

    // As IQ has commented, this is a regular Voronoi pass, so it should be
    // pretty self explanatory.
    //
    // First pass: Regular Voronoi.
	vec3 mo, o;

    // Minimum distance, "smooth" distance to the nearest cell edge, regular
    // distance to the nearest cell edge, and a line distance place holder.
    float md = 8., lMd = 8., lMd2 = 8., lnDist, d;

    // Note the ugly "gIFrame" hack. The idea is to force the compiler not
    // to unroll the loops, thus keep the program size down... or something.
    // GPU compilation is not my area... Come to think of it, none of this
    // is my area. :D
    for( int k=min(-2, gIFrame); k<=2; k++ ){
    for( int j=min(-2, gIFrame); j<=2; j++ ){
    for( int i=min(-2, gIFrame); i<=2; i++ ){

        o = vec3(i, j, k);
        o += hash33(n + o) - p;
        // Saving the hash values for reuse in the next pass. I don't know for sure,
        // but I've been informed that it's faster to recalculate the had values in
        // the following pass.
        //svO[j*3 + i] = o;

        // Regular squared cell point to nearest node point.
        d = dot(o, o);

        if( d<md ){

            md = d;  // Update the minimum distance.
            // Keep note of the position of the nearest cell point - with respect
            // to "p," of course. It will be used in the second pass.
            mo = o;
            cellID = vec3(i, j, k) + n; // Record the cell ID also.
        }

    }
    }
    }

    // Second pass: Distance to closest border edge. The closest edge will be one of the edges of
    // the cell containing the closest cell point, so you need to check all surrounding edges of
    // that cell, hence the second pass... It'd be nice if there were a faster way.
    for( int k=min(-3, gIFrame); k<=3; k++ ){
    for( int j=min(-3, gIFrame); j<=3; j++ ){
    for( int i=min(-3, gIFrame); i<=3; i++ ){

        // I've been informed that it's faster to recalculate the hash values, rather than
        // access an array of saved values.
        o = vec3(i, j, k);
        o += hash33(n + o) - p;
        // I went through the trouble to save all sixteen expensive hash values in the first
        // pass in the hope that it'd speed thing up, but due to the evolving nature of
        // modern architecture that likes everything to be declared locally, I might be making
        // things worse. Who knows? I miss the times when lookup tables were a good thing. :)
        //
        //o = svO[j*3 + i];

        // Skip the same cell... I found that out the hard way. :D
        if( dot(o - mo, o - mo)>.00001 ){

            // This tiny line is the crux of the whole example, believe it or not. Basically, it's
            // a bit of simple trigonometry to determine the distance from the cell point to the
            // cell border line. See IQ's article for a visual representation.
            lnDist = dot(0.5*(o + mo), normalize(o - mo));

            // Abje's addition. Border distance using a smooth minimum. Insightful, and simple.
            //
            // On a side note, IQ reminded me that the order in which the polynomial-based smooth
            // minimum is applied effects the result. However, the exponentional-based smooth
            // minimum is associative and commutative, so is more correct. In this particular case,
            // the effects appear to be negligible, so I'm sticking with the cheaper polynomial-based
            // smooth minimum, but it's something you should keep in mind. By the way, feel free to
            // uncomment the exponential one and try it out to see if you notice a difference.
            //
            // // Polynomial-based smooth minimum.
            //lMd = smin(lMd, lnDist, lnDist*.75); //lnDist*.75
            //
            // Exponential-based smooth minimum. By the way, this is here to provide a visual reference
            // only, and is definitely not the most efficient way to apply it. To see the minor
            // adjustments necessary, refer to Tomkh's example here: Rounded Voronoi Edges Analysis -
            // https://www.shadertoy.com/view/MdSfzD
            lMd = sminExp(lMd, lnDist, 16.);

            // Minimum regular straight-edged border distance. If you only used this distance,
            // the web lattice would have sharp edges.
            lMd2 = min(lMd2, lnDist);
        }

    }
    }
    }

    // Return the smoothed and unsmoothed distance. I think they need capping at zero... but
    // I'm not positive.
    return max(vec3(lMd, lMd2, md), 0.);
}


/*

	Polar Grid Motion
	-----------------

	Utilizing a repeat polar grid to plot objects moving along a transcendental rose
    curve path.

	I wanted to post something nice and simple that the average GPU could handle. This
    was inspired in part by Vovosunt's "Dots and Spirals," but was based on Fabrice's
    original "rosace 3" example. They say small things amuse small minds, and at the
    time, that particular example kept me busy for ages. :)

	The particle movement is pretty standard. Each object moves in a circular path with a
    varying radius and Z position based on transcendental functions according to angular
    position. For anyone not quite familiar with it, I've provided some links below.

	The only mildly tricky bit was the rendering process. That involved keeping track of
	potentially overlapping objects, then sorting prior to rendering. There are three
	polar partitioned grids overlapping one another, but arranged radially and depthwise
	to look like they form a continuous closed curve. Potential overlapping objects
	consist of three overlapping polar cells (one for each loop) and each of their adjacent
    polar neighbors (the two on either side). That's nine altogether, which is not very
    taxing on the average GPU. Therefore, I wouldn't expect any frame rate problems.

	I've tried to keep things relatively simple, but there's still a bit of esoteric
	dressing up code in there. For anyone who'd like to make something similar, you'd be
    better off starting from scratch, then using this, Vovosunt's or Fabrices example as
    a guide.

	Anyway, I'd call this a pseudo 3D example. At some stage, I'll try to come up with a
	raymarched variation. By the way, there's a "TRANSPARENT" and "RANDOM_VARIATION"
	define below, for anyone bored enough to try them. :)


    Examples:

	// Really nice.
	Dots and Spirals - Vovosunt
	https://www.shadertoy.com/view/MltyzN

    // The original that I based this off of.
    rosace 3 ( 215 chars) - FabriceNeyret2
	https://www.shadertoy.com/view/ls3XWM

    // A much simpler demonstration of the key motion concept.
	Linear motion - ABizard
	https://www.shadertoy.com/view/4lKyzd

	Links:

    Rose (mathematics)
    https://en.wikipedia.org/wiki/Rose_(mathematics)

	ROSE
    https://www.mathcurve.com/courbes2d.gb/rosace/rosace.shtml

*/



// A custom transparent overlay effect. The thing I like about pixel shaders is that
// Photoshop layer effects are almost rudimentary.
// #define TRANSPARENT


// Making use of the individual particle ID to produce a random variation.
//#define RANDOM_VARIATION



// Cheap and nasty 2D smooth noise function with inbuilt hash function -- based on IQ's
// original. Very trimmed down. In fact, I probably went a little overboard. I think it
// might also degrade with large time values, but that's not an issue here.
float n2D(vec2 p) {

	vec2 i = floor(p); p -= i; p *= p*(3. - p*2.);

	return dot(mat2(fract(sin(vec4(0, 1, 113, 114) + dot(i, vec2(1, 113)))*43758.5453))*
                vec2(1. - p.y, p.y), vec2(1. - p.x, p.x) );

}

// FBM -- 4 accumulated noise layers of modulated amplitudes and frequencies.
float fbm(vec2 p){ return n2D(p)*.533 + n2D(p*2.)*.267 + n2D(p*4.)*.133 + n2D(p*8.)*.067; }


vec3 getLight(vec3 p, vec3 n, vec3 lp){

    vec3 ld = lp - p;
    float lDist = length(ld);
    ld /= lDist;
    float diff = max(dot(ld, n), 0.);
    float atten = 1.5/(1. + lDist*lDist);
    vec3 light = vec3(1)*(diff + .5)*atten;

    return light;
}


vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 U = _xy;


	vec2 R = RENDERSIZE.xy;

    // Setting a minimum resolution, since fullscreen looks too bloated... Of course,
    // that would ruin a mobile phone fullscreen settings... Too many systems, so it's
    // impossible to win without seriously ruining your code. I miss the days when we were
    // all on 17 inch screens... but not the grainy PPI. Definitely don't miss that. :)
    float yRes = min(R.y, 800. + 400*sin(morph_time));

    // Screen coordinates. This started with a discussion on Fabrice's rosace example
    // (See the link above). Hence, the confusing minimal variable names. :)
    U = (2.*U - R)/yRes;


    //U *= 1. + dot(U, U)*.025; // Makeshift fisheye, if that's your thing.


    // Three lines of 21, so 63 objects in all. There's nothing special about 21. It's
    // just the number I settled on. Higher numbers work, but object size needs to be
    // reduced, since overlap becomes a problem.
    float num = 21. + int(num_circles);


    // The scene light. Placed just above the scene.
    vec3 lp = vec3(0, 0, 1); // Moving light: vec3(.3*cos(t), .2*sin(t), 1).
    // The object normals. Trivial, in this case, since all are facing the same way.
    vec3 n = vec3(0, 0, -1);






    // SCENE BACKGROUND.
    //
    // Just some noise, lines, and square geometry. Hopefully, self explanatory.

    // Initialize the background do a brownish gradient.
    //vec3 bg = vec3(.5, .45, .4);
    vec3 bg = vec3(0); // nope, we want it black so we can automask it

    //#ifdef TRANSPARENT
    //vec3 bgTransparent = bg / 2.;
    //bg = mix(bg, bgTransparent, transparency);
    //#endif

    // Apply some light to the background.
    vec3 light = getLight( vec3(U, 2.6), n,lp);
    bg *= light;

    // Apply some subtle marbly noise.
    float ns = fbm(U*5. + 17.3);
    ns = mix(ns, sin(ns*32. - cos(ns*34.)), .125);
    bg *= max(ns, 0.)*.4 + .8;

    // Apply a grainy randomized diagonal pattern. It's subtle, but I prefer it.
    // Without it, the background seems a little too clean.
    float pat = clamp(sin((U.x - U.y)*min(R.y, 800.)/1.5)*1. + .5, 0., 1.);
    float rnd = fract(sin(dot(U, vec2(11.27, 112.43)))*43758.5453);
    if(rnd>.5) pat *= .6;
    else pat *= 1.4;
    bg *= pat*.3 + .75;



    // Initiate the scene color to the background.
    vec3 col = bg;

    // Render some border objects to frame things a little bit.
    //
    // Border sights: The background corners looked a little empty, so I threw
    // these in to balance things out... Not sure if it worked, but it's done now. :)
    vec2 b = vec2(RENDERSIZE.x/RENDERSIZE.y, 1) - .1;
    vec2 q = abs(U);
    float bord = max(q.x - b.x, q.y - b.y);
    bord = max(bord, -(bord + .11));
    bord = max(bord, -min(q.x - b.x + .22, q.y - b.y + .22));
    //bord = max(bord, -(bord + .02));


    // Render the border sight... edge things, or whatever they are.
    float falloff = 1./min(R.y, 800.);
    col = mix(col, vec3(0), (1. - smoothstep(0., falloff*12., bord ))*.35);
    col = mix(col, vec3(0), (1. - smoothstep(0., falloff, bord))*.7);
    col = mix(col, bg*2.2, (1. - smoothstep(0., falloff, bord + .01)));
    col = mix(col, vec3(0), (1. - smoothstep(0., falloff, bord + .035)));

    col = mix(col, bg*1.2, (1. - smoothstep(0., falloff, bord + .044)));


    // OBJECT MOVEMENT.
    //
    // Moving the cell objects around the rose curve path whilst storing their
    // positions for rendering.

    // Nine storage vectors: XYZ hold the 3D positions of the center object and each
    // of its 8 polar neighbors. The W position holds the object ID.
    //
    // When taking polar coordinates into account, there are a potential 9 objects that
    // could possibly overlap. Therefore, Z distances on all 9 need to be sorted to get
    // the rendering order correct. If you spaced out the objects more, then you could
    // probably get away with 3.
    vec4 c[9];

    vec4 p; // Utility storage for the object position and ID.

    // A bit of global rotation -- Rotating the collection of objects as a whole.
    float t = fly_time/4.;
    // Alternatively, you could take "t" out of the loop below, and globally rotate
    // "U" itself.
    //float t = TIME/4., cs = cos(t), sn = sin(t);
    //U *= mat2(cs, sn, -sn, cs);

    // Storing "atan" to save a couple of extra calls. Not overly necessary, but my
    // oldschool brain still thinks of it as an expensive function. :)
    float a0 = atan(U.y, U.x) - t;

    // Due to overlap, polar neighbors need to be considered. If the cell objects are
    // smaller, then it's not a problem, but I wanted a bit of object density.
    //
    for(int i=0; i<3; i++){ // Adjacent polar angle cells.

        float a = a0; // Polar angle.

        for(int j=0; j<3; j++) {  // Three intertwining overlapping revolutions.

            // Current cell angle.
            //
            // Note the "i - 1" figure. That's because we're considering overlapping
            // cells to the left and right of "ia," and not two to the right... It took
            // me a while to realize that oversight. :)
            //
            // Cell index.
        	float ia = mod(floor(a*num/6.283) + float(i - 1), num*3.);
            // Set the object ID to the cell index.
            p.w = ia;
            // Covert cell index to a polar angle.
            ia = (ia/num + .5/num)*6.283;

            // Move X and Y along a rose curve path... or a rosace path, as Fabrice calls
            // it, which I'll assume is the French rosette discription. Without going into
            // detail, it's a circular path with a varying sinusoidal radius that gives it
            // that interesting overlapping look.
            //
            // By the way, figures of 2./3., 4./3, 7./3. will also produce patterns. For
            // other patterns, an adjustment of the "j" loop here and below and the "c"
            // array size would be necessary... I'm sure you'll figure it out. :)
            float off = ia*5./3. + fly_time;
        	// The X and Y positions. Basically, a circle of varying radius.
            p.xy = (.55 - .25*sin(off))*vec2(cos(ia + t), sin(ia + t));

            // By varying the Z component with a complementing offset, the objects move
            // along an interwoven closed path. Obviously, if you set Z to a constant, all
            // objects would be coplanar and things wouldn't work.
        	p.z = 2. + cos(off - 3.14159/5.)*.35;
            //p.z = (2. + cos(a*5./3. + TIME - 3.14159/5.)*.35);


            // Store the current cell postion and ID for usage in the rendering loop.
            c[i*3 + j] = p;

            // Increase the polar angle. Altogether, there'll be three whole revolutions.
            a += 6.283;


    	}

    }

    // OVERLAPPING NEIGHBORING OBJECT SORTING.

    // Super lazy distance ordering: Since there are only 45 ((9 + 1)*9/2) iterations --
    // or something along those lines -- performing a quick swap, the GPU shouldn't really
    // notice. Also, I've heard that for small datasets, keeping the algorithm simple
    // (branchless, etc) is more important than iteration count, but I don't know for sure.
    //
    // By the way, I think there's a quick vector swap somewhere, so I should probably
    // track that down.
    //
    // On a side note, you could probably get away with a Z-buffer test and do away with the
    // ordering, but it might disturb the smooth rendering.
    //
	for(int i=0; i<9; i++){
        for(int j=i + 1; j<9; j++){
            // Branchless swaps are possible (see below), but I'm a little paranoid
            // regarding precision issues, so I'm sticking to what I know. :)
            if(c[i].z<c[j].z){
                vec4 temp = c[i]; c[i] = c[j]; c[j] = temp;

                // Branchless swap: It works fine on my machine, but I can't
                // guarantee it'll work in all situations, or that it's faster.
                // If someone knows one way or the other, I'd love to know.
                //c[i] = c[i] + c[j]; c[j] = c[i] - c[j]; c[i] = c[i] - c[j];
            }
        }
	}

    // OBJECT RENDERING.

    for(int i=0; i<3; i++){ // Adjacent polar angle cells.

        for (int j=0; j<3; j++) {  // Three intertwining overlapping revolutions.

            // Obtaining the position and ID for the current cell.
            p = c[i*3 + j];

            // Using the Z coordi
            float sz = .175/p.z;
            float d = length(U - p.xy) - sz;

            // Main object color
            //#ifdef RANDOM_VARIATION
                // Object ID based random value.
            	float rnd = fract(sin(p.w + 37.)*43758.5453);
                // Annulus. Equivalent to: max(d, -(d + sz*.75)).
            	if(rnd>.5) d = (abs(d + sz*.375) - sz*.375);

            	// Random color variations.
            	//#ifdef TRANSPARENT
              vec3 pColTransparent = mod(p.w, 3.)==0.?  vec3(1.5, .9, .3) : vec3(1.5, .24, .52);
                //#else
            	//vec3 pCol = rnd>.25? vec3(1, .75, .25) : vec3(.6, .9, .25);
            	vec3 pColRand = mix(mod(p.w, 3.)==0.? vec3(.35) : vec3(1, .22, .45), pColTransparent, transparency);

            	vec3 pCol = mix(vec3(1, .75, .25), pColRand, random_colors);


            // Lighting the object. Very simple.
            light = getLight(p.xyz, n, lp);

            // Circular smoothstep falloff, based in the radial inverse.
            falloff = .0005/sz;

            // Rendering the simple pattern on the discs. By the way, you could use some
            // repeat trickery, and cut these steps down, but this isn't a GPU intensive
            // example, and I wanted to try different variations, and so forth. Also, for
            // readability, I wanted "col" written on the left, so these could be trimmed
            // down further.
            //
            // Shadow, edges, color, etc.
            //

            	pCol = mix(pCol * light,(col + .1)*pCol*light*3., transparency);
                // Alternate: Fade between transparent and opaque.
                //pCol = mix(pCol*light, (col + .1)*pCol*light*3., smoothstep(-.1, .1, sin(TIME/4.)));
            //	pCol *= light;

            col = mix(col, vec3(0), (1. - smoothstep(0., falloff*10., d - .0035))*.5);
            col = mix(col, vec3(0), 1. - smoothstep(0., falloff, d));
            col = mix(col, pCol, 1. - smoothstep(0., falloff, d + .01));
            col = mix(col, vec3(0), 1. - smoothstep(0., falloff, d + sz - sz*.4));
            col = mix(col, vec3(light), 1. - smoothstep(0., falloff, d + sz - sz*.4 + .01));
            col = mix(col, vec3(0), 1. - smoothstep(0., falloff, d + sz - sz*.2 + .01));


        }

    }


	// POSTPROCESSING.

	// A bit of color mixing, based on the canvas Y coordinate.
    //col = mix(col.xzy, col, .75);
    col = mix(col.xzy, col, U.y*.3 + .65);

    #ifndef RANDOM_VARIATION

    col = mix(col, col.zyx, transparency); // transparency check
    //col = mix(); // random_colors check
    #endif

    // Rough gamma correction.
    O = vec4(sqrt(max(col, 0.)), 1);
	return O;
 }




vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}
