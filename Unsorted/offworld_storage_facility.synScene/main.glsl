

			//******** BuffA Code Begins ********

/*

    Offworld Storage Facility
    -------------------------
    
    Using a mixture of raymarching and cell by cell traversal techniques 
    to raymarch a pseudo random box subdivision in order to render a
    basic sci-fi looking scene in realtime. The methods used here aren't 
    what I'd call common, but most are not new either.
    
    I love perusing static geometry-based sci-fi scenes rendered with 
    applications like Blender, Cinema 4D, etc. Of course, this scene
    doesn't compare to the beauty and sheer complexity of some of the
    imagery out there, but hopefully it conveys that it's possible to do 
    more in a pixelshader than rendering a bunch of perfectly aligned 
    flat cubes.
    
    The methods used have been described before, but for those not
    familiar, this is a rendering of a basic 3D subdivided grid, with some
    extra XYZ axes shuffling to give the appearance of haphazardness.
    This is definitely not what I'd call a decent random 3D packing, but
    it was simple to make and efficient enough for the purposes of this
    realtime demonstration. By the way, with some basic tweaks, this 
    particular packing can look quite random in appearance.
    
    A standard raymarching algorithm has been used with the addition of 
    ray-to-cell wall boundary collisions performed in order to advance 
    the ray from cell to cell. This way, you can render artifact free 
    scenes with most of the benefits or raymarching. More importantly, 
    only one cell per pass need be rendered, which means it's much faster. 
    It also means that you can render in more detail without slowing down 
    your GPU.
    
    Anyway, the purpose of this example was to demonstrate a traversal of
    a quasi random box scene in realtime, so not a lot of effort was put
    into the design -- I'm not sure where I was going with this, but it 
    has subtle Quake 3 overtones. :)
    
    
    
    Related examples:
    
    // This is more of a sliced layer example, but it's still a 3D traversal 
    // of what I'd call an irregular grid and it uses a similar technique.
    asymmetric blocks layers tower - jt
    https://www.shadertoy.com/view/DdKBDh
    
    // IQ was using a similar technique to this before it was cool. :D
    Cubescape  - iq
	https://www.shadertoy.com/view/Msl3Rr 
    

*/

// Maximum scene distance.
#define FAR  15.

// Color scheme - Copper with lights: 0, Titanium with purple: 1
#define SCHEME 0

// Light type - Direct: 0, Point: 1
#define LIGHT_TYPE 1

// Object ID.
int gOID;


// Standard 2D rotation formula.
mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, s, -s, c); }

 
// Tri-Planar blending function: Based on an old Nvidia writeup:
// GPU Gems 3 - Ryan Geiss: http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
vec3 tex3D(sampler2D tex, in vec3 p, in vec3 n){    
    
    // Ryan Geiss effectively multiplies the first line by 7. It took me a while to realize that 
    // it's largely redundant, due to the division process that follows. I'd never noticed 
    // on account of the fact that I'm not in the habit of questioning stuff written by 
    // Ryan Geiss. :)
    n = max(n*n - .2, .001); // max(abs(n), 0.001), etc.
    n /= dot(n, vec3(1)); 
    //n /= length(n); 
    
    // Texure samples. One for each plane.
    vec3 tx = texture(tex, p.yz).xyz;
    vec3 ty = texture(tex, p.zx).xyz;
    vec3 tz = texture(tex, p.xy).xyz;
    
    // Multiply each texture plane by its normal dominance factor.... or however you wish
    // to describe it. For instance, if the normal faces up or down, the "ty" texture 
    // sample, represnting the XZ plane, will be used, which makes sense.
    
    // Textures are stored in sRGB (I think), so you have to convert them to linear space 
    // (squaring is a rough approximation) prior to working with them... or something like 
    // that. :) Once the final color value is gamma corrected, you should see correct 
    // looking colors.
    return mat3(tx*tx, ty*ty, tz*tz)*n; 
}


float ubox(vec3 p, vec3 b) {
  p = abs(p) - b;
  return length(max(p, 0.)) + min(max(p.x, max(p.y, p.z)), 0.);
}

 
// IQ's vec2 to float hash.
float hash31(vec3 p){ 
    //return texture(iChannel2, p).x;
    return fract(sin(dot(p, vec3(113.619, 57.583, 27.897)))*43758.5453); 
}


float hashV(vec2 a, float b){

    vec3 p = vec3(a.x, b, a.y);
    return fract(sin(dot(p, vec3(113.619, 57.583, 27.897)))*43758.5453); 

}
/*
// This is IQ's WebGL 2 hash formula: I've modified it slightly to 
// take in the normal decimal floats that we're used to passing. It 
// works here, I think, but I'd consult the experts before using it.
//
// I remember reading through a detailed explanation of the C++ hash 
// we all used to use many years ago (which the following would be
// similar to), but have long since forgotten why it works. By the 
// way Nimitz, and Dave Hoskins have formulae on Shadertoy that's worth
// looking up.
//
// Integer Hash - III - IQ
// https://www.shadertoy.com/view/4tXyWN
float hash21(vec2 p){
    
    uvec2 q = floatBitsToUint(p);
	q = 1103515245U*((q>>1U)^q.yx);
    uint n = 1103515245U*(q.x^(q.y>>3U));
    return float(n)/float(0xffffffffU);
}
*/


 

vec4 uTone(vec4 x){
    return ((x*(x*.6 + .1) + .004)/(x*(x*.6 + 1.)  + .06) - .0667)*1.933526;    
}


/*
// Commutative smooth minimum function. Provided by Tomkh, and taken 
// from Alex Evans's (aka Statix) talk: 
// http://media.lolrus.mediamolecule.com/AlexEvans_SIGGRAPH-2015.pdf
// Credited to Dave Smith @media molecule.
float smin(float a, float b, float k){

   float f = max(0., 1. - abs(b - a)/k);
   return min(a, b) - k*.25*f*f;
}

vec3 smin(vec3 a, vec3 b, float k){

   vec3 f = max(vec3(0), 1. - abs(b - a)/k);
   return min(a, b) - k*.25*f*f;
}

// Commutative smooth maximum function. Provided by Tomkh, and taken 
// from Alex Evans's (aka Statix) talk: 
// http://media.lolrus.mediamolecule.com/AlexEvans_SIGGRAPH-2015.pdf
// Credited to Dave Smith @media molecule.
float smax(float a, float b, float k){
    
   float f = max(0., 1. - abs(b - a)/k);
   return max(a, b) + k*.25*f*f;
}

vec3 smax(vec3 a, vec3 b, float k){
    
   vec3 f = max(vec3(0), 1. - abs(b - a)/k);
   return max(a, b) + k*.25*f*f;
}
*/


// Texture sample.
//
vec3 getTex(sampler2D iCh, vec2 p){
    
    // Strething things out so that the image fills up the window. You don't need to,
    // but this looks better. I think the original video is in the oldschool 4 to 3
    // format, whereas the canvas is along the order of 16 to 9, which we're used to.
    // If using repeat textures, you'd comment the first line out.
    //p *= vec2(RENDERSIZE.y/RENDERSIZE.x, 1);
    vec3 tx = texture(iCh, p).xyz;
    return tx*tx; // Rough sRGB to linear conversion.
}

// Camera path. Arranged to coincide with the frequency of the lattice.
vec3 path(float t){
  
    //return vec3(0, 0, t); // Straight path.
    //return vec3(-sin(t/2.), sin(t/2.)*.5 + 1.57, t); // Windy path.
    
    //float s = sin(t/24.)*cos(t/12.);
    //return vec3(s*12., 1., t);
     
    // Transcendental path moving around the XY plane.
    float a = sin(t*.22);
    float b = cos(t*.28);
    return vec3(a*4. - b*1.5, (b*1.2 + a*1.), t);
     
    
}


/*
// Height map value, which is just the pixel's greyscale value.
float hm(in vec2 p){ 

    //float h = getTex(BuffA, p/32.).x; 
    float h = dot(getTex(BuffA, p/32.), vec3(.299, .587, .114)); 
    float hv = min(h, abs(p.x)/8.);
    return hv*8. + h*1. + .05;//h*.975 + .025;
 }
 */
 
 // The 3D distance field function to test voxels against. This is
 // just some 3D noise with a path guided tunnel running through it.
 float hm(vec3 p3){

    // Path.
    vec2 pth = path(p3.z).xy;
    vec3 p = p3;
    p.xy -= pth.xy;
    
    // Noise with a cylinder (wrapped around a path) taken out of it.
    float tx = texture(iChannel2, p3/48.).x;
    return max((tx - .5), -(length(p.xy/vec2(1, 1.5)) -  .85));

}



// IQ's 3D signed box formula: I tried saving calculations by using the unsigned one, and
// couldn't figure out why the edges and a few other things weren't working. It was because
// functions that rely on signs require signed distance fields... Who would have guessed? :D
float sBoxS(vec3 p, vec3 b, float sf){

  p = abs(p) - b + sf;
  return min(max(p.x, max(p.y, p.z)), 0.) + length(max(p, 0.)) - sf;
}

// IQ's 2D signed box formula with some added rounding.
float sBoxS(vec2 p, vec2 b, float sf){

  p = abs(p) - b + sf;
  return min(max(p.x, p.y), 0.) + length(max(p, 0.)) - sf;
}

// A slight variation on a function from Nimitz's hash collection, here: 
// Quality hashes collection WebGL2 - https://www.shadertoy.com/view/Xt3cDn
vec3 hash23(vec3 f){
    
    uvec3 p = floatBitsToUint(f);
    p = 1103515245U*((p >> 2U)^(p.yzx>>1U)^p.zxy);
    uint h32 = 1103515245U*(((p.x)^(p.y>>3U))^(p.z>>6U));

    uint n = h32^(h32>>16);
    
    // See: http://random.mat.sbg.ac.at/results/karl/server/node4.html
    uvec3 rz = uvec3(n, n*16807U, n*48271U); 
    return vec3((rz >> 1) & uvec3(0x7fffffffU))/float(0x7fffffff);
}

// Subdivided rectangle grid.
vec3 getGrid(vec3 p, inout vec3 sc, inout vec3 id){
    
 
    // Block offsets.
    vec3 ipOffs = vec3(0);
    /*
    // Alternate Y floor bricks.
    if(mod(floor(p.y/sc.y) + floor(p.z/sc.z), 2.)<.5){
        p.x -= sc.x/2.; // Row offset. //*hash21(vec2(0, id.y))
        ipOffs.x += .5; //*hash21(vec2(0, id.y))
    }*/
    
    // X, Y and Z offsets.
    vec3 ii = floor(p/sc);
    //vec3 h3 = vec3(hash21(ii.xy), hash21(ii.yz), hash21(ii.zx))*.5 + .5;//
    vec3 h3 = vec3(2./3.);//
    
    vec3 mp = mod(ii, 2.);//mod(floor(p/sc), 2.);
    if(mp.x<.5 && mp.y>.5){
        p.z -= sc.z/2.*h3.x; // Row offset. //*hash21(vec2(0, id.y))
        ipOffs.z += .5*h3.x; //*hash21(vec2(0, id.y))
    } 
    if(mp.y<.5 && mp.z>.5){
        p.x -= sc.x/2.*h3.y; // Row offset. //*hash21(vec2(0, id.y))
        ipOffs.x += .5*h3.y; //*hash21(vec2(0, id.y))
    }
    if(mp.z<.5 && mp.x>.5){
        p.y -= sc.y/2.*h3.z; // Row offset. //*hash21(vec2(0, id.y))
        ipOffs.y += .5*h3.z; //*hash21(vec2(0, id.y))
    }    

 
    // Original position.
    vec3 oP = p;
    
    
    // Block ID.
    vec3 ip;
    
    //ip = floor(p/sc) + .5;
    
    //#define EQUAL_SIDES
    
    // Subdivide.
    for(int i = 0; i<2; i++){
        
        // Current block ID.
        ip = floor(p/sc) + .5;
        float fi = float(i)*.0617; // Unique loop number.
        #ifdef EQUAL_SIDES        
        // Squares.
        
        // Random split.
        if(hash31(ip*sc + .253 + fi)<.333){
           sc /= 2.;
           p = oP;
           ip = floor(p/sc) + .5; 
        }
        
        #else
        
        // Powers of two rectangles.
        
        //vec3 h23 = hash23(ip*sc + .253 + fi);
        // h42 = texture(iChannel2, ip*sc*113.619 + .253 + fi);
        
        // Random X-split.
        if(hash31(ip*sc + .253 + fi)<.333){//3 && sc.x>1./8.
        //if(h23.x<.333){
           sc.x /= 2.;
           p.x = oP.x;
           ip.x = floor(p.x/sc.x) + .5;
        }
        // Random Y-split.
        if(hash31(ip*sc + .453 + fi)<.333){ // && sc.y>1./8.
        //if(h23.y<.333){
           sc.y /= 2.;
           p.y = oP.y;
           ip.y = floor(p.y/sc.y) + .5;
        }
        // Random Z-split.
        if(hash31(ip*sc + .653 + fi)<.333){ // && sc.y>1./8.
        //if(h23.z<.333){
           sc.z /= 2.;
           p.z = oP.z;
           ip.z = floor(p.z/sc.z) + .5;
        }        
        #endif
         
    }
    
    // Cell ID (id is an "inout" variable).
    id = (ip + ipOffs)*sc;    
     
    
    // Return the local coordinates.
    return p - ip*sc;

}



// Global cell boundary distance variables.
vec3 gDir; // Cell traversing direction.
vec3 gRd; // Ray direction.
float gCD; // Cell boundary distance.
// Box dimension and local XY coordinates.
vec3 gSc; 
vec3 gID;
vec3 gP;


vec3 glow; // Glow.

float map(vec3 q3){
  


    // Floor. You can barely see it, but it's down there.
    float fl = abs(q3.y + .1 + 4.) - .1;
  
     
    vec3 sc = vec3(.5, 1, .5); // Scale.
 
     // Local coordinates and cell ID.
    vec3 id; // The cell ID is an "inout" variable.
    vec3 p3 = getGrid(q3, sc, id); 
    // 3D coordinates.
    vec3 p = p3;
    
   

    // The distance from the current ray position to the cell boundary
    // wall in the direction of the unit direction ray. This is different
    // to the minimum wall distance, so you need to trace out instead
    // of merely doing a box calculation. Anyway, the following are pretty 
    // standard cell by cell traversal calculations. The resultant cell
    // distance, "gCD", is used by the "trace" and "shadow" functions to 
    // restrict the ray from overshooting, which in turn restricts artifacts.
    vec3 rC = (gDir*sc - p)/gRd;
    //vec2 rC = (gDir.xz*sc.xz - p)/gRd.xz; // For 2D, this will work too.
    
    // Minimum of all distances, plus not allowing negative distances, which
    // stops the ray from tracing backwards... I'm not entirely sure it's
    // necessary here, but it stops artifacts from appearing with other 
    // non-rectangular grids.
    gCD = max(min(min(rC.x, rC.y), rC.z), 0.) + .0001;
    //gCD = max(min(rC.x, rC.y), 0.) + .001; // Adding a touch to advance to the next cell.


    // The 3D distance field for that particular cell. 
    // See the "hm" function, above.
    float h = hm(id);
 

    // Change the prism rectangle scale just a touch to create some subtle
    // visual randomness. You could comment this out if you prefer more order.
    vec3 oSc = sc;
    sc -= .03;//.1*(hash31(id + .09)*.8 + .2);
    
    // Saving the global scale.
    gSc = sc;
    
    
    
    // Object distance, inner window object and window light.
    float d = 1e5;
    float innerD = 1e5;
    float winLight = 1e5;
    
    
    // If we're under the 3D distance field zero mark, render an object,
    // and whatever else, in that cell.
    if(h<0.){
   

        // Edge factor.
        float ef = sqrt(min(min(oSc.x, oSc.y), oSc.z));
        
        d = sBoxS(p, vec3(sc.x, sc.y, sc.z)/2., ef*.04);
        float oD = d;
        innerD = d + .08;
        vec3 win;
        win.x = sBoxS(p, vec3(oSc)/vec3(3.5, 3.5, 2), ef*.04);
        win.y = sBoxS(p, vec3(oSc)/vec3(2, 3.5, 3.5), ef*.04);
        win.z = sBoxS(p, vec3(oSc)/vec3(3.5, 2, 3.5), ef*.04);
        //win.z = sBoxS(p.xz, sc.xz/2., 0.);

        

        win = max(win, -innerD);
        ////float tx = texture(iChannel2, q3*4.).x;
        //win += tx*.002;

        // Window frames.
        d = min(d, min(min(win.x, win.y), win.z));
        //d = min(d, win.y);
        //d = min(d, win.z);



        // Empty out the space in the windowed rooms.
        d = max(d, -innerD);
        innerD += .005;

        
        // Adding some very basic striated greeble lines. You could
        // do much more interesting stuff than this.
        vec3 q = p3;
        q = mod(q, 1./16.) - .5/16.;
         
        // X-direction. Left and right box striations.
        if(hash31(id + .11)<.35){
            float lx = sBoxS(p.xy, oSc.xy/vec2(3.5) - ef*.06, 0.);
            d = max(d, -lx);
            winLight = min(winLight, lx);
            /*vec2 q = p.xy;
            d = max(d, -sBoxS(q, sc.xy/8., 0.));        
            q = mod(q- sc.yz/12., sc.xy/3.) - sc.xy/3.;
            d = max(d, -sBoxS(q , sc.xy/8., 0.));*/
        }
        else d = max(d, -max(abs(q.x) - .25/16., -oD));

        // Z-direction box striations.
        if(hash31(id + .21)<.35){
           float lx = sBoxS(p.yz, oSc.yz/vec2(3.5) - ef*.06, 0.);
            d = max(d, -lx);
            winLight = min(winLight, lx);
            /*vec2 q = p.yz;
            d = max(d, -sBoxS(q, sc.yz/8., 0.));        
            q = mod(q - sc.yz/12., sc.yz/3.) - sc.xy/3.;
            d = max(d, -sBoxS(q , sc.yz/8., 0.));*/
        }
        else d = max(d, -max(abs(q.z) - .25/16., -oD));
    
    }
    
 
    // Object ID.
    gOID = d<innerD && d<fl? 0 : innerD<fl? 1 : 2;
    

    // Minimum distance.
    d = min(fl, d);


    // If we've hit the light object, add some glow.
    if(max(winLight, innerD - .1)<d){
         float rnd = hash31(id + .06);
         vec3 gCol = .5 + .47*cos(6.2831*rnd/12. + vec3(0, 1.2, 2) + 1.);//vec3(1, .2, .1)
         glow += gCol/(.001 + dot(p, p)*256.); // Truchet cable lights.
    }


    // Saving the box dimensions and local coordinates.
    //gSc = vec3(sc.x, h, sc.z);
    gSc = vec3(sc.x, sc.y, sc.z);
    gID = id;
    gP = p;
    
 
    // Return the minimum distance.
    return d;
}




float rayMarch(vec3 ro, vec3 rd, out int iNum) {
    
    float d, t = hash31(ro + rd + fract(TIME))*.05;
    //const float tol = TOLERANCE;
    vec2 dt = vec2(1e8, 0);
    int i = 0;
    
    glow = vec3(0);


    // Set the global ray direction varibles -- Used to calculate
    // the cell boundary distance inside the "map" function.
    gDir = sign(rd)*.5;
    gRd = rd; 

    const int iter = 160;
    for (i = 0; i<iter; i++) {
       
        d = map(ro + rd*t);
        // IQ's clever edge desparkle trick. :)
        if (d<dt.x) { dt = vec2(d, t); } 

        if (d < .001*(1. + t*.05) || t > FAR) { break; }

        t += min(d*.7, gCD);
        //t += min(min(d*.7, gCD), .1);
    }
    
    
    //t = min(t, FAR);
    
    if(i == iter - 1) { t = dt.y; }

    iNum = i;

    return min(t, FAR);
}

float softShadow(in vec3 ro, in vec3 ld, float lDist, in float k) {
    
    float sh = 1.;
    float t = 0.;
    
    /*
    #if LIGHT_TYPE == 1
    ld = (ld - ro);
    float lDist = length(ld);
    ld /= max(lDist, .0001);
    #endif
    */

    // Set the global ray direction varibles -- Used to calculate
    // the cell boundary distance inside the "map" function.
    gDir = sign(ld)*.5;
    gRd = ld; 

    for (int i=0; i<48; ++i) {
        float d = map(ro + ld*t);
        sh = min(sh, k*d/t);
        if (sh<0. || t>lDist) break;
        t += clamp(min(d, gCD), .01, .5);
    }
    return max(sh, 0.);
}

// Texture bump mapping. Four tri-planar lookups, or 12 texture lookups in total. 
// I tried to make it as concise as possible. Whether that translates to speed, 
// or not, I couldn't say.
vec3 texBump( sampler2D tx, in vec3 p, in vec3 n, float bf){
   
    const vec2 e = vec2(.001, 0);
    
    // Three gradient vectors rolled into a matrix, constructed with offset greyscale 
    // texture values.    
    mat3 m = mat3(tex3D(tx, p - e.xyy, n), tex3D(tx, p - e.yxy, n), 
                  tex3D(tx, p - e.yyx, n));
    
    vec3 g = vec3(.299, .587, .114)*m; // Converting to greyscale.
    g = (g - dot(tex3D(tx,  p , n), vec3(.299, .587, .114)))/e.x; 
    
    // Adjusting the tangent vector so that it's perpendicular to the normal -- Thanks 
    // to EvilRyu for reminding me why we perform this step. It's been a while, but I 
    // vaguely recall that it's some kind of orthogonal space fix using the Gram-Schmidt 
    // process. However, all you need to know is that it works. :)
    g -= n*dot(n, g);
                      
    return normalize( n + g*bf ); // Bumped normal. "bf" - bump factor.
	
} 


// Normal function. It's not as fast as the tetrahedral calculation, but more symmetrical.
vec3 normal(in vec3 p) {
	
    //return normalize(vec3(map(p + e.xyy) - map(p - e.xyy), map
    //                      map(p + e.yxy) - map(p - e.yxy),	
    //                      map(p + e.yyx) - map(p - e.yyx)));
    
    // This mess is an attempt to speed up compiler time by contriving a break... It's 
    // based on a suggestion by IQ. I think it works, but I really couldn't say for sure.
    float sgn = 1.;
    vec3 e = vec3(.002, 0, 0), mp = e.zzz; // Spalmer's clever zeroing.
    for(int i = min(FRAMECOUNT, 0); i<6; i++){
		mp.x += map(p + sgn*e)*sgn;
        sgn = -sgn;
        if((i&1)==1){ mp = mp.yzx; e = e.zxy; }
    }
    
    return normalize(mp);
}



// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)
// Anyway, I like this one. I'm assuming it's based on IQ's original.
float calcAO(in vec3 p, in vec3 n){
    

	float sca = 2., occ = 0.;
    for(int i = 0; i<5; i++){
    
        float hr = float(i + 1)*.15/5.;        
        float d = map(p + n*hr);
        occ += max(hr - d, 0.)*sca;
        sca *= .75;
    }
    
    return clamp(1. - occ, 0., 1.);    
    
}

///////////////////////////
const float PI = 3.14159265;

// Microfaceted normal distribution function.
float D_GGX(float NoH, float roughness) {
    float alpha = pow(roughness, 4.);
    float b = (NoH*NoH*(alpha - 1.) + 1.);
    return alpha/(PI*b*b);
}

// Surface geometry function.
float G1_GGX_Schlick(float NoV, float roughness) {
    //float r = roughness; // original
    float r = .5 + .5*roughness; // Disney remapping.
    float k = (r*r)/2.;
    float denom = NoV*(1. - k) + k;
    return max(NoV, .001)/denom;
}

float G_Smith(float NoV, float NoL, float roughness) {
    float g1_l = G1_GGX_Schlick(NoL, roughness);
    float g1_v = G1_GGX_Schlick(NoV, roughness);
    return g1_l*g1_v;
}

// Bidirectional Reflectance Distribution Function (BRDF). 
//
// If you want a quick crash course in BRDF, see the following:
// Microfacet BRDF: Theory and Implementation of Basic PBR Materials
// https://www.youtube.com/watch?v=gya7x9H3mV0&t=730s
//
vec3 BRDF(vec3 col, vec3 n, vec3 l, vec3 v, 
          float type, float rough, float fresRef){
     
  vec3 h = normalize(v + l); // Half vector.

  // Standard BRDF dot product calculations.
  float nv = clamp(dot(n, v), 0., 1.);
  float nl = clamp(dot(n, l), 0., 1.);
  float nh = clamp(dot(n, h), 0., 1.);
  float vh = clamp(dot(v, h), 0., 1.);   

  
  // Specular microfacet (Cook- Torrance) BRDF.
  //
  // F0 for dielectics in range [0., .16] 
  // Default FO is (.16 * .5^2) = .04
  // Common Fresnel values, F(0), or F0 here.
  // Water: .02, Plastic: .05, Glass: .08, Diamond: .17
  // Copper: vec3(.95, .64, .54), Aluminium: vec3(.91, .92, .92), Gold: vec3(1, .88, .6).
  vec3 f0 = vec3(.16*(fresRef*fresRef)); 
  // For metals, the base color is used for F0.
  f0 = mix(f0, col, type);
  vec3 F = f0 + (1. - f0)*pow(1. - vh, 5.);  // Fresnel-Schlick reflected light term.
  // Microfacet distribution... Roughness.
  float D = D_GGX(nh, rough); 
  // Geometry self shadowing term.
  float G = G_Smith(nv, nl, rough); 
  // Combining the terms above.
  vec3 spec = F*D*G/(4.*max(nv, .001));
  
  
  // Diffuse calculations.
  vec3 diff = vec3(nl);
  diff *= 1. - F; // If not specular, use as diffuse (optional).
  diff *= (1. - type); // No diffuse for metals.

  
  // Combining diffuse and specular.
  // You could specify a specular color, multiply it by the base
  // color, or multiply by a constant. It's up to you.
  return col*(diff + vec3(.2, .4, 1)*spec*PI);
  
}
////////////////////


vec3 render(vec3 ro, vec3 rd, vec3 light, inout float tt){
 
  
  int iter;
  float t = rayMarch(ro, rd, iter);
  tt = t;
  

  // Saving the cell ID, cell scale, and local coordinates.
  vec3 svID = gID;
  vec3 svSc = gSc;
  vec3 svP = gP;
  
  // Glow.
  vec3 svGlow = glow;
  
  // Object ID.
  int svOID = gOID;
  
  // Dark fog color.
  vec3 fog = mix(vec3(.7, 1, .5), vec3(.5, .7, 1), rd.y*.5 + .5)/64.;
  
  // Diferent color scheme fog.
  #if SCHEME == 1
  fog = fog.zxy*4.;
  #endif
  
  // Initialize the scene background to the fog.
  vec3 col = fog;
   
  
  // If we've hit an object, color it in.
  if (t < FAR) {
  
    // Hit point and normal.
    vec3 p = ro + rd*t;
    vec3 n = normal(p);
    
    
    // Texture base bump mapping.
    n = texBump(image3, p*2. + svID/4.*0., n, .005);///(1. + t/FAR)
    
    
    // Light type: Point light or directional light.
    #if LIGHT_TYPE == 1
    vec3 ld = light - p; // Point light.
    float lDist = length(ld);
    ld /= max(lDist, .001);
    float atten = 1./(1. + lDist*lDist*.5);
    #else
    float lDist = FAR; 
    vec3 ld = normalize(light); // Directional light.
    float sDist = length(ld*FAR - t*rd)/FAR;
    float atten = 1./(1. + sDist*sDist*.5);
    #endif
    
    
    // Shadow and AO.
    float shd = softShadow(p + n*.0015, ld, lDist, 8.);
    float ao = calcAO(p, n);
    
    // Texturing the surface.
    float txSc = svOID==0? 1. : 1./3.;
    vec3 tx = tex3D(image3, (p + svID/4.*0.)*txSc, n);
  
   
    // Material properties.
    // Adding in some very artificial microfaceted surface roughness.
    float matType = 1.; // Metallic.
    float roughness = min(dot(tx, vec3(.299, .587, .114))*3., 1.);
    float reflectance = .5;
 
    
    // Main object color. Some kind of copper, I guess.
    vec3 oCol = vec3(.8, .6, .4)/3.;
    
    
    #if SCHEME == 1
    oCol = mix(oCol.yzx, oCol.yyy, .8); // Different color scheme.
    #endif
    
    // Floor color.
    if(svOID==2) oCol /= 2.;
    
    // Window light color. Obviously brighter.
    if(svOID==1) oCol = vec3(4);
  
   
    // Adding some texture to the color.
    oCol *= tx*2. + .05;
    
    // Randomly shading various objects.
    float h = hash31(svID + .12);
    oCol *= h*.5 + .5;
    
 
    // Specular reflection. 
    vec3 hv = normalize(ld - rd); // Half vector.
    vec3 ref = reflect(rd, n); // Surface reflection.
    vec3 refTx = texture(iChannel3, ref).xyz; refTx *= refTx; // Cube map.
    float spRef = pow(max(dot(hv, n), 0.), 5.); // Specular reflection.
    float rf = 1.;//mix(.5, 1., 1. - smoothstep(0., .01, d + .02));
    rf *= (svOID == 0)? 32. : .1;
    oCol += oCol*spRef*refTx.xxx*rf*1.;
    
    
    // BRDF lighting.
    vec3 ct = BRDF(oCol, n, ld, -rd, matType, roughness, reflectance);
    
    // Ambient light.
    #if LIGHT_TYPE == 1
    float amb = .4;
    #else 
    float amb = .3; // A little less ambient light for direct lighting.
    #endif
    
    // Applying lighting and shadows.
    col = (oCol*amb*(shd*.5 + .5) + ct*shd); 
    
     
    #if SCHEME == 1
    svGlow = svGlow.yzx; // Alternate color scheme.
    #endif
    
    // Adding the window glow.
    col += oCol*mix(svGlow.zyx, svGlow.yxz, min(length(svP/gSc), 1.))*4.;
    
    // AO and attenuation.
    col *= ao*atten; 
    

  }

  // Applying the fog.
  col = mix(col, fog, smoothstep(.3, .99, t/FAR));
  
  // Scene color.
  return col;
  
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

 
  vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;
 
  
  	// Camera Setup.
    float speed = 1.;
    float tm = TIME*speed;
    vec3 lk = path(tm + .25) + vec3(0, -.05, 0);  // "Look At" position. Look down a bit.
    vec3 ro = path(tm);// + vec3(0, 1.5, 0); // Camera position, doubling as the ray origin.
    
    #if LIGHT_TYPE == 1
    vec3 lp = path(tm + 3.25); // Light position, somewhere near the moving camera.
   // Light position offset. Redundant here.
    vec3 loffs =  vec3(0, 0, 0);
    //vec2 a = sin(vec2(1.57, 0) - lp.z*1.57/10.);
    //loffs.xy = mat2(a, -a.y, a.x)*loffs.xy; 
    lp += loffs;
    #else
    vec3 lp = -vec3(1, -4, -2); // Directional light option.
    #endif
    

    // Using the above to produce the unit ray-direction vector.
    float FOV = 1.; // FOV - Field of view.
    vec3 fwd = normalize(lk - ro);
    vec3 rgt = normalize(vec3(fwd.z, 0., -fwd.x )); 
    vec3 up = cross(fwd, rgt);

    // Unit direction ray.
    vec3 rd = normalize(fwd/FOV + uv.x*rgt + uv.y*up);
    // Lens distortion.
    //vec3 r = fwd + FOV*(u.x*rgt + u.y*up);
    //r = normalize(vec3(r.xy, (r.z - length(r.xy)*.25)));
    
    // Swiveling the camera from left to right when turning corners.
    rd.xy = rot2(path(lk.z).x/24.)*rd.xy; 

    float t;
    vec3 col = render(ro, rd, lp, t);
    
    // Applying a bit of last minute Uncharted 2 toning.
    col = uTone(col.xyzx).xyz;
    //col = col/(1. + col/2.5); // Rough extended Reinhard toning.

    
    
    #if 0
    // Mix the previous frames in with no camera reprojection.
    // It's OK, but full temporal blur will be experienced.
    vec4 preCol = texelFetch(BuffA, ivec2(fragCoord), 0);
    float blend = (FRAMECOUNT < 2) ? 1. : 1./2.; 
    fragColor = mix(preCol, vec4(clamp(col, 0., 1.), t), blend);
    #else
    fragColor = vec4(max(col, 0.), t);
    #endif
	return fragColor; 
 } 



/*

    Offworld Storage Facility
    -------------------------

    See "Buffer A" for an explanation.

*/


/*

// Just a very basic depth of field routine -- I find a lot of it is
// common sense. Basically, you store the scene distance from the camera 
// in the fourth channel, then use it to determine how blurry you want
// your image to be at that particular distance.
//
// For instance, in this case, I want pixels that are 2.25 units away from 
// the camera to be in focus (not blurred) and for things to get more 
// blurry as you move away from that point -- aptly named the focal point 
// for non camera people. :)
//
// I based this on old code of mine, but adopted things that I found in 
// IQ and Nesvi7's examples, which you can find here:
//
// Ladybug - IQ
// https://www.shadertoy.com/view/4tByz3
//
// Cube surface II - Nesvi7
// https://www.shadertoy.com/view/Mty3DV
//
vec3 DpthFld(sampler2D iCh, vec2 uv){
	
    // Focal point and circle of confusion.
    const float focD = 3., coc = 2.5;
    // Linear distance from either side of the focal point.
    float l = abs(texture(iCh, uv).w - focD - coc) - coc;
    // Using it to calculate the DOF.
    float dof = clamp(l/coc, 0., 2.); 
    
    // Combine samples. Samples with a larger DOF value are taken further 
    // away from the original point, and as such appear blurrier.
    vec3 acc = vec3(0);

    for(int i = 0; i<25; i++){
        // Accumulate samples.
        acc += texture(iCh, uv + (vec2(i/5, i%5) - 2.)/vec2(800, 450)*dof).xyz;
        //acc.x *= dof/2.;
    }

    // Return the new variably blurred value.
    return acc /= 25.;
    // Visual debug representation of DOF value.
    //return vec3(length(dof)*450./2.5);
}

*/


// IQ's float to float hash. I've added an extra sine wrapping modulo to
// cater for my annoying AMD based system, which can't wrap sine with a 
// proper degree of accuracy.
float hash11B(float x){ return fract(sin(mod(x, 6.2831853))*43758.5453); }

// This is an amalgamation of old blur and DOF functions with a heap of borrowed 
// lines from Dave Hoskins's much nicer Fibonacci based "Bokeh disc" function, which 
// you can find here: https://www.shadertoy.com/view/4d2Xzw
//
// If you're interested in bokeh, Dave's function above and some of Shadertoy user, 
// Hornet's, are probably the one's you should be looking at. Xor has some cool simple 
// ones on here too.
//
vec4 bloom(sampler2D iCh, vec2 uv){


    // UV based DOF. Focused on the horizontal line, then blurring further away.
    //float r = smoothstep(0., 1., abs(uv.y - .57)/.57)*2.;
    // Focal point and circle of confusion.
    const float focD = 3., coc = 2.5;
    // Linear distance from either side of the focal point.
    float l = abs(texture(iCh, uv).w - focD - coc) - coc;
    // Using it to calculate the DOF.
    float ra = clamp(l/coc, 0., 2.);
    //float ra = smoothstep(.1, .9, abs(uv.y - .5)*2.)*2.;
    //float ra = mix(clamp(l/coc, 0., 2.), smoothstep(.3, 1., abs(uv.y - .5)*2.), .25);
    //float ra = (smoothstep(.2, 1., length(uv - .5)));

 

    // Standard Fibonacci distribution calculations, compliments of Dave Hoskins.
    const int iter = 96;
    float rad = 2.5;//max(2.*ra, .5); // Bokeh radius.
    float r = 1.;
	vec4 tot = vec4(0), div = tot;
    vec2 vangle = vec2(0., rad*.01/sqrt(float(iter)));
    #define GA 2.3999632 // Golden angle.
    const mat2 rot = mat2(cos(GA), sin(GA), -sin(GA), cos(GA));

    // Aspect ratio.
    vec2 aspect = vec2(RENDERSIZE.y/RENDERSIZE.x, 1);
    
    
	for (int i = 0; i<iter; i++){
        
        #if 1
        
        // Dave Hoskin's Fibonacci based scattering. Cheaper and much nicer, so
        // it's set as the default.
        // The approx increase in the scale of sqrt(0, 1, 2, 3...).
        r += 1./r;
	    vangle = rot*vangle;
        vec4 col = texture(iCh, uv - (r - 1.)*vangle*aspect, RENDERSIZE.y/450.*1.); 
        
        #else
        
        // A hash based random distribution, for anyone who wants to try it.
        //int ii = i%10; // For square bokeh.
        //int jj = i/10;
    
        // Random offset contained within a disk or radius n.
        float fi = float(i) + fract(TIME);
        //vec2 fi = vec2(ii, jj) - 5. + fract(TIME);
        vec2 rnd2 = vec2(hash11B(fi), hash11B(fi + .1)*6.2831);
        vec2 offs = 6.*sqrt(rnd2.x)*vec2(cos(rnd2.y), sin(rnd2.y));
        offs *= rad;
        //offs = rad*(offs + (vec2(hash11B(fi), hash11B(fi + .21)) - .5));
        vec4 col = texture(iCh, uv - offs/RENDERSIZE.xy, RENDERSIZE.y/450.*1.);  
  
        #endif
         
        // Thanks to Dave for figuring out how to tweak the colors to produce brighter 
        // contrast. It's common sense... once someone figures it out for you. :D 
        //
        // Linear falloff -- Not always necessary, but I need things to fade out
        // toward the edges. Nonlinear falloff is possible too.
        col *= 1. - float(i)/float(iter); 
        vec4 bokeh = col*col;
		tot += bokeh*bokeh;
		div += bokeh;
        
	}
    
    
    // Mixing the original value with the bokeh tweaked value according
    // to the depth of field and a hand tweaked brightness factor.
    vec4 colOrig = texture(iCh, uv);
	return colOrig + tot/div*2.*ra; //mix(colOrig, colOrig + tot/div*2., ra);///
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord{ = _xy;



    // Screen oordinates.
    vec2 uv = fragCoord/RENDERSIZE.xy;

    // Uncomment the "DpthFld" function first.
    //vec4 col = DpthFld(BuffA, uv).xyzz;

    
    // A bloom function. A lot of it comes from Dave Hoskins's bokeh function.
    // See the function above for a link to the original.
    vec4 col = bloom(BuffA, uv);
    
    // Retrieving the stored color only.
    //vec4 col = texture(BuffA, uv);


    // Subtle vignette.
    col *= pow(16.*uv.x*uv.y*(1. - uv.x)*(1. - uv.y) , 1./32.);

    // Rough gamma correction and screen presentation.
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