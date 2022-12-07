

/*

    Extruded Semi-regular Tiling
    ----------------------------    

    A while ago, after looking at Fizzer's "Wythoff Uniform Tilings + Duals" example, I 
    got curious as to what one of the semi-regular tilings would look like in extruded
    form. The variation I chose was a 3.4.6.4 semi-regular tessellation, which is formed 
    using a triangle, square and regular hexagon. I like it because it resembles 
    overlapping dodecahedrons. I'm not sure if the configuration has a special name, but 
    it's very common.
    
    Raymarching extruded grid tiles in realtime requires more planning and trickery than 
    their 2D counterparts. Luckily, the variation I chose came together fairly easily. 
    The extruded heightfield in planar form was kind of interesting, but a little 
    underwhelming, so I decided to jazz it up a bit by using a few cheap demoscene 
    cliches -- Warped space, cheap transparent tubes to emulate light rays, etc.
    
    I had the cathedral scene from Farbrausch's "fr-08, .the .product" 64k demo in the 
    back of my mind, when making this, for anyone who still remembers that. This is the 
    geometric alien chamber version... It's a very loose connection. :)
    
    As an aside, I went searching for examples of extruded semi-regular pattern imagery 
    and found none, which surprised me, since extruded regular prismatic square, hexagon, 
    etc, and simple variations are commonplace amongst the Blender, stock imagery crowd, 
    and so forth. Semi-regular variations are a natural extension of the aforementioned,
    so I figured it'd be a common thing. There are probably examples out there, but it
    was clear that the graphics crowd haven't adopted the more interesting variations yet.
    
    Anyway, for anyone ever in need of one of these, here's the code. I've explained the 
    general process inside the distance function. I'll leave the cell by cell traversal 
    as an exercise for the reader. :)


    
    Related examples:
    
    // There are a heap of standard semi-regular tilings in this example, and
    // their duals, which is very handy when trying to visualize things.
    Wythoff Uniform Tilings + Duals  - fizzer
    https://www.shadertoy.com/view/3tyXWw
    
    // Hyperbolic semi-regular tilings. A really nice example. As an aside,
    // I have an extruded hyperbolic regular tiled example somewhere.
    Wythoffian Tiling Generator - mla
    https://www.shadertoy.com/view/wlGSWc

    
*/

  
// Scene theme -- Timber: 0, Metallic\Green: 1, Metallic\Purple: 2.
#define SCENE 0

// I put this in as an afterthought, just to show that it could be done. At 
// present, it doesn't look very interesting, but I intend to produce a tailored
// example with a different tiling arrangement later.
//#define CYLINDRICAL

// Global tile scale -- Left over from an old shader, so it needs to
// be left alone... I'll tidy it up later.
const vec2 scale = vec2(1./8.);

// Max ray distance.
#define FAR 30.


// Scene object ID.
float objID;


// Standard 2D rotation formula.
mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }


// IQ's vec2 to float hash.
float hash21(vec2 p){  
    #ifdef CYLINDRICAL
    // Wrapping for a cylinder.
    p.x = mod(p.x, 1./scale.x);
    #endif
    return fract(sin(dot(p, vec2(27.619, 57.583)))*43758.5453); 
}

// IQ's vec2 to float hash.
float hash31(vec3 p){  
    return fract(sin(dot(p, vec3(113.619, 57.583, 27.897)))*43758.5453); 
}


// Tri-Planar blending function. Based on an old Nvidia tutorial by Ryan Geiss.
vec3 tex3D(sampler2D t, in vec3 p, in vec3 n){ 
    
    n = max(abs(n) - .2, .001); // max(abs(n), 0.001), etc.
    //n /= dot(n, vec3(.8)); 
    n /= length(n);
    
    // Texure samples. One for each plane.
    vec3 tx = texture(t, p.yz).xyz;
    vec3 ty = texture(t, p.zx).xyz;
    vec3 tz = texture(t, p.xy).xyz;
    
    // Multiply each texture plane by its normal dominance factor.... or however you wish
    // to describe it. For instance, if the normal faces up or down, the "ty" texture sample,
    // represnting the XZ plane, will be used, which makes sense.
    
    // Textures are stored in sRGB (I think), so you have to convert them to linear space 
    // (squaring is a rough approximation) prior to working with them... or something like that. :)
    // Once the final color value is gamma corrected, you should see correct looking colors.
    return mat3(tx*tx, ty*ty, tz*tz)*n; // Equivalent to: tx*tx*n.x + ty*ty*n.y + tz*tz*n.z;

}

// Height map value.
float hm(in vec2 p){ 

    // Reading into "Buffer A".
    // Stretching to account for the varying buffer size.
    //p *= vec2(RENDERSIZE.y/RENDERSIZE.x, 1);
    return texture(image6, p + .5).x;
    
}


// IQ's extrusion formula.
float opExtrusion(in float sdf, in float pz, in float h){
    
    vec2 w = vec2( sdf, abs(pz) - h );
  	return min(max(w.x, w.y), 0.) + length(max(w, 0.));

    /*
    // Slight rounding. A little nicer, but slower.
    const float sf = .015;
    vec2 w = vec2( sdf, abs(pz) - h - sf/2.);
  	return min(max(w.x, w.y), 0.) + length(max(w + sf, 0.)) - sf;
    */
}

//////////////////

// Standard polar partitioning.
vec2 polRot(vec2 p, inout float na, float aN){

    float a = atan(p.y, p.x);
    na = mod(floor(a/6.2831*aN), aN);
    float ia = (na + .5)/aN;
    p *= rot2(-ia*6.2831);

    return p;
}


// Flat top hexagon scaling.
const vec2 s = vec2(1.7320508, 1);

// Hexagon edge and vertex IDs. They're useful for neighboring edge comparisons,
// etc. Multiplying them by "s" gives the actual vertex postion.
//
// Vertices and edges: Clockwise from the left.
//
// Note that these are all six times larger than usual. We're doing this to 
// get rid of decimal places, especially those that involve division by three.
// I't a common accuracy hack. Unfortunately, "1. - 1./3." is not always the 
// same as "2./3." on a GPU.

// Multiplied by 12 to give integer entries only.
const vec2[6] vID = vec2[6](vec2(-4, 0), vec2(-2, 6), vec2(2, 6), 
                      vec2(4, 0), vec2(2, -6), vec2(-2, -6)); 

const vec2[6] eID = vec2[6](vec2(-3, 3), vec2(0, 6), vec2(3), 
                      vec2(3, -3), vec2(0, -6), vec2(-3));


// Signed distance to a regular hexagon, with a hacky smoothing variable thrown
// in. -- It's based off of IQ's more exact pentagon method.
float getHex(in vec2 p, float r, in float sf){
    
      // Flat top.
      const vec3 k = vec3(-.8660254, .5, .57735); // pi/6: cos, sin, tan.
      // Pointed top.
      //const vec3 k = vec3(.5, -.8660254, .57735); // pi/6: cos, sin, tan.
     
      // X and Y reflection.  
      p = abs(p); 
      p -= 2.*min(dot(k.xy, p), 0.)*k.xy;

      r -= sf;
      // Polygon side.
      // Flat top.
      return length(p - vec2(clamp(p.x, -k.z*r, k.z*r), r))*sign(p.y - r) - sf;
      // Pointed top.
      //return length(p - vec2(r, clamp(p.y, -k.z*r, k.z*r)))*sign(p.x - r) - sf;
    
}

// IQ;s signed distance to an equilateral triangle.
// https://www.shadertoy.com/view/Xl2yDW
float getTri(in vec2 p, in float r){

    const float k = sqrt(3.0);
    p.x = abs(p.x) - r;
    p.y = p.y + r/k;
    if(p.x + k*p.y>0.) p = vec2(p.x - k*p.y, -k*p.x - p.y)/2.;
    p.x -= clamp(p.x, -2.*r, 0.);
    return -length(p)*sign(p.y);
}

// IQ's signed box formula.
float getSq(in vec2 p, in vec2 b, in float sf){

  vec2 d = abs(p) - b + sf;
  return min(max(d.x, d.y), 0.) + length(max(d, 0.)) - sf;
}


// Hexagonal grid coordinates. This returns the local coordinates and the cell's center.
// The process is explained in more detail here:
//
// Minimal Hexagon Grid - Shane
// https://www.shadertoy.com/view/Xljczw
//
vec4 hexGrid(vec2 p){

    // Extra scaling to wrap the pattern around a cylinder along the longer X-axis.
    // If we were using the shorter one, it wouldn't matter.
    //p *= s.x; 
   //p = p.yx;
    vec4 ip = floor(vec4(p/s, p/s - .5)) + .5;
    vec4 q = p.xyxy - vec4(ip.xy, ip.zw + .5)*s.xyxy;
    return dot(q.xy, q.xy)<dot(q.zw, q.zw)? vec4(q.xy, ip.xy) : vec4(q.zw, ip.zw + .5);
   
} 

// This mess is to account for cylindrical mapping, when used, but
// it's not really needed here.
const float gridScale = 3./3.14159*1.*s.x;//3./3.14159*2.;

float getHeight(vec3 c){ return dot(c, vec3(.299, .587, .114))*.3; }

vec3 getColor(vec2 p){
    
    // The ID is in index form (multiplied by 12), so
    // needs to be multiplied by "s". Plus there's the global scale, "gridScale" to consider.
    // Viewport stretching, if needed.
    // p *= vec2(RENDERSIZE.y/RENDERSIZE.x, 1);

    p.x = mod(p.x, 1./scale.x);
    
    #if SCENE == 0
    vec3 tx = texture(image6, (p*s/gridScale/12./4. + .5), -100.).xyz; // Timber.
    #else
    vec3 tx = texture(image3, (p*s/gridScale/12./4. + .5), -100.).xyz; // Metallic.
    #endif
    return tx*tx;
    
    //float rnd = hash21(p);
    //return .5 + .45*cos(6.2831*rnd + vec3(0, 1, 2)*1.5);
}


// Z-axis twist. It's another demoscene cliche. Return mat2(1, 0, 0, 1) to see the
// non-warped planes.
mat2 mTwist(float z){ return rot2(z*.35); }

// Space warping.
vec3 coordTrans(vec3 p){
    
    #ifdef CYLINDRICAL
    // Cylindrical transform.
    const float aN = 1.;
    //p.xy *= rot2(3.14159/aN);
    float a = atan(p.x, p.y); // XY-plane pixel angle.
    //float ia = floor(a/6.2831*aN);
    //ia = (ia + .5)/aN;
    //p.xy *= rot2(ia*6.2831);
    // XY-plane polar coordinate.
    p.xy = vec2(mod(a, 6.2831/aN) - 6.2831/aN/2., length(p.xy));
    // Moving the radial surface out.
    p.y = .5 - p.y;
    #else
    // Twisted planes transform. It's a demoscene favorite,
    // since it's simple cheap and cool looking.
    p.xy *= mTwist(p.z);
    // Breaking up the symmetry between the top and bottom planes.
    //p.xz += p.y<.0? vec2(0) : vec2(0, 4.5);
    // Splitting into two planes.
    p.y = .5 - abs(p.y);
    #endif
    
    return p;

}

//////////////////


// A global to stored the pylon's 2D face field value, the pylon height,
// and the face's central position based ID.
vec4 gVal;
vec3 glow, lCol; // Global glow and glow color variables.

// The scene.
float map(vec3 p){

    
    // Apply the space transform.
    p = coordTrans(p);
    
    // Floor -- I'm leaving it in, but we don't need it here.
    float fl = 1e5;//(p.y + 32.) - 32.005; // Adding a touch of thickness to alleviate artifacts.
 
 
    // The extruded 3.4.6.4 semi-regular tiling. First things first, pick any vertex on
    // the plane, find the face with the lowest number of vertices (in this case a triangle),
    // then proceed a full revolution counting the face vertices as you go along. You'll see
    // that it's triangle, square, hexagon, square, or 3.4.6.4. That's how mathematicians
    // decided to name these things.
    //
    // Anyway, in order to raymarch any extruded tiling pattern, you need to find a way to
    // render overlaying grids of repeat objects that don't touch. It can be confusing at first,
    // but you get used to it. For instance, in this pattern, the repeat hexagons don't touch
    // one another, so you can render those. Note that the triangles don't touch either, so you
    // can render those in polar form. That leaves the squares, which you'll see do touch in
    // the corners. The trick there is to split them into two groups of three that are spaced
    // apart... Yeah, as mentioned it's confusing at first, but you get used to it.
    //
    // The following is just an application of the above. I've renderd four interwoven repeat
    // grids... which I complicated by including glow. However, for anyone wishing to do this,
    // start from scratch, and use this as a guide, if need be.
    
    
    // Hexagon grid. Returns the local coordinate and hexagon ID.
    vec4 h = hexGrid(p.xz*gridScale);
    
    // Local coordinates.
    vec2 q = h.xy;
    // Face ID. It's multiplied by 12 to avoid GPU accuracy errors when dividing by 3...
    // I've explained it before, but basically, "1./3." and "1. - 2./3." are not the same
    // on a GPU due to rounding errors. However, "4." and "12. - 8." are equal.
    vec2 id = h.zw*12.;
    
    // Face heights, 2D distance field values and associated 3D prism field values. 
    vec4 hgt4, d2D4, d3D4;
    mat4x2 mID; // Four face IDs.
    
    mat4x2 svP; // Saving local coordinates.
    
    float r = sqrt(1./7.)*.8660254;
    float ew = .01*gridScale;
    

    // Central hexagon values.
    d2D4.x = getHex(q, r - ew, 0.);  
    d2D4.x /= gridScale;
    svP[0] = q;
    mID[0] = id;
    vec3 hCol = getColor(mID[0]);
     // Extruding the 2D field.
    hgt4.x = getHeight(hCol);
    
    
    // Surrounding triangle values.
    float na;
    q = rot2(3.14159/6.)*h.xy;
    q = polRot(q, na, 6.);
    q.x -= .5/.8660254;
    q *= rot2(3.14159/6.);
    d2D4.y = getTri(q, r*.57735 - ew); // Triangle.
    d2D4.y /= gridScale;
    svP[1] = q;
    //
    // Moving the polar index to the correct starting position, since the polar index
    // represents a different part of the hexagon. Thankfully, that didn't take me
    // too long to figure out. :)
    int ind = int(mod(8. - na, 6.));
    mID[1] = id + vID[(ind)%6];
    vec3 tCol = getColor(mID[1]);//vec3(.25, .6, 1)
    hgt4.y = getHeight(tCol);
 
    
    // Three spaced out surrounding square values.
    q = rot2(-3.14159/6.)*h.xy;
    q = polRot(q, na, 3.);
    q.x -= .5;
    //d2D4.z = max(abs(q.x), abs(q.y)) - r*.57735 + ew; // Square one.
    d2D4.z = getSq(q, vec2(r*.57735 - ew), 0.); // Square one.
    d2D4.z /= gridScale;
    svP[2] = q;
    //
    ind = int(mod(8. - na*2., 6.));
    mID[2] = id + eID[(ind)%6];
    vec3 sCol1 = getColor(mID[2]);//vec3(1, .9, .1)
    //
    hgt4.z =  getHeight(sCol1);


    // The other three spaced out surrounding square values.
    q = rot2(3.14159/6.)*h.xy;
    q = polRot(q, na, 3.);
    q.x -= .5;
    //
    d2D4.w = getSq(q, vec2(r*.57735 - ew), 0.); // Square two.
    d2D4.w /= gridScale;
    svP[3] = q;
    //
    ind = int(mod(7. - na*2., 6.)); // One less rotation, so one less index.
    mID[3] = id + eID[(ind)%6];
    vec3 sCol2 = getColor(mID[3]);//vec3(.5, 1, .2)  
    //
    hgt4.w = getHeight(sCol2);

    // Glow colum variables.
    float dd = 1e5, cLight = 1e5;
    vec2 lgtID = vec2(0);
    
    
    // Tapering ratio, relative to shape area.
    //vec4 sL = vec4(1, 1./4.5, 1./18.,  1./18.);
    //vec4 sL = vec4(1, .57735, .57735/2., .57735/2.);
   

    float d = 1e5;
    for(int i = 0; i<4; i++){
         
        // Random face value.
        float rndI = hash21(mID[i] + .05);
        // Column type.
        int type = 0;
        float sv2D = d2D4.x; // Face distance value.
        if(rndI<.35){ 
           
            type = 1; // Non glowing, but flat.
            // Random columns not running through the center.
            if(rndI<.25 && abs(mID[i].x)>1.){
                type = 2; // Flag it to glow.
            }
            
            // For glowing or flat faced column, we only want the face edge.
            d2D4.x = abs(d2D4.x + .02) - .02;
         
        }
        
        // The 3D extruded face. 
        d3D4.x = opExtrusion(d2D4.x, p.y - (hgt4.x/2. - 8.), hgt4.x/2. + 8.);
        d3D4.x += d2D4.x*.25; // Raised tops.
        // Not glowing, so cap off the hole... I was going to open and close
        // the top to stream light through, but I ran out of steam. I might
        // arrange for that later though.
        if(type==1) d3D4.x = min(d3D4.x, max(sv2D, abs(p.y + .5) - .5)); 
        
        // Minimum pylon distance.
        if(d3D4.x<d){

           d = d3D4.x;
           // Saving global pylon values for later use.
           gVal.x = d2D4.x; // 2D face field.
           gVal.y = hgt4.x; // Pylon height.
           gVal.zw = mID[i]; // Pylong ID.

        }
        
        
        // Constructing the light beams eminating from various grid shapes.
        if(type == 2){
            float cL = sv2D + .02; // Light colum, the same shape as the face. 
            //float cL = length(svP[i]) - pw.x*2.;//length(svP[i]) - .1;//
            cL /= gridScale;
            // Taper off the light rays as they get further from the source. The smaller windows
            // are tapered off by a factor relative to the square of the side length.
            cL += p.y*p.y*.2;//*sL[i]; 
            //cL *= 1./(1. + p.y*p.y*8.);    // Alternative way to do it.
            // Save the closest beam and beam ID.
            if(cL<dd){
                dd = cL;
                lgtID = mID[i];
            }
        }
        
        // Cyberjax's loop shuffling trick, which saves considerably on compiler time.
        d3D4 = d3D4.yzwx;
        d2D4 = d2D4.yzwx;
        hgt4 = hgt4.yzwx;
        
    }
    
    // Overall object ID.
    objID = fl<d? 1. : 0.;
    
    // Combining the floor with the extruded object. Redundant here.
    d = min(fl, d);


    // The glow color, which is just a smoothstepped 2D distance field
    // column. I don't see glow calculated this way, but it's the way
    // I prefer to do it. The glow calculation doesn't effect the scene 
    // geometry, but is passed on to the distance function for accumulation.
    //
    // No glow.
    lCol = vec3(0);
    // If the pylon is flagged to pass a glow colum through it, then 
    // calculate the glow field, and color.
    if(dd<d){
        dd = max(-dd, 0.);
        vec3 cCol = .5 + .45*cos(6.2831*hash21(lgtID + .12)/4. + vec3(0, 1, 2));
        //lCol = cCol.yxz*clamp(.0001/(.001 + dd*dd*4.), 0., 1.);
        lCol = cCol*smoothstep(.0, .25, dd); //max(dd, 0.)*8.;
    }
        
    
    
    // Return the minimum scene distance.
    return  d;
 
}

// 3D noise texture. It's way faster than a handcoded one.
vec3 n3DT(vec3 p){ return texture(iChannel3, p).xyz; }

// Basic raymarcher.
float trace(in vec3 ro, in vec3 rd){

    // Overall ray distance and scene distance.
    float t = hash31(fract(ro/7.319) + rd)*.1, d;
    //float t = 0., d;
    
    glow = vec3(0);
    
    for(int i = min(0, FRAMECOUNT); i<128; i++){
    
        vec3 p = ro + rd*t;
        d = map(p);

        // Note the "t*b + a" addition. Basically, we're putting less emphasis on accuracy, as
        // "t" increases. It's a cheap trick that works in most situations... Not all, though.
        if((abs(d)<.001) || t>FAR) break; // Alternative: 0.001*max(t*.25, 1.), etc.
       
        // Accumulate the glow color.
        glow += lCol;///(1. + t);
        
        // Note that the ray is capped (to .1). It's slower, but is necessary for the
        // glow to work. I guess it could also help with overstepping the mark a bit.
        t += min(d*.7, .1); 
    }

    // Minimum distance. Since FAR is theoretically the furthest the ray can do, 
    // technically, this should be done. Depending on the situation, it can sometimes 
    // help avoid far plane artifacts as well.
    return min(t, FAR);
}


// Normal function. It's not as fast as the tetrahedral calculation, but more symmetrical.
vec3 getNormal(in vec3 p, float t) {
	
    const vec2 e = vec2(.001, 0);
    
    //return normalize(vec3(m(p + e.xyy) - m(p - e.xyy), m(p + e.yxy) - m(p - e.yxy),	
    //                      m(p + e.yyx) - m(p - e.yyx)));
    
    // This mess is an attempt to speed up compiler time by contriving a break... It's 
    // based on a suggestion by IQ. I think it works, but I really couldn't say for sure.
    float sgn = 1.;
    float mp[6];
    vec3[3] e6 = vec3[3](e.xyy, e.yxy, e.yyx);
    for(int i = min(FRAMECOUNT, 0); i<6; i++){
		mp[i] = map(p + sgn*e6[i/2]);
        sgn = -sgn;
        if(sgn>2.) break; // Fake conditional break;
    }
    
    return normalize(vec3(mp[0] - mp[1], mp[2] - mp[3], mp[4] - mp[5]));
}



// Cheap shadows are hard. In fact, I'd almost say, shadowing particular scenes with limited 
// iterations is impossible... However, I'd be very grateful if someone could prove me wrong. :)
float softShadow(vec3 ro, vec3 lp, vec3 n, float k){

    // More would be nicer. More is always nicer, but not always affordable. :)
    const int maxIterationsShad = 32; 
    
    ro += n*.0015; // Coincides with the hit condition in the "trace" function.  
    vec3 rd = lp - ro; // Unnormalized direction ray.

    float shade = 1.;
    float t = 0.; 
    float end = max(length(rd), .0001);
    //float stepDist = end/float(maxIterationsShad);
    rd /= end;

    // Max shadow iterations - More iterations make nicer shadows, but slow things down. Obviously, 
    // the lowest number to give a decent shadow is the best one to choose. 
    for (int i = min(FRAMECOUNT, 0); i<maxIterationsShad; i++){

        float d = map(ro + rd*t);
        shade = min(shade, k*d/t);
        //shade = min(shade, smoothstep(0., 1., k*d/t)); // Thanks to IQ for this tidbit.
        // So many options here, and none are perfect: dist += clamp(d, .01, stepDist), etc.
        t += clamp(d, .005, .15); 
        
        
        // Early exits from accumulative distance function calls tend to be a good thing.
        if (d<0. || t>end) break; 
    }

    // Sometimes, I'll add a constant to the final shade value, which lightens the shadow a bit --
    // It's a preference thing. Really dark shadows look too brutal to me. Sometimes, I'll add 
    // AO also just for kicks. :)
    return max(shade, 0.); 
}


// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)
// Anyway, I like this one. I'm assuming it's based on IQ's original.
float calcAO(in vec3 p, in vec3 n){

	float sca = 1., occ = 0.;
    for( int i = 0; i<5; i++ ){
    
        float hr = float(i + 1)*.15/5.;        
        float d = map(p + n*hr);
        occ += (hr - d)*sca;
        sca *= .7;
    }
    
    return clamp(1. - occ, 0., 1.);  
}



vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    
    // Screen coordinates.
	vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;

    
	// Camera Setup.
	vec3 lk = vec3(0, 0, TIME/2.);//vec3(0, -.25, TIME);  // "Look At" position.
    vec3 camH = vec3(0, 0, -1);
    camH.xy *= mTwist(-(lk.z + camH.z));
    vec3 ro = lk + camH; // Camera position, doubling as the ray origin.
 
    // Light positioning. One is just in front of the camera, and the other is in front of that.
    vec3 lpH = vec3(0, 0, .75);
    #ifdef CYLINDRICAL
    lpH = vec3(0, 0, 1.);
    #endif
    lpH.xy *= mTwist(-(lk.z + lpH.z));
 	vec3 lp = lk + lpH;// Put it a bit in front of the camera.
	
    
    // Using the above to produce the unit ray-direction vector.
    //
    float FOV = 1.; // FOV - Field of view.
    vec3 fwd = normalize(lk - ro); // Forward.
    vec3 rgt = normalize(cross(vec3(0, 1, 0), fwd));// Right. 
    // "right" and "forward" are perpendicular normals, so the result is normalized.
    vec3 up = cross(fwd, rgt); // Up.
    // Camera.
    //mat3 mCam = mat3(rgt, up, fwd);
    // rd - Ray direction.
    //vec3 rd = mCam*normalize(vec3(uv, 1./FOV));//
    vec3 rd = normalize(uv.x*rgt + uv.y*up + fwd/FOV);
    
    
    // Substle camera direction movement.
    rd.xy *= rot2(-TIME/32.);
    rd.yz *= rot2(sin(TIME/8.)*.05);
 
    
    // Raymarch to the scene.
    float t = trace(ro, rd);
    
    // Save the object ID.
    float svObjID = objID;
    
    vec4 svVal = gVal;
  
	
    // Initiate the scene color to black.
	vec3 col = vec3(0);
	
	// The ray has effectively hit the surface, so light it up.
	if(t < FAR){
        
  	
    	// Surface position and surface normal.
	    vec3 sp = ro + rd*t;
        vec3 sn = getNormal(sp, t);
        
            	// Light direction vector.
	    vec3 ld = lp - sp;

        // Distance from respective light to the surface point.
	    float lDist = max(length(ld), .001);
    	
    	// Normalize the light direction vector.
	    ld /= lDist;

        
        
        // Shadows and ambient self shadowing.
    	float sh = softShadow(sp, lp, sn, 16.);
    	float ao = calcAO(sp, sn); // Ambient occlusion.
        //sh = min(sh + ao*.25, 1.);
	    
	    // Light attenuation, based on the distances above.
	    float atten = 1./(1. + lDist*lDist*.1);

    	
    	// Diffuse lighting.
	    float diff = max( dot(sn, ld), 0.);
        //diff = pow(diff, 4.)*2.; // Ramping up the diffuse.
    	
    	// Specular lighting.
	    float spec = pow(max(dot(reflect(ld, sn), rd ), 0.), 32.); 
	    
	    // Fresnel term. Good for giving a surface a bit of a reflective glow.
        //float fre = pow(clamp(1. - abs(dot(sn, rd))*.5, 0., 1.), 2.);
        
		// Schlick approximation. I use it to tone down the specular term. It's pretty subtle,
        // so could almost be aproximated by a constant, but I prefer it. Here, it's being
        // used to give a hard clay consistency... It "kind of" works.
		float Schlick = pow( 1. - max(dot(rd, normalize(rd + ld)), 0.), 5.);
		float freS = mix(.15, 1., Schlick);  //F0 = .2 - Glass... or close enough. 
        
          
        // Texture coordinates.
        vec3 txP = sp;
       
        // Matching the coordinate space warping in the 
        // distance function with the texture coordinates.
        //
        // Apply the space transform.
        txP = coordTrans(txP);
        
        float sf = 1.5/RENDERSIZE.y; // Smoothing factor.
        const float ew = .001; // Edge width.
        
         // Texel color. 
	    vec3 texCol = vec3(0); 
       

        
        if(svObjID<.5){
            
            
            // The extruded pattern.
            
            // Texture color lookup, based on the central ID of the
            // extruded prism we've hit. In this case, it'll be a hexagon,
            // square or triangle.
            texCol = getColor(svVal.zw)*2.;
            
            // Apply some edging.
            float d = abs(svVal.x) - ew;
            float h = svVal.y; // Object height.
            d = max(d, abs(txP.y - h) - ew);
            
            d *= 1./(1. + t*.25); // Taper the edging effect further away to reduce banding.
            
            // Apply the edge.
            texCol = mix(texCol, texCol*.1, 1. - smoothstep(0., sf, d)); 
            
            // Tri-planar texture lookup.
            vec3 tx = tex3D(image3, sp, sn);
            
            // Combine for the final surface color.
            texCol *= tx*2.;
           
            // Ramping up the diffuse for a more metallic look.
            //diff = pow(diff, 4.)*2.;
      
 
        }
        else {
            
            // The floor. Redundant here, but sometimes the 
            // background floor can be visible.
             
            
            // Background.
            texCol = vec3(0);
            
            /*
            vec3 tx = tex3D(image3, sp, sn);
            texCol *= tx*2.;
            // Ground rim.
            float d = abs(svVal.x) - ew*1.5;
            texCol = mix(texCol, vec3(0), 1. - smoothstep(0., sf, d));
            */
          
        }
        
        // Specular reflection.
        vec3 hv = normalize(-rd + ld); // Half vector.
        vec3 ref = reflect(rd, sn); // Surface reflection.
        vec3 refTx = texture(iChannel2, ref).xyz; refTx *= refTx;
        refTx = (texCol*1.5 + .66)*refTx;//smoothstep(.2, .5, refTx);
        float spRef = pow(max(dot(hv, sn), 0.), 8.); // Specular reflection.
        float rf = (svObjID == 1.)? .25 : 1.;//mix(.5, 1., 1. - smoothstep(0., .01, d + .08));
        texCol += spRef*refTx*rf*2.; //smoothstep(.03, 1., spRef) 


        
        // Combining the above terms to produce the final color.
        col = texCol*(diff*sh + .2 + vec3(1, .97, .92)*spec*freS*2.*sh);
        
        //col += col.xxx*spRef*refTx*rf*16.; //smoothstep(.03, 1., spRef) 
        
        // Arranging for some fine particles to stream into the chamber by
        // applying some fine dust to the rays... Technically, this should be
        // processed inside the raymarching loop, but this is way cheaper and
        // still conveys the general idea.
        vec3 nSp = vec3(8, 4, 8)*txP + vec3(1, -4, 1)*TIME/8.;
        vec3 ns = mix(n3DT(nSp), n3DT(nSp*2.), 1./3.);
        mix(ns, n3DT(nSp*4.), 1./7.);
        glow = vec3(1)*glow*(ns*.75 + .25);//*hash31(sp/vec3(1, 8, 1) + tm);// dot(glow, vec3(.299, .587, .114))
        
        // Apply the glow to the warped planes.
        // It doesn't quite look right with a straight cylinder, but I think I'll
        // make a specific tunnel example later.
        #ifndef CYLINDRICAL
        
            #if SCENE == 0
            // Warm hues.
            glow = mix(glow, glow.yxz, max(txP.y, 0.));
            glow = mix(glow, glow.xxx, .5);
            col = col/2. + col*glow*32.;
            #else
            // Cliche Borg green.
            glow = mix(glow.yxz, glow.zyx, min(txP.y*txP.y*12., 1.)); 
                #if SCENE == 2
                // Purple.
                glow = glow.xzy; 
                //glow = mix(glow.xzy, glow.yzx, floor(hash21(svVal.zw)*4.999)/4.);
                #endif
            col = col/3. + col*glow*32.;
            #endif
        
        #endif
 
      
        // Shading.
        col *= ao*atten;
        
        
        
        // It's sometimes helpful to check things like shadows and AO by themselves.
        //col = vec3(ao);
          
	
	}
    
    // Horizon fog. Not visible here, but provided for completeness.
    col = mix(col, vec3(.5, .7, 1)*0., smoothstep(0., .99, t/FAR));
          
    
    // Rought gamma correction.
	fragColor = vec4(sqrt(max(col, 0.)), 1);
	
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}