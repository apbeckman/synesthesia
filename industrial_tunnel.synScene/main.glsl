/*

	Industrial Tunnel
	-----------------

	Using repeat cylindrical coordinates (polar plus Z) to construct a steel-plated tunnel
	with an industrial feel. I think I've mentioned it before, but I never quite appreciate
	just how well Dr2's examples run until I try to construct a scene with a fraction of the
	detail... then watch the FPS counter plummet. :)

	My fastest machine pushes this scene out with ease. However, my slowest machine didn't
	fare as well. I'd set myself a limit of 40 FPS, but could only manage 30 FPS... for now,
	anyway. I could get to 40 FPS, but would really have to mess the code up, so I've left
	it alone. Obviously, the key to faster scenes is insuring that the distance function is
	easy to hone-in on and involves fewer instructions. Unfortunately, that's not so easy to
	achieve when adding more detail. On top of that, my slow machine isn't a fan of the "atan"
	function, which is pretty difficult to avoid when doing angular-based calculations.

	The scene itself is very loosely based on a section of the Greenwich pedestrian tunnel.
	The lighting is not very realistic, but works well enough. The scene relies on a heavy
	curvature setting to give it a quasi hand painted look. Anyway, I'd like to put up a
	simpler abstract geometric version at some stage that'll be easier to read, and
	hopefully, faster.


	// Other examples:

	// Cylindrical-coordinate-based tunnel. I've seen it pop up in many forms on the net.
	Tunnel #1 -  WAHa_06x36
	https://www.shadertoy.com/view/4dfGDr

	// Great example making use of polar coordinates. Really stylish and great atmosphere.
	Metro Tunnel - fb39ca4
	https://www.shadertoy.com/view/ldsGRS

	// Very cool. Amazingly detailed for a shader.
	Gotthard Tunnel - dr2
    https://www.shadertoy.com/view/MlSXRR



*/

// Maximum ray distance.
#define FAR 50. // MS: Can't do fog calculations here because FAR is also used for the ray tracer

// Comment this out to omit the detailing. Basically, the bump mapping won't be included.
#define SHOW_DETAILS

// Object ID, used for the gold trimming in the bump mapping section.
float svObjID;
float objID;

// MS: Object IDs
#define TUN 0. // Tunnel: Basically, the metal plates.
#define FLR 1. // Floor:  Concrete curbs and underflooring.
#define MTL 2. // Metal:  The metallic mesh and the light casings.
#define LGT 3. // Lights: The lights attached to the top mesh.
#define BLT 4. // Bolts:  The hexagonal bolts.
#define PIP 5. // Pipes:  The cyclinders beside the lights.

// 2D rotation. Always handy. Angle vector, courtesy of Fabrice.
mat2 rot( float th ){ vec2 a = sin(vec2(1.5707963, 0) + th); return mat2(a, -a.y, a.x); }


// Compact, self-contained version of IQ's 3D value noise function.
float n3D(vec3 p){

	const vec3 s = vec3(7, 157, 113);
	vec3 ip = floor(p); p -= ip;
    vec4 h = vec4(0., s.yz, s.y + s.z) + dot(ip, s);
    p = p*p*(3. - 2.*p); //p *= p*p*(p*(p * 6. - 15.) + 10.);
    h = mix(fract(sin(h)*43758.5453), fract(sin(h + s.x)*43758.5453), p.x);
    h.xy = mix(h.xz, h.yw, p.y);
    return mix(h.x, h.y, p.z); // Range: [0, 1].
}

// Layered noise.
float fBm(vec3 p){
    return n3D(p)*.57 + n3D(p*2.)*.28 + n3D(p*4.)*.15;
    //return n3D(p)*.5333 + n3D(p*2.)*.2667 + n3D(p*4.)*.1333 + n3D(p*8.)*.0667;
}

// Camera path. Arranged to coincide with the frequency of the lattice.
vec3 camPath(float t){

    //return vec3(0, 0, t); // Straight path.

    // Curvy path. Weaving around the columns.
    float a = sin(t*3.14159265/16. + 1.5707963*1.);
    float b = cos(t*3.14159265/16.);

    // MS: Mix these together so we can slide in between straight and curved path
    float mixA = mix(0, 2*a, Curve_Path);
    float mixB = mix(0, b*a, Curve_Path);

    return vec3(mixA, mixB, t);
}

// Minimum - with corresponding object ID.
vec2 objMin(vec2 a, vec2 b){

    // Returning the minimum distance along with the ID of the
    // object. This is one way to do it. There are others.

    //Equivalent to: return a.x < b.x ? a: b;
    float s = step(a.x, b.x);
    return s*a + (1. - s)*b;
}

// The tunnel scene. There's a bit of code, but it's nothing more than a bunch of boxes
// and cylinders spread out around some repeat cylindrical coordinates.
float map(vec3 p){

    float depth = .1 * syn_BassPresence; // Depth of the rounded metal plates.

    // Mold the scene around the path.
    p.xy -= camPath(p.z).xy;

    p.y += .1*Wave_Height*sin(p.z + syn_Time*.5);

    // The edge of the tunnel. Set at a radius of one, plus the depth.
    float tun = (1. + depth) - length(p.xy);

    //////////////

    // The concrete floor. Just a plane, with a circular chunk cut out. It gives it
    // a curb-like appearance.
    float flr = p.y + .695;
    flr = max(flr, tun - depth - .1);

    ///////////////

    // The tunnel walls. Due to the amount of detailing and the polar conversion, it
    // looks more complicated than it actually is. To repeat across the XY plane we
    // do something along the lines of "p.xz = mod(p.xz, c) - c/2." To repeat around
    // a circle in the XY plane, we convert to polar coordinates, "p.xy = rot(angle),"
    // (angle is based on "atan(p.y, p.x)," then do the same thing. The rest is basic
    // unit circle trigonometry, etc. By the way, this is a rough description, so if
    // something doesn't quite make sense, it probably doesn't. :)

    // Converting the XY plane to polar coordinates. I'm handling the panels (six per
    // circle) and the bolts (18 per circle) at the same time to share some calculations.
    // I'd love to use the geometry of one to constuct the other - in order to save
    // some instructions, but I'm leaving it alone for now.
    vec3 q = p;
    vec3 q2 = p;

    float a = atan(q.y, q.x)/6.2831853; // Polar angle of "p.xy" coordinate.
    float ia = (floor(a*6.) + .5)/6.*6.2831853; // Angle between "PI/6" intervals.
    float ia2 = (floor(a*18.) + .5)/18.*6.2831853; // Angle between "PI/18" intervals.

     // Polar conversion for 6 segments, but offset every second panel... and shifted
    // to the center-cell position (see the Z-repetition).
    q.xy *= rot(ia + sign(mod(q.z + .25, 1.) - .5)*3.14159/18.);
    q2.xy *= rot(ia2); // Polar conversion for 18 segments (for the bolts).

    q2.xy *= Remove_Bolts;

    // The X-coordinate is now the radial coordinate, which radiates from the center
    // to infinity. We want to break it into cells that are 2 units wide, but centered
    // in the middle. The result is that the panels will start at radius one.
    q.x = mod(q.x, 2.) - 1.;
    // Plain old linear Z repetion. We want the panels and bolts to be repeated in the
    // Z-direction (down the tunnel) every half unit.
    q.z = mod(q.z, .5) - .25;

    // Moving the bolts out to a distance of 2.1.
    q2.x = mod(q2.x, (2. + depth)) - (2. + depth)/2.;

    // Now, it's just a case of drawing an positioning some basic shapes. Boxes and
    // tubes with a hexagonal cross-section.
    q = abs(q);
    q2 = abs(q2);

    // Bolts. Hexagon shapes spaced out eighteen times around the tunnel walls. The
    // panels are spaced out in sixths, so that means three per panel.
    float blt = max(max(q2.x*.866025 + q2.y*.5, q2.y) - .02, q.z - .08);


    // Putting in some extra rails where the mesh and curb meets the tunnel. The extra
    // code is fiddly (not to mention, slows things down), but it makes the structure
    // look a little neater.
    q2 = p;
    q2.xy *= rot(ia - sign(p.x)*3.14159/18.);
    q2 = abs(q2);

    // Lines and gaps on the tunnel to give the illusion of metal plating.
    float tunDetail = max(min(min(q.y - .06, q.z - .06), max(q2.y - .06, p.y)),
                          -min(min(q.y - .01, q.z - .01), max(q2.y - .01, p.y)));

    // Adding the tunnel details (with a circular center taken out) to the tunnel.
    tun = min(tun, max(tunDetail, tun-depth));

    ///////////////

    // The metalic mesh elements and light casings. The lights are calculated in this
    // block too.
        // The metalic mesh elements and light casings. The lights are calculated in this
    // block too.
/*
    q = abs(p);
    float mtl = max(q.x - .14, abs(p.y - .88) - .02);  // Top mesh.
    mtl = min(mtl, max(q.x - .396, abs(p.y + .82) - .02)); // Bottom mesh.//.81

    q.z = abs(mod(p.z, 2.) - 1.);

    float lgt = max(max(q.x - .07, q.z - .07), abs(p.y - 1.) - .255);
    float casings = max(max(q.x - .1, q.z - .1), abs(p.y - 1.) - .23);

    q.xz = abs(mod(q.xz, 1./8.) - .5/8.);

    mtl = max(mtl, -max(max(q.x - .045, q.z - .045), abs(abs(p.x) - .19) - .14)); // Holes in the mesh.
    mtl = min(mtl, casings ); // Add the light casings to the top mesh.
*/

    // Alternative mesh setup with smaller holes. I like it more, but Moire patterns are a problem
    // with smaller window sizes.
    q = abs(p);
    float mtl = max(q.x - .13, abs(p.y - .88) - .02);  // Top mesh.
    mtl = min(mtl, max(q.x - .396, abs(p.y + .82) - .02)); // Bottom mesh.//.81

    q.z = abs(mod(p.z, 2.) - 1.);

    float lgt = max(max(q.x - .07, q.z - .07), abs(p.y - 1.) - .255);
    float casings = max(max(q.x - .1, q.z - .1), abs(p.y - 1.) - .23);

    q.xz = abs(mod(q.xz, 1./16.) - .5/16.);

    mtl = max(mtl, -max(max(q.x - .025, q.z - .025), abs(abs(p.x) - .19) - .155)); // Holes in the mesh.
    mtl = min(mtl, casings ); // Add the light casings to the top mesh.

    ///////////////

    // Pipes. Electricity... water? Not sure what their function is, but I figured I
    // should slow the distance function down even more, so put some in. :)
    q = p;
    const float th = 6.283/18.;
    float sx = sign(p.x);
    float pip = length(q.xy - vec2(sx*sin(th*1.4), cos(th*1.4))*1.05) - .015;
    pip = min(pip, length(q.xy - vec2(sx*sin(th*1.6), cos(th*1.6))*1.05) - .015);

    ///////////////

    // Determine the overall closest object and its corresponding object ID. There's a way
    // to save some cycles and take the object-ID calculations out of the distance function,
    // but I'm leaving them here for simplicity.
    vec2 d = objMin(vec2(tun, TUN), vec2(blt, BLT));
    d = objMin(d, vec2(mtl, MTL));
    d = objMin(d, vec2(lgt, LGT));
    d = objMin(d, vec2(flr, FLR));
	d = objMin(d, vec2(pip, PIP));

    ///////////////


    objID = d.y; // Set the global object ID.

    return d.x; // Return the closest distance.


}


// Raymarching.
float trace(vec3 ro, vec3 rd){

    float t = 0., d;

    int numSteps = int( 96 - 92*pow(Ray_Steps,.25) );
    for (int i=0; i<numSteps; i++){

        d = map(ro + rd*t);
        if(abs(d)<.001*(t*.125 + 1.) || t>FAR) break;
        t += d;
    }
    return min(t, FAR);
}

// Tri-Planar blending function. Based on an old Nvidia writeup:
// GPU Gems 3 - Ryan Geiss: http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
vec3 tex3D(sampler2D channel, vec3 p, vec3 n){

    n = max(abs(n) - .2, 0.001);
    n /= dot(n, vec3(1));
	vec3 tx = texture(channel, p.yz).xyz;
    vec3 ty = texture(channel, p.xz).xyz;
    vec3 tz = texture(channel, p.xy).xyz;

    // Textures are stored in sRGB (I think), so you have to convert them to linear space
    // (squaring is a rough approximation) prior to working with them... or something like that. :)
    // Once the final color value is gamma corrected, you should see correct looking colors.
    return tx*tx*n.x + ty*ty*n.y + tz*tz*n.z;
}



// The bump mapping function.
float bumpFunction(in vec3 p){

    p.xy -= camPath(p.z).xy;

    // Adding a bit of functional noise variations to the
    // concrete floor and the tunnel. Not a great deal of
    // effort went into it.
    float res = 0.;
    if(svObjID==FLR) {
        p.xy *= 16.;
        res = n3D(p*4.)*.66 + n3D(p*8.)*.34;
        res = 1.-abs(res - .75)/.75;
    }
    else if(svObjID==TUN){
        //res = fBm(p*16.);
        res = n3D(p*16.)*.66 + n3D(p*32.)*.34;
    }


    // Subtle metal bump. More thought needs to be put into it. :)
    if(svObjID==MTL){

        p.xz = abs(mod(p.xz + 1./8.,  1./8.) - .5/8.);
        res = max(p.x, p.z) - .25/8.;

        res = max(res*8., 0.);

        //res = 1. - smoothstep(0., .5, res);

    }



    return res; // Range: [0, 1].

}


// Standard function-based bump mapping function.
vec3 doBumpMap(in vec3 p, in vec3 nor, float bumpfactor){

    const vec2 e = vec2(0.001, 0);
    float ref = bumpFunction(p + Surface_Time*.01*sin(syn_HighTime*.05));
    vec3 grad = (vec3(bumpFunction(p - e.xyy),
                      bumpFunction(p - e.yxy),
                      bumpFunction(p - e.yyx) )-ref)/e.x;

    grad -= nor*dot(nor, grad);

    return normalize( nor + grad*bumpfactor );

}


// Texture bump mapping. Four tri-planar lookups, or 12 texture lookups in total. I tried to
// make it as concise as possible. Whether that translates to speed, or not, I couldn't say.
vec3 texBump( sampler2D tx, in vec3 p, in vec3 n, float bf){

    const vec2 e = vec2(0.001, 0);

    // Three gradient vectors rolled into a matrix, constructed with offset greyscale texture values.
    mat3 m = mat3( tex3D(tx, p - e.xyy, n), tex3D(tx, p - e.yxy, n), tex3D(tx, p - e.yyx, n));

    vec3 g = vec3(0.299, 0.587, 0.114)*m; // Converting to greyscale.
    g = (g - dot(tex3D(tx,  p , n), vec3(0.299, 0.587, 0.114)) )/e.x; g -= n*dot(n, g);

    return normalize( n + g*bf ); // Bumped normal. "bf" - bump factor.

}


// The normal function with some curvature rolled into it. Sometimes, it's possible to get away
// with six taps, but we need a bit of epsilon value variance here, so there's an extra six.
vec3 nr(vec3 p, inout float crv, float ef){
	//ef/RENDERSIZE.y
    vec2 e = vec2(ef/450., 0); // Larger epsilon for greater sample spread, thus thicker edges.

    // Take some distance function measurements from either side of the hit point on all three axes.
	float d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	float d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	float d5 = map(p + e.yyx), d6 = map(p - e.yyx);
	float d = map(p);	// The hit point itself - Doubled to cut down on calculations. See below.

    // Seven-tap curvature calculation. You can get away with four taps, but this is a little
    // more accurate.
    crv = clamp((d1 + d2 + d3 + d4 + d5 + d6 - d*6.)*32. + .5, 0., 1.);

    // Redoing the calculations for the normal with a more precise epsilon value.
    e = vec2(.002, 0);
	d1 = map(p + e.xyy), d2 = map(p - e.xyy);
	d3 = map(p + e.yxy), d4 = map(p - e.yxy);
	d5 = map(p + e.yyx), d6 = map(p - e.yyx);

    // Return the normal.
    // Standard, normalized gradient mearsurement.
    return normalize(vec3(d1 - d2, d3 - d4, d5 - d6));
}


// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)
// Anyway, I like this one. I'm assuming it's based on IQ's original.
float cao(in vec3 p, in vec3 n){

    float sca = 1., occ = 0.;
    for(float i=0.; i<5.; i++){

        float hr = .01 + i*.5/4.;
        float dd = map(n * hr + p);
        occ += (hr - dd)*sca;
        sca *= 0.7;
    }
    return clamp(1.0 - occ - pow(Darkness,.3), 0., 1.);
}


// Cheap shadows are hard. In fact, I'd almost say, shadowing particular scenes with limited
// iterations is impossible... However, I'd be very grateful if someone could prove me wrong. :)
float softShadow(vec3 ro, vec3 lp, float k){

    // More would be nicer. More is always nicer, but not really affordable... Not on my slow test machine, anyway.
    const int maxIterationsShad = 60;

    vec3 rd = (lp-ro); // Unnormalized direction ray.

    float shade = 1.0;
    float dist = 0.05;
    float end = max(length(rd), 0.001);
    //float stepDist = end/float(maxIterationsShad);

    rd /= end;

    // Max shadow iterations - More iterations make nicer shadows, but slow things down. Obviously, the lowest
    // number to give a decent shadow is the best one to choose.
    for (int i=0; i<maxIterationsShad; i++){

        float h = map(ro + rd*dist);
        shade = min(shade, k*h/dist);
        //shade = min(shade, smoothstep(0.0, 1.0, k*h/dist)); // Subtle difference. Thanks to IQ for this tidbit.
        //dist += min( h, stepDist ); // So many options here: dist += clamp( h, 0.0005, 0.2 ), etc.
        dist += clamp(h, 0.01, 0.25);

        // Early exits from accumulative distance function calls tend to be a good thing.
        if (h<0.001 || dist > end) break;
    }

    // I've added 0.5 to the final shade value, which lightens the shadow a bit. It's a preference thing.
    return min(max(shade, 0.) + 0.2, 1.0);
}


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    // Screen coordinates.
	vec2 u = (fragCoord - RENDERSIZE.xy*0.5)/RENDERSIZE.y;

	// Camera Setup.
    const float speed = 2.;
    float camTime = bass_time;
    vec3 ro = camPath(camTime*speed); // Camera position, doubling as the ray origin.
    ro.xy += LookXY * -.35;


    vec3 lk = camPath(camTime*speed + .5);  // "Look At" position.
    lk.xy += LookXY;
    //ro.y -= .1; // Hack to lower the camera.
    //lk.y -= .1;

    vec3 lp = camPath(camTime*speed + 2.25); // Light position, somewhere near the moving camera.
    lp.y += .6;
    lp.z -= 5 * bass_lights * syn_BassLevel;
    // Alternative. Moving the light to the right a bit.
    //float th = 6.283*1./12.;
    //lp += vec3(sin(th)*.6, cos(th)*.6, 0);


    // Using the above to produce the unit ray-direction vector.
    float FOV = .75 * FOV; // FOV - Field of view.
    vec3 fwd = normalize(lk-ro);
    vec3 rgt = normalize(vec3(fwd.z, 0., -fwd.x ));
    vec3 up = cross(fwd, rgt);

    // Unit direction ray.
    //vec3 rd = normalize(fwd + FOV*(u.x*rgt + u.y*up));

    // Mild lens distortion to fit more of the scene in, and to mix things up a little.
    vec3 rd = fwd + FOV*(u.x*rgt + u.y*up);
    rd = normalize(vec3(rd.xy, (rd.z - length(rd.xy)*.3)*.7));

    // Swiveling the camera from left to right when turning corners.
    float swivel = camPath(lk.z).x;
    rd.xy = rot(swivel/48. )*rd.xy;
    rd.xz = rot(swivel/16. )*rd.xz;


    // Raymarch.
    float t = trace(ro, rd);
    svObjID = objID;

    // Surface hit point.
    vec3 sp = ro + rd*t;

    // Offset path. Only used for one thing here, but it's handy to have an
    // adjusted hit point that aligns with the path contour.
    vec3 pathSp = sp - camPath(sp.z);


    // Normal with curvature component.
    float crv = 1., ef = 12.; // ef - Edge and curvature factor.
    vec3 sn = nr(sp, crv, ef);

    // Shadows and ambient self shadowing.
    float sh = softShadow(sp, lp, 16.); // Soft shadows.
    float ao = cao(sp, sn); // Ambient occlusion.

    // Light direction vector setup and light to surface distance.
    lp -= sp;
    float lDist = max(length(lp), .0001);
    lp /= lDist;

    // Attenuation.
    float atten = 1./(1.0 + lDist*.25 + lDist*lDist*.025);

    // Texturing the object.
    const float tSize0 = 1./2.;
    vec3 texture = tex3D(image3, sp*tSize0, sn);
    texture = smoothstep(0., .5, texture);
    vec3 sexyColor = vec3(1, .2, 0);
    vec3 tx = vec3(0);
    tx = mix(sexyColor, texture, Texture);

    // Ugly "if" statements for object coloring. They do the job though.
	if(svObjID==BLT || svObjID==PIP) tx *= vec3(1.25, 1, .75);
    else if(svObjID==TUN) tx *= vec3(1., .4, .2);
    else if(svObjID==MTL) tx *= vec3(1.1, .8, .7);
    else if(svObjID==FLR) tx *= vec3(1.5, .78, .62);
    else if(svObjID==LGT) tx *= vec3(5, 4, 3); // Really bright for fake glow.


    // More fake lighting. This was just a bit of trial-and-error to produce some repetitive,
    // slightly overhead, spotlights throughout the space. Cylinder in XY, sine repeat
    // in the Z direction over three rows... Something like that. I tried to code in some
    // flickering, but in the end, decided working lights were better. :)

    float spot = max(2. - length(pathSp.xy - vec2(0, 1.)), 0.)*(cos((sp.z + 1.)*3.14159)*.5+.5);
    spot = smoothstep(0.5, 1., spot);
    //float flicker = smoothstep(-1., -.9, sin(TIME*1. + hash(floor(sp.z/2.))*6.283));
    //float flicker = sin(TIME*8. + hash(floor(sp.z/2.))*6.283)*.2 + 1.;
    //spot *= flicker;
    //tx += (tx)*spot*2.;



    // Function bump.
    #ifdef SHOW_DETAILS
    float bf =.005;
    if(svObjID==FLR || svObjID==MTL) bf = .02;
    sn = doBumpMap(sp, sn, bf/(1. + t/FAR));
    // tx *= bumpFunction(sp)*.75 + .25; // Accentuates the bump, but not needed here.
    #endif

    // Texture-based bump mapping.
    float tbf = .03;
    if(svObjID==LGT) tbf = .0;
    sn = texBump(image3, sp*tSize0*2., sn, tbf);

    // Adding a bit of mold build-up to the scene. Very hacky, and in need of a tweak, but it works well enough.
    float slm = fBm(sp*4.);
    float slf = 1.;  // Slime factor.
    if(svObjID!=TUN && svObjID!=FLR) slf = .75;
    tx = mix(vec3(1), vec3(.25, .6, .2), (slm*.75 + .25)*slf)*tx; //
    tx = mix(vec3(1), vec3(.25, .4, .15)*.1, (1.-abs(slm - .5)*2.)*.75*slf)*tx;


    // Diffuse, specular and Fresnel.
    float dif = max(dot(lp, sn), 0.);
    float spe = pow(max(dot(reflect(-lp, sn), -rd), 0.), 32.);
    float fre = pow(clamp(dot(rd, sn) + 1., 0., 1.), 2.);



    // Combining the terms above to produce the final color.
    vec3 fc = tx*(dif + .125 + vec3(1, .8, .5)*fre*8. + vec3(1, .9, .7)*spot*4.) + vec3(1, .7, .5)*spe*1.;
    fc *= atten*sh*ao*clamp(crv*1.5, 0., 1.);


    //fc = mix(vec3(1), vec3(.25, 1, .15)*.2, slm*.75)*fc; // More mold.


    // Mixing in some fog.
    vec3 bg = vec3(.4, .35, .3);
    fc = mix(fc, bg, smoothstep(0., .95, t/(FAR-48.5*pow(Fog,.1))));


    // Post processing.
    //float gr = dot(fc, vec3(.299, .587, .114));
    //fc = fc*.5 + pow(min(vec3(1.5, 1, 1)*gr, 1.), vec3(1, 3, 16))*.5;

     // Approximate gamma correction.
	fragColor = vec4(sqrt(clamp(fc, 0., 1.)), 1.0);
	return fragColor;
 }


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}
