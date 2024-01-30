

/*

    Loxodrome Spiral Projection
    ---------------------------    
    
    I don't have a lot to do with math these days, but I still get a kick out 
    of visual art demonstrations that convey mathematical ideas. Casting 
    shadows onto planes from spherical objects to illustrate various 
    projection patterns is one such example.
    
    In this case, an extruded dashed spherical loxodrome spiral path has been 
    light projected to a back plane to form a Mobius spiral shadow -- I think a 
    loxodrome is a form of Mobius transformation, so I guess that's not 
    surprising. At any rate, in various circles, it's a fairly common 
    demonstration, and there's already at least one such example on Shadertoy, 
    for whch I've provided a link to below.
    
    My first attempts some time ago involved mapping out a spiral onto a sphere 
    using standard spherical coordinates, which kind of produced a Mobius spiral 
    looking shadow, but wasn't what I was after. After inquiring into MLA's 
    mercator coordinate examples, I realized that I wasn't using an angle 
    preserving formula, so once I adjusted for that, everything came together. 
    By the way, if you're interested in spherical spiral mapping, I've provided 
    a link to one of MLA's examples below, which is written in a much nicer way.
    
    I did originally have a very simple spiral projected onto a plain background,
    but since there were already a few of those on the internet, I applied some
    CSG to create dashed lines, then got carried away with glow, etc. The idea to 
    give the inner sphere a different color to the outside came from Kamoshika's 
    "Cubic Shell with Slits" example, which is a interesting in its own right.
    
    I've seen some really nice spherical spiral imagery on the internet that 
    doesn't involve projection, and I'd like to make one of those at some point.
    In regard to shadow projection, I intend to post one of the common polyhedral 
    to Poincare disk presentations. By the way Neozhaoliang has a beautiful 
    example along those lines that I'll provide a link to below.
    
    
    

    // Other examples:
    
    // A very cool shader. Unfortunately, I wrote most of the code to this
    // before looking at MLA's example, which is written in a much nicer way. 
    // The one line I "borrowed" from it is the reason I got my example to 
    // work properly, so I'm thankful anyway. :)
    Indexed Mercator Spiral - mla
    https://www.shadertoy.com/view/NsdBzn
    
    // More of a classy rendering, and reminiscent of Paul Nylander's 
    // classic loxodrome projection image.
    Loxodrome - fb39ca4
    https://www.shadertoy.com/view/MsX3D2

    // Beautiful example.
    Hyperbolic Group Limit Set - Neozhaoliang
    https://www.shadertoy.com/view/NstSDs
    
    // Nice and concise. Kamoshika's other examples are worth
    // looking at too.
    Cubic Shell with Slits - Kamoshika
    https://www.shadertoy.com/view/ftj3WG
    

*/

// Dashed lines. Commenting this out will display the more
// traditional connected spiral arms.
#define DASHED

// Outer and inner glow colors.
#define glowCol1 vec3(3, 1, .2)
#define glowCol2 vec3(1, 3, .2) // vec3(.3, 2.5, 5)

// Maximum ray distance.
#define FAR 27.

// Trying to stop the loop from unrolling to reduce compile times.
#define ZERO min(int(FRAMECOUNT), 0)

// Standard 2D rotation formula.
mat2 rot2(in float a) {
    float c = cos(a), s = sin(a);
    return mat2(c, -s, s, c);
}

// IQ's vec2 to float hash.
float hash21(vec2 p) {
    return fract(sin(dot(p, vec2(27.617, 57.643))) * 43758.5453);
}

// IQ's vec2 to float hash.
float hash31(vec3 p) {
    return fract(sin(dot(p, vec3(113.619, 57.583, 27.897))) * 43758.5453);
}

// Commutative smooth minimum function. Provided by Tomkh, and taken 
// from Alex Evans's (aka Statix) talk: 
// http://media.lolrus.mediamolecule.com/AlexEvans_SIGGRAPH-2015.pdf
// Credited to Dave Smith @media molecule.
float smin(float a, float b, float k) {

    float f = max(0., 1. - abs(b - a) / k);
    return min(a, b) - k * .25 * f * f;
}

// The smooth maximum equivalent.
float smax(float a, float b, float k) {

    float f = max(0., 1. - abs(b - a) / k);
    return max(a, b) + k * .25 * f * f;
}

/*
// Tri-Planar blending function. Based on an old Nvidia tutorial by Ryan Geiss.
vec3 tex3D(sampler2D t, in vec3 p, in vec3 n){ 
    
    n = max(abs(n) - .2, .001); // max(abs(n), 0.001), etc.
    n /= dot(n, vec3(.8)); 
    //n /= length(n);
    
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
*/

// Global object ID and 2D line distance.
float gOID, gLn;

// Global lattitudinal and logitudinal (sphere surface) texture coordinates.
vec2 gSUV;

// Global glow and glow color variables.
vec3 glow, lCol; 

// Loxodromic sphere scene.
float map(vec3 p) {

    // Back plane.
    float pln = -p.z + 2.;//abs(p.z - 2. + .01) - .01;

    // Sphere rotation.
    p.xy *= rot2(PI / 3. + spin_time);
    p.xz *= rot2(-bass_time);
    //p.yz *= rot2(-TIME/4.);
    p.yz *= rot2(mid_time);

    // Sphere coordinates. 
    vec3 sphP = p;
    // Coordinates for the polar caps. Note the abs(Y) to draw both poles with one call. 
    vec3 q = vec3(sphP.x, abs(sphP.y), sphP.z) - vec3(0, 1, 0);

    // The X value controls the number of dashes per revolution and the Y value
    // dictates the number of spiral arms.
    vec2 sSc = vec2(density.x, density.y);
    float ax = atan(sphP.z, sphP.x) / 6.2831; // Longitudinal coordinate.
    float ay = asin(sphP.y / length(sphP));// Or atan(sphP.y, length(sphP.xz)); // Latitude. 
    // The above are only spherical coordinates, which you can arrange to spiral around
    // a sphere. However, to spiral along a spherical rhumb line in an angle preserving
    // manner, the following adjustment needs to be applied. I knew something was necessary,
    // but didn't know what it was. Thankfully, MLA did, so I was able to use his calculation. 
    // I usually prefer to derive these things myself, but I'm going to take it at face value. 
    // I'm guessing it'd be an application of differential surface geometry, or something 
    // along those lines.
    // 
    // Indexed Mercator Spiral - mla
    // https://www.shadertoy.com/view/NsdBzn
    //
    // The following are equivalent. I'm not sure which one the GPU likes best, 
    // but I'm going to go with the most concise one. :)
    ay = asinh(tan(ay)) / (PI*2.0);
    //ay = log(tan(3.14159265/4. + ay/2.))/6.2831; // Mercator
    //ay = sign(ay)*acosh(1./cos(ay))/6.2831;// + .5;

    // The number of spirals per arm.
    float spirals = Spirals; 
    // Spiral lines running perpendicular to one another. You can combine
    // these to form objects, provide texture coordinates, etc.
    vec2 spir = vec2(ax * spirals + ay, ax - ay * spirals);
    gSUV = spir * sSc;

    // Unique dash line UV. Not used here.
    //vec2 ia = floor(gSUV)*sSc + .5;

    // Line thickness. The first line arranges for increased thickness near the poles.
    float th = 1. / (.001 + pow(length(sphP - sign(sphP) * vec3(0, 1, 0)), .5) * 256.);
    th += .04+thickness*0.025; // Constant thickness.

    // Divisor much be a factor of the scale (sSc) above.
    vec2 ln = abs(mod(spir, 1. / sSc) - .5 / sSc) - th;

    // spir = _rotate(spir, TIME);
    #ifdef DASHED
    ln.x = smax(ln.y, -(ln.x - .01), .02); // Using the other spiral to cut out holes.
    #else 
    ln.x = smax(ln.y, -(ln.x + .044), .015); // Closing the gaps.
    //ln.x = ln.y; // The less interesting continuous spiral arms.
    #endif

    // Smoothly blending in the polar caps.
    ln.x = smin(ln.x, length(q) - .05, .1);

    // Saving the 2D dashed line (the rounded rectangle shape) object to a global
    // to be used later.
    gLn = ln.x;

    // Unit sphere.
    float sph = length(p) - 1.;
    // Extruded dashed line distance one.
    float d = max(ln.x, abs(sph) - .01);// + (abs(fract(oAx*64. + TIME*1.) - .5) + .1)*.01*0.
    // Putting a frame around the dashed line.
    float d2 = max(abs(ln.x) - .01, abs(sph) - .02);

    // Object ID: Plane, extruded inner rectangle, extruded rectangle frame.
    gOID = pln < d && pln < d2 ? 0. : d < d2 ? 1. : 2.;

    // Glow tubes.
    lCol = vec3(0);
    if((d2 < .2 + hash31(p) * .1) && gOID > .5) {
       // Glow color. Reddish on the outside, and aqua on the inside.
        vec3 gCol = sph > .0 ? glowCol1 : glowCol2;

       // Record the glow distance for this iteration. "ln" is just a
       // tube running through the yellow inner rectangle in the direction
       // of the normal.
        lCol = gCol * clamp(-(ln.x + .01), 0., .5);
    }

    // Return the minimum distance.
    return min(pln, min(d, d2));
}

// Basic raymarcher.
float trace(in vec3 ro, in vec3 rd) {

    // Overall ray distance and scene distance.
    // Adding some jitter to the jump off point to alleviate banding.
    float t = hash31(fract(ro / 7.319 + TIME) + rd) * .1, d;

    glow = vec3(0);

    for(int i = ZERO; i < 96; i++) {

        d = map(ro + rd * t);
        // Note the "t*b + a" addition. Basically, we're putting less emphasis on accuracy, as
        // "t" increases. It's a cheap trick that works in most situations... Not all, though.
        if(abs(d) < .00025 || t > FAR)
            break; // Alternative: 0.001*max(t*.25, 1.), etc.

        // Accumulate the glow color.
        glow += lCol;///(1. + t);

        // Note that the ray is capped (to .1). It's slower, but is necessary for the
        // glow to work. I guess it could also help with overstepping the mark a bit.
        t += min(d * .5, .1); //d*.5;// 
    }

    return min(t, FAR);
}

// fb39ca4's inverse mix function.
// Inverse mix takes a value between "a" and "b" and maps it to zero to one range.
float invMix(float a, float b, float x) {
    x = (x - a) / (b - a);
    return x * x; // Returning the square for darker tones... My tweak, and not correct.
}

float softShadow(vec3 ro, vec3 lp, vec3 n, float k) {

	// More would be nicer. More is always nicer, but not always affordable. :)
    const int maxIterationsShad = 64;

    ro += n * .0015; // Coincides with the hit condition in the "trace" function.  
    vec3 rd = lp - ro; // Unnormalized direction ray.

    float shade = 1.; // Initialize the shadow to 1., or no shadow.
    float t = 0.;//hash31(fract(ro/7.319) + n)*.01; // Scene distance.
    float maxD = max(length(rd), .0001); // Max light distance.
    //float stepDist = end/float(maxIterationsShad);
    rd /= maxD; // Normalize.

    // Max shadow iterations - More iterations make nicer shadows, but slow things down. Obviously, 
    // the lowest number to give a decent shadow is the best one to choose. 
    for(int i = ZERO; i < 64; i++) {

        float d = map(ro + rd * t); // Distance to the scene.
        #if 1
        // This is a tweak I found in fb39ca4's Loxodrom example. It makes sense,
        // but I'd need to investigate further. The shadows are more succint, but lighter.
        // https://www.shadertoy.com/view/MsX3D2
        float penumbraDist = t / k;
        shade = min(shade, invMix(-penumbraDist, penumbraDist, d));
        t += min((d + penumbraDist) * .5, .2);
        #else
        // IQ's simpler calculation. If feel the shade itself is more constistant, but
        // the shape isn't perfect. Emulating soft shadows isn't easy, if not impossible.
        shade = min(shade, k * d / t);
        //shade = min(shade, smoothstep(0., 1., k*d/t)); // Thanks to IQ for this tidbit.
        t += clamp(d, .005, .2);
        #endif

        // Early exit, and not exceeding the maximum light distance.
        if(d < 0. || t > maxD)
            break;
    }

    shade = max(shade, 0.); // Capping the shadow value above zero.
    //return shade;
    // Another one of fb39ca4's additions. Penumbra stuff. :) 
    shade = shade * 2. - 1.;
    return ((sqrt(1. - shade * shade) * shade + asin(shade)) + 3.14159265 / 2.) / 3.14159265;

}

// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)
// Anyway, I like this one. I'm assuming it's based on IQ's original.
float calcAO(in vec3 p, in vec3 n) {

    float sca = 3., occ = 0.;
    for(int i = ZERO; i < 5; i++) {

        float hr = float(i + 1) * .15 / 5.;
        float d = map(p + n * hr);
        occ += (hr - d) * sca;
        sca *= .7;
        if(d > 1e5)
            break; // Fake break.
    }

    return clamp(1. - occ, 0., 1.);
}

// Normal function. It's not as fast as the tetrahedral calculation, but more symmetrical.
vec3 norm(in vec3 p) {

    const vec2 e = vec2(.001, 0);

    //return normalize(vec3(m(p + e.xyy) - m(p - e.xyy), m(p + e.yxy) - m(p - e.yxy),	
    //                      m(p + e.yyx) - m(p - e.yyx)));

    // This mess is an attempt to speed up compiler time by contriving a break... It's 
    // based on a suggestion by IQ. I think it works, but I really couldn't say for sure.
    float sgn = 1.;
    float mp[6];
    vec3[3] e6 = vec3[3](e.xyy, e.yxy, e.yyx);
    for(int i = ZERO; i < 6; i++) {
        mp[i] = map(p + sgn * e6[i / 2]);
        sgn = -sgn;
        if(sgn > 2.)
            break; // Fake conditional break;
    }

    return normalize(vec3(mp[0] - mp[1], mp[2] - mp[3], mp[4] - mp[5]));
}

vec4 renderMainImage() {
    vec4 fragColor = vec4(0.0);
    vec2 fragCoord = _xy;

    // Pixel coordinates.
    vec2 uv = (fragCoord - RENDERSIZE.xy * .5) / RENDERSIZE.y;
    uv = mix(uv, uv + _uvc, fov);
    // Unit direction ray and ray origin.
    vec3 rd = normalize(vec3(uv, .5));
    vec3 ro = vec3(sin(TIME / 2.) * .1, 0, -1.7);
    float mirror_x = smin(rd.x*-1, rd.x, 0.5);
    float mirror_y = smin(rd.y*-1, rd.y, 0.5);
    vec2 mirror = vec2(x_mirror, y_mirror);
    vec2 rd_mirrored = vec2(mix(mirror_x, mirror_x*-1, invert), mix(mirror_y, mirror_y*-1, invert));

    rd.xy = mix(rd_mirrored , rd.xy, 1.0-mirror);

    // The inner light. The shadow from this creates the Mobius shadow pattern.
    // Where you place it depends on the shadow pattern you're after.
    vec3 lp = vec3(0, 0, -.85);
    // Outer light, just to light the spherical object a bit. Without this
    // some things would be dark.
    vec3 lp2 = vec3(0, .5, -1.5);

    // Rotate the camera slightly.
    rd.yz *= rot2(sin(smoothTime*.1) * .05);

    // Raymarch.
    float t = trace(ro, rd);

    // Object ID. 
    float oID = gOID;

    // Saving the spherical texture coordinates, and 2D dashed line field.
    vec2 sUV = gSUV;
    float ln = gLn;

    // Scene color. Initialized to zero.
    vec3 col = vec3(0);

    if(t < FAR) {

        // Surface hit point and normal.
        vec3 p = ro + rd * t;
        vec3 n = norm(ro);

        // Light distances.
        float lDst = length(lp - p);
        float lDst2 = length(lp2 - p);

        // Unit direction lights -- Inner and outer.
        vec3 ld = (lp - p) / lDst;
        vec3 ld2 = (lp2 - p) / lDst2;

        float dif = max(dot(ld, n), 0.);
        float dif2 = max(dot(lp2, n), 0.);

        // Object color.
        vec3 oCol;

        #ifndef DASHED
        // If not dashed, tone the glow down a little.
        glow *= .7;
        #endif

        if(oID > 0.) {

            // Sphere object. 

            if(oID == 1.) { 

                // Inner face.

                // Inside and outside sphere colors.
                if(length(p) < 1.) {
                    oCol = glowCol2 + glow / 2.;
                } else
                    oCol = glowCol1 + glow / 2.;

                // Dotted LED pattern on the face.
                const float sc = 1. / 16.;
                vec2 q = sUV + sc / 2.;
                if(mod(floor(q.y / sc), 2.) < .5)
                    q.x += .5 * sc;
                vec2 iq = floor(q / sc) + .5;
                q -= iq * sc;

                // Slight random colorization to each cell.
                oCol = mix(oCol, oCol * vec3(2, .5, .5), hash21(iq + .2));

                // Cell edging and application.
                float d = -(smax(abs(q.x), abs(q.y), .01) - .4 * sc);
                oCol = mix(oCol, oCol * .05, 1. - smoothstep(0., 2. / RENDERSIZE.y, d)); 

                // Dark edges on the frame using the saved 2D rectangle distance
                // from the distance function.
                oCol = mix(oCol, oCol * .1, 1. - smoothstep(0., .001, abs(ln + .01) - .0025));

            } else {

                // Outer frame.

                // The purple outer frames with edging.
                oCol = vec3(3, 1, .2).yzx / 2.5 + glow;
                oCol = mix(oCol, oCol * .1, 1. - smoothstep(0., .003, abs(ln - .0075) - .0025));
                oCol = mix(oCol, oCol * .1, 1. - smoothstep(0., .003, abs(ln + .01) - .0025));

            }

            // Texturing. Not necessary here.
            //vec3 tx = tex3D(iChannel1, vec3(sUV + .5, length(p)), n);
            //oCol *= tx*2. + .5;

        } else {

            // Backgournd plane.

            // Dark purple.
            oCol = vec3(1, .2, 3) / 3.5;

            // Apply a repeat hexagonal dot pattern.
            vec2 sc = vec2(1. / .8660254, 1) * 1. / 7.;
            vec2 q = p.xy;
            if(mod(floor(q.y / sc.y), 2.) < .5)
                q.x -= .5 * sc.x;
            vec2 iq = floor(q / sc) + .5;
            q -= iq * sc;

            // Subtle random coloring from cell to cell.
            oCol = mix(oCol, vec3(1, .3, .05) / 3.5, hash21(iq + .06) * .35);

            // The dot background edges and application.
            float d = length(q) - .3 * sc.y*background_on_off;
            oCol = mix(oCol * .2, oCol, 1. - smoothstep(0., 1. / sc.y / RENDERSIZE.y, d))*background_on_off;
            //oCol = mix(oCol, oCol*.5, 1. - smoothstep(0., 1./sc.y/RENDERSIZE.y, d + .15*sc.y));

            // Background texture. Not used.
            //vec3 tx = texture(iChannel1, p.xy/8.).xyz; tx *= tx;
            //oCol *= tx*3. + .2;

        }

        // Shadows and ambient self shadowing.
        float sh = softShadow(p, lp, n, 128.);
        float ao = calcAO(p, n);

        // Inner and outer light attenuation.
        float att = 1. / (1. + lDst * lDst * .05);
        float att2 = 1. / (1. + lDst2 * lDst2 * .05);

        // Specular lighting.
        float spe = pow(max(dot(reflect(ld, n), rd), 0.), 16.);
        float spe2 = pow(max(dot(reflect(ld2, n), rd), 0.), 16.);

        // Specular reflection.
        vec3 hv = normalize(-rd + ld2); // Half vector.
        vec3 ref = reflect(rd, n); // Surface reflection.
        //vec3 refTx = texture(image8, ref).xyz; 
        //refTx *= refTx;
        //refTx = (col*1.5 + .66)*refTx;//smoothstep(.2, .5, refTx);
        vec3 refTx = (col * 1.5 + .66);//smoothstep(.2, .5, refTx);

        float spRef = pow(max(dot(hv, n), 0.), 8.); // Specular reflection.
        float rf = oID == 0. ? oCol.z * .5 : oID == 1. ? oCol.x : oCol.z;
        float highs = pow(syn_HighLevel * 0.5 + 0.5 * syn_MidHighLevel, 2.0);
        // Adding the specular reflection and glow for the inner light.
        oCol += spRef * refTx * rf * (.25 * (1.0 + highs));
        oCol += glow / (1.5 - highs);         
        oCol *= 1.0 + highs;
        // Combining all terms for the inner light.
        col = oCol * (dif * sh + vec3(1, .5, .8) * spe * sh + .2) * att * ao;

        // Placing another light on the outside of the sphere. Technically, there should
        // be another shadow calculation, but it'd doesn't have much effect, so we're
        // saving the extra calculations.
        oCol += spRef * refTx * rf * .25; // Add more reflection to the outside.
        col += oCol * ((dif2 * .15 + .05 + glow * 2.) + vec3(1, .5, .8) * spe2 * 2.) * att2 * ao;

    }

    fragColor = vec4(sqrt(max(col, 0.)), 1.);
    return fragColor;
}

vec4 renderMain() {
    if(PASSINDEX == 0) {
        return renderMainImage();
    }
}