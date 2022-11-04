vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********

/*

    The scene itself: It's just standard single pass raymarching code. The only 
    interesting function is the "block" function that renders the extruded 
    quadtree structure.
    
    
*/

// Block scale.
#define GSCALE vec2(1./2.)

// Pylon face shape: Each one has its own appeal.
// Circles: 0, Square: 1, Rounded Square: 2.
#define SHAPE 2

// Give the pylon sides the same color as the trim.
//#define LIGHT_SIDES

// Grayscale, for that artsy look.
//#define GRAYSCALE

// Max ray distance.
#define FAR 20.

highp int FC = int(FRAMECOUNT);

// Scene object ID to separate the mesh object from the terrain.
float objID;


// Standard 2D rotation formula.
mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }


// IQ's vec2 to float hash.
float hash21(vec2 p){  return fract(sin(dot(p, vec2(27.609, 57.583)))*43758.5453); }

// IQ's texure hash -- I'm taking Flockaroo's word for it that the textureLod
// function will help avoid artifacts.
vec3 hash23T(vec2 p){  return textureLod(image30, p, 0.).xyz; }



// Tri-Planar blending function. Based on an old Nvidia tutorial by Ryan Geiss.
vec3 tex3D(sampler2D t, in vec3 p, in vec3 n){ 
    
    n = max(abs(n) - .2, .001); // max(abs(n), 0.001), etc.
    n /= dot(n, vec3(1)); 
    //n /= length(n);
    
	vec3 tx = texture(t, p.yz).xyz;
    vec3 ty = texture(t, p.zx).xyz;
    vec3 tz = texture(t, p.xy).xyz;
    
    // Textures are stored in sRGB (I think), so you have to convert them to linear space 
    // (squaring is a rough approximation) prior to working with them... or something like that. :)
    // Once the final color value is gamma corrected, you should see correct looking colors.
    return (tx*tx*n.x + ty*ty*n.y + tz*tz*n.z);
}

// Getting the video texture. I've deliberately stretched it out to fit across the screen,
// which means messing with the natural aspect ratio.
//
// By the way, it'd be nice to have a couple of naturally wider ratio videos to choose from. :)
//
vec3 getTex(sampler2D ch, vec2 p){
    
    // Strething things out so that the image fills up the window. You don't need to,
    // but this looks better. I think the original video is in the oldschool 4 to 3
    // format, whereas the canvas is along the order of 16 to 9, which we're used to.
    // If using repeat textures, you'd comment the first line out.
    //p *= vec2(RENDERSIZE.y/RENDERSIZE.x, 1);
    vec3 tx = texture(ch, fract(p)).xyz;
    return tx*tx; // Rough sRGB to linear conversion.
}

float hm(in vec2 p){ 
    return dot(getTex(image2, p/4.), vec3(.299, .587, .114));
}



// IQ's extrusion formula.
float opExtrusion(in float sdf, in float pz, in float h){
    
    vec2 w = vec2( sdf, abs(pz) - h );
  	return min(max(w.x, w.y), 0.) + length(max(w, 0.));

    
    // Slight rounding. A little nicer, but slower.
    //const float sf = .025;
    //vec2 w = vec2( sdf, abs(pz) - h) + sf;
  	//return min(max(w.x, w.y), 0.) + length(max(w, 0.)) - sf;
    
}


// IQ's signed box formula.
float sBoxS(in vec2 p, in vec2 b, in float sf){
   

  vec2 d = abs(p) - b + sf;
  return min(max(d.x, d.y), 0.) + length(max(d, 0.)) - sf;
}

// Pylon or moving object ID.
vec2 oID;

// The extruded quadtree.
vec4 blocks(vec3 q3){
    
    // Block dimension.
    vec2 oDim = GSCALE;
    vec2 dim = GSCALE; 
    
    // Distance.
    float d = 1e5;
    
    // Final entry needs to fill in the rest, so you give it a 100% chance of success.
    // I'd rather not say how long it took me to figure that out. :D
    //float rndTh[4] = float[4](.25, .333, .5, 1.01);
    float rndTh[4] = float[4](.333, .5, 1.01, 1.01);    
    
    vec2 id = vec2(0);
    float boxID = 0.;
    oID = vec2(0);
    
     
    // Three quadtree levels.
    for(int k=0; k<3; k++){
 
        // Four block neighbors.
        for(int j = 0; j<=1; j++){
            for(int i = 0; i<=1; i++){
                
                // For the 12 tap (4-tap IJ loop), we need the "-.5" figure to center things.
        		// For the 27 tap (9-tap IJ loop), take it off the expression.
        	    vec2 ip0 = q3.xy/dim + vec2(i, j) - .5;

               
                float rndIJ[4];
                rndIJ[0] = hash21(floor(ip0));
                if(rndIJ[0]>=rndTh[k]) continue;

                rndIJ[1] = hash21(floor(ip0/2.));
                if(k==1 && rndIJ[1]<rndTh[0]) continue;

                rndIJ[2] = hash21(floor(ip0/4.));
                if(k==2 && (rndIJ[1]<rndTh[1] || rndIJ[2]<rndTh[0])) continue;

                //rndIJ[3] = hash21(floor(ip0/8.)); 
                //if(k==3 && (rndIJ[1]<rndTh[2] || rndIJ[2]<rndTh[1] || rndIJ[3]<rndTh[0])) continue; 


                //if(rndIJ[0]<rndTh[k])
                {

                    //vec2 p = mod(oP, dim) - dim/2.; // Last term for 8 iterations.
                    vec2 p = q3.xy;
                    vec2 ip = floor(ip0) + .5;
                    p -= ip*dim; // Last term for 8 iterations.

                    vec2 idi = ip*dim;

                    // The extruded block height. See the height map function, above.
                    //float h = hm2(idi/256.);//15 + float(j)*.15;//
                    float h = hm(idi);//15 + float(j)*.15;//

                    //h = floor(h*15.999)/15.; // Or just, "h *= .15," for nondiscreet heights.
                    h *= .2;//.2;
                   
                    // Offset. Not used here, but it can be.
                    vec2 off = vec2(0);//(hash22(idi + .37*dim) - .5)*dim*.35;

                    // The 2D blocks.
                    #if SHAPE == 0
                    // Circle.
                    float di2D = length(p - off) - dim.x/2.;
                    #elif SHAPE == 1
                    // Rounded squares.
                    float di2D = sBoxS(p, dim/2. - .01*oDim.x, .0);
                    #else
                    // Rounded squares.
                    float di2D = sBoxS(p, dim/2. - .01*oDim.x, .3*(dim.x));
                    #endif
                    
                    // Some random number.
                    vec3 r3 = hash23T(idi*7.27183 + .36);
                    float rnd = r3.x;
                    if(rnd<.5) {
                        // Hollow out random pylons.
                        di2D = abs(di2D + .04*oDim.x) - .04*oDim.x;
                    
                    }


                    // The extruded distance function value.
                    float di = opExtrusion(di2D, (q3.z + h - 32.125), h + 32.125);
                    
                    // Adding a bit of the face value to put a bit of variation on the top.
                    di += di2D/8.;
                    

                    // The moving ball bearings. I made this up quickly, so there'd be
                    // better ways to go about it.
                    vec2 toID = vec2(0);
                    // Only do it for hollowed out pylons for the first two levels.
                    if(k<2 && rnd<.5) {
                    

                        float bh = sin(TIME*(r3.y + 5.)/4. + dot(idi, vec2(27.163, 113.457)) + r3.z*6.2831)*(dim.x + 1.);
                        float qz = q3.z + h - 1. + bh;
                        float ball = length(vec3(p, qz)) - .3*dim.x*(r3.y*.5 + .5);
                        //float ball = max(oDi2D + oDim.x*.2, abs(qz) - dim.x*.2);
                        //float ball2 = max(length(p) - .3*dim.x*(rnd2*.5 + .5), abs(qz) - .0*dim.x);
                        //ball = max(ball, -(abs(ball2) - .02*dim.x));
                        
                        toID.x = di<ball? 0. : 1.;
                        toID.y = bh;
                        di = min(di, ball);
                    
                    }
                    
                    // Random sprinkles.
                    //di += hash23T((q3.xy + q3.z)/8.).z*.003;
                    


                    // If applicable, update the overall minimum distance value,
                    // ID, and box ID. 
                    if(di<d){
                        
                        d = di;
                        id = idi;
                        // Not used in this example, so we're saving the calulation.
                        boxID = di2D;//float(k*4 + j*2 + (1 - i));
                        oID = toID;
                    }


                }
                


            }
        }

        // Halve the dimension for the next iteration.
        dim /= 2.; 
        
    }  
    
    return vec4(d, id, boxID);
    
    
}



// Block ID -- It's a bit lazy putting it here, but it works. :)
vec3 gID;

// The extruded image.
float map(vec3 p){
    
    // Floor.
    float fl = -p.z + 64.251;

    // The extruded blocks.
    vec4 d4 = blocks(p);
    gID = d4.yzw; // Individual block ID.
    
 
    // Overall object ID.
    objID = fl<d4.x? 1. : 0.;
    
    // Combining the floor with the extruded image
    return  min(fl, d4.x);
 
}

 
// Basic raymarcher.
float trace(in vec3 ro, in vec3 rd){

    // Overall ray distance and scene distance.
    float t = 0., d;
    
    for(int i = min((FC), 0); i<96; i++){
    
        d = map(ro + rd*t);
        // Note the "t*b + a" addition. Basically, we're putting less emphasis on accuracy, as
        // "t" increases. It's a cheap trick that works in most situations... Not all, though.
        if(abs(d)<.001 || t>FAR) break; // Alternative: 0.001*max(t*.25, 1.), etc.
        
        //t += i<32? d*.5 : d*.85; 
        t += d*.7; 
    }

    return min(t, FAR);
}

/*
// Standard normal function. It's not as fast as the tetrahedral calculation, but more symmetrical.
vec3 getNormal(in vec3 p, float t) {
	const vec2 e = vec2(.001, 0);
	return normalize(vec3(map(p + e.xyy) - map(p - e.xyy), map(p + e.yxy) - map(p - e.yxy),	map(p + e.yyx) - map(p - e.yyx)));
}
*/

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
    for(int i = min(FC, 0); i<6; i++){
		mp[i] = map(p + sgn*e6[i/2]);
        sgn = -sgn;
        if(sgn>2.) break; // Fake conditional break;
    }
    
    return normalize(vec3(mp[0] - mp[1], mp[2] - mp[3], mp[4] - mp[5]));
}

// Cheap shadows are hard. In fact, I'd almost say, shadowing particular scenes with limited 
// iterations is impossible... However, I'd be very grateful if someone could prove me wrong. :)
float softShadow(vec3 ro, vec3 lp, vec3 n, float k){

    // More would be nicer. More is always nicer, but not really affordable... Not on my slow test machine, anyway.
    const int maxIterationsShad = 24; 
    
    ro += n*.0015;
    vec3 rd = lp - ro; // Unnormalized direction ray.
    

    float shade = 1.;
    float t = 0.;//.0015; // Coincides with the hit condition in the "trace" function.  
    float end = max(length(rd), 0.0001);
    //float stepDist = end/float(maxIterationsShad);
    rd /= end;

    // Max shadow iterations - More iterations make nicer shadows, but slow things down. Obviously, the lowest 
    // number to give a decent shadow is the best one to choose. 
    for (int i = min(FC, 0); i<maxIterationsShad; i++){

        float d = map(ro + rd*t)*.7;
        shade = min(shade, k*d/t);
        //shade = min(shade, smoothstep(0., 1., k*h/dist)); // Subtle difference. Thanks to IQ for this tidbit.
        // So many options here, and none are perfect: dist += min(h, .2), dist += clamp(h, .01, stepDist), etc.
        t += clamp(d, .01, .2); 
        
        
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

	float sca = 2., occ = 0.;
    for( int i = min(FC, 0); i<5; i++ ){
    
        float hr = float(i + 1)*.15/5.;        
        float d = map(p + n*hr);
        occ += (hr - d)*sca;
        sca *= .7;
        
        // Deliberately redundant line that may or may not stop the 
        // compiler from unrolling.
        if(sca>1e5) break;
    }
    
    return clamp(1. - occ, 0., 1.);
}


// Compact, self-contained version of IQ's 3D value noise function. I have a transparent noise
// example that explains it, if you require it.
float n3D(in vec3 p){
    
	const vec3 s = vec3(7, 157, 113);
	vec3 ip = floor(p); p -= ip; 
    vec4 h = vec4(0., s.yz, s.y + s.z) + dot(ip, s);
    p = p*p*(3. - 2.*p); //p *= p*p*(p*(p * 6. - 15.) + 10.);
    h = mix(fract(sin(h)*43758.5453), fract(sin(h + s.x)*43758.5453), p.x);
    h.xy = mix(h.xz, h.yw, p.y);
    return mix(h.x, h.y, p.z); // Range: [0, 1].
}

// Very basic pseudo environment mapping... and by that, I mean it's fake. :) However, it 
// does give the impression that the surface is reflecting the surrounds in some way.
//
// More sophisticated environment mapping:
// UI easy to integrate - XT95    
// https://www.shadertoy.com/view/ldKSDm
vec3 envMap(vec3 p){
    
    p *= 4.;
    p.y += TIME;
    
    float n3D2 = n3D(p*2.);
   
    // A bit of fBm.
    float c = n3D(p)*.57 + n3D2*.28 + n3D(p*4.)*.15;
    c = smoothstep(.45, 1., c); // Putting in some dark space.
    
    p = vec3(c*c*c, c*c, c); // Blueish tinge.
    
    return mix(p.zxy, p, n3D2); // Mixing in a bit of red.

}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


 
    // Screen coordinates.
	vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;
	
	// Camera Setup.
	vec3 ro = vec3(0., TIME/2., -1.8); // Camera position, doubling as the ray origin.
	vec3 lk = ro + vec3(0, 0, .25);//vec3(0, -.25, TIME);  // "Look At" position.
 
    // Light positioning. One is just in front of the camera, and the other is in front of that.
 	vec3 lp = ro + vec3(-.5, 4, .9);// Put it a bit in front of the camera.
	

    // Using the above to produce the unit ray-direction vector.
    float FOV = 1.; // FOV - Field of view.
    vec3 fwd = normalize(lk - ro);
    vec3 rgt = normalize(vec3(fwd.z, 0., -fwd.x )); 
    // "right" and "forward" are perpendicular, due to the dot product being zero. Therefore, I'm 
    // assuming no normalization is necessary? The only reason I ask is that lots of people do 
    // normalize, so perhaps I'm overlooking something?
    vec3 up = cross(fwd, rgt); 

    // rd - Ray direction.
    vec3 rd = normalize(uv.x*rgt + uv.y*up + fwd/FOV);
    
    // Swiveling the camera about the XY-plane.
	//rd.xy *= rot2( sin(TIME)/32. );
	rd.xy *= rot2(3.14159/5.); // + sin(TIME/6.)/8. + .125
	rd.yz *= rot2(-3.14159/5.);

	 
    
    // Raymarch to the scene.
    float t = trace(ro, rd);
    
    // Global distance.
    float gDist = t;
    
    // Save the block ID, object ID, etc.
    vec3 svGID = gID;
    vec2 svoID = oID;
    float svObjID = objID;
  
	
    // Initiate the scene color to black.
	vec3 col = vec3(0);
	
	// The ray has effectively hit the surface, so light it up.
	if(t < FAR){
        
  	
    	// Surface position and surface normal.
	    vec3 sp = ro + rd*t;
	    //vec3 sn = getNormal(sp, edge, crv, ef, t);
        vec3 sn = getNormal(sp, t);
        
        // Light direction vector.
	    vec3 ld = lp - sp;

        // Distance from respective light to the surface point.
	    float lDist = max(length(ld), .001);
    	
    	// Normalize the light direction vector.
	    ld /= lDist;

        
        
        // Shadows and ambient self shadowing.
    	float sh = softShadow(sp, lp, sn, 12.);
    	float ao = calcAO(sp, sn); // Ambient occlusion.
 	    
	    // Light attenuation, based on the distances above.
	    float atten = 1./(1. + lDist*.05);

    	
    	// Diffuse lighting.
	    float diff = max( dot(sn, ld), 0.);
        //diff = pow(diff, 4.)*2.; // Ramping up the diffuse.
    	
    	// Specular lighting.
	    float spec = pow(max(dot(reflect(ld, sn), rd ), 0.), 16.); 
	    
	    // Fresnel term. Good for giving a surface a bit of a reflective glow.
        float fre = pow(clamp(1. - abs(dot(sn, rd))*.5, 0., 1.), 2.);
        
		// Schlick approximation. I use it to tone down the specular term. It's pretty subtle,
        // so could almost be aproximated by a constant, but I prefer it. Here, it's being
        // used to give a hard clay consistency... It "kind of" works.
		float Schlick = pow( 1. - max(dot(rd, normalize(rd + ld)), 0.), 5.);
		float freS = mix(.15, 1., Schlick);  //F0 = .2 - Glass... or close enough.  
        
          
        // Obtaining the texel color. 
	    vec3 texCol;   

        // The extruded grid.
        if(svObjID<.5){
        
            // Coloring the pylons and moving objects.
            
            // Coloring the individual blocks with the saved ID.
            vec3 tx = getTex(image3, svGID.xy/8. + .12);
            
            // Texture coordinates.
            vec3 txP = sp;
           
            if(svoID.x==1.){
            
                // Moving bearing object texture coordinates.
                txP += vec3(0, 0, svoID.y);

                // Evening out the moving bearing object color. 
                tx = mix(tx, min(1./vec3(dot(tx, vec3(.299, .587, .114))), 1.), .2);
            }
            
            tx = min(tx*1.5, 1.);
            
            // Multicolored.
            //tx = .58 + .42*cos(6.2831*n3D(vec3(svGID.xy, 1))/2. + vec3(0, 1, 2));
            
            
            // Tri-planar texturing.
            vec3 tx2 = tex3D(image3, txP/2. + svoID.x/2.*0., sn);
            tx2 = smoothstep(0., .5, tx2);
            
            // Pylon color... I need to give the names more meaning, I think. :)
            texCol = mix(tx, tx*tx2*4., 1.);
            
           
             
            //diff = pow(diff, 4.)*2.;
            // Central dots.
            //vec2 svP = sp.xy - svGID.xy;
            //texCol = mix(texCol, vec3(0), 1. - smoothstep(0., .005, length(svP) - .02));
            
            // Coloring the pylons only.
            if(svoID.x==0.){
            
                // Face height and value.
                float ht = hm(svGID.xy)*.2;
                float face = svGID.z;

                // Top pylon face coloring.
                float face2 = face;
                face = max(abs(face), abs(sp.z + ht*2.)) - .02; // Face border.
                //face = min(face, abs(face2 + .01) - .00125); // Extra border.
                
                
                // Line pattern.
                float lSc = 30.;
                float pat = (abs(fract((sp.x - sp.y)*lSc - .5) - .5)*2. - .5)/lSc;
                //float pat = (abs(fract(hex2*lSc - .5) - .5)*2. - .25)/lSc;
                
                // Top color.
                vec3 tCol = texCol/4.;
               
                // Top face color... I hacked this together in a hurry, but I'll tidy it later.
                vec3 faceCol = mix(texCol, vec3(0), (1. - smoothstep(0., .01, pat))*.25);
                float ew = .03;
                tCol = mix(tCol, tx2*vec3(1.3, .9, .3)*4., (1. - smoothstep(0., .002, face2 + .01)));
                tCol = mix(tCol, texCol/4., (1. - smoothstep(0., .002, face2 + ew)));
                tCol = mix(tCol, faceCol, (1. - smoothstep(0., .002, face2 + ew + .01)));
                
                #ifdef LIGHT_SIDES
                // Colored sides.
                texCol = tx2*vec3(1.3, .9, .3)*3.;
                #endif
                
                // Applying the top face color.
                texCol = mix(texCol, tCol, (1. - smoothstep(0., .002, sp.z + ht*2.)));
 
            
                
            }
            //else texCol = (tx2*vec3(1.3, .9, .3)*4.);
           
            #ifdef GRAYSCALE
            // Grayscale, for that artistic feel.
            texCol = vec3(dot(texCol, vec3(.299, .587, .114)));    
            #endif
            
            // Darkening the edges of the open silos, or whatever they are. :)
            texCol = mix(texCol, vec3(0), smoothstep(0., 1., sp.z));
            
        
 
 
        }
        else {
            
            // The dark floor in the background. Hidden behind the pylons, but
            // you still need it.
            texCol = vec3(0);
        }
       
    	
        
        // Combining the above terms to procude the final color.
        col = texCol*(diff*sh + vec3(.04, .08, .12) + vec3(.2, .4, 1)*fre*sh*0. + vec3(1, .7, .4)*spec*freS*sh*8.);


        // Fake environment mapping.
        vec3 cTex = envMap(reflect(rd, sn));
        col += col*cTex*15.; 
        
        
        // Shading.
        col *= ao*atten;
          
	
	}
    

    // Mild temporal blur... Not really need here.
    //vec4 oCol = texelFetch(BuffA, ivec2(fragCoord), 0);
    //col = mix(oCol.xyz, max(col, 0.), 1./2.);
    
    // Rought gamma correction.
	fragColor = vec4(max(col, 0.), gDist);
	
	return fragColor; 
 } 


/*


    Realtime Extruded Quadtree
    --------------------------
    
    FMS_Cat posted a beautifully rendered, static extruded quadtree example the 
    other day that appealed to me on many levels. Since I already had a realtime  
    quadtree example tucked away in my account, it inspired me to put on some
    finishing touches and post it. Designwise, I "borrowed" two aspects from 
    FMS_Cat's scene. One was the camera angle, and the other was the hollowing 
    out of random blocks. I added some moving objects to the hollowed blocks to 
    provide a personal touch and to further illustrate the idea that this is a 
    realtime example. 
    
    The quadtree code itself was based on a 2D quadtree demonstration I posted 
    on Shadertoy a few months ago, which ironically was much harder to produce 
    than this.
	
	Coding up an extruded quadtree in realtime is still a bit of an ask. 
    Thankfully, unlike my coding skills, machines have improved considerably 
    over the past few years, so it runs reasonably efficiently... Not fantastic, 
    but not too bad on decent machines. Apologies to anyone with a slower system,
    but even with cost cutting, there's a fair bit of processing going on.

    If you were to code an extruded tri-level quadtree using nested brute force 
    neighboring column repeat methods, you'd need a GPU-burning 64 taps, which 
    isn't satisfactory. However, using a mixture of simultaneous scaling and 
    repeat neighboring methods can get it right down to just 12, which 
    conceptually is about as good as you're going to get... Having said that, 
    there are some freakishly good coders on this site, so it wouldn't shock me 
    if someone got the number down. :)
    
    FMS_Cat's scene also included water, which would have been a great addition, 
    but I had to let it go, as it would have required extra passes -- Fake 
    environmental lighting was all I could afford this time around. However, I 
    have a regular extruded block example with water that I'll post later.
	
    

	Inspired by:
    
	// I love static renderings like this.
    "Wooden Structure" - FMS_Cat
	https://www.shadertoy.com/view/sdjXWy
    
    // I based the code on this particular example.
    Sorted Overlapping Quadtree - Shane
    https://www.shadertoy.com/view/wtjfDy


*/


// Just a very basic depth of field routine -- I find a lot of it is
// common sense. Basically, you store the scene distance from the camera 
// in the fourth channel, then use it to determine how blurry you want
// your image to be at that particular distance.
//
// For instance, in this case, I want pixels that are 3.5 units away from 
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
    const float focD = 3.5, coc = 1.;
    // Linear distance from either side of the focal point.
    float l = abs(focD - texture(iCh, uv).w);
    // Using it to calculate the DOF.
    vec2 dof = clamp((l - coc)/(1.*coc), 0., 1.)/vec2(800, 450); 
    
    // Combine samples. Samples with a larger DOF value are taken further 
    // away from the original point, and as such appear blurrier.
    vec3 acc = vec3(0);

    for(int i = 0; i<25; i++){
        // Accumulate samples.
        acc += texture(iCh, uv + (vec2(i/5, i%5) - 2.)*dof).xyz;
    }

    // Return the new variably blurred value.
    return acc /= 25.;
    // Visual debug representation of DOF value.
    //return vec3(length(dof)*450./2.5);
}


// This would normally be a very quick routine that displays
// the scene and gives it a distance of field effect, but I 
// wanted to put in a little loading bar graphic just to let
// people know that some precalculation is happening in the 
// background... and to give impatient people like me a visual 
// representation of the time it's going to take. :D
//
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

     
    
    // Apply some depth of field, then present to the screen.
    vec3 col = DpthFld(BuffA, fragCoord/RENDERSIZE.xy);
    
    // Rough gamma correction.
	fragColor = vec4(sqrt(max(col, 0.)), 1);
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