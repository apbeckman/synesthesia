#include "hg_sdf.glsl"
#include "lygia/sdf.glsl"

// using the raymarching loop and lighting from shane's abstract corridor (https://www.shadertoy.com/view/MlXSWX)

// Grey scale.
float getGrey(vec3 p){ return p.x*0.299 + p.y*0.587 + p.z*0.114; }

// Non-standard vec3-to-vec3 hash function.
vec3 hash33(vec3 p){ 
    
    float n = sin(dot(p, vec3(7, 157, 113)));    
    return fract(vec3(2097152, 262144, 32768)*n); 
}

// 2x2 matrix rotation.
mat2 rot2(float a){
    float c = cos(a); float s = sin(a);
	return mat2(c, s, -s, c);
}

// Tri-Planar blending function. Based on an old Nvidia tutorial.
vec3 tex3D( sampler2D tex, in vec3 p, in vec3 n ){
  
    n = max((abs(n) - 0.2)*7., 0.001); // max(abs(n), 0.001), etc.
    n /= (n.x + n.y + n.z );  
    
	return (texture(tex, p.yz)*n.x + texture(tex, p.zx)*n.y + texture(tex, p.xy)*n.z).xyz;
}

// The triangle function that Shadertoy user Nimitz has used in various triangle noise demonstrations.
// See Xyptonjtroz - Very cool. Anyway, it's not really being used to its full potential here.
vec3 tri(in vec3 x){return abs(x-floor(x)-.5);} // Triangle function.


// The path is a 2D sinusoid that varies over time, depending upon the frequencies, and amplitudes.
vec2 path(in float z){ float s = sin(z/24.)*cos(z/12.); return vec2(s*12., 0.); }




float map(vec3 p){
    vec2 hexSize = vec2(.40, .25);
    vec3 bfSize = vec3(01.5, 01.5, 01.5);
    float bfThiccness = thiccness;
    float bfMod = boxMod;
    vec2 cSize = vec2(.5, .5);
    vec3 x = p;
    vec3 h = p;
    vec3 offset = vec3(3.0, 0.0, 0.0);
    h.z -= 5.0;
    // h.yz = _rotate(h.yz, PI*0.5);
    // h *= .650;
    h.xy = _rotate(h.xy, TIME*0.05);
	pR(h.xy,pow(sin(h.z*-.05), 2.));
    h.x = opUnion(h.x + offset.x, -h.x + offset.x);
    h.y = opUnion(h.y + offset.x, -h.y + offset.x);
    vec3 s = h;
    vec3 h1 = h;
    
	pR(h.xy, h.z*.2);
	pR(h1.xy, h1.z*-.2);
    // h.x = opUnion(h.x + offset.x*2., -h.x + offset.x*2.);
    // h.y = opUnion(h.y + offset.x*2., -h.y + offset.x*2.);
    h.xy = _rotate(h.xy, TIME*-0.2);
    h1.xy = _rotate(h1.xy, TIME*0.2);
    // h.y = pModPolar(h.yz, 4);
    // s.z += fract(TIME);
    // s.z = mod(s.z, 4.);
	// pR(s.xy, s.z*-.);
    s.z += bass_time*6.;
    pR(s.xy, s.z*0.13);
    s = opRepeat(s, vec3(0.0, 0.0, 9.0));
    // s.xz = _rotate(s.xz, TIME*01.2);
    // s.xy = _rotate(s.xy, smoothTimeC*0.2);
    // s.x += 3.0;
	// pR(h.xy, h.z*-1.2);
    h.z = _smin(h.z, -h.z, 0.02);
    h1.z = _smin(h.z, -h.z, 0.02);
    vec3 c = h;
    c.xy = abs(c.xy);
    float cMod = 9.0;
    c.z -= bass_time;
    c.z = mod(c.z, cMod);
    c.z -= cMod*0.5;
    c.xz = _rotate(_rotate(c.xz-.1, .5), smoothTimeC*0.1);
    h.z = mod(h1.z, bfMod);
    h1.z = mod(h1.z, bfMod);
    h.xz = _rotate(h.xz, 0.013*bfThiccness);
    h.yz = _rotate(h.yz, -0.013*bfThiccness);
    h1.xy *= 1.25;
    h.z -= bfMod*0.5;
    h1.z -= bfMod*0.5;
    c.yz = _rotate(c.yz, 0.5*PI);
    // h.xy = _rotate(h.xy,TIME);
    // h.z -= 1.2520;
    // x.z -= 5.0;
    float mSize = .5+ objSize*0.5;
    float result = 0.;
    // float t = dodecahedronSDF(x, 1.20);
    float cy = icosahedronSDF(c, mSize);
    float dec = dodecahedronSDF(c, mSize);
    float tet = tetrahedronSDF(c, mSize);
    float first_part = mix(cy, dec, .35*sin(TIME*0.5));
    float second_part = mix(first_part, tet, .35*cos(TIME*.5));
    // float r = fOpUnionChamfer(first_part, second_part, .05);
    
    
    
    // p.z = mod(p.z, 1.5);
    // float test = hexPrismSDF(p, vec2(0.19,1.0));
    // float hP = hexPrismSDF(h, hexSize);
    float bF = boxFrameSDF(h, bfSize, bfThiccness);
    // float bfBig = boxFrameSDF()
    float bF1 = boxFrameSDF(h1, bfSize, bfThiccness);
    // hP = fOpDifferenceRound(hP, cy, .12);
    float sP = hexPrismSDF(s, vec2(0.735, 21.5));
    result = opUnion(bF, second_part, .01);
    result = opUnion(result, bF1, .01);
    // result = fOpUnionRound(result, cy, .02);
    // float n= 1.0 + (0.25 - test);
    float n = result;
    // n = fOpUnionStairs(n, t, 01.2, 14.0);
    // return min(n, 1.0);
    return n;
 
}

// Texture bump mapping. Four tri-planar lookups, or 12 texture lookups in total.
vec3 doBumpMap( sampler2D tex, in vec3 p, in vec3 nor, float bumpfactor){
   
    const float eps = 0.001;
    float ref = getGrey(tex3D(tex,  p , nor));                 
    vec3 grad = vec3( getGrey(tex3D(tex, vec3(p.x - eps, p.y, p.z), nor)) - ref,
                      getGrey(tex3D(tex, vec3(p.x, p.y - eps, p.z), nor)) - ref,
                      getGrey(tex3D(tex, vec3(p.x, p.y, p.z - eps), nor)) - ref )/eps;
             
    grad -= nor*dot(nor, grad);          
                      
    return normalize( nor + grad*bumpfactor );
	
}

// Surface normal.
// vec3 getNormal(in vec3 p) {
	
// 	const float eps = 0.001;
// 	return normalize(vec3(
// 		map(vec3(p.x + eps, p.y, p.z)) - map(vec3(p.x - eps, p.y, p.z)),
// 		map(vec3(p.x, p.y + eps, p.z)) - map(vec3(p.x, p.y - eps, p.z)),
// 		map(vec3(p.x, p.y, p.z + eps)) - map(vec3(p.x, p.y, p.z - eps))
// 	));

// }
vec3 getNormal(vec3 p) {

    const vec2 e = vec2(0.002, 0);
    return normalize(vec3(map(p + e.xyy) - map(p - e.xyy), map(p + e.yxy) - map(p - e.yxy), map(p + e.yyx) - map(p - e.yyx)));
}
// Based on original by IQ.
float calculateAO(vec3 p, vec3 n){

    const float AO_SAMPLES = 6.0;
    float r = 0.0, w = 1.0, d;
    
    for (float i=1.0; i<AO_SAMPLES+1.1; i++){
        d = i/AO_SAMPLES;
        r += w*(d - map(p + n*d));
        w *= 0.5;
    }
    
    return 1.0-clamp(r,0.0,1.0);
}

// Cool curve function, by Shadertoy user, Nimitz.
//
// From an intuitive sense, the function returns a weighted difference between a surface 
// value and some surrounding values. Almost common sense... almost. :) If anyone 
// could provide links to some useful articles on the function, I'd be greatful.
//
// Original usage (I think?) - Cheap curvature: https://www.shadertoy.com/view/Xts3WM
// Other usage: Xyptonjtroz: https://www.shadertoy.com/view/4ts3z2
float curve(in vec3 p, in float w){

    vec2 e = vec2(-1., 1.)*w;
    
    float t1 = map(p + e.yxx), t2 = map(p + e.xxy);
    float t3 = map(p + e.xyx), t4 = map(p + e.yyy);
    
    return 0.01/(w*w) *(t1 + t2 + t3 + t4 - 4.*map(p));
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	
	// Screen coordinates.
	vec2 uv = (fragCoord - RENDERSIZE.xy*0.5)/RENDERSIZE.y;
    uv.xy = mix(uv.xy, uv.xy +_uvc, fov);
	float travel = bass_time*4.;
	// Camera Setup.
	vec3 camPos = vec3(0.0, 0.0, 0.0 + travel); // Camera position, doubling as the ray origin.
	vec3 lookAt = camPos + vec3(0.0, sin(TIME)*0.01, 0.5);  // "Look At" position.
    // Light positioning. One is a little behind the camera, and the other is further down the tunnel.
 	vec3 light_pos = camPos + vec3(0.0, 0.125, -0.125);// Put it a bit in front of the camera.
	vec3 light_pos2 = camPos + vec3(0.0, 0.0, 7.0);// Put it a bit in front of the camera.

	// Using the Z-value to perturb the XY-plane.
	// Sending the camera, "look at," and two light vectors down the tunnel. The "path" function is 
	// synchronized with the distance function.
	// lookAt.xy += path(lookAt.z);
	// camPos.xy += path(camPos.z);
	// light_pos.xy += path(light_pos.z);
	// light_pos2.xy += path(light_pos2.z);

    // Using the above to produce the unit ray-direction vector.
    float FOV = PI/3.; // FOV - Field of view.
    vec3 forward = normalize(lookAt-camPos);
    vec3 right = normalize(vec3(forward.z, 0., -forward.x )); 
    vec3 up = cross(forward, right);

    // rd - Ray direction.
    vec3 rd = normalize(forward + FOV*uv.x*right + FOV*uv.y*up);
    float perspectiveMod = _uvc.x*PI*Perspective*fov;

    rd.xz = _rotate(rd.xz, lookxy.x*PI + perspectiveMod);
    rd.yz = _rotate(rd.yz, lookxy.y*PI);
    rd.xy = _rotate(rd.xy, rotate*PI);
    
    // Swiveling the camera from left to right when turning corners.
    // rd.xy = rot2( path(lookAt.z).x/32. )*rd.xy;
		
    // Standard ray marching routine. I find that some system setups don't like anything other than
    // a "break" statement (by itself) to exit. 
	float t = 0.0, dt;
	for(int i=0; i<128; i++){
		dt = map(camPos + rd*t);
		if(dt<0.005 || t>150.){ break; } 
		t += dt*0.75;
	}
	
    // The final scene color. Initated to black.
	vec3 sceneCol = vec3(0.);
	
	// The ray has effectively hit the surface, so light it up.
	if(dt<0.005){
    	
    	// Surface position and surface normal.
	    vec3 sp = t * rd+camPos;
	    vec3 sn = getNormal(sp);
        
        // Texture scale factor.
        const float tSize0 = 1./.4; 
        const float tSize1 = 1./4.;
    	
    	// Texture-based bump mapping.
	    // if (sp.y<-(FH-0.005)) sn = doBumpMap(image6, sp*tSize1, sn, 0.025); // Floor.
	    // else sn = doBumpMap(image9, sp*tSize0, sn, 0.025); // Walls.
	    // sn = doBumpMap(image9, sp*tSize0, sn, 0.025); // Walls.
	    
	    // Ambient occlusion.
	    float ao = calculateAO(sp, sn);
    	
    	// Light direction vectors.
	    vec3 ld = light_pos-sp;
	    vec3 ld2 = light_pos2-sp;

        // Distance from respective lights to the surface point.
	    float distlpsp = max(length(ld), 0.001);
	    float distlpsp2 = max(length(ld2), 0.001);
    	
    	// Normalize the light direction vectors.
	    ld /= distlpsp;
	    ld2 /= distlpsp2;
	    
	    // Light attenuation, based on the distances above. In case it isn't obvious, this
        // is a cheap fudge to save a few extra lines. Normally, the individual light
        // attenuations would be handled separately... No one will notice, nor care. :)
	    float atten = min(2./(distlpsp) + 1./(distlpsp2), 1.);
    	
    	// Ambient light.
	    float ambience = 0.25;
    	
    	// Diffuse lighting.
	    float diff = max( dot(sn, ld), 0.0);
	    float diff2 = max( dot(sn, ld2), 0.0);
    	
    	// Specular lighting.
	    float spec = pow(max( dot( reflect(-ld, sn), -rd ), 0.0 ), 18.);
	    float spec2 = pow(max( dot( reflect(-ld2, sn), -rd ), 0.0 ), 18.);
    	
    	// Curvature.
	    float crv = clamp(curve(sp, 0.0125)*0.5 + 0.5, .0, 1.);
	    
	    // Fresnel term. Good for giving a surface a bit of a reflective glow.
        float fre = pow( clamp(dot(sn, rd) + 1., .0, 1.), 1.);
        
        // Obtaining the texel color. If the surface point is above the floor
        // height use the wall texture, otherwise use the floor texture.
        vec3 texCol;
        // if (sp.y<-(FH - 0.005)) texCol = tex3D(image6, sp*tSize1, sn); // Floor.
 	    // else texCol = tex3D(image9, sp*tSize0, sn); // Walls.
 	    // texCol = tex3D(image9, sp*tSize0, sn); // Walls.
        texCol = vec3(1.0);
        // Shadertoy doesn't appear to have anisotropic filtering turned on... although,
        // I could be wrong. Texture-bumped objects don't appear to look as crisp. Anyway, 
        // this is just a very lame, and not particularly well though out, way to sparkle 
        // up the blurry bits. It's not really that necessary.
        //vec3 aniso = (0.5 - hash33(sp))*fre*0.35;
	    //texCol = clamp(texCol + aniso, 0., 1.);
    	
    	// Darkening the crevices. Otherwise known as cheap, scientifically-incorrect shadowing.	
	    float shading =  crv*0.55 + 0.95; 
    	
    	// Combining the above terms to produce the final color. It was based more on acheiving a
        // certain aesthetic than science.
        //
        // Glow.
        // sceneCol = getGrey(texCol)*((diff + diff2)*0.75 + ambience*0.25) + (spec + spec2)*texCol*2. + fre*crv*texCol.zyx*2.;
        //
        // Other combinations:
        //
        // Shiny.
        sceneCol = texCol*((diff + diff2)*vec3(1.0, 0.95, 0.9) + ambience + fre*fre*texCol) + (spec + spec2);
        // Abstract pen and ink?
        //float c = getGrey(texCol)*((diff + diff2)*1.75 + ambience + fre*fre) + (spec + spec2)*0.75;
        //sceneCol = vec3(c*c*c, c*c, c);

	    
        // Shading.
        sceneCol *= atten*shading*ao;
        
        // Drawing the lines on the walls. Comment this out and change the first texture to
        // granite for a granite corridor effect.
        // sceneCol *= clamp(1.-abs(curve(sp, 0.0125)), .0, 1.);        
	
	}
	float background;
    // Edit: No gamma correction -- I can't remember whether it was a style choice, or whether I forgot at
    // the time, but you should always gamma correct. In this case, just think of it as rough gamma correction 
    // on a postprocessed color: sceneCol = sqrt(sceneCol*sceneCol); :D
	fragColor = vec4(clamp(sceneCol, 0., 1.), 1.0);
	fragColor = sqrt(pow(fragColor, vec4(2.0)));
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}