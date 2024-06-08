#include "hg_sdf.glsl"
#include "lygia/sdf.glsl"
#define FAR 90.0
// using the raymarching loop and lighting from shane's abstract corridor (https://www.shadertoy.com/view/MlXSWX)
mat3 rotationMatrix(vec3 axis, float angle) {
//  from 'hex array' 
    float s = sin(angle);
    float c = cos(angle);
    float oc = 2.0 - c;
    
    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}

// Grey scale.
float getGrey(vec3 p){ return p.x*0.299 + p.y*0.587 + p.z*0.114; }

// Non-standard vec3-to-vec3 hash function.
vec3 hash33(vec3 p){ 
    
    float n = sin(dot(p, vec3(7, 157, 113)));    
    return fract(vec3(2097152, 262144, 32768)*n); 
}
vec2 _smin2(vec2 a, vec2 b, float c) {
    vec2 r;
    r.x = _smin(a.x, b.x, c);
    r.y = _smin(a.x, b.x, c);
    return r;
}
vec3 _smin3(vec3 a, vec3 b, float c) {
    vec3 r;
    r.x = _smin(a.x, b.x, c);
    r.y = _smin(a.x, b.x, c);
    r.z = _smin(a.z, b.z, c);
    return r;
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

/*
}{-|-}Raymarch Mapping{-|-}{
*/
vec3 g;

float map(vec3 p){
    float result;

    vec2 xyoffset = vec2(2.50, 2.5*(1.0))*0.;
    vec3 bfRepeat = vec3(4.0, 4.0, 4.0); 
    // coordinates
    vec3 bF;
    vec3 bF1;
    p.z -= 1.0;
    bF = p;

    pR(bF.xy, bF.z * 0.035);
    bF.xy = min(bF.xy + xyoffset.xy, -bF.xy + xyoffset.xy);
    bF.xy += xyoffset.x;
    // bF.x -= 2.0;
    // bF.x += .50;
    // bF.y += .50;
    bF.z -=5.0;
    vec3 s = bF;
    float bFMod = 4.0;
    float sMod = bFMod*(pow(1.1618, 4.0));
   
    s.z -= mod(shape_time*1.1, sMod);
    s.xy = mod(s.xy, bFMod);
    bF.xy = mod(bF.xy, bFMod);
    s.z = mod(s.z, sMod);
    bF.z = mod(bF.z, bFMod*0.5);
    s.z -= sMod*0.5;
    s.xy -= bFMod*0.5;
    bF.z -= bFMod*0.25;
    bF.xy -= bFMod*0.5;
    bF = abs(bF);
    vec3 c = bF;
    
    // bF.xy += dot(bF.x, bF.y)*0.1;
    pR(bF.xy, bF.z * (0.125 + _nsin(-smoothTimeC*0.075)*.25) );
    bF.xy = _rotate(bF.xy-0.125, smoothTimeC*0.01);
    bF = abs(bF);
    bF.xy = min(bF.xy, _rotate(bF.xy, .25));
    bF1 = bF;
    s.xz = abs(s.xz);
    s.xy = _rotate(s.xy, smoothTimeC*0.04);
    s.xz = _rotate(s.xz, smoothTimeC*0.05);
    // defining SDF params
    vec4 bfSize = vec4(.750, .750, 1.0, .05); //x, y, z, w
    // SDFs
    float shapeSize = .60;
    float shape = icosahedronSDF(s, shapeSize);
    float shape2 = tetrahedronSDF(s, shapeSize);
    float cube = boxFrameSDF(bF, bfSize.xyz, bfSize.w);
    float tubes = cylinderSDF(c.yzx, vec2(.1, 1.0));
    shape = mix(shape, shape2, .5);
    result = min(shape, cube);
    g += vec3(.52, .4, .5) * .025 / (.0125 + pow(result, 2.0)) * 0.95;

    // result = min(result, tubes);s
    return result; 
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

//curve function by Shadertoy user Nimitz.
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

	// Camera Setup.
	float travel = bass_time*4.;
	vec3 camPos = vec3(0.0, 0.0, 0.0 + travel); // Camera position, doubling as the ray origin.

	vec3 lookAt = camPos + vec3(cos(TIME*0.6)*0.0125, sin(TIME*0.4)*0.0125, 0.5 + travel);  // "Look At" position.
    // Light positioning. One is a little behind the camera, and the other is further down the tunnel.
 	vec3 light_pos = camPos + vec3(0.0, 0.125, -0.125);// Put it a bit in front of the camera.
	vec3 light_pos2 = camPos + vec3(0.0, 0.0, 7.0);// Put it a bit in front of the camera.

    
    float FOV = PI/4. + .5* length(_uvc)*perspective; // FOV - Field of view.
    vec3 forward = normalize(lookAt-camPos);
    vec3 right = normalize(vec3(forward.z, 0., -forward.x )); 
    vec3 up = cross(forward, right);

    // rd - Ray direction.
    vec3 rd = normalize(forward + FOV*uv.x*right + FOV*uv.y*up);
    float perspectiveMod = _uvc.x*PI*perspective*fov;
    rd.xz = _rotate(rd.xz, lookxy.x*PI + perspectiveMod);
    rd.yz = _rotate(rd.yz, lookxy.y*PI);
    rd.xy = _rotate(rd.xy, rotate*PI);
    
	float t = 0.0, dt;
	for(int i=0; i<128; i++){
		dt = map(camPos + rd*t);
		if(dt<0.005 || t> FAR){ break; } 
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
        const float tSize0 = 1./.6; 
        const float tSize1 = 1./4.;
    	
    	
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
	    float ambience = 0.125;
    	// Diffuse lighting.
	    float diff = max( dot(sn, ld), 0.0);
	    float diff2 = max( dot(sn, ld2), 0.0);
    	
    	// Specular lighting.
	    float spec = pow(max( dot( reflect(-ld, sn), -rd ), 0.0 ), 32.)*0.8;
	    float spec2 = pow(max( dot( reflect(-ld2, sn), -rd ), 0.0 ), 32.);
    	
    	// Curvature.
	    float crv = clamp(curve(sp, 0.0125)*0.5 + 0.5, .0, 1.);
	    
	    // Fresnel term. Good for giving a surface a bit of a reflective glow.
        float fre = pow( clamp(dot(sn, rd) + 1., .0, 1.4), 1.);
        
        // Obtaining the texel color. 
        vec3 texCol;
        texCol = vec3(1.0);
        texCol = _hueRotate(tex3D(image6, sp*tSize0, sn), -40);
    	// Darkening the crevices. Otherwise known as cheap, scientifically-incorrect shadowing.	
	    float shading =  crv*.25 + 0.75; 
    	
    	// Combining the above terms to produce the final color
        //
        // Glow.
        sceneCol = _grayscale(texCol, 1.0)*((diff + diff2)*0.75 + ambience*0.25) + (spec + spec2)*texCol*2. + fre*crv*texCol.zyx*2.;
        //
        // Other combinations:
        //
        // Shiny.
        // sceneCol = texCol*((diff + diff2)*vec3(1.0, 0.95, 0.9) + ambience + fre*fre*texCol) + (spec + spec2);
	    
        // Shading.
        sceneCol *= atten*shading*ao;
	    sceneCol += g * (.0125 * (1.0 + 3 * pow(syn_HighHits * 0.5 + syn_HighLevel * 0.5, 2.) * 0.1));	
	}
	float background;
    vec3 bg = vec3(0.0);
    sceneCol = sqrt(sceneCol * sceneCol); //gamma correction
    sceneCol = mix(sceneCol, bg, smoothstep(0., FAR-25., t));//min(bg.zyx*vec3(1.3, .6, .2)*1.5, 1.)
	fragColor = vec4(clamp(sceneCol, 0., 1.), 1.0);
	// fragColor = sqrt(pow(fragColor, vec4(1.990)));
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}