

// Author: zackpudil
// My First attempt at using Physical Based Rendering (pbr)
//  	and Image Based Lighting (ibl).
// I feel like the look could be coerced from simpler, more efficient techniques,
// but I wanna build a memeory of the lighting to have a general go to technique.

//--------------------GEOMETRY functions-------------------
mat2 rot(float a) {
    float c = cos(a), s = sin(a);
    return mat2(c, s, -s, c);
}

float box(vec3 p, vec3 b) {
    vec3 q = abs(p) - b;
    return max(q.x, max(q.y, q.z));
}

float box(vec2 p, vec2 b) {
    vec2 q = abs(p) - b;
    return max(q.x, q.y);
}

// A plastic box with metal edges.
vec2 shape(vec3 p) {
    float a = box(p, vec3(1)); // plastic box.
    
    p = -abs(p) + vec3(1); // mirror space to save on distance evaluations
    float ce = 0.2;
	// edges.
    float b = box(p - vec3(1, 0, 0), vec3(1.0 + ce, ce, ce));
    float c = box(p - vec3(0, 0, 1), vec3(ce, ce, 1.0 + ce));
    float d = box(p - vec3(0, 1, 0), vec3(ce, 1.0 + ce, ce));
    
    vec2 s = vec2(a, 2.0); // 2.0 == material id for white plastic.
    vec2 t = vec2(min(min(b, c), d), 1.0); // 1.0 == material id for metal.
    
    return s.x < t.x ? s : t;
}

// path function is used to carve out the path the camera and light takes.
vec2 path(float z) {
    return vec2(0, 2.0*cos(0.13*z));
}

// Yet Another low iteration Kali (AbsBox) fractal from me.
vec2 de(vec3 p) {
    vec3 op = p;
    p = mod(p + 2.0, 4.0) - 2.0; // repeat space infinitly on all axises.
    vec4 q = vec4(p, 1);
    
    // very basic absbox, or "kali" fractal.
    // See this to get better ideas of it's potential: 
	// https://www.shadertoy.com/user/Kali - THE GOAT
    q.xyz -= 1.0;
    for(int i = 0; i < 4+Iters; i++) {
        q.xyz = abs(q.xyz + 1.+Open) - 1.0+Open;
        q /= clamp(dot(q.xyz, q.xyz), 0.5+0.25*sin(smoothTimeC*0.0125), 1.0+0.0125*cos(smoothTimeC*0.05));
        q.xz *= rot(0.9);
        q.xy *= rot(-0.3);
        q.yz *= rot(0.1);
        q *= 1.045;
    }
    
    vec2 s = shape(q.xyz)/vec2(q.w, 1);
    s.x = max(s.x, -box(op.xy - path(op.z), vec2(0.1)));
    
    return s;
}

// straight forward sphere tracing, got from IQ (the legend).
vec2 trace(vec3 ro, vec3 rd, float mx) {
    float t = 0.0, m = -1.0;
    for(int i = 0; i < 120; i++) {
        vec2 d = de(ro + rd*t);
        if(d.x < 0.0001 || t >= mx) break;
        t += d.x*0.85;
        m = d.y;
    }
    return vec2(t, t < mx ? m : -1.0);
}

// central difference based normal, got from again IQ.
vec3 normal(vec3 p) {
    vec2 h = vec2(0.00001, 0.0);
    vec3 n = vec3(
        de(p + h.xyy).x - de(p - h.xyy).x,
        de(p + h.yxy).x - de(p - h.yxy).x,
        de(p + h.yyx).x - de(p - h.yyx).x);
    
    return normalize(n);
}

// pbr == Physical Based Rendering.  This better mimics how real world light
// interacts with surfaces. See full explanation of algorithm here:
// https://learnopengl.com/PBR/Lighting
vec3 pbr(vec3 p, vec3 n, vec3 l, vec3 rd, vec3 a, float r, float m,
         inout vec3 f0, inout float hov, inout float nov) {
    
    // set up view and halfway vector.
    vec3 v = normalize(-rd);
    vec3 h = normalize(l + v);
    
    // get angles needed for calculation.
    float nol = clamp(dot(n, l), 0.0, 1.0);
    float noh = clamp(dot(n, h), 0.0, 1.0);
    nov = clamp(dot(n, v), 0.0, 1.0);
    hov = clamp(dot(h, v), 0.0, 1.0);
    
    // F == The fresnel term (pronounced Freh-nel).  
    // This is the Schlick approximation.
    f0 = mix(vec3(0.04), a, m);
    vec3 F = f0 + (1.0 - f0)*pow(1.0 - hov, 5.0);
    
    // The "Normal Distrubution Function".
    // approximates the "microfacets" of the surface.
    float a2 = pow(r, 4.0);
    float D = a2/(3.141*pow(noh*noh*(a2 - 1.0) + 1.0, 2.0));
    
    // The "Geometry Function" approximates self shadowing of microfaces on material.
    // This is the "Smith" approximation.
    float k = 0.5*pow(0.5*r + 0.5, 2.0);
    float kl = nol*(1.0 - k) + k;
    float kv = nov*(1.0 - k) + k;
    float V = 1.0/(4.0*kl*kv);
    
    // Normal Distrubution Function * Geometry Function * Fresnel term = specular term.
    vec3 spe = F*D*V; // spe aka BRDF (bidirectional reflective distribution function)
    
    // diffuse (metal meterials don't have a diffuse term).
    vec3 dif = (1.0 - F)*(1.0 - m);

    return (a*dif/3.141 + spe)*nol; // all times the normal dot light
}

// Image Based Lighting.  This models the indirect light coming from iChannel0
vec3 ibl(vec3 p, vec3 n, vec3 r, vec3 a, float ro, float m,
         vec3 f0, float hov, float nov) {
    
    // Fresnel term modification based on Sebastien Lagarde.  Includes roughness to fix edge problem.
    vec3 F = f0 + (max(vec3(1.0 - ro), f0) - f0)*pow(1.0 - hov, 5.0);
    //vec3 irr = texture(iChannel0, n, 10.0).rgb; // irradiance from environment.
    //vec3 dif = (1.0 - F)*(1.0 - m)*irr*a; // indirect diffuse lighting.
    vec3 dif = (1.0 - F)*(1.0 - m)*a; // indirect diffuse lighting.
    
    // prefiltered specular color.
    //vec3 pc = textureLod(iChannel0, r, ro*15.0).rgb;
    
    // === this code is an approximation to the environment BRDF function describe in learnopengl.com.
    // find explanation here: https://www.unrealengine.com/en-US/blog/physically-based-shading-on-mobile
    vec4 c0 = vec4(-1, -0.0275, -0.572, 0.022);
    vec4 c1 = vec4(1, 0.0425, 1.04, -0.04);
    vec4 cr = ro*c0 + c1;
    float a004 = min(cr.x*cr.x, exp2(-9.28*nov))*cr.x + cr.y;
    vec2 ab = vec2(-1.04, 1.04)*a004 + cr.zw;
    
    // indirect brdf is prefilteredColor*(Fresnel*brdf.x + brdf.y)
    //vec3 spe = pc*(F*ab.x + ab.y);
    vec3 spe = (F*ab.x + ab.y);
    
    // very quick and dirty ambient occlusion.
    // Got this from evvvil_ (https://www.shadertoy.com/user/evvvvil) What up Broski.
    float occ = exp2(-pow(max(0.0, 1.0 - de(p + n*0.05).x/0.05), 2.0));
    
    // full indirect lighting based on environment (or image from cube map in Channel0).
    return (dif + spe)*occ;
}

// Triplaniar blending of 2d materials to 3d surfaces.
// Got this from Shane (THE GOAT).
vec3 tex3D(vec3 p, vec3 n, sampler2D s) {
    vec3 m = pow(abs(n), vec3(10.0));
    m /= dot(m, vec3(1));
    
    vec3 x = texture(s, p.yz).rgb;
    vec3 y = texture(s, p.xz).rgb;
    vec3 z = texture(s, p.xy).rgb;
    
    return m.x*x*x + m.y*y*y + m.z*z*z;
}

// Normal bumping based on material. Also got this from Shane.
vec3 bump(vec3 p, vec3 n, float bf, float f) {
    p *= f;
    vec2 h = vec2(0.001, 0.0);
    vec3 g = mat3(
        tex3D(p + h.xyy, n, image3) - tex3D(p - h.xyy, n, image3),
        tex3D(p + h.yxy, n, image3) - tex3D(p - h.yxy, n, image3),
        tex3D(p + h.yyx, n, image3) - tex3D(p - h.yyx, n, image3))
        *vec3(0.299, 0.584, 0.114); // grey scale.
    
    g -= n*dot(n, g);
    
    return normalize(n + bf*g);
}

// calculate the material for each material id in the scene.
// very simple since we have 1 white plastic, 1 white metal.
void material(float mid, vec3 p, inout vec3 n,
              inout vec3 a, inout float r, inout float m) {
    if(mid == 1.0) { // metal material.
        a = vec3(0.7, 1, 1); // white.
        n = bump(p, n, 03.4, 03.4); // minimal bumping.
        r = 0.13; // low roughness.
        m = 1.0; // metal.
    } else if(mid == 2.0) { // plastic material.
        a = vec3(0.7, 1, 1); // all white.
        n = bump(p, n, 6.0, 6.0); // lots of bumping. 
        r = 0.1; // as smooth as you can get.
        m = 0.0; // not metal.
    }   
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    uv.xy = _rotate(uv.xy, Rotate*PI);

    vec2 mo = (3.0*_mouse.xy - 1.5*RENDERSIZE.xy)/RENDERSIZE.y;
    if(_mouse.z <= 0.0) mo = vec2(0);
    // camera animation.
    float at = bass_time;
    vec3 ro = vec3(path(at), at);
    vec3 la = vec3(path(at + 1.0) + mo, at + 1.0);
    la.xy += _uvc*PI*FOV;
    
    // ray direction setup.
    vec3 ww = normalize(la-ro);
    vec3 uu = normalize(cross(vec3(0, 1, 0), ww));
    vec3 vv = normalize(cross(ww, uu));
    vec3 rd = normalize(mat3(uu, vv, ww)*vec3(uv, 01.));
    rd.xz = _rotate(rd.xz, (LookXY.x*PI)+PI*_uvc.x*Perspective*FOV);
    rd.yz = _rotate(rd.yz, (LookXY.y)*PI);
   // rd.xy = _rotate(rd.xy, Rotate*PI);

    vec3 col = texture(image3, uv).rgb;
    
    vec3 lp = la; // light position is the same as lookat vector.
    lp.xy -= _uvc*PI*FOV;
    
    vec2 t = trace(ro, rd, 50.0);
    if(t.y > 0.0) {
        // geometry.
        vec3 p = ro + rd*t.x;
        vec3 n = normal(p);
        vec3 l = normalize(lp - p); // light direction.

        // initialize material props.
        vec3 a = vec3(1);
        float ro = 0.1, m = 0.0;
        
        material(t.y, p, n, a, ro, m); // assign material props.
        
        // some calculated variables that we can reuse for ibl.
        vec3 f0;
        float hov, nov;
        col = pbr(p, n, l, rd, a, ro, m, f0, hov, nov); // direct lighting.
        
        vec3 r = reflect(rd, n); 
        col += ibl(p, n, r, a, ro, m, f0, hov, nov); // indirect lighting.
        col *= .24 +(syn_HighLevel*0.25+0.25*syn_MidHighLevel)*(.750+syn_Intensity);
    }
    
    col = mix(col, vec3(0), 1.0 - exp(-0.1*t.x)); // some black fog.  lazy attenuation.
    col = 1.0 - exp(-0.9*col); // contrast/tone mapping.
    fragColor = vec4(pow(col, vec3(0.454545)), 1); // some gamma correction, and that's it.
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}