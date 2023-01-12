

////////////////////////////////////////////////////////////
// Equirec_MengerTunnel  by mojovideotech
//
// based on : shadertoy/lsjcWV
//
// Creative Commons Attribution-NonCommercial-ShareAlike 3.0
////////////////////////////////////////////////////////////

//Modifications by Meebs
//Additional mods by Aidan
float random (in float x) {
	return fract(sin(x) * 1e4);
}
vec4 texture2DAA(sampler2D tex, in vec3 p, vec2 uv) {
    vec2 texsize = vec2(textureSize(tex,0));
    vec2 uv_texspace = uv*texsize;
    vec2 seam = floor(uv_texspace+.5);
    uv_texspace = (uv_texspace-seam)/fwidth(uv_texspace)+seam;
    uv_texspace = clamp(uv_texspace, seam-.5, seam+.5);
    return texture(tex, uv_texspace/texsize);
}


vec4 texture2DAA1(sampler2D tex, in vec3 p, in vec3 n) {
    vec2 uv = _xy;
    vec2 texsize = vec2(textureSize(tex,0));
    vec2 uv_texspace = uv*texsize;
    vec2 seam = floor(uv_texspace+.5);
    uv_texspace = (uv_texspace-seam)/fwidth(uv_texspace)+seam;
    uv_texspace = clamp(uv_texspace, seam-.5, seam+.5);
    return texture(tex, uv_texspace/texsize);
}

// Based on Morgan McGuire @morgan3d
// https://www.shadertoy.com/view/4dS3Wd
float n3D (in vec3 p) {
	const vec3 step = vec3(110.0, 241.0, 171.0);

	vec3 i = floor(p);
	vec3 f = fract(p);

    // For performance, compute the base input to a
    // 1D random from the integer part of the
    // argument and the incremental change to the
    // 1D based on the 3D -> 1D wrapping
	float n = dot(i, step);

	vec3 u = f * f * (3.0 - 2.0 * f);

	return mix( mix(mix(random(n + dot(step, vec3(0,0,0))),
                        random(n + dot(step, vec3(1,0,0))),
                        u.x),
                    mix(random(n + dot(step, vec3(0,1,0))),
                        random(n + dot(step, vec3(1,1,0))),
                        u.x),
                u.y),
                mix(mix(random(n + dot(step, vec3(0,0,1))),
                        random(n + dot(step, vec3(1,0,1))),
                        u.x),
                    mix(random(n + dot(step, vec3(0,1,1))),
                        random(n + dot(step, vec3(1,1,1))),
                        u.x),
                u.y),
                u.z);
}

float Noise = 0.065*n3D(vec3(smoothTime*0.025));
#define 	twpi  	PI*2.0  	// two pi, 2*pi
#define 	pi   	PI 	// pi

float mid(vec3 p) { p = min(p, p.yzx); return max( max(p.x, p.y), p.z); }

vec3 tpl( sampler2D t, in vec3 p, in vec3 n ){
   
    n = max(abs(n), 0.001);
    n /= (n.x + n.y + n.z );  
    vec3 tx = (texture(t, p.yz)*n.x + texture(t, p.zx*.5)*n.y + texture(t, p.xy)*n.z).xyz;
    
    return tx*tx;
}

vec3 colGet(float t){
    return _palette(_triWave(t, 1.0), vec3(0.340, 0.380-0.12*syn_HighLevel, -0.180), vec3(0.339, 0.040+0.2*syn_MidHighLevel, 0.670), vec3(0.2555, 0.900, 0.315), vec3(0.500, 0.890, 0.000));
}

vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);

	float T = smoothTime;
	vec2 v = (_xy.xy / RENDERSIZE.xy) + RENDERSIZE.y;
	v.y -= 0.5;
	float th =  v.y * pi, ph = v.x * twpi;
    vec3 sp = vec3( sin(ph) * cos(th), sin(th), cos(ph) * cos(th) );
    // sp.xy += lookXY;
    sp = mix(sp, normalize(vec3(_uvc*PI, 1.0)), perspective);

    sp.yz = _rotate(sp.yz, lookXY.y*PI+Noise+_uvc.y*PI*perspective*perspective2);
    sp.xz = _rotate(sp.xz, lookXY.x*PI+Noise);

    sp.xy = _rotate(sp.xy, Rotate*PI+0.025*n3D(vec3(smoothTime*0.01)));
    sp.xy += sp.xy*_uvc*PI*QuadMirror*PI;

    vec3 camPos = vec3( pi, pi, bass_time*2.);
    vec3 pos = camPos;
    vec3 dir = normalize(sp);
    vec3 signdir = sign(dir);
    float stepsize = 1.0;
    float dist = 0.0;
    vec3 normal;
    for (int i = 0; i < 45; i++) {
        vec3 p = mod(pos, twpi * stepsize) - pi * stepsize;
        p.xy -= 0.25*MengerMod*p.xy;
    
        p.z = mix(p.z, 0.0, mix(thin_structure, 0.78, fully_open));
        vec3 num = (stepsize - p * signdir) * step( abs(p), vec3(stepsize) ) / dir * signdir;
        float len = mid(num);
        if (len < 0.01) {
            if (stepsize < 0.0125) break;
            stepsize /= noise;
            //stepsize += sin(smoothTimeC*0.0125)*0.001+cos(smoothTimeC*0.0125)*0.001;
        } else normal = vec3( equal( vec3(len), num) );
        pos += dir*len;
        // pos.xy = _rotate(pos.xy, i*sin(TIME)*0.000);
        dist += len;
    }

    vec3 finalCol = vec3(0.0);
    vec3 lp = camPos;
    vec3  li = normalize( lp - pos ); // Point light.
    float spe = pow(max(dot(reflect(-li, normal), -dir), 0.), 32.); // Object specular.
    float dif = clamp(dot(normal, li), 0.0, 1.0); // Diffuse.
    finalCol = vec3(0.0);

    finalCol += colGet(dot(pos, normal*10.0)*20.01)*( (.50+syn_Intensity*0.4+syn_HighLevel*0.2) / dist);

    //cyan circuits
   // finalCol.rgb += tpl(circuits2, pos*0.149, normal)*vec3(0.2,0.8,1.0)*(syn_HighLevel+syn_MidHighLevel)*sin(pos.z+0.125*smoothTimeB);
    finalCol.rgb += tpl(circuits2, pos*0.24, normal)*vec3(01.2,01.75,01.50)*(syn_HighLevel+syn_HighHits+syn_MidHighLevel)*sin(pos.z+0.125*smoothTimeB);

    finalCol.rgb += tpl(circuits, pos*0.17+vec3(0.0, 0.0, -syn_BPMTwitcher*0.05), normal)*vec3(1.40,0.764,01.27370)*pow(_fbm(pos+smoothTimeB*0.06),2.50)*(0.850+syn_Intensity+.85*syn_HighLevel+syn_MidHighHits);

    if (syn_MediaType > 0.5){
        vec3 mediaCol = tpl(syn_UserImage, pos*media_scale, normal*10.0);
        finalCol.rgb = mix(finalCol.rgb + finalCol.rgb*mediaCol*1.5, finalCol.rgb + mediaCol, media_mult_or_add);
    }
    finalCol = clamp(finalCol, 0.0, 1.0);
    
    finalCol += 0.2*pow(mix(0.0, dist, smoothstep(20.0,200.0,dist)),0.5)*(1.0+syn_Intensity+syn_HighLevel);

    return vec4(finalCol,1.0); 
} 
