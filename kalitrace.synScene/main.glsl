

/** KaliTrace
	https://www.shadertoy.com/view/4sKXWG

	(cc) 2016, stefan berke

	Another attempt to ray-march the kali-set. 
	Quite happy this time, not so much artefacts and good framerate. 
	
	It's quite inspired by Shane's general visual style,
	especially bump-mapping and evironment-mapping.
	Both drawn from the kaliset as well.
	

*/


// minimum distance to axis-aligned planes in kali-space
// uses eiffie's mod (/p.w) https://www.shadertoy.com/view/XtlGRj
// to keep the result close to a true distance function
vec3 kali_set(in vec3 pos, in vec3 param)
{
    vec4 p = vec4(pos, 1.);
    vec3 d = vec3(100.);
    for (int i=0; i<9; ++i)
    {
        p = abs(p) / dot(p.xyz,p.xyz);
        d = min(d, p.xyz/p.w);
        p.xyz -= param;
    }
    return d;
}


float DE(in vec3 p, in vec3 param)
{
    // floor and ceiling
    float d = min(p.y, -p.y+.2);

    // displaced by kaliset
	d -= kali_set(p*vec3(1,2,1), param).x;
    
    return d;
}

vec3 DE_norm(in vec3 p, in vec3 param)
{
    vec2 e = vec2(0.0001, 0.);
    return normalize(vec3(
        DE(p+e.xyy, param) - DE(p-e.xyy, param),
        DE(p+e.yxy, param) - DE(p-e.yxy, param),
        DE(p+e.yyx, param) - DE(p-e.yyx, param)));
}


const float max_t = 1.;

// common sphere tracing
// note the check against abs(d) to get closer to surface
// in case of overstepping
float trace(in vec3 ro, in vec3 rd, in vec3 param)
{
    float t = 0.001, d = max_t;
    for (int i=0; i<50; ++i)
    {
        vec3 p = ro + rd * t;
        d = DE(p, param);
        if (abs(d) <= 0.0001 || t >= max_t)
            break;
        t += d * .5; // above kali-distance still needs a lot of fudging
    }
    return t;
}

// "Enhanced Sphere Tracing"
// Benjamin Keinert(1) Henry SchÃ¤fer(1) Johann KorndÃ¶rfer Urs Ganse(2) Marc Stamminger(1)
// 1 University of Erlangen-Nuremberg, 2 University of Helsinki
// 
// It was a try... disabled by default (see rayColor() below)
// Obviously the algorithm does not like "fudging" which is needed for my distance field..
// It renders more stuff close to edges but creates a lot of artifacts elsewhere
float trace_enhanced(in vec3 ro, in vec3 rd, in vec3 param)
{
    float omega = 1.2; // overstepping
    float t = 0.001;
    float candidate_error = 100000.;
    float candidate_t = t;
    float previousRadius = 0.;
    float stepLength = .0;
    float signedRadius;
    float pixelRadius = .012;
    float fudge = 0.6;
    for (int i = 0; i < 50; ++i) 
    {
        signedRadius = DE(rd*t + ro, param);
        float radius = abs(signedRadius);
        bool sorFail = omega > 1. && (radius + previousRadius) < stepLength;
        if (sorFail) 
        {
        	stepLength -= omega * stepLength;
        	omega = 1.;
        } 
        else 
        {
        	stepLength = signedRadius * omega;
        }
        previousRadius = radius;
        float error = radius / t;
        if (!sorFail && error < candidate_error) 
        {
        	candidate_t = t;
        	candidate_error = error;
    	}
    	if (!sorFail && error < pixelRadius || t > max_t)
    		break;
    	t += stepLength * fudge;
    }
    return (t > max_t || candidate_error > pixelRadius)
        ? max_t : candidate_t;
}

// common ambient occlusion
float traceAO(in vec3 ro, in vec3 rd, in vec3 param)
{
    float a = 0., t = 0.01;
    for (int i=0; i<5; ++i)
    {
        float d = DE(ro+t*rd, param);
        a += d / t;
        t += d;
    }
    return min(1., a / 5.);
}

// environment map, also drawn from kaliset
vec3 skyColor(in vec3 rd)
{
    //vec3 par = vec3(0.075, 0.565, .03);
    vec3 par = vec3(1.2, 1.01, .71);
    
    vec3 c = kali_set(rd*2., par);
    c = vec3(.9*c.x,.7,1.)*pow(vec3(c.x),vec3(.7,.5,.5));
    
    return clamp(c, 0., 1.);
}

// trace and color
vec3 rayColor(in vec3 ro, in vec3 rd)
{
    // magic params for kali-set
    vec3 par1 = vec3(1.),				// scene geometry 
         par2 = vec3(.63, .55, .73);	// normal/bump map
    
#if 1
    float t = trace(ro, rd, par1);
#else    
    float t = trace_enhanced(ro, rd, par1);
#endif    
    vec3 p = ro + t * rd;
    float d = DE(p, par1);
    
    vec3 col = vec3(0.);

    // did ray hit?
    if (d < 0.03) 
    // note, we always find a surface in this scene except for rays parallel to the 
	// two enclosing planes. The 0.03 is quite large, just to remove the blackness
    // close to edges
    {
        // surface normal
        vec3 n = DE_norm(p, par1);
        // normal displacement
        n = normalize(n + min(p.y+0.05,.14)*DE_norm(p+.1*n, par2));
        n = normalize(n + 0.04*DE_norm(sin(p*30.+n*10.), par2)); // micro-bumps
        // reflected ray
        vec3 rrd = reflect(rd,n);
		// normal towards light
        vec3 ln = normalize(vec3(0.7,0.2,0) - p);
		// 1. - occlusion
        float ao = traceAO(p, n, par1);
        // surface color
        vec3 surf = .1*mix(vec3(1,1.4,1), vec3(3,3,3), ao);

        // lighting
        surf += .25 * ao * max(0., dot(n, ln));
        float d = max(0., dot(rrd, ln));
        surf += ao * (.5 * d + .7 * pow(d, 8.));

        // environment map
        surf += .5 * ao * skyColor(rrd);
    
        // distance fog
    	col = surf * (1.-t / max_t);
    }
    
    return col;
}

// camera path
vec3 path(in float ti)
{
    ti /= 7.;
    vec3 p = vec3(sin(ti)*.5+.5, 
                  .05+.03*sin(ti*3.16), 
                  -.5*cos(ti));
    return p;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 suv = fragCoord.xy / RENDERSIZE.xy;
	vec2 uv = (fragCoord.xy - RENDERSIZE.xy*.5) / RENDERSIZE.y * 2.;
        
    float ti = smoothTime*0.25+23.;
    
    vec3 ro = path(ti);
    vec3 look = path(ti+1.+.5*sin(ti/2.3))+vec3(0,.02+.03*sin(ti/5.3),0);
    float turn = sin(ti/6.1); 
        
    // lazily copied from Shane
    // (except the hacky turn param)
    float FOV = .7; // FOV - Field of view.
    vec3 fwd = normalize(look-ro);
    vec3 rgt = normalize(vec3(fwd.z, turn, -fwd.x));
    vec3 up = cross(fwd, rgt);
    
    vec3 rd = normalize(fwd + FOV*(uv.x*rgt + uv.y*up));
    
    
    
    //vec3 col = kali_set(vec3(uv, 0.), vec3(1.));
    vec3 col = rayColor(ro, rd);
    //col = skyColor(rd);
    
    
    col *= pow(1.-dot(suv-.5,suv-.5)/.5, .6);
    
	fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}