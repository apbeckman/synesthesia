

////////////////////////////////////////////////////////////
// Equirec_MengerTunnel  by mojovideotech
//
// based on : shadertoy/lsjcWV
//
// Creative Commons Attribution-NonCommercial-ShareAlike 3.0
////////////////////////////////////////////////////////////


#define 	twpi  	6.283185307179586  	// two pi, 2*pi
#define 	pi   	3.141592653589793 	// pi

float mid(vec3 p) { p = min(p, p.yzx); return max( max(p.x, p.y), p.z); }

vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);

	float T = TIME * rate;
	vec2 v = (_xy.xy / RENDERSIZE.xy) + RENDERSIZE.y;
	v.y -= 0.5;
	float th =  v.y * pi, ph = v.x * twpi;
    vec3 sp = vec3( sin(ph) * cos(th), sin(th), cos(ph) * cos(th) );
    vec3 pos = vec3( pi, pi, T);
    vec3 dir = normalize(sp);
    vec3 signdir = sign(dir);
    float stepsize = 1.0;
    float dist = 0.0;
    vec3 normal;
    for (int i = 0; i < 20; i++) {
        vec3 p = mod(pos, twpi * stepsize) - pi * stepsize;
        vec3 num = (stepsize - p * signdir) * step( abs(p), vec3(stepsize) ) / dir * signdir;
        float len = mid(num);
        if (len < 0.01) {
            if (stepsize < 0.05) break;
            stepsize /= depth;
        } else normal = vec3( equal( vec3(len), num) );
        pos += dir*len;
        dist += len;
    }
	out_FragColor  = vec4((( sin(pos - T * colorCycle) * 0.5 + 0.5) + 0.25 * normal) * 1.0 / dist, 1.0);

return out_FragColor; 
 } 







