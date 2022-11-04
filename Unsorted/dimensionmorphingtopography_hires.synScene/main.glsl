


////////////////////////////////////////////////////////////
// DimensionMorphingTopography  by mojovideotech
//
// based on :
// shadertoy.com/\lsBBD1  by bal-khan
//
// Creative Commons Attribution-NonCommercial-ShareAlike 3.0
////////////////////////////////////////////////////////////


vec3	lux;
float tc = ((smoothTimeB)*0.25+cycle);
float distanceToL,ii,m;

float sdTorus( vec3 p, vec2 tx ) {
	vec2 q = vec2(length(p.zy)-tx.x,p.x);
	return length(q)-tx.y;
}

float sdHexPrism( vec3 p, vec2 hx ) {
    vec3 q = abs(p);
    return max(q.z-hx.y,max((q.x*0.866025+q.y*0.5),q.y)-hx.x);
}

float sdBox( vec3 p, vec3 b ) { return length(max(abs(p)-b,0.0)); }

void rotate(inout vec2 v, float angle) { v = vec2(cos(angle)*v.x+sin(angle)*v.y,-sin(angle)*v.x+cos(angle)*v.y); }

float	scene(vec3 p) {
    distanceToL = 1e3;
    mat2 ma;
    float r2 = 1e5, k = 1.0;
    ii=0.0;
    m = r2;
    //float aa = ((TIME+syn_Time*0.3)*0.5)*rate*0.025;
   // float aa = rate*(-0.5)+(script_time*0.06+script_bass_time*0.015);
    float aa = ((smoothTime)+rate);
    aa *= 0.025;
    p.z+=5.0;
    float tr = (smoothTimeC*0.1)+roto;
    rotate(p.zx, rXY.x+cos(tr)/4.);
    rotate(p.zy, rXY.y+sin(tr)/4.);
    rotate(p.xy, rZ-sin(tr)/4.+tr/4.);
    for(float	i = -1.0; i < 36.0; ++i) {
        ++ii;
        if (i > floor(subdivisions*3.0)) { break; }
		r2= min(r2, sdHexPrism(p, vec2(0.3,0.3)) );
		distanceToL = sdHexPrism(p, vec2(0.3, 0.0))*60.0;
		aa=aa+0.5/(i+2.0);
        if (mod(i, 3.0) == 0.0) {
            ma = mat2(cos(aa+1.0*ii*0.25),sin(aa+1.0*ii*0.25), -sin(aa+1.0*ii*0.25), cos(aa+1.0*ii*0.25) );
	        p.xy*=ma;
	        p.xy = abs(p.xy)-0.125;
			p.z -= 0.2;
        }
        else if (mod(i, 3.0) == 1.0) {
            ma = mat2(cos(aa*3.0+1.04+1.0*ii*0.1),sin(aa*3.0+1.04+1.0*ii*0.1), -sin(aa*3.0+1.04+1.0*ii*0.1), cos(aa*3.0+1.04+1.0*ii*0.1) );
	        p.yz*=ma;
	        p.zy = abs(p.zy)-0.125;
            p.x -= 0.2;
        }
        else if (mod(i, 3.0) == 2.0) {
            ma = mat2(cos(aa*2.0+2.08+1.0*ii*0.5),sin(aa*2.0+2.08+1.0*ii*0.5), -sin(aa*2.0+2.08+1.0*ii*0.5), cos(aa*2.0+2.08+1.0*ii*0.5) );
            p.zx*=ma;
	        p.xz = abs(p.xz)-0.125-(expand-(subdivisions*0.0125));
	       	p.y -= 0.2;
        }
	m = min(m, log(sdBox(p,vec3(0.0612510))/(k*k) ) );
    k *= 1.125;
    }
    return r2;
}

vec3 evaluateLight(in vec3 pos) {
    vec3 lightCol = vec3(blend.x,length(blend),blend.y)+sin(tc);
    return vec3(lightCol * 1.0/(distanceToL*distanceToL)/(31.0-(light+(12.*highhits))))*mix(lightCol+pos+highhits,lightCol-pos-highhits,sin(tc))*1.0;
}

vec2	march(vec3 pos, vec3 dir) {
    vec3 p ;
    vec2 s , dist ;
    for (int i = 0; i < 64; ++i) {
    	p = pos + dir * dist.y;
        dist.x = scene(p);
        dist.y += dist.x;
        lux += evaluateLight(p);
        if (dist.x < 0.001 || dist.y > 30.50) { break; }
        s.x++;
    }
    s.y = dist.y;
    return (s);
}

vec3	camera(vec2 uv) {
    float   fov = (zoom)/floor(subdivisions*1.50);
	vec3    forw  = vec3(0.0, 0.0, -1.0);
	vec3    right = vec3(1.0, 0.0, 0.0);
	vec3    up    = vec3(0.0, 1.0, 0.0);
    return (normalize((uv.x) * right + (uv.y) * up + fov * forw));
}

vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);

    vec4 o = vec4(0.0,0.0,0.0,1.0);
    vec2 uv  = vec2(_xy.xy - RENDERSIZE.xy/2.0) / RENDERSIZE.y;
	vec3 dir = camera(uv);
    vec3 pos = vec3(0.0, 0.0, 0.0);
    vec2 inter = (march(pos, dir));
    if (inter.y < 20.0) {
    	o.xyz += vec3( abs(sin(tc*1.0+ii*0.1+m+1.04)-blend.x), abs(sin(tc*1.0+ii*0.1+m+2.09+blend.x)), abs(sin(tc*1.0+ii*0.1+m+PI+blend.y)))*(1.0-inter.x*0.05);
   		o.xyz += lux;
    }
    out_FragColor = vec4(o);

return out_FragColor; 
 } 
