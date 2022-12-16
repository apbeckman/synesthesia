

// Created by Stephane Cuillerdier - Aiekick/2017 (twitter:@aiekick)
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Tuned via XShade (http://www.funparadigm.com/xshade/)

//#define MOUSE

const vec3 ld = vec3(0.,1., .5);

vec2 m;

float fullAtan(vec2 p)
{
    return step(0.0,-p.x)*3.14159 + sign(p.x) * atan(p.x, sign(p.x) * p.y);
}

float shape(vec2 p)
{
    return length(min(abs(p.xy) - 1., - 0.1)) - 0.2;
}

float df(vec3 p)
{
	float a = fullAtan(p.xz) * floor(m.x*5.) + TIME;
    vec2 rev = vec2(length(p.xz) - 3.5, p.y);
    mat2 rot = mat2(cos(a),-sin(a),sin(a),cos(a));
	return abs(shape(rev * rot)) - 0.03; 
}

vec3 nor( vec3 p, float prec )
{
    vec2 e = vec2( prec, 0. );
    vec3 n = vec3(
		df(p+e.xyy) - df(p-e.xyy),
		df(p+e.yxy) - df(p-e.yxy),
		df(p+e.yyx) - df(p-e.yyx));
    return normalize(n);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    m = _mouse.xy/RENDERSIZE.xy;
    if (m.x == 0.) m.x = .1;
	
#ifndef MOUSE
	m.x = sin(TIME * .2) * .5;
#endif
    
	fragColor = vec4(1);
	
	vec2 g = fragCoord.xy;
	vec2 si = RENDERSIZE.xy;
	vec2 uv = (g+g-si)/min(si.x, si.y);
	
	vec3 ro = vec3(5,0.,0);
	vec3 cu = vec3(0,1,0);
	vec3 co = vec3(0,0,0);
	
	float fov = .5;
	vec3 z = normalize(co - ro);
	vec3 x = normalize(cross(cu, z));
	vec3 y = normalize(cross(z, x));
	vec3 rd = normalize(z + fov * uv.x * x + fov * uv.y * y);
	
	float s = 1.;
	float d = 2.;
	float dm = 20.;
	
	for (float i=0.; i<150.; i++)
	{
		if(d*d/s > 1e6 || d>dm) break;
        s = df(ro + rd * d);
		d += s * 0.2;
	}
	
   	fragColor.rgb = vec3(0.2,0.5,0.8);
         
    if (d<dm)
	{
		vec3 p = ro + rd * d;	
		vec3 n = nor(p, 0.001);
			       
        fragColor.rgb = mix( vec3(
			df(p-n*1.7),
			df(p-n*1.9),
			df(p-n*2.1)) * (dot(n, ld)+2.)/2., 
			fragColor.rgb, 1.0-exp( -0.01*d*d ) ); 
        
	}
    
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}