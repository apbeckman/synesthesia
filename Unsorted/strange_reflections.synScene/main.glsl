//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


//#define PI 3.14159265359
#define rot(a) mat2(cos(a+PI*vec4(0,1,0.5,0)))

// dave hoskins hash function
vec4 hash( in vec2 p ) {
	vec4 p4 = fract(vec4(p.xyxy) * vec4(.1031, .1030, .0973, .1099));
    p4 += dot(p4, p4.wzxy+19.19);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);
}

// value noise
vec4 noise( in vec2 p ) {
    p*=250.0;
    vec2 i = floor( p );
    vec2 f = fract( p );
	vec2 u = f*f*(3.0-2.0*f);
    return mix( mix( hash( i + vec2(0.0,0.0) ), 
                     hash( i + vec2(1.0,0.0) ), u.x),
                mix( hash( i + vec2(0.0,1.0) ), 
                     hash( i + vec2(1.0,1.0) ), u.x), u.y);
}

// domain warped noise
float liquid( in vec2 p ) {
    p += noise(vec2(0+smoothTime*0.00025, smoothTimeC*0.00025+(smoothTime*0.0003))-(p*=rot(0.1)*vec2(2.5, 0.5))).rg*0.01;
    p += noise(vec2(0+smoothTime*0.00065, smoothTime*0.000038+(smoothTimeC*0.0002))+(p*=rot(0.2)*2.5)).ba*0.01;
    p += noise(p*=6.5+sin(smoothTimeC*0.01)).rg*0.005+sin(smoothTimeC*0.00051);
    return noise(p*0.1+cos(smoothTimeC*0.0002)).a+basshits;
}

// used for normal calculation
float height( in vec3 p ) {
    return p.z-liquid(p.xy)*0.000825;
}

// normal from central differences
vec3 normal( in vec2 uv ) {
    const vec2 e = vec2(0.0, 0.0000125);
    vec3 p = vec3(uv, 0);
    return normalize(vec3(height(p-e.yxx)-height(p+e.yxx),
                          height(p-e.xyx)-height(p+e.xyx),
                          height(p-e.xxy)-height(p+e.xxy)));
}

// custom cubemap
vec3 cubemap( in vec3 dir ) {
    vec3 color = cos(dir*vec3(1, 9, 2)+vec3(2, 3, 1))*0.5+0.5;
    color = (color * vec3(0.8, 0.3, 0.7)) + vec3(0.2);
    color *= dir.y*0.5+0.5;
    color += exp(4.0*dir.y-4.0)*0.05;
    color = pow(color, vec3(1.0/2.2));
    return color*(1.0+highhits*Flash*0.75);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    
    vec2 uv = (fragCoord-RENDERSIZE.xy*0.5)/RENDERSIZE.x;
    vec3 dir = normalize(vec3(uv, 0.2));
    
    vec3 norm = normal(uv*0.02/Zoom);
    dir = reflect(dir, norm);
    
	dir.xz *= rot(cos(smoothTimeB*0.125));
    dir.yz *= rot(sin(smoothTimeB*0.125)*0.7);
    
    fragColor.rgb = cubemap(dir);
    
    fragColor.rgb = clamp(fragColor.rgb, vec3(0), vec3(1));
    fragColor.rgb = mix(fragColor.rgb, vec3(0), dot(uv, uv)*1.0);
    
    fragColor.a = 1.0;
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}