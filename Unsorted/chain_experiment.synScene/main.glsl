

// All the distance functions from:https://iquilezles.org/articles/distfunctions
// Edge detection code from:https://www.shadertoy.com/view/MsSGD1
#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))
#define matRotateX(rad) mat3(1,0,0,0,cos(rad),-sin(rad),0,sin(rad),cos(rad))
#define matRotateY(rad) mat3(cos(rad),0,-sin(rad),0,1,0,sin(rad),0,cos(rad))
#define matRotateZ(rad) mat3(cos(rad),-sin(rad),0,sin(rad),cos(rad),0,0,0,1)
#define hash(h) fract(sin(h) * 43758.5453123)
#define EDGE_WIDTH 0.025

// Tweaked this based on Iq's suggestion.
float sdChain(vec3 p, float speed) {
    p.y -= smoothTimeC*speed;
    
    const float le = 0.2, r1 = 0.3, r2 = 0.1;
    
    float ya = max(abs(mod(p.y,     1.5)-0.75)-le,0.0);
    float yb = max(abs(mod(p.y+0.75,1.5)-0.75)-le,0.0);

    float la = ya*ya - 2.0*r1*sqrt(p.x*p.x+ya*ya);
    float lb = yb*yb - 2.0*r1*sqrt(p.z*p.z+yb*yb);
    
    return sqrt(dot(p.xz,p.xz) + r1*r1 + min(la,lb)) - r2;
}

vec4 combine(vec4 val1, vec4 val2 ){
    if ( val1.w < val2.w ) {
        return val1;
    }
    return val2;
}

vec4 map(vec3 p){    
    vec3 pref = p;
    vec2 uv = p.xy;
    vec4 f = vec4(vec3(1.0,0.0,0.0),p.y+1.0);
    
    float c = sdChain((p+vec3(-1.0,0.0,0.0))*matRotateX(radians(30.0*Test)),-2.0);
    float c2 = sdChain((p+vec3(1.0,0.0,0.0))*matRotateZ(radians(-20.0*Test)),-3.5);
    float c3 = sdChain(p*matRotateZ(radians(20.0))*matRotateX(radians(-30.0*Test)),-2.7);
    
    float res = min(min(c,c2),c3);
    return vec4(vec3(0.9,0.3,0.0),res);
}

vec3 normalMap(vec3 p){
    float d = 0.0001;
    return normalize(vec3(
        map(p + vec3(  d, 0.0, 0.0)).w - map(p + vec3( -d, 0.0, 0.0)).w,
        map(p + vec3(0.0,   d, 0.0)).w - map(p + vec3(0.0,  -d, 0.0)).w,
        map(p + vec3(0.0, 0.0,   d)).w - map(p + vec3(0.0, 0.0,  -d)).w
    ));
}

float shadowMap(vec3 ro, vec3 rd){
    float h = 0.0;
    float c = 0.001;
    float r = 1.0;
    float shadow = 0.5;
    for(float t = 0.0; t < 30.0; t++){
        h = map(ro + rd * c).w;
        if(h < 0.001){
            return shadow;
        }
        r = min(r, h * 16.0 / c);
        c += h;
    }
    return 1.0 - shadow + r * shadow;
}

// from simon green and others
float ambientOcclusion(vec3 p, vec3 n)
{
    const int steps = 4;
    const float delta = 0.15;

    float a = 0.0;
    float weight = 4.;
    for(int i=1; i<=steps; i++) {
        float d = (float(i) / float(steps)) * delta; 
        a += weight*(d - map(p + n*d).w);
        weight *= 0.5;
    }
    return clamp(1.0 - a, 0.0, 1.0);
}

mat3 setCamera( in vec3 ro, in vec3 ta, float cr )
{
    vec3 cw = normalize(ta-ro);
    vec3 cp = vec3(sin(cr), cos(cr),0.0);
    vec3 cu = normalize( cross(cw,cp) );
    vec3 cv =          ( cross(cu,cw) );
    return mat3( cu, cv, cw );
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 p = (fragCoord.xy * 2.0 - RENDERSIZE.xy) / min(RENDERSIZE.x, RENDERSIZE.y);
    vec2 uv = p;
    
    float time = TIME*2.0;
    
    vec3 ro = vec3( 3.5*cos(0.1*smoothTime + 6.0), 0.0, -0.5+5.5*sin(0.1*smoothTime + 6.0) );
    vec3 ta = vec3( 0.0, -0.4, -0.7 );
    mat3 ca = setCamera( ro, ta, 0.0 );
    float zoom = 1.5;
    vec3 rd = ca * normalize( vec3(p.xy,zoom) );
    
    float t, dist;
    float lastDistEval = 1e10;
	float edge = 0.0;
    t = 0.0;
    vec3 distPos = ro+rd;
    vec4 distCl = vec4(0.0);
    for(int i = 0; i < 240; i++){
        distCl = map(distPos);
        dist = distCl.w;
        t += dist;
        distPos = ro+rd*t;
        
		if (lastDistEval < EDGE_WIDTH && dist > lastDistEval + 0.001) {
			edge = 1.0;
		}
        if (dist < lastDistEval) lastDistEval = dist;
        if(dist < 0.01 || dist > 30.0) break;
    }

    vec3 color;
    float shadow = 1.0;
    if(dist < 1.0){
        // lighting
        vec3 lightDir = vec3(0.0, 1.0, 0.0);
        vec3 light = normalize(lightDir + vec3(0.5, 0.0, 0.9));
        vec3 normal = normalMap(distPos);

        // difuse color
        float diffuse = clamp(dot(light, normal), 0.5, 1.0);
        float lambert = max(.0, dot( normal, light));
        
        // ambient occlusion
        float ao = ambientOcclusion(distPos,normal);
        
        // shadow
        shadow = shadowMap(distPos + normal * 0.001, light);

        // result
        color += vec3(lambert);
        color = ao*diffuse*(distCl.xyz+(.1-length(p.xy)/3.))*vec3(1.0, 1.0, 1.0);
        
    }else{
        color =.0*max(mix(vec3(1.1,0.9,0.0)+(.1-length(p.xy)/3.),vec3(1),.1),0.);
    }

    // rendering result
    float brightness = 1.5;
    vec3 dst = (color * max(0.8, shadow))*brightness;
    
    // add edge detection result
    dst = mix(dst,vec3(0.1,0.1,0.1),edge);
    
    fragColor = vec4(dst, 1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}