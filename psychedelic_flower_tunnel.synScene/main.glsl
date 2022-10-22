vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


//------------------------------------------------------------------------
// # refs
// ## repeat, mod, pmod
// - https://qiita.com/kaneta1992/items/21149c78159bd27e0860#ifs%E3%81%AB%E3%82%88%E3%82%8B%E8%A4%87%E9%9B%91%E3%81%AA%E5%BD%A2%E7%8A%B6
// - https://qiita.com/edo_m18/items/8c20c9c80d44e8b1dfe4#_reference-7c25132332209c2cd12e
// - https://gam0022.net/blog/2017/03/02/raymarching-fold/
// ## collor pallete
// - https://www.shadertoy.com/view/ll2GD3
// ## material mechanism
// - https://www.shadertoy.com/view/ldlcRf
// ## fake motion blur ?
// - https://www.shadertoy.com/view/WllSDM
// ## distance functions, fog ... 
// - https://iquilezles.org/www/index.htm
//------------------------------------------------------------------------

precision highp float;

#define EPS .0001
#define PI 3.1415
#define repeat(p, o) mod(p, o) - o * .5

// ref: https://www.shadertoy.com/view/ldlcRf
#define TEST_MAT_LESS (a, b) a < (b + .1)
#define TEST_MAT_GREATER (a, b) a > (b - .1)

const float stopThreshold = .0001;

struct Light {
	vec3 position;
  	float intensity;
  	vec3 color;
  	vec3 diffuse;
  	vec3 specular;
  	float attenuation;
};
    
struct Surface {
  	float depth;
  	float dist;
  	vec3 position;
  	vec3 baseColor;
  	vec3 normal;
  	vec3 emissiveColor;
  	float material;
};  

// ref: https://www.shadertoy.com/view/ldlcRf
vec2 minMat(vec2 d1, vec2 d2) {
	return (d1.x < d2.x) ? d1 : d2;
}

// ref: https://www.shadertoy.com/view/WllSDM
float n3(vec3 p) {
	vec3 r = vec3(1, 99, 999);
	vec4 s = dot(floor(p), r) + vec4(0., r.yz, r.y + r.z);
	p = smoothstep(0., 1., fract(p));
	vec4 a = mix(fract(sin(s) * 5555.), fract(sin(s + 1.) * 5555.), p.x);
	vec2 b = mix(a.xz, a.yw, p.y);
	return mix(b.x, b.y, p.z);
}
    
// ref: https://www.shadertoy.com/view/ll2GD3
// t: 0-1, a: contrast, b: brightness, c: times, d: offset
vec3 palette(float t, vec3 a, vec3 b, vec3 c, vec3 d) {
    return a + b * cos(6.28318 * (c * t + d));
}

mat2 rot(float a) {
    float c = cos(a), s = sin(a);
    return mat2(c,s,-s,c);
}

float sdBox(vec3 p, vec3 b) {
  vec3 d = abs(p) - b;
  return length(max(d,0.0))
         + min(max(d.x,max(d.y,d.z)),0.0);
}

// d1 subtract to d2
float opSub(float d1, float d2) {
	return max(-d1, d2);
}

float sdWireBox(vec3 p, float s, float b) {
    float d = 0.;
    float o = s + .01;
    float i = s - b;
    d = opSub(
        sdBox(p, vec3(o, i, i)),
        sdBox(p, vec3(s))
    );
    d = opSub(
        sdBox(p, vec3(i, o, i)),
		d
    );
    d = opSub(
        sdBox(p, vec3(i, i, o)),
		d
    );
    return d;
}

vec2 pmod(vec2 p, float r) {
	float a = atan(p.x, p.y) + PI / r;
    float n = PI * 2. / r;
    a = floor(a / n) * n;
    return p * rot(-a);
}

vec2 scene(vec3 p) {
    float d = 0.;
    vec3 _p = p;
    vec3 q = p;
    float m = 0.;
    
    float scale = .8;
    vec3 size = vec3(1.);
    
    float repeatNum = 12.;
    
    vec3 id = floor(_p / repeatNum);
    
    _p = repeat(_p, repeatNum);  
    
    d = sdWireBox(_p, .22, .07);    
    
    vec2 obj = vec2(d, 1.);
    
	_p.xy = pmod(_p.xy, repeatNum);
	size = vec3(1.);
    
    for(int i = 0; i < 7; i++) {
        _p.y -= size.x * .32;
        _p.y -= size.y;
        _p.xz *= rot(sin(TIME * .4+ float(i) * .4) * PI * 2.);

        float currentDist = sdWireBox(p, .22 * size.y, .07 * size.y);
		vec2 currentObj = vec2(currentDist, float(i) + 1.);
        m = minMat(obj, currentObj).y;
        obj.x = min(obj.x, currentObj.x);
        obj.y = m;
        
        p = _p;
        size *= scale;
    }
    
    return vec2(obj.x, m);
}

mat3 camera(vec3 o, vec3 t, vec3 u) {
	vec3 forward = normalize(t - o);
    vec3 right = cross(forward, u);
    vec3 up = cross(right, forward);
    return mat3(right, up, forward);
}

vec3 getNormal(vec3 p, float eps) {
    vec2 e = vec2(eps, 0);
    return normalize(
    	vec3(
            scene(p + e.xyy).x - scene(p - e.xyy).x,
            scene(p + e.yxy).x - scene(p - e.yxy).x,
            scene(p + e.yyx).x - scene(p - e.yyx).x
        )
    );
}

float getSpecular(vec3 position, vec3 normal, Light light, float diffuse, vec3 cameraPos) {
  	vec3 lightDir = light.position - position;
  	vec3 ref = reflect(-normalize(lightDir), normal);
  	float specular = 0.;
  	if(diffuse > 0.) {
    	specular = max(0., dot(ref, normalize(cameraPos - normal)));
    	float specularPower = 64.;
    	specular = pow(specular, specularPower);
  	}
  	return specular;
}

void calcDirectionalLight(inout Light light, Surface surface, vec3 cameraPos) {
  	// diffuse
  	float diffuseCoef = max(0., dot(surface.normal, normalize(light.position)));
  	vec3 diffuse = diffuseCoef * light.attenuation * light.color * light.intensity;
  	// specular
  	float specularCoef = getSpecular(surface.position, surface.normal, light, diffuseCoef, cameraPos);
  	vec3 specular = vec3(specularCoef * light.attenuation * light.color * light.intensity);  

  	light.diffuse = diffuse;
  	light.specular = specular;
}

vec3 lighting(Surface surface, vec3 cameraPos) {
  	vec3 position = surface.position;

  	vec3 color = vec3(0.);
  	vec3 normal = surface.normal;

  	Light directionalLight;
  	directionalLight.position = vec3(-.2, .2, -1.);
  	directionalLight.intensity = .8;
  	directionalLight.color = vec3(1., 1., 1.);
  	directionalLight.attenuation = 1.;
  	calcDirectionalLight(directionalLight, surface, cameraPos);

  	vec3 diffuse = directionalLight.diffuse;
  	vec3 specular = directionalLight.specular;
    
  	vec3 ambient = vec3(.2);
  
  	color = surface.baseColor * diffuse + ambient;  
  
  	return color;
}

vec3 fog(vec3 color, float distance, vec3 fogColor, float b) {
  	float fogAmount = 1. - exp(-distance * b);
  	return mix(color, fogColor, fogAmount);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = (fragCoord.xy - RENDERSIZE.xy * .5) / min(RENDERSIZE.x, RENDERSIZE.y);
    vec2 mouse = (iMouse.xy - RENDERSIZE.xy * .5) / min(RENDERSIZE.x, RENDERSIZE.y);
    
    // ref: https://www.shadertoy.com/view/WllSDM
    float nois = (
        n3(vec3(normalize(uv) * 20., 1.)) +
        n3(vec3(normalize(uv) * 40., 2.)) * 1.5  +
        n3(vec3(normalize(uv) * 100., 3.)) * 1.5
    ) / 4.;
    nois = pow(nois, 4.) * 2. - 1.;    
    float tt = TIME * .11 + (nois * .01 + n3(vec3(uv * 600., 0.)) * .003) * pow(dot(uv, uv), 0.3);    
    float speed = tt * -50.;    
    
    vec3 movement = vec3(0., 0., - speed);    
    vec3 ro = vec3(6., 6., 4.) + movement;
    vec3 target = ro + vec3(0., 0., 0.) + movement;
    float fov = .7;
    
    vec3 up = vec3(0., 1., 0.);
    up.xy *= rot(TIME * .12);
    vec3 rd = camera(ro, target, up) * normalize(vec3(uv, fov));
    
    // raymarching
    float depth = 0.;
    float dist = 0.;
    vec2 result = vec2(0.);
    for(int i = 0; i < 64; i++) {
    	result = scene(ro + rd * depth);
        dist = result.x;
        if(dist < stopThreshold) {
            break;
        }
        depth += result.x;
    }   
    
    vec3 color = vec3(0.);
   
    vec3 position = ro + rd * depth;
    vec3 normal = getNormal(position, .01);
    
    Surface surface;
    surface.depth = depth;
    surface.dist = dist;
    surface.position = position;
    surface.normal = normal;
    surface.material = result.y;
    
	vec3 sceneColor = vec3(0.);
    
    vec3 bgColor = vec3(.05, .05, .15);
       
    if(dist >= stopThreshold) {
	    // no hit
        sceneColor = bgColor;
    } else {
	    // hit
	    surface.baseColor = palette(
    		surface.material / 9.,
        	vec3(.7),
        	vec3(.75),
        	vec3(1.),
        	vec3(.5 + TIME * .14, TIME * .18, TIME * .24)    
    	);
		sceneColor = lighting(surface, ro);           
    }

    // fog
    sceneColor = fog(sceneColor, depth, bgColor, .045);
	// vignet
    sceneColor *= smoothstep(1.3, .6, length(uv.xy));
 
    fragColor = vec4(sceneColor, 1.);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}