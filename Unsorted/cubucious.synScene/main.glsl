

// ponk (Leon Denise) 19/07/2018
// most lines below are from the shadertoy community
// licensed under hippie love conspiracy
// happy tweaking

// Geometry
float range = .8;
float radius = .4;
float blend = .3;
float count = 8.+(Count);

// Light
vec3 lightPos = vec3(1, 1, 1);
float specularSharpness = 10.;
float glowSharpness = 1.;

// Colors
vec3 ambient = vec3(.1);
vec3 light = vec3(0);
vec3 specular = vec3(1);
vec3 glow = vec3(1);

// Raymarching
const float epsilon = .0001;
const float steps = 25.;
const float far = 20.;
#define repeat(p,r) (mod(p,r)-r/2.)
#define sdist(p,r) (length(p)-r)
float box (vec3 p, vec3 b) { vec3 d = abs(p) - b; return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0)); }
float torus (vec3 p, vec2 t) { vec2 q = vec2(length(p.xz)-t.x,p.y); return length(q)-t.y; }
float smoothmin (float a, float b, float r) { float h = clamp(.5+.5*(b-a)/r, 0., 1.); return mix(b, a, h)-r*h*(1.-h); }
mat2 rot (float a) { float c=cos(a),s=sin(a); return mat2(c,-s,s,c); }
vec3 look (vec3 eye, vec3 target, vec2 anchor) {
	vec3 forward = normalize(target-eye);
	vec3 right = normalize(cross(forward, vec3(0,1,0)));
	vec3 up = normalize(cross(right, forward));
	return normalize(forward + right * anchor.x + up * anchor.y);
}

// Miscellaneous
#define time TIME*.4
//#define PI 3.141592
#define TAU PI *2.0
#define PIHALF PI/2.0
#define PIQUART PI/4.0
#define saturate(p) clamp(p,0.,1.)
float random (in vec2 st) { return fract(sin(dot(st.xy, vec2(12.9898+sin(smoothTime),78.233)))* 43758.5453123); }

float geometry (vec3 pos)
{
    float scene = 10.;
    vec3 p = pos;
	for (float index = count; index > 0.; --index) {
		float ratio = index / count;
        
        // easing
		ratio *= ratio;
        
        // domain reptition and translation offset
		p.xz = abs(p.xz) - range * ratio*(Expand);
        p.xy *= -rot(smoothTimeC*0.0475+sin(smoothTimeC*0.0475));
        // rotations
		p.xz *= rot(PIQUART+smoothTime*0.0675);
		p.yz *= -rot(smoothTimeC*0.0375+TAU);
		p.yx *= rot(PIHALF);

		scene = smoothmin(scene, box(p, vec3(radius * ratio)), blend * ratio);
	}
    return scene;
}

vec3 getNormal (vec3 p) {
    vec2 e = vec2(epsilon,0);
    return normalize(vec3(geometry(p+e.xyy)-geometry(p-e.xyy),
                          geometry(p+e.yxy)-geometry(p-e.yxy),
                          geometry(p+e.yyx)-geometry(p-e.yyx)));
}

void raymarching (vec3 pos, vec3 ray, inout vec4 hit)
{
	float total = 0.;
	for (float i = steps; i >= 0.; --i) {
		float dist = geometry(pos);
		if (dist < epsilon * total || total > far) {
			hit.xyz = pos;
			hit.w = i/steps;
			break;
		}
		total += dist;
		pos += ray * dist;
	}
}

vec4 renderMainImage() {
	vec4 color = vec4(0.0);
	vec2 coordinate = _xy;
    
    vec2 mouse = _mouse.xy / RENDERSIZE.xy;
    mouse.x = (mouse.x * 2. - 1.) * PI;
    mouse.y *= 1.5;
    
	vec2 uv = (coordinate.xy-.5*RENDERSIZE.xy)/RENDERSIZE.y;
    vec3 eye = vec3(0,2,2) * (2. - Zoom);
    vec3 target = vec3(0);
	vec4 hit;
    eye.xz *= rot(Rotate.x);
	eye.yz *= rot(Rotate.y); //adding rotate control AB
	eye.xy += (_rotate(eye.xy, smoothTime*0.01)*PI/10.);

    lightPos.xz *= rot(smoothTimeB*0.5);
    
	vec3 ray = look(eye, target, uv);
    raymarching(eye, ray, hit);
    
    vec3 pos = hit.xyz;
	vec3 normal = getNormal(pos);
	vec3 lightDir = normalize(lightPos);
	float lightIntensity = clamp(dot(lightDir, normal),0.,1.);
	float specularIntensity = saturate(pow(max(0., dot(reflect(lightDir, normal), ray)), specularSharpness));
	float glowIntensity = saturate(pow(abs(1.-abs(dot(normal, ray))), glowSharpness));

	color.rgb = ambient + light * lightIntensity + specular * (specularIntensity*(1.0+highhits)) + glow * glowIntensity;
	color.rgb *= hit.w;
    color.rgb *= step(length(eye-pos), far);
	//color.rgb = normal * .5 + .5;
	return color; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}