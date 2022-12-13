vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


#define PI 3.14159

#define FIXED_UP vec3(0.,1.,0.)
#define TAN_HFOVY .5773502691896257

#define GE1 vec3(.04,0.,0.)
#define GE2 vec3(0.,.04,0.)
#define GE3 vec3(0.,0.,.04)
#define STEP_D .06

#define RES RENDERSIZE
#define MS iMouse
#define PHASE smoothstep(.2, .8, MS.z > 0. ? MS.y/RES.y : .5+.25*cos(smoothTimeC))
#define ANG (MS.z > 0. ? 2.*PI*MS.x/RES.x : PI + (smoothTime))
float smoothTime = (smooth_basstime * 0.75 + TIME * 0.125 + syn_Time * 0.55)*0.25;
float smoothTimeB = (smooth_hightime * 0.75 + TIME * 0.125 + syn_Time * 0.55)*0.25;
float smoothTimeC = (smooth_midtime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.25;

// ================
// Generic helpers
// ================

float square(float x) { return x*x; }

float rand(vec2 p) {
    return fract(sin(dot(p,vec2(12.9898,78.233))) * 43758.5453);
}

vec3 nvCamDirFromClip(vec3 nvFw, vec2 c) {
    vec3 nvRt = normalize(cross(nvFw, FIXED_UP));
    vec3 nvUp = cross(nvRt, nvFw);
    return normalize(TAN_HFOVY*(c.x*(RES.x/RES.y)*nvRt + c.y*nvUp) + nvFw);
}

// Rodrigues' formula: v -> (v.n)n + (v-(v.n)n)c - (vxn)s
mat3 oProd(vec3 n) {
    float xy = n.x*n.y, xz = n.x*n.z, yz = n.y*n.z;
    return mat3(n.x*n.x,xy,xz,  xy,n.y*n.y,yz,  xz,yz,n.z*n.z);
}
mat3 glRotate(vec3 axis, float angle) {
    float c = cos(angle), s = sin(angle);
    vec3 n = normalize(axis);
    return (
        (1.-c)*oProd(n) +
        mat3(c,s*n.z,-s*n.y,  -s*n.z,c,s*n.x,  s*n.y,-s*n.x,c)
	);
}

// Colormap adapted from: iq - Palettes (https://www.shadertoy.com/view/ll2GD3)
vec3 colormap(float t) {
    return .5 + .5*cos(2.*PI*( t + vec3(0.,.1,.2) ));
}

// ==================
// Volume definition
// ==================

vec3 hopf(vec3 p) {
	// inverse stereographic R^3 -> S^3
    float psq = dot(p, p);
    vec4 q = vec4(2.*p.xyz, -1. + psq) / (1. + psq);

    // hopf map S^3 -> S^2 (map versor q -> qjq*, i.e., rotate j using q)
    return vec3(
        2.*(q.y*q.z - q.x*q.w),
        1. - 2.*(q.y*q.y + q.w*q.w), // == q.x^2 - q.y^2 + q.z^2 - q.w^2
        2.*(q.x*q.y + q.z*q.w)
    );
}

vec4 getV(vec3 p, mat3 m) {
    // Get point on S^2, and rotate it for visual effect
    vec3 rp = m * hopf(p);

    // Fixed y == union of linked circles of same radii (interp. by theta)
	float theta = atan(-rp.z, rp.x);
    float sdy = square(.5 + .5*cos(rp.y*5.*PI));
    float sdt = .5 + .5*sin(theta*20.);
    float d1 = max(sdt*sdy, .55*sdy);
    d1 = mix(d1, .8*sdy, PHASE); // Sampling from a cubemap here is interesting

    return vec4(.5 + .5*rp.y, 0., 0., d1);
}

vec4 getC(vec3 p, vec3 camPos, mat3 m) {
    // Quick exit if volume has low alpha
    vec4 data0 = getV(p, m);
    if (data0.a < .4) { return vec4(0.); }

    // Normal
    vec3 gradA = vec3(
        getV(p + GE1, m).a - getV(p - GE1, m).a,
        getV(p + GE2, m).a - getV(p - GE2, m).a,
        getV(p + GE3, m).a - getV(p - GE3, m).a
    );
    vec3 normal = (-gradA/(length(gradA) + 1e-5));

    // Material color with lighting
    vec3 matColor = mix(colormap(data0.r), .5+.5*normal.xyz, .25);
    vec3 lightPos = camPos + vec3(0., 1., 0.);
    vec3 nvFragToLight = normalize(lightPos - p);
    vec3 nvFragToCam = normalize(camPos - p);

	vec3 diffuse = clamp(dot(normal, nvFragToLight), 0., 1.) * matColor;
    vec3 blinnH = normalize(nvFragToLight + nvFragToCam);
    vec3 specular = pow(clamp(dot(normal, blinnH), 0., 1. ), 60.) * vec3(1.);

    return vec4(.4*matColor+.6*(diffuse+specular), pow(data0.w, 3.));
}

// =============
// Render scene
// =============

void march(in vec3 p, in vec3 nv, out vec4 color) {
    nv.xy += _uvc*PI*nv.xy*MouseXY*PI;
    p.xy += _uvc*PI*p.xy*MouseXY*PI;

    color = vec4(0.);
    vec2 tRange = vec2(mix(1., 2.75, PHASE), 6.);

	mat3 m = glRotate( normalize(vec3(1.,sin(smoothTime*0.5),1.)), -ANG );

    float t = tRange.s;
    for (int i=0; i<200; i++) {
        vec4 cHit = getC(p + t*nv, p, m);

        vec4 ci = vec4(cHit.rgb, 1.)*( (STEP_D/.08)*cHit.a ); // Premul alpha
        color += (1.-color.a)*ci;

        t += STEP_D;
        if (t > tRange.t || color.a > .95) { return; }
    }
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RES.xy;
    vec3 camPos = vec3(0., 1., 3.5);
    vec3 nvCamDir = nvCamDirFromClip(normalize(-camPos), uv*2. - 1.);

    vec4 objColor;
    march(camPos + rand(fragCoord)*nvCamDir*STEP_D, nvCamDir, objColor);
    vec3 finalColor = objColor.rgb + (1. - objColor.a)*vec3(.0);


    fragColor = vec4(finalColor, 1.);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}