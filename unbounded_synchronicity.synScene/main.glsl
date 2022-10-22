vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


#define fdist 0.5
#define iters 90
#define tol 0.003
//#define maxdist 80.
#define eps 0.01
#define tmin 0.1
#define bevel 0.1
#define innerradius 0.5
#define sphereradius 0.45
#define sidelength 0.55
float gridlength;
#define speed 0.5

#define PI 3.1415926

float box(vec3 ro, vec3 rd) {
    ro = mod(ro, gridlength) - gridlength/2.;
    vec3 t = -ro/rd + abs(gridlength/rd/2.);
    return min(t.x, min(t.y, t.z));
}

void rotX(inout vec3 ro, float a) {
    float c = cos(a);
    float s = sin(a);
    ro.yz = mat2(c, -s, s, c) * ro.yz;
}
void rotY(inout vec3 ro, float a) {
    float c = cos(a);
    float s = sin(a);
    ro.zx = mat2(c, -s, s, c) * ro.zx;
}
void rotZ(inout vec3 ro, float a) {
    float c = cos(a);
    float s = sin(a);
    ro.xy = mat2(c, -s, s, c) * ro.xy;
}

vec2 map(vec3 ro) {
    vec3 roo = ro/gridlength;
    ivec3 b = ivec3(floor(roo.x), floor(roo.y), floor(roo.z));
    vec3 sgn = vec3(2*(abs(b) % 2)-1);
    vec3 offset = (cos(min(mod(speed*smoothTimeC + vec3(2., 3., 4.) * PI, 3.*PI),PI))+1.) * sgn.zxy;
    vec3 displacement = offset * gridlength/2.;
    vec3 rotation = offset * PI/4.;
    roo = ro;
    ro = mod(roo + displacement, gridlength) - gridlength/2.;
    rotX(ro, rotation.z);
    rotY(ro, rotation.x);
    rotZ(ro, rotation.y);
    //beveled cube
    vec3 rad = clamp(ro, -sidelength, sidelength);
    float d = length(ro-rad) - bevel;
    
    d = max(d, -length(ro.xy)+innerradius);
    d = max(d, -length(ro.yz)+innerradius);
    d = max(d, -length(ro.zx)+innerradius);
    
    ro = mod(roo+displacement.zxy+displacement.xyz, gridlength) - gridlength/2.;
    float d2 = length(ro) - sphereradius;
    if (d2 < d) {
        return vec2(d2, 2.0);
    } else {
        return vec2(d, 1.0);
    }
}

vec2 map(vec3 ro, vec3 rd) {
    float t = box(ro, rd);
    vec2 m = map(ro);
    if (m.x < t) {
        return m;
    }
    return vec2(t + tol, 0.0);
}

vec3 raytrace(vec3 ro, vec3 rd) {
    float t = 0.;
    vec2 m = map(ro, rd);
    int i;
    for (i=0; i<iters; i++) {
        t += m.x;
        m = map(ro + rd * t, rd);
        if (abs(m.x) < tol) {
            return vec3(t, float(i)/float(iters), m.y);
        } //else if (m.x > maxdist) break;
    }
    return vec3(t, float(i)/float(iters), 0.0);
}

vec4 getnormal(vec3 ro) {
    vec2 d = vec2(eps, 0.0);
    float x1 = map(ro+d.xyy).x;
    float x2 = map(ro-d.xyy).x;
    float y1 = map(ro+d.yxy).x;
    float y2 = map(ro-d.yxy).x;
    float z1 = map(ro+d.yyx).x;
    float z2 = map(ro-d.yyx).x;
    return vec4(normalize(vec3(
        x1-x2,
        y1-y2,
        z1-z2)),
        x1+x2+y1+y2+z1+z2-6.*map(ro).x);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    gridlength = 2.1 + 2.5*(-cos(smoothTime/25.)+1.);
    //vec3 ro = vec3(4.*cos(TIME), (iMouse.y-RENDERSIZE.y/2.)/10., 4.*sin(TIME));
    vec3 ro = vec3(gridlength/2., gridlength/2., 0.);
    vec3 w = vec3(-1., 0., 0.);
    vec3 u = normalize(cross(w, vec3(0., 1.0, 0.)));
    vec3 v = cross(u, w);
    vec3 rd = normalize(w*fdist+(fragCoord.x/RENDERSIZE.x-0.5)*u+(fragCoord.y-RENDERSIZE.y/2.0)/RENDERSIZE.x*v);
    vec3 d = raytrace(ro, rd);
    vec4 n = getnormal(ro + rd*d.x);
	vec3 col = clamp(500.*n.w, 0., 0.25) + d.y*(n.xyz+1.0)/2.0;
    if (d.z > 1.5) {
        vec3 refl = reflect(rd, n.xyz);
        d = raytrace(ro + rd*d.x + refl*tmin, refl);
        vec4 n2 = getnormal(ro + rd*d.x);
        col = mix(col, d.y*(n2.xyz+1.0)/2.0, pow(1.-max(0., dot(n.xyz, -rd)), 1.4));
    }
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}