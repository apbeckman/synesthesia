			//******** Common Code Begins ********
float dist(vec2 p0, vec2 pf){return sqrt((pf.x-p0.x)*(pf.x-p0.x)+(pf.y-p0.y)*(pf.y-p0.y));}
float d = dist(RENDERSIZE.xy*0.5,_xy.xy)*(_mouse.x/RENDERSIZE.x+0.1)*0.005;
float d2 = dist(RENDERSIZE.xy*0.5,_xy.xy)*0.001;
float shock = shock_back - shock_forward;

#define aT(p) texelFetch(BuffA, ivec2(mod(p,R)), 0)
#define bT(p) texelFetch(BuffB, ivec2(mod(p,R)), 0)
#define A(p) texture(BuffA, mod(p,R)/R)
#define B(p) texture(BuffB, mod(p,R)/R)
#define C(p) texture(BuffC, mod(p,R)/R)
#define D(p) texture(BuffD, mod(p,R)/R)
#define PI 3.14159265
vec4 media_mixed = texture(media_pass_fx, _uv) * media_impact;
#define dt 1. 
#define R RENDERSIZE.xy

const vec2 dx = vec2(0, 1);

float hash11(float p) {
    p = fract(p * .1031);
    p *= p + 33.33;
    p *= p + p;
    return fract(p);
}

#define rand_interval 250
#define random_gen(a, b, seed) ((a) + ((b)-(a))*hash11(seed + float(FRAMECOUNT/rand_interval)))

#define distribution_size 1.4 +length(media_mixed)

/* FIRE
//mold stuff 
#define sense_num 6
#define sense_ang 1.
#define sense_dis 420.
#define sense_oscil 0.1
#define oscil_scale 1.
#define oscil_pow 2.
#define sense_force -0.35
#define force_scale 1.
#define trailing 0.
#define acceleration 0.
*/

#define sense_num 12
#define sense_ang 0.2
#define sense_dis 11.
#define sense_oscil 0.0
#define oscil_scale 1.
#define oscil_pow 1.
#define sense_force 0.2
#define distance_scale 2.
#define force_scale 1.
#define trailing 0.
#define acceleration 0.09

//useful functions
#define GS(x) exp(-dot(x,x))
#define GSS(x) exp(-dot(x,x))
#define GS0(x) exp(-length(x))
#define Dir(ang) vec2(cos(ang), sin(ang))
#define Rot(ang) mat2(cos(ang), sin(ang), -sin(ang), cos(ang))
#define loop(i,x) for(int i = 0; i < x; i++)
#define range(i,a,b) for(int i = a; i <= b; i++)

//SPH pressure
#define Pressure(rho) 0.3*rho.z +length(media_mixed)
#define fluid_rho 0.2

//data packing
#define PACK(X) ( uint(round(65534.0*clamp(0.5*X.x+0.5, 0., 1.))) + \
65535u * uint(round(65534.0 * clamp(0.5 * X.y + 0.5, 0., 1.))))   

#define UNPACK(X) (clamp(vec2(X%65535u, X/65535u)/65534.0, 0.,1.)*2.0 - 1.0)              

#define DECODE(X) UNPACK(floatBitsToUint(X))
#define ENCODE(X) uintBitsToFloat(PACK(X))

			//******** BuffA Code Begins ********

// const int KEY_SPACE = 32;
// bool isKeyPressed(int KEY)
// {
// 	return texelFetch( iChannel3, ivec2(KEY,0), 0 ).x > 0.5;
// }

vec4 renderPassA() {
vec4 U = vec4(0.0);
vec2 pos = _xy;

ivec2 p = ivec2(pos);

vec2 X = vec2(0);
vec2 V = vec2(0);
float M = 0.;

    //basically integral over all updated neighbor distributions
    //that fall inside of this pixel
    //this makes the tracking conservative
range(i, - 2, 2) range(j, - 2, 2) {
vec2 tpos = pos + vec2(i, j);
vec4 data = bT(tpos);
tpos -= _uvc*(Zoom+shock);
tpos += d2*_uvc*(Fisheye+0.5*shock)* (1.0 + 2.0 * low + syn_Intensity * 1.5);

vec2 X0 = DECODE(data.x) + tpos;
vec2 V0 = DECODE(data.y);
vec2 M0 = data.zw;

X0 += V0 * dt; //integrate position

        //particle distribution size
float K = distribution_size;

vec4 aabbX = vec4(max(pos - 0.5, X0 - K * 0.5), min(pos + 0.5, X0 + K * 0.5)); //overlap aabb
vec2 center = 0.5 * (aabbX.xy + aabbX.zw); //center of mass
vec2 size = max(aabbX.zw - aabbX.xy, 0.); //only positive

        //the deposited mass into this cell
float m = M0.x * size.x * size.y / (K * K); 
        //add weighted by mass
X += center * m;
V += V0 * m;

        //add mass
M += m;
}

    //normalization
if(M != 0.) {
X /= M;
V /= M;
}

    //mass renormalization
float prevM = M;
M = mix(M, 0.2, 0.05);
V = V * prevM / M;

    //initial condition
if(FRAMECOUNT <= 1 || reset != 0.) {
X = pos;
vec2 dx0 = (pos - R * 0.3);
vec2 dx1 = (pos - R * 0.7);
V = 0.5 * Rot(PI * 0.5) * dx0 * GS(dx0 / 30.) - 0.5 * Rot(PI * 0.5) * dx1 * GS(dx1 / 30.);
V += 0.2 * Dir(2. * PI * hash11(floor(pos.x / 20.) + R.x * floor(pos.y / 20.)));
M = 0.1 + pos.x / R.x * 0.01 + pos.y / R.x * 0.01;
}

X = clamp(X - pos, vec2(- 0.5), vec2(0.5));
U = vec4(ENCODE(X), ENCODE(V), M, 0.);
return U;
} 

			//******** BuffB Code Begins ********

vec4 renderPassB() {
vec4 U = vec4(0.0);
vec2 pos = _xy;

vec2 uv = pos / R;
ivec2 p = ivec2(pos);

vec4 data = aT(pos);
vec2 X = DECODE(data.x) + pos;
vec2 V = DECODE(data.y);
float M = data.z;

if(M != 0.) //not vacuum
{
        //Compute the force
vec2 F = vec2(0.);

        //get neighbor data
vec4 d_u = aT(pos + dx.xy), d_d = aT(pos - dx.xy);
vec4 d_r = aT(pos + dx.yx), d_l = aT(pos - dx.yx);

        //position deltas
vec2 p_u = DECODE(d_u.x), p_d = DECODE(d_d.x);
vec2 p_r = DECODE(d_r.x), p_l = DECODE(d_l.x);

        //velocities
vec2 v_u = DECODE(d_u.y), v_d = DECODE(d_d.y);
vec2 v_r = DECODE(d_r.y), v_l = DECODE(d_l.y);

        //pressure gradient
vec2 p = vec2(Pressure(d_r) - Pressure(d_l), Pressure(d_u) - Pressure(d_d));

        //density gradient
vec2 dgrad = vec2(d_r.z - d_l.z, d_u.z - d_d.z);

        //velocity operators
float div = v_r.x - v_l.x + v_u.y - v_d.y;
div += 0.01 * length(media_mixed) * 0.1;

float curl = v_r.y - v_l.y - v_u.x + v_d.x;
        //vec2 laplacian = 

F -= M * p;

float ang = atan(V.y, V.x);
float dang = sense_ang * PI / float(sense_num);
vec2 slimeF = vec2(0.);
        //slime mold sensors
range(i, - sense_num, sense_num) {
float cang = ang + float(i) * dang;
vec2 dir = (1. + sense_dis * pow(M, distance_scale)) * Dir(cang);
vec3 s0 = C(X + dir).xyz;
float fs = pow(s0.z, force_scale);
float os = oscil_scale * pow(s0.z - M, oscil_pow);
os += length(media_mixed);
slimeF += sense_oscil * Rot(os) * s0.xy + sense_force * Dir(ang + sign(float(i)) * PI * 0.5) * fs;
}

        //remove acceleration component and leave rotation
       // slimeF -= dot(slimeF, normalize(V))*normalize(V);
F += slimeF / float(2 * sense_num);

if(_mouse.z > 0.) {
vec2 dx = pos - _mouse.xy;
F += 0.1 * Rot(PI * 0.5) * dx * GS(dx / 30.);
}

        //integrate velocity
V += Rot(- 0. * curl) * F * dt / M;

        //acceleration for fun effects
V *= 1. + acceleration;

        //velocity limit
float v = length(V);
V /= (v > 1.) ? v / 1. : 1.;
}

    //mass decay
   // M *= 0.999;

    //input
    //if(_mouse.z > 0.)
    //\\	M = mix(M, 0.5, GS((pos - _mouse.xy)/13.));
    //else
     //   M = mix(M, 0.5, GS((pos - R*0.5)/13.));

    //save
X = clamp(X - pos, vec2(- 0.5), vec2(0.5));
U = vec4(ENCODE(X), ENCODE(V), M, 0.);
return U;
} 

			//******** BuffC Code Begins ********

#define radius 2.

vec4 renderPassC() {
vec4 fragColor = vec4(0.0);
vec2 pos = _xy;
float rho = 0.001;
vec2 vel = vec2(0., 0.);

    //compute the smoothed density and velocity
range(i, - 2, 2) range(j, - 2, 2) {
vec2 tpos = pos + vec2(i, j);
vec4 data = bT(tpos);

vec2 X0 = DECODE(data.x) + tpos;
vec2 V0 = DECODE(data.y);
float M0 = data.z;
vec2 dx = X0 - pos;

float K = GS(dx / radius) / (radius);
rho += M0 * K;
vel += M0 * K * V0;
}

vel /= rho;

fragColor = vec4(vel, rho, 1.0);
return fragColor;
} 

			//******** BuffD Code Begins ********

//normal advection
const int KEY_SPACE = 32;
// bool isKeyPressed(int KEY)
// {
	// return bool(reset);
// }

vec4 renderPassD() {
vec4 fragColor = vec4(0.0);
vec2 pos = _xy;

vec2 V0 = vec2(0.);
if(FRAMECOUNT <= 1 || reset != 0.) {
vec4 data = bT(pos);
V0 = 1. * DECODE(data.y);
float M0 = data.z;
} else {

}

fragColor.xy = C(pos - V0 * dt).xy - V0 * dt / R;

    //initial condition
if(FRAMECOUNT <= 1 || reset != 0.) {
fragColor.xy = vec2(0.);
}
return fragColor;
} 

// Fork of "Fireballs" by michael0884. https://shadertoy.com/view/tlfBDX
// 2020-08-20 00:44:41

// Fork of "Random slime mold generator" by michael0884. https://shadertoy.com/view/ttsfWn
// 2020-08-19 23:28:40

// Fork of "Everflow" by michael0884. https://shadertoy.com/view/ttBcWm
// 2020-07-19 18:18:22

// Fork of "Paint streams" by michael0884. https://shadertoy.com/view/WtfyDj
// 2020-07-11 22:38:47

//3d mode
//#define heightmap

vec3 hsv2rgb(in vec3 c) {
vec3 rgb = clamp(abs(mod(c.x * 6.0 + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0, 0.0, 1.0);

rgb = rgb * rgb * (3.0 - 2.0 * rgb); // cubic smoothing	

return c.z * mix(vec3(1.0), rgb, c.y);
}

#define FOV 1.56+fov
#define RAD R.x*0.7

float gauss(float x, float r) {
x /= r;
return exp(- x * x);
}

float sdSegment(in vec2 p, in vec2 a, in vec2 b) {
vec2 pa = p - a, ba = b - a;
float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
return length(pa - ba * h);
}

float sdSphere(vec3 p, float s) {
return length(p) - s;
}

float sdBox(vec3 p, vec3 b) {
vec3 q = abs(p) - b;
return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float rho(vec3 pos) {
pos.xy += R * 0.5;
pos.xy = mod(pos.xy, R - 1.);
vec4 v = B(pos.xy);
return v.z;
}

vec3 color(vec3 pos) {
pos.xy += R * 0.5;
pos.xy = mod(pos.xy, R - 1.);

vec4 v = C(pos.xy);
v += texture(media_pass, _uv);
return v.xyz;
}

float DE(vec3 pos) {
float y = 1. * rho(pos);

pos.xy += R * 0.5;
pos.xy = mod(pos.xy, R - 1.);
float de = 1e10;
de = min(de, 0.7 * sdBox(pos - vec3(R, 4. * y) * 0.5, vec3(R * 0.51, 3.)));
return de;
}

vec4 calcNormal(vec3 p, float dx) {
const vec3 K = vec3(1, - 1, 0);
return (K.xyyx * DE(p + K.xyy * dx) +
    K.yyxx * DE(p + K.yyx * dx) +
    K.yxyx * DE(p + K.yxy * dx) +
    K.xxxx * DE(p + K.xxx * dx)) / vec4(4. * dx, 4. * dx, 4. * dx, 4.);
}

#define marches 70.
#define min_d 1.
vec4 ray_march(vec3 p, vec3 r) {
float d;
for(float i = 0.;
i < marches;
i ++) {
d = DE(p);
p += r * d;
if(d < min_d || d > R.x) break;
}
return vec4(p, d);
}

vec4 renderMainImage() {
vec4 col = vec4(0.0);
vec2 pos = _xy;

    #ifdef heightmap
        // Normalized pixel coordinates 
pos = (pos - R * 0.5) / max(R.x, R.y);
pos = vec2(pos.x, pos.y);
vec2 uv = _mouse.xy / R;
vec2 angles = vec2(- 1.5 + 0.0 * smoothTimeC, - 0.5) * PI;

vec3 camera_z = vec3(cos(angles.x) * cos(angles.y), sin(angles.x) * cos(angles.y), sin(angles.y));
vec3 camera_x = normalize(vec3(cos(angles.x + PI * 0.5), sin(angles.x + PI * 0.5), 0.));
vec3 camera_y = - normalize(cross(camera_x, camera_z));

        //tracking particle
vec4 fp = vec4(R * 0.5 + 150. * vec2(smoothTimeC, 0.), 0., 0.);

vec3 ray = normalize(camera_z + FOV * (pos.x * camera_x + pos.y * camera_y));
vec3 cam_pos = vec3(fp.xy - R * 0.5, 0.) - RAD * vec3(cos(angles.x) * cos(angles.y), sin(angles.x) * cos(angles.y), sin(angles.y));

vec4 X = ray_march(cam_pos, ray);

if(X.w < min_d) {

float D = rho(X.xyz);
vec3 c = color(X.xyz);
if(_exists(syn_UserImage)) {
vec3 albedo = 0.5 * (D + 0.07) * mix(texture(image5, c.xy).xyz, _grayscale(texture(media_pass_fx, _uv)).xyz, media_impact).xyz;
                // vec3 colA = mix(texture(image5,  rd.yzx).xyz, texture(media_pass_fx, rd.yxz).xyz, media_impact);
} else {
vec3 albedo = 0.5 * (D + 0.07) * texture(image5, c.xy).xyz;
}

float rough = 1. - 0.1 * distance(albedo, vec3(1.));

vec4 N0 = calcNormal(X.xyz, 2. * X.w) * vec4(4., 4., 1., 1.);
vec3 n = normalize(N0.xyz);
vec3 rd = reflect(ray, n);
if(_exists(syn_UserImage)) {
vec3 colA = mix(texture(image5, rd.yzx).xyz, texture(media_pass_fx, rd.yzx).xyz, media_impact).xyz;
} else {
vec3 colA = texture(image5, rd.yzx).xyz;
}
vec3 colB = (vec3(0.5) + 0.5 * dot(rd, normalize(vec3(1.))));
colB += 3. * rough * pow(max(dot(rd, normalize(vec3(1.))), 0.), 10.);
colB += 3. * rough * pow(max(dot(rd, normalize(vec3(- 0.5, - 0.9, 0.8))), 0.), 10.);
float b = clamp(0.5 + 0.5 * dot(n, normalize(vec3(1, 1, 1))), 0., 1.);
float K = 1. - pow(max(dot(n, rd), 0.), 2.);
col.xyz = 1. * albedo * colB + 0. * rough * colA * K;
} else {
if(_exists(syn_UserImage)) {
vec3 cc = mix(texture(image5, ray.yzx).xyz, texture(media_pass_fx, _uv).xyz, media_impact);
col = 1. * cc;
} else {
col = 1. * texture(image5, ray.yzx);
}

            //background
}
col = tanh(8. * col);
    #else
float r = B(pos.xy).z;
vec4 c = D(pos.xy);

    	//get neighbor data
vec4 d_u = bT(pos + dx.xy), d_d = bT(pos - dx.xy);
vec4 d_r = bT(pos + dx.yx), d_l = bT(pos - dx.yx);

        //position deltas
vec2 p_u = DECODE(d_u.x), p_d = DECODE(d_d.x);
vec2 p_r = DECODE(d_r.x), p_l = DECODE(d_l.x);

        //velocities
vec2 v_u = DECODE(d_u.y), v_d = DECODE(d_d.y);
vec2 v_r = DECODE(d_r.y), v_l = DECODE(d_l.y);

        //pressure gradient
vec2 p = vec2(Pressure(d_r) - Pressure(d_l), Pressure(d_u) - Pressure(d_d));

        //velocity operators
float div = (v_r.x - v_l.x + v_u.y - v_d.y);
float curl = (v_r.y - v_l.y - v_u.x + v_d.x);

if(_exists(syn_UserImage)) {
col = 1.7 * r * mix(texture(image5, c.xy + pos / R), texture(media_pass_fx,_uv), media_impact);
} else {
col = 1.7 * r * texture(image5, c.xy + pos / R);
}
    	//col.xyz += vec3(1,0.1,0.1)*max(curl,0.) + vec3(0.1,0.1,1.)*max(-curl,0.);

    #endif
return col;
}
vec4 mediaPass() {
vec4 media = _loadMedia();
return media;
}
vec4 mediaPassFX() {
vec4 mediaFX = mix(texture(media_pass, _uv), _edgeDetectSobel(media_pass, _uv), edge_mix);
mediaFX.xyz = _pixelate(_grayscale(mediaFX.xyz), 0.25);
return mediaFX;

}

vec4 renderMain() {
if(PASSINDEX == 0) {
return renderPassA();
}
if(PASSINDEX == 1) {
return renderPassB();
}
if(PASSINDEX == 2) {
return renderPassC();
}
if(PASSINDEX == 3) {
return renderPassD();
}
if(PASSINDEX == 4) {
return mediaPass();
}
if(PASSINDEX == 5) {
return mediaPassFX();
}
if(PASSINDEX == 6) {
return renderMainImage();
}
}