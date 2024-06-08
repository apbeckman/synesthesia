#include "hg_sdf.glsl"
#include "lygia/sdf.glsl"
#include "lygia/math/const.glsl"
#define CAMERA_POSITION     cam
#define SAMPLE2DCUBE_FLIP_Y
uniform sampler2D u_tex0;

#define RESOLUTION          RENDERSIZE
#define RAYMARCH_SAMPLES    100
#define RAYMARCH_MULTISAMPLE 4
#define RAYMARCH_BACKGROUND ( vec3(0.1, 0.1, .10) + ray.y * 0.8 )
#define RAYMARCH_AMBIENT    vec3(0.7, 0.9, 1.0)
#define RAYMARCH_MAX_DIST 55

#define SAMPLE2DCUBE_FLIP_Y
// #define SAMPLE2DCUBE_CELL_SIZE 64.0
// #define SAMPLE2DCUBE_CELLS_PER_SIDE 8.0
#define SAMPLE2DCUBE_FNC(TEX, UV) sampleBicubic(TEX, UV, vec2(512.0))
// #define SAMPLE2DCUBE_FNC(TEX, UV) sampleSmooth(TEX, UV, vec2(512.0))
// #define SAMPLE2DCUBE_FNC(TEX, UV) sampleNearest(TEX, UV, vec2(512.0))

#include "lygia/math/cubic.glsl"
#include "lygia/math/quartic.glsl"
#include "lygia/math/quintic.glsl"
#include "lygia/sdf.glsl"
#include "lygia/space/ratio.glsl"
#include "lygia/space/scale.glsl"
#include "lygia/sample/bicubic.glsl"
#include "lygia/sample/smooth.glsl"
#include "lygia/sample/nearest.glsl"
#include "lygia/sample/3DSdf.glsl"
#include "lygia/math/saturate.glsl"
#include "lygia/lighting/raymarch.glsl"
#include "lygia/color/space/linear2gamma.glsl"

vec3 cam = vec3(vec2(1., -1.), -100);
float smin(float c, float b, float k) {
    float h = max(k - abs(c - b), 0.) / k;
    return min(c, b) - h * h * h * k * 1. / 6.;
}

//3D Rotation
vec3 rot3d(in vec3 x, in vec3 r) {
    x.xy = _rotate(x.xy, r.z);
    x.yz = _rotate(x.yz, r.x);
    x.zx = _rotate(x.zx, r.y);
    return x;
}
#include "lygia/lighting/raymarch.glsl"

float checkBoard(vec2 uv, vec2 _scale) {
for (int i = 0; i <=2;i++){
    uv += _noise(abs(_uv)+smoothTimeC*0.01);
}
uv = floor(fract(uv * _scale) * 12.0);
return min(1.0, uv.x + uv.y) - (uv.x * uv.y);
}

vec4 raymarchMap( in vec3 pos ) {
    vec3 ones = vec3(1.0);
    vec3 pyr_pos = pos;
    pyr_pos += vec3(0., 1.0, 0.)*-1;
    for (int i = 0; i <=4; i++){
        if (i % 2 == 0.){
            // pyr_pos.
            pyr_pos= rot3d(pyr_pos, vec3(0.0, 0.0, -smoothTime*0.01));
            pyr_pos.y += 0.5;
            
            // pyr_pos.yz = -pyr_pos.yz;
        }
        else{
            pyr_pos= rot3d(pyr_pos, vec3(0.0125*smoothTimeC, 0., 0.));
        }
        pyr_pos = abs(pyr_pos);
        pyr_pos.x += -.125;
        pyr_pos.z += -.125;
        pyr_pos.y += -.125;
        pyr_pos = abs(pyr_pos);
    }
    vec3 pos_mod = vec3(0., .75, 0.)*-1;
    vec3 c_pos = pos;
    for (int i = 0; i <= 4; i++){
        // c_pos.x = smin(-c_pos.x, c_pos.x, 0.02);
        if (i % 2 == 0.){
            // c_pos.xy = _rotate(c_pos.xy, 0.2*TIME);
            c_pos *= 1.52;
            c_pos= rot3d(c_pos, vec3(0.50, 0.05*smoothTime, -smoothTime*0.01));
            c_pos.x -= 0.5;
            c_pos.y -= 0.5;
            c_pos.z -= 0.45;
            c_pos *= 0.7;
        }
        c_pos= rot3d(c_pos, vec3(0.0, 0.1, smoothTime*0.01));
        c_pos.x = abs(c_pos.x);
            c_pos.x -= 0.5;
            c_pos.y -= 0.5;
            c_pos.z -= 0.45;
        c_pos.y = abs(c_pos.y);
        c_pos.y -= 0.3;
        c_pos= rot3d(c_pos, vec3(0.0, smoothTimeC*0.02, 0.1));
    }
    // c_pos+=pos_mod;
    vec4 res = vec4(1.);
    float tetrahedron = tetrahedronSDF(pyr_pos, size);
    // res = opUnion(res, vec4(ones, tetrahedron), 0.1);
    res = opUnion(res, vec4( 1.0, 1.0, 1.0, sample3DSdf(u_tex0, c_pos) ));
    
    float check = checkBoard(vec2(_noise(pos.xz)), vec2(5.0));
    res = opUnion( res, vec4( vec3( 0.5 + check * 0.5), planeSDF(pos + _noise(_uv+TIME*0.01) + 4.0) ), 0.25 );
        
    return res;
}


vec4 renderMain() {
vec4 fragColor = vec4(0.);
vec3 color = vec3(0.0);
vec2 pixel = 1.0 / RENDERSIZE;
vec2 st = gl_FragCoord.xy * pixel;
vec2 uv = ratio(st, RENDERSIZE);

vec2 mo = _mouse.xy * pixel;
float time = 32.0 + TIME * 1.5;
vec3 cam = vec3(6.5 , 8.2, 5.5 ) * 15.*fov;
cam = rot3d(cam, vec3(0., orbit.y, 0.));
cam.xz = _rotate(cam.xz, TIME*0.05 + orbit.x*PI);
// cam.yz = _rotate(cam.yz, 0.5);
// vec3 cam = vec3(4.5 * cos(0.1 * time + 7.0 * mo.x), 2.2, 4.5 * sin(0.1 * time + 7.0 * mo.x)) * 10.;

color = raymarch(cam, uv).rgb;
color = linear2gamma(color);

fragColor = vec4(color, 1.0);
return fragColor;
}
    // vec3 n = normal(p);
    // vec3 l = normalize(lightPos - p);
    // float diff = clamp(dot(n, l), 0.0, 1.0);
    // vec4 fragColor = vec4(vec3(diff), 1.0);
