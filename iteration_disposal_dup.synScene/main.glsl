// lygia includes
#include "lygia/space/ratio.glsl"
#include "lygia/sdf.glsl"
#include "lygia/lighting/ray.glsl"
#include "lygia/color/saturationMatrix.glsl"
#include "lygia/color/levels.glsl"
#include "lygia/color/blend.glsl"
#define RESOLUTION RENDERSIZE
//Raymarch Definitions before calling #include "lygia/lighting/raymarch.glsl"

//Raymarch Light
vec3 light_dir = vec3(0.0,0.0,50.0);
vec3 light_pos = vec3(0.0,10.0,50.0);
vec3 cam = vec3(-cam_offset*vec2(1.,-1.),-100+cam_zoom);
vec3 surface = vec3(0.);

#define LIGHT_POSITION light_pos
#define LIGHT_DIRECTION     light_dir
#define RESOLUTION          RENDERSIZE
#define LIGHT_COLOR vec3(0.1, 0.1, 0.1)
#define RAYMARCHSOFTSHADOW_ITERATIONS 4
#define SURFACE_POSITION    surface
#define CAMERA_POSITION     cam
vec3 raymarchPbrMaterial(vec3 ray, vec3 pos, vec3 nor, vec3 map);

//Raymarch Step Count
#define RAYMARCH_SAMPLES 40
#define RAYMARCH_MULTISAMPLE 1
#define RAYMARCH_MAX_DIST 20
#define RAYMARCH_SAMPLES_AO ambient_occlusion

//Background
#define RAYMARCH_BACKGROUND vec3(0.0)
#define RAYMARCH_AMBIENT vec3(0.1)
#define RAYMARCH_MAP_MATERIAL rgb

//Material
#define RAYMARCH_MATERIAL_FNC   raymarchPbrMaterial

//more lygia includes
#include "lygia/lighting/raymarch.glsl"
#include "lygia/color/space/linear2gamma.glsl"
#include "lygia/lighting/envMap.glsl"
#include "lygia/lighting/reflection.glsl"
#include "lygia/lighting/pbr.glsl"
#include "lygia/lighting/material/new.glsl"
#include "lygia/lighting/raymarch/glass.glsl"

//Mercury shader library
#include "hg_sdf.glsl"
float smin(float a, float b, float k) {
    float h = max(k - abs(a-b), 0.) / k;
    return min(a, b) - h*h*h*k*1./6.;
}
vec3 adjustSaturation(vec3 color, float saturationAmount) {
    mat4 satMat = saturationMatrix(saturationAmount);
    return vec3(satMat * vec4(color, 1.0));
}

vec3 adjustLevels(vec3 color, vec3 iMin, vec3 gamma, vec3 iMax, vec3 oMin, vec3 oMax) {
    return levels(color, iMin, gamma, iMax, oMin, oMax);
}

//3D Rotation
vec3 rot3d(in vec3 x, in vec3 r) {
    x.xy = _rotate(x.xy, r.z);
    x.yz = _rotate(x.yz, r.x);
    x.zx = _rotate(x.zx, r.y);
    return x;
}

//Raymarch Map
vec4 raymarchMap( in vec3 pos ) {

    vec3 pos2 = pos;
    vec3 pos3 = pos;
    vec3 pos4 = pos;

    vec4 res = vec4(1.0);
    
    //click to move
    pos.xy += _muvc.xy*_click.x;
    pos.xy += vec2(0.0,0.0);
    pos2.xy += _muvc.xy*_click.x;
    pos2.xy += vec2(0.0,0.0);
    pos3.xy += _muvc.xy*_click.x;
    pos3.xy += vec2(0.0,0.0);
    
    pos.z += 0.5 + -_ncos(syn_BassLevel*2.5*kick);
    pos2.z += 0.8 + -_ncos(inner_kick);
    pos.z += 1.0 + -_ncos(middle_kick);
    pos3.z += 1.0 + -_ncos(outer_kick);
    pos.z += hitvals;
    pos2.z += hitvals;
    pos3.z += hitvals;

    //repeat positive direction z axis
    float repite = pModSingle1(pos.z, 6.5);
    float repite2 = pModSingle1(pos2.z, 6.5);
    float repite3 = pModSingle1(pos3.z, 6.5);

    //ring rotations
    float clckwse =synvals;
    float outward =bassvals;
    float inner =midvals;
    pos = rot3d(pos, vec3(0.0,0.0,clckwse));
    pos2 = rot3d(pos2, vec3(0.0,0.0,-clckwse));
    pos3 = rot3d(pos3, vec3(0.0,0.0,-clckwse));    

    //polar repeat
    float cell = pModPolar(pos.xy, more_its);
    float cell2 = pModPolar(pos2.xy, more_its);
    float cell3 = pModPolar(pos3.xy, more_its);

    //rotations
    vec3 boxrot = vec3(0.0,-bassvals,0.0);
    vec3 boxrot2 = vec3(0.0,bassvals,0.0);
    vec3 trirot = vec3(0.0,-midvals,0.0);
    
    //SDF
    float hole_size = 1. -cos(syn_Presence+hole);
    
    float o1 = fOctahedron(rot3d(pos2-vec3( 1.0+hole_size,0.0,0.0), trirot), 0.23);
    float b1 = boxSDF(rot3d(pos2-vec3( 1.0+hole_size,0.0, 0.0), trirot), vec3(0.26) );
    float b2 = boxSDF(rot3d(pos-vec3( 1.94+hole_size,0.0, 0.0), boxrot2), vec3(0.3) );    
    float o2 = fOctahedron(rot3d(pos-vec3( 1.94+hole_size,0.0,0.0), boxrot2), 0.3);
    float b3 = boxSDF(rot3d(pos3-vec3( 3.05+hole_size,0.0, 0.1), trirot), vec3(0.40) );
    float o3 = fOctahedron(rot3d(pos3-vec3( 3.05+hole_size,0.0,0.1), trirot), 0.4);
    float b4 = boxSDF(rot3d(pos-vec3( 4.61+hole_size,0.0, 0.1), boxrot2), vec3(0.55) );    
    float o4 = fOctahedron(rot3d(pos-vec3( 4.61+hole_size,0.0,0.1), boxrot2), 0.6);
    float b5 = boxSDF(rot3d(pos3-vec3( 6.75+hole_size,0.0, 0.1), trirot), vec3(0.75) );    
    float o5 = fOctahedron(rot3d(pos3-vec3( 6.75+hole_size,0.0,0.1), trirot), 0.8);
    float b6 = boxSDF(rot3d(pos-vec3( 9.55+hole_size,0.0, 0.1), boxrot2), vec3(0.85) );    
    float o6 = fOctahedron(rot3d(pos-vec3( 9.55+hole_size,0.0,0.1), boxrot2), 0.8);
    
    float sCam = fSphere(pos4-cam/10, 1.0);
    
    res = opUnion( res, vec4( box_color, b1 ) );
    res = opUnion( res, vec4( tri_color, o1 ) );
    res = opUnion( res, vec4( box_color, b2 ) );
    res = opUnion( res, vec4( tri_color, o2 ) );
    res = opUnion( res, vec4( box_color, b3 ) );
    res = opUnion( res, vec4( tri_color, o3 ) );
    res = opUnion( res, vec4( box_color, b4 ) );
    res = opUnion( res, vec4( tri_color, o4 ) );
    res = opUnion( res, vec4( box_color, b5 ) );
    res = opUnion( res, vec4( tri_color, o5 ) );
    res = opUnion( res, vec4( box_color, b6 ) );
    res = opUnion( res, vec4( tri_color, o6 ) );
    
    float stairs = fOpUnionStairs(b1, o1, 0.2, 8.);
    float stairs2 = fOpUnionStairs(b2, o2, 0.30, 8.);
    float stairs3 = fOpUnionStairs(b3, o3, 0.30, 8.);
    float stairs4 = fOpUnionStairs(b4, o4, 0.45, 8.);
    float stairs5 = fOpUnionStairs(b5, o5, 0.50, 8.);
    float stairs6 = fOpUnionStairs(b6, o6, 0.50, 8.);
    
    res.w = min(stairs, res.w);
    res.w = min(stairs2, res.w);
    res.w = min(stairs3, res.w);
    res.w = min(stairs4, res.w);
    res.w = min(stairs5, res.w);
    res.w = min(stairs6, res.w);
    res.w = opSubtraction(sCam, res.w);

    return res;
}

vec3 raymarchPbrMaterial(vec3 ray, vec3 pos, vec3 nor, vec3 map) {
    // Define a background color for the scene
    vec3 backgroundColor = vec3(0.0, 0.0, 0.0);

    // If no object is hit, return the background color
    if (sum(map) <= 0.0) {
        return backgroundColor;
    }
    // Initialize material properties
    Material material = materialNew();
    material.albedo.rgb = vec3(map); // Color of the material
    material.roughness = .92;       // Controls the spread of the specular highlight
    material.metallic = 1.0;        // Partially metallic material
    material.normal = nor;          // Surface normal
    material.f0 = vec3(0.04);       // Reflectivity at normal incidence
    material.ior = vec3(0.0);       // Index of refraction

    // Calculate color using PBR glass material
    vec3 color = pbr(material).rgb;
    color.rgb *= tonemap(sphericalHarmonics(nor));

    // Calculate ambient occlusion
    float ao = raymarchAO(pos, nor);
    float aoStrength = 1.0; // Value between 0 and 1 to scale the AO effect
    color *= mix(1.0, ao, aoStrength); // Apply the ambient occlusion to the color

    // Apply level adjustments

    color.rgb = adjustLevels(color.rgb, vec3(.04), vec3(0.95), vec3(1.0), vec3(0.0), vec3(1.0));

    // Calculate glow based on the Fresnel effect
    float viewAngle = max(dot(nor, normalize(CAMERA_POSITION - pos)), 0.0);
    float fresnelEffect = pow(1.0 - viewAngle, 1.77); // Adjust the exponent as needed
    vec3 glowColor = glow_color; // glow color

    // Blend the glow effect on top of the material color
    vec3 finalColor = blendGlow(color, glowColor * fresnelEffect, fresnelEffect);

    return linear2gamma(finalColor); // Convert to gamma space before returning
}


vec4 renderMain(void) {
    vec2 pixel = _uvc + vec2(0.5);
    vec3 color = raymarch(CAMERA_POSITION, pixel).rgb;

    // Apply brightness before bloom effect
    color = _brightness(color, 1.0);

    // Extract bright areas for bloom
    vec3 brightAreas = max(color - 10., .61); // Threshold can be adjusted

    // Sample the final pass texture to simulate a bloom effect
    vec3 bloom = vec3(0.0);
    float bloomRadius = 0.00001; // The radius of the bloom effect
    int samples = 1; // How many samples in each direction to take
    int totalSamples = samples * samples;
    for (int i = -samples / 2; i <= samples / 2; i++) {
        for (int j = -samples / 2; j <= samples / 2; j++) {
            vec2 offset = bloomRadius * vec2(i, j);
            bloom += texture(syn_FinalPass, _uv + offset).rgb;
        }
    }
    bloom /= float(totalSamples); // Average the sampled values

    // Combine original color with "bloomed" bright areas
    vec3 finalColor = color + bloom * brightAreas;

    return vec4(finalColor, 1.0);
}
