

			//******** Common Code Begins ********


#define R RENDERSIZE.xy

// Dave Hoskins https://www.shadertoy.com/view/4djSRW
vec2 hash23(vec3 p3)
{
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+33.33);
    return fract((p3.xx+p3.yz)*p3.zy);
}

float gyroid (vec3 seed)
{
    return dot(sin(seed),cos(seed.yzx));
}

float fbm (vec3 seed)
{
    float result = 0.;
    float a = .5;
    for (int i = 0; i < 3; ++i) {
        seed += result / 2.;
        result += gyroid(seed/a)*a;
        a /= 2.;
    }
    return result;
}

// Fork of "generative art deco" by morisil. https://shadertoy.com/view/7sKfDd
// 2022-09-28 11:25:15

// Copyright Kazimierz Pogoda, 2022 - https://xemantic.com/
// I am the sole copyright owner of this Work.
// You cannot host, display, distribute or share this Work in any form,
// including physical and digital. You cannot use this Work in any
// commercial or non-commercial product, website or project. You cannot
// sell this Work and you cannot mint an NFTs of it.
// I share this Work for educational purposes, and you can link to it,
// through an URL, proper attribution and unmodified screenshot, as part
// of your educational material. If these conditions are too restrictive
// please contact me and we'll definitely work it out.

// copyright statement borrowed from Inigo Quilez

// Music by Giovanni Sollima, L'invenzione del nero:
// https://soundcloud.com/giovanni-sollima/linvenzione-del-nero

// See also The Mathematics of Perception to check the ideas behind:
// https://www.shadertoy.com/view/7sVBzK

const float SHAPE_SIZE = .618;
const float CHROMATIC_ABBERATION = .02;
const float ITERATIONS = 10.;
const float INITIAL_LUMA = .4;

//const float PI = 3.14159265359;
const float TWO_PI = 6.28318530718;

mat2 rotate2d(float _angle){
    return mat2(cos(_angle),-sin(_angle),
                sin(_angle),cos(_angle));
}

float sdPolygon(in float angle, in float distance) {
  float segment = TWO_PI / 4.0;
  return cos(floor(.5 + angle / segment) * segment - angle) * distance;
}

float getColorComponent(in vec2 st, in float modScale, in float blur) {
    vec2 modSt = mod(st, 1. / modScale) * modScale * 2. - 1.;
    float dist = length(modSt);
    float angle = atan(modSt.x, modSt.y) + sin(smoothTimeC * .08) * 9.0;
    dist = sdPolygon(angle, dist);
    float shapeMap = smoothstep(SHAPE_SIZE + blur, SHAPE_SIZE - blur, dist);
    return shapeMap;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 st =
        (2.* fragCoord - RENDERSIZE.xy)
        / min(RENDERSIZE.x, RENDERSIZE.y);
    vec2 origSt = st;
    st *= rotate2d(sin(TIME * .14) * .3);
    st *= (sin(TIME * .15) + 2.) * .3;
    st *= log(length(st * .28)) * .8;

    float modScale = 1.;

    vec3 color = vec3(0);
    float luma = INITIAL_LUMA;
    float blur = .5 + sin(smoothTimeC * .5) * .3;
    for (float i = 0.; i < ITERATIONS; i++) {
        vec2 center = st + vec2(sin(smoothTime * .12), cos(smoothTime * .13));
        //center += pow(length(center), 1.);
        vec3 shapeColor = vec3(
            getColorComponent(center - st * CHROMATIC_ABBERATION, modScale, blur),
            getColorComponent(center, modScale, blur),
            getColorComponent(center + st * CHROMATIC_ABBERATION, modScale, blur)        
        ) * luma;
        st *= 1.1 + getColorComponent(center, modScale, .04) * 1.;
        st *= rotate2d(sin(TIME  * .05) * .3);
        color += shapeColor;
        color = clamp(color, 0., 1.);
        if (color == vec3(1)) break;
        luma *= .6;
        blur *= .63;
    }
    const float GRADING_INTENSITY = .4;
    vec3 topGrading = vec3(
        1. - sin(TIME * 1.1 * .3) * GRADING_INTENSITY,
        1. - sin(TIME * 1.2 * .3) * GRADING_INTENSITY,
        1. - sin(TIME * 1.3 * .3) * GRADING_INTENSITY
    );
    vec3 bottomGrading = vec3(
        1. + cos(TIME * 1.4 * .3) * GRADING_INTENSITY,
        1. - cos(TIME * 1.5 * .3) * GRADING_INTENSITY,
        1. - cos(TIME * 1.6 * .3) * GRADING_INTENSITY
    );
    vec3 colorGrading = mix(topGrading, bottomGrading, length(origSt));
    fragColor = vec4(pow(color.rgb, colorGrading), 1.);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}