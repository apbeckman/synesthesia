

// "Fork" of "Geodesic tiling (abs position)" by tdhooper. https://shadertoy.com/view/XtKSWc
// Most of the geodesic stuff is gone, but the march/render code is still in place
// various sdf functions via iq and mercury
// 2021-01-01 23:50:17


// --------------------------------------------------------
// Modelling
// --------------------------------------------------------

struct Model {
    float dist;
    vec3 color;
};



// --------------------------------------------------------
// Camera
// https://www.shadertoy.com/view/Xl2XWt
// --------------------------------------------------------

mat3 calcLookAtMatrix( in vec3 ro, in vec3 ta, in float roll )
{
    vec3 ww = normalize( ta - ro );
    vec3 uu = normalize( cross(ww,vec3(sin(roll),cos(roll),0.0) ) );
    vec3 vv = normalize( cross(uu,ww));
    return mat3( uu, vv, ww );
}

float pModPolar(inout vec2 p, float repetitions, float phase) {
	float angle = 2.*3.14159/repetitions;
	float a = atan(p.y, p.x) + angle/2.+phase/repetitions;
	float r = length(p);
	float c = floor(a/angle);
	a = mod(a,angle) - angle*.5;
    a = abs(a);
	p = vec2(cos(a), sin(a))*r;
	// For an odd number of repetitions, fix cell index of the cell in -x direction
	// (cell index would be e.g. -5 and 5 in the two halves of the cell):
	if (abs(c) >= (repetitions/2.)) c = abs(c);
	return c;
}
float pMod1(inout float p, float size) {
	float halfsize = size*0.5;
	float c = floor((p + halfsize)/size);
	p = mod(p + halfsize, size) - halfsize;
	return c;
}


vec4 fOpUnionRound(vec4 a, vec4 b, float r) {
	float h = clamp( 0.5 + 0.5*(b.x-a.x)/r, 0.0, 1.0 );
    return mix( b, a, h ) - r*h*(1.0-h);
}
    
vec4 opSmoothSubtraction( vec4 d1, vec4 d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2.x+d1.x)/k, 0.0, 1.0 );
    return mix( d2, vec4(-d1.x, d1.yzw), h ) + k*h*(1.0-h); 
}
vec4 opUnion( vec4 d1, vec4 d2 ) {  return mix(d1,d2, clamp(ceil(d1.x-d2.x), 0., 1.)); }

void pR(inout vec2 p, float a) {
	p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}
float vmax(vec2 v) {
	return max(v.x, v.y);
}

float vmax(vec3 v) {
	return max(max(v.x, v.y), v.z);
}

float vmax(vec4 v) {
	return max(max(v.x, v.y), max(v.z, v.w));
}

float vmin(vec2 v) {
	return min(v.x, v.y);
}

float vmin(vec3 v) {
	return min(min(v.x, v.y), v.z);
}

float vmin(vec4 v) {
	return min(min(v.x, v.y), min(v.z, v.w));
}

float fBox(vec3 p, vec3 b) {
	vec3 d = abs(p) - b;
	return length(max(d, vec3(0))) + vmax(min(d, vec3(0)));
}

float fTorus(vec3 p, float smallRadius, float largeRadius) {
	return length(vec2(length(p.xz) - largeRadius, p.y)) - smallRadius;
}


// The actual model
Model map(vec3 p) {
    
    float sineTime = TIME*6.28*.5;

	float sphere = length(p) - 2.; 
    sphere = max(-(length(p)-1.8), sphere);
    vec2 uv = p.xy*2.;
   
    float f = floor(uv.x+uv.y) + floor(uv.y-uv.x);
    uv = vec2(fract(uv.x+uv.y),fract(uv.y-uv.x));
    uv -= .5;
    
    vec3 pos = p;
    //pos.y-=.08;
    float i = pMod1(pos.y, .25);
    float j = pModPolar(pos.xz,30.,i*3.14159+(sineTime+sin(sineTime+3.1415)*1.5)*(mod(i, 2.)-.5));
    pos -= vec3(.8-cos((i)*.2+3.14159), 0., 0.);
    float eye = length(pos+vec3(0., .04, 0.))-.1;
    eye = max(eye, length(pos-vec3(0., .04, 0.))-.1);
    vec4 s = fOpUnionRound(vec4(sphere, 0., 0., 0.), vec4(eye, 2., 2., 2.), 0.05);
    s = opSmoothSubtraction(vec4(eye-.005, 0., 0., 0.), s, .01);
    pR(pos.xy, p.y*.5+sin(i+j+sineTime)*.05);
    float a =atan(pos.z, pos.y);
    vec3 eyeColor = sin(vec3((p.x+p.y*2.+vec3(0., .3, .6))*7.+sineTime))*.4+.4;
    eyeColor+=cos(a*20.+a*3.)*.2+.2;
    float pupilSize = 0.025+sin(sineTime-(p.y-p.x)*2.)*0.005;
    eyeColor*=smoothstep(pupilSize, pupilSize+.01, length(pos.zy));
    s = opUnion(s, vec4(length(pos)-.06, eyeColor));
    p*=.8;
    //p.y -= .01;
    pR(p.yz, 1.);
    pos = p;
    i = pModPolar(pos.xz, 10., sineTime);
    pR(pos.xy, -sineTime*.5);
    j = pModPolar(pos.yx, 5., sineTime);
    float box = fBox(pos, vec3(0., .1*smoothstep(0., .1,length(p-vec3(0., .1, 0.))), 0.));
    vec4 b = vec4(box, sin((pos.y*8.+vec3(0., .66, .33))*7.-sineTime+atan(p.x, p.z))*.5+.5);
    b = fOpUnionRound(b, vec4(length(p)-.05, 1, 1, 1), .1);
    s = opUnion(s, b); 
	return Model(s.x, s.yzw);
}


// --------------------------------------------------------
// Ray Marching
// Adapted from cabbibo https://www.shadertoy.com/view/Xl2XWt
// --------------------------------------------------------

const float MAX_TRACE_DISTANCE = 10.;
const float INTERSECTION_PRECISION = .00001;
const int NUM_OF_TRACE_STEPS = 100;

struct CastRay {
    vec3 origin;
    vec3 direction;
};

struct Ray {
    vec3 origin;
    vec3 direction;
    float len;
};

struct Hit {
    Ray ray;
    Model model;
    vec3 pos;
    bool isBackground;
    vec3 normal;
    vec3 color;
};

vec3 calcNormal( in vec3 pos ){
    vec3 eps = vec3( 0.001, 0.0, 0.0 );
    vec3 nor = vec3(
        map(pos+eps.xyy).dist - map(pos-eps.xyy).dist,
        map(pos+eps.yxy).dist - map(pos-eps.yxy).dist,
        map(pos+eps.yyx).dist - map(pos-eps.yyx).dist );
    return normalize(nor);
}
    
Hit raymarch(CastRay castRay){

    float currentDist = INTERSECTION_PRECISION * 2.0;
    Model model;
    
    Ray ray = Ray(castRay.origin, castRay.direction, 0.);

    for( int i=0; i< NUM_OF_TRACE_STEPS ; i++ ){
        //if (currentDist < INTERSECTION_PRECISION || ray.len > MAX_TRACE_DISTANCE) {
        //    break;
        //}
        model = map(ray.origin + ray.direction * ray.len);
        currentDist = model.dist;
        ray.len += currentDist;
    }
    
    bool isBackground = false;
    vec3 pos = vec3(0);
    vec3 normal = vec3(0);
    vec3 color = vec3(0);
    
    if (ray.len > MAX_TRACE_DISTANCE) {
        isBackground = true;
    } else {
        pos = ray.origin + ray.direction * ray.len;
        normal = calcNormal(pos);
    }

    return Hit(ray, model, pos, isBackground, normal, color);
}


// --------------------------------------------------------
// Rendering
// --------------------------------------------------------

vec3 render(Hit hit){
    if (hit.isBackground) {
        return vec3(0);
    }
    vec3 color = hit.model.color;
    color += sin(dot(hit.normal, vec3(0,1,0))) * .2; // lighting
    color *= 1. - clamp(hit.ray.len * .4 - .8, 0., 1.); // fog
    return color;
}


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    vec2 p = (-RENDERSIZE.xy + 2. * fragCoord.xy) / RENDERSIZE.y;

    vec3 camPos = vec3(0, 0, .4);
    vec3 camTar = vec3(0);
    float camRoll = 0.;
    mat3 camMat = calcLookAtMatrix(camPos, camTar, camRoll);
    
    vec3 rd = normalize(camMat * vec3(p.xy, 2.));
    Hit hit = raymarch(CastRay(camPos, rd));

    vec3 color = render(hit);
    float f = .85;
    //color *= smoothstep(f+.05, f, abs(p.y-.08)-abs(p.x)*.1);
    fragColor = vec4(color,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}