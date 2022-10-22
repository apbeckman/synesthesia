vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********


#define DARK_MODE

const float MAX_DISPERSE = 5.;
const float MAX_BOUNCE = 6.;

mat4 cameraMatrix = mat4(
    0.9780874848365784,
    -0.07870610803365707,
    -0.19274398684501648,
    0,
    0.20812205970287323,
    0.3452317416667938,
    0.9151507019996643,
    0,
    -0.005486522801220417,
    -0.935211718082428,
    0.35404732823371887,
    0,
    -4.901630745735019e-7,
    0.00000874467059475137,
    -9.502106666564941,
    1
);


// HG_SDF
// https://www.shadertoy.com/view/Xs3GRB
//--------------------------------------------------------

#define PI 3.14159265359
#define TAU 6.28318530718

#define saturate(x) clamp(x, 0., 1.)

void pR(inout vec2 p, float a) {
    p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

float smax(float a, float b, float r) {
    vec2 u = max(vec2(r + a,r + b), vec2(0));
    return min(-r, max (a, b)) + length(u);
}

float vmax(vec2 v) {
	return max(v.x, v.y);
}

float fBox(vec2 p, vec2 b) {
	vec2 d = abs(p) - b;
	return length(max(d, vec2(0))) + vmax(min(d, vec2(0)));
}

float range(float vmin, float vmax, float value) {
  return clamp((value - vmin) / (vmax - vmin), 0., 1.);
}

mat4 rotX(float a) {
    return mat4(1,0,0,0, 0,cos(a),-sin(a),0, 0,sin(a),cos(a),0, 0,0,0,1);
}

mat4 rotZ(float a) {
    return mat4(cos(a),-sin(a),0,0, sin(a),cos(a),0,0, 0,0,1,0, 0,0,0,1);
}

// Spectrum palette
// IQ https://www.shadertoy.com/view/ll2GD3
//--------------------------------------------------------

vec3 pal( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d ) {
    return a + b*cos( 6.28318*(c*t+d) );
}

vec3 spectrum(float n) {
    return pal( n, vec3(0.5,0.5,0.5),vec3(0.5,0.5,0.5),vec3(1.0,1.0,1.0),vec3(0.0,0.33,0.67) );
}

// Polyhedra
// Knighty https://www.shadertoy.com/view/MsKGzw
//--------------------------------------------------------

int Type=3;

vec3 nc,pab,pbc,pca;
void initPoly() {//setup folding planes and vertex
    float cospin=cos(PI/float(Type)), scospin=sqrt(0.75-cospin*cospin);
	nc=vec3(-0.5,-cospin,scospin);//3rd folding plane. The two others are xz and yz planes
	pab=vec3(0.,0.,1.);
	pbc=vec3(scospin,0.,0.5);//No normalization in order to have 'barycentric' coordinates work evenly
	pca=vec3(0.,scospin,cospin);
	pbc=normalize(pbc);	pca=normalize(pca);//for slightly better DE. In reality it's not necesary to apply normalization :) 
}

vec3 fold(vec3 pos) {
	for(int i=0;i<3 /*Type*/;i++){
		pos.xy=abs(pos.xy);//fold about xz and yz planes
		pos-=2.*min(0.,dot(pos,nc))*nc;//fold about nc plane
	}
	return pos;
}


//========================================================
// Modelling
//========================================================

float time;

float tweenOffset(float t, float start, float duration) {
    t = range(start, start + duration, t);
    t = pow(t, 2.);
    return t;
}

float tweenBlend(float t, float start, float duration) {
    t = range(start, start + duration, t);
    t = pow(t, .5);
    t = smoothstep(0., 1., t);
    return t;
}

const float STEP_SCALE = 1./3.;

float tetAnim(vec3 p, float time) {
    
    p = fold(p);

    // config
    float sz = .3;
    float rBase = .04;
    float rInner = rBase * STEP_SCALE;
    float blendDuration = .75;
    float offsetDuration = .75;
    float t = time * (blendDuration + offsetDuration);
    offsetDuration *= 2.; // extend animation beyond loop section
    float offsetDistance = .6;

    // animation
    float blend = tweenBlend(t, .0, blendDuration);
    float offsetT = tweenOffset(t, blendDuration, offsetDuration);
    float offset = offsetT * offsetDistance;

    // skip if animation hasn't started or is complete
    if (t < 0. || offsetT >= 1.) {
        return 1e12;
    }

    // tetrahedron planes
    vec3 n1 = pca;
    vec3 n2 = normalize(pca * vec3(-1,-1,1));
    vec3 n3 = normalize(pbc * vec3(1,-1,-1));
    vec3 n4 = normalize(pbc * vec3(-1,-1,-1));

    float sep = .001 * (1. - offsetT);
    float scale = 1. - offsetT;
    
    float bound = (dot((p + (n4 + n3) * offset) / scale, n1) - sz) * scale;
    if (bound > .004) {
        return bound;
    }

    vec3 pp = p;
    
    // octahedrons
    p = pp + n4 * offset;
    p /= scale;
    float oct = dot(p, n1) - sz;
    oct = smax(oct, (dot(p, n2) - sz), rBase);
    oct = smax(oct, -(dot(p, n4) + .5 - sep), rInner);
    oct = smax(oct, -(dot(p, n3) + .1 - sep), rInner);
    oct = smax(oct, (dot(p, n4) + .1 + sep), rInner);
    oct *= scale;

    // shared edge
	p = pp + (n4 + n3) * offset;
    p /= scale;
    float edge = dot(p, n1) - sz;
    edge = smax(edge, (dot(p, n2) - sz), rBase);    
    
    // side tets
    p = pp + (n4 + n3) * offset;
    p /= scale;
    float side = edge;
    side = smax(side, (dot(p, n3) + .1 + sep), rInner);
    side = smax(side, (dot(p, n4) + .1 + sep), rInner);
    side *= scale;

    // vertex tets
    p = pp + n4 * (offset + offset);
    p /= scale;
    float vert = edge; vert = smax(vert, (dot(p, n3) - sz), rBase);
    vert = smax(vert, (dot(p, n4) + .5 + sep), rInner);
    vert *= scale;

    float sliced = 1e12;
    sliced = min(sliced, oct);
    sliced = min(sliced, vert);
    sliced = min(sliced, side);

    // inner tet gets replaced with the next iteration
    if (time < 1.) {
        p = pp;
        float inner = -(dot(p, n4) + .1 - sep);
        inner = smax(inner, -(dot(p, n3) + .1 - sep), rInner);
        inner = smax(inner, -(dot(p, n2) + .1 - sep), rInner);
        sliced = min(sliced, inner);
    }

    if (blend >= 1.) {
        return sliced;
    }

    // base tet
    float base = dot(p, n1) - sz;
    base = smax(base, (dot(p, n2) - sz), rBase);
    base = smax(base, (dot(p, n3) - sz), rBase);

    float surface = 1. - saturate(-base / sz); // 1 at surface, 0 at center
    
    // blend indentations into the surface
    float surfaceBlend = saturate(blend * .66 * range(.9, 1., surface));
    base = mix(base, sliced, surfaceBlend);
    
    // grow the sliced tet from the center of the unsliced tet
    float slicedS = min(sliced, -base - (.3 - .3 * blend));
    float d = max(base, slicedS);
    d = mix(d, sliced, smoothstep(.9, 1., blend));

    return d;
}

float tetLoop(vec3 p) {
    pR(p.xy, PI/2. * -time + PI/2.);

    float t = time;
    float scale = pow(STEP_SCALE, t);
    float d = tetAnim(p * scale, time) / scale;

    scale *= STEP_SCALE;
    pR(p.xy, PI/2. * -1.);
    d = min(d, tetAnim(p * scale, time + 1.) / scale);

    return d;
}

vec2 map(vec3 p) {
    float d = tetLoop(p);
    return vec2(d, 1.);
}


//========================================================
// Lighting
//========================================================

float intersectPlane(vec3 rOrigin, vec3 rayDir, vec3 origin, vec3 normal, vec3 up, out vec2 uv) {
    float d = dot(normal, (origin - rOrigin)) / dot(rayDir, normal);
  	vec3 point = rOrigin + d * rayDir;
	vec3 tangent = cross(normal, up);
	vec3 bitangent = cross(normal, tangent);
    point -= origin;
    uv = vec2(dot(tangent, point), dot(bitangent, point));
    return max(sign(d), 0.);
}

mat3 envOrientation;

vec3 light(vec3 origin, vec3 rayDir) {
    origin = -(cameraMatrix * vec4(origin, 1)).xyz;
    rayDir = -(cameraMatrix * vec4(rayDir, 0)).xyz;

    origin *= envOrientation;
    rayDir *= envOrientation;

    vec2 uv;
    float hit = intersectPlane(origin, rayDir, vec3(5,-2,-8), normalize(vec3(1,-.5,-.1)), normalize(vec3(0,1,0)), uv);
    float l = smoothstep(.75, .0, fBox(uv, vec2(.4,1.2) * 6.));
	return vec3(l) * hit;
}

vec3 env(vec3 origin, vec3 rayDir) {    
    origin = -(cameraMatrix * vec4(origin, 1)).xyz;
    rayDir = -(cameraMatrix * vec4(rayDir, 0)).xyz;

    origin *= envOrientation;
    rayDir *= envOrientation;

    float l = smoothstep(.0, 1.7, dot(rayDir, vec3(.5,-.3,1))) * .4;
   	return vec3(l) * vec3(1,1,1);
}



//========================================================
// Marching
//========================================================

highp int ZERO = int(min(FRAMECOUNT,0));

// https://iquilezles.org/articles/normalsSDF
vec3 normal( in vec3 pos )
{
    vec3 n = vec3(0.0);
    for( int i=(min(ZERO, 0)); i<4; i++ )
    {
        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
        n += e*map(pos+0.001*e).x;
    }
    return normalize(n);
}

struct Hit {
    vec2 res;
    vec3 p;
    float len;
    float steps;
};

Hit march(vec3 origin, vec3 rayDir, float invert, float maxDist, float understep) {
    vec3 p;
    float len = 0.;
    float dist = 0.;
    vec2 res = vec2(0.);
    vec2 candidate = vec2(0.);
    float steps = 0.;

    for (float i = 0.; i < 60.; i++) {
        len += dist * understep;
        p = origin + len * rayDir;
        candidate = map(p);
        dist = candidate.x * invert;
        steps += 1.;
        if (dist < .001) {
            res = candidate;
            break;
        }
        if (len >= maxDist) {
            len = maxDist;
            break;
        }
    }   

    return Hit(res, p, len, steps);
}

mat3 sphericalMatrix(vec2 tp) {
    float theta = tp.x;
    float phi = tp.y;
    float cx = cos(theta);
    float cy = cos(phi);
    float sx = sin(theta);
    float sy = sin(phi);
    return mat3(
        cy, -sy * -sx, -sy * cx,
        0, cx, sx,
        sy, cy * -sx, cy * cx
    );
}

mat3 calcLookAtMatrix(vec3 ro, vec3 ta, vec3 up) {
    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(ww,up));
    vec3 vv = normalize(cross(uu,ww));
    return mat3(uu, vv, ww);
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    initPoly();

    time = smoothTime / 5.;
    //time = 0.5;
    //time = 0.1;
    //time = 1.;
    time = fract(time + .4);
    
    #ifndef DARK_MODE
        envOrientation = sphericalMatrix(((vec2(81.5, 119) / vec2(187)) * 2. - 1.) * 2.);
    #else
        envOrientation = sphericalMatrix((vec2(0.7299465240641712,0.3048128342245989) * 2. - 1.) * 2.);
    #endif

    vec2 uv = (2. * fragCoord - RENDERSIZE.xy) / RENDERSIZE.y;

    Hit hit, firstHit;
    vec2 res;
    vec3 p, rayDir, origin, sam, ref, raf, nor, camOrigin, camDir;
    float invert, ior, offset, extinctionDist, maxDist, firstLen, bounceCount, wavelength;
    
    vec3 col = vec3(0);
    float focal = 3.8;
    bool refracted;

    vec3 bgCol = vec3(.22);

    invert = 1.;
    maxDist = 15.; 
    
        cameraMatrix *= rotX(((iMouse.y / RENDERSIZE.y) * 2. - 1.) * 2.);
	    cameraMatrix *= rotZ(((iMouse.x / RENDERSIZE.x) * 2. - 1.) * 2.);
    

	camOrigin = -(cameraMatrix[3].xyz) * mat3(cameraMatrix);
   	camDir = normalize(vec3(uv * .168, -1.) * mat3(cameraMatrix));


    firstHit = march(camOrigin, camDir, invert, maxDist, .9);
    firstLen = firstHit.len;

    float steps = 0.;
    
    for (float disperse = 0.; disperse < MAX_DISPERSE; disperse++) {
        invert = 1.;
    	sam = vec3(0);

        origin = camOrigin;
        rayDir = camDir;

        extinctionDist = 0.;
        wavelength = disperse / MAX_DISPERSE;
   		float rand = texture(Image30, fragCoord.xy / RENDERSIZE.xy).r;
        rand = fract(rand + float(FRAMECOUNT) * 1.61803398875);

        wavelength += (rand * 2. - 1.) * (.5 / MAX_DISPERSE);
        
		bounceCount = 0.;

        for (float bounce = 0.; bounce < MAX_BOUNCE; bounce++) {

            if (bounce == 0.) {
                hit = firstHit;
            } else {
                hit = march(origin, rayDir, invert, maxDist / 2., 1.);
            }
            
            steps += hit.steps;
            
            res = hit.res;
            p = hit.p;
            
            if (invert < 0.) {
	            extinctionDist += hit.len;
            }

            // hit background
            if ( res.y == 0.) {
                break;
            }

            vec3 nor = normal(p) * invert;            
            ref = reflect(rayDir, nor);
            
            // shade
            sam += light(p, ref) * .5;
            sam += pow(1. - abs(dot(rayDir, nor)), 5.) * .1;
            sam *= vec3(.85,.85,.98);

            // refract
            float ior = mix(1.3, 1.6, wavelength);
            ior = invert < 0. ? ior : 1. / ior;
            raf = refract(rayDir, nor, ior);
            bool tif = raf == vec3(0); // total internal reflection
            rayDir = tif ? ref : raf;
            offset = .01 / abs(dot(rayDir, nor));
            origin = p + offset * rayDir;
            //invert = tif ? invert : invert * -1.;
            invert *= -1.; // not correct but gives more interesting results

            bounceCount = bounce;
        }

        #ifndef DARK_MODE
            sam += bounceCount == 0. ? bgCol : env(p, rayDir);	
        #endif

        if (bounceCount == 0.) {
            // didn't bounce, so don't bother calculating dispersion
            col += sam * MAX_DISPERSE / 2.;
            break;
        } else {
            vec3 extinction = vec3(.3,.3,1.) * .5;
            extinction = 1. / (1. + (extinction * extinctionDist));	
            col += sam * extinction * spectrum(-wavelength+.2);
        }
	}
    
    // debug
 	//fragColor = vec4(spectrum(steps / 2000.), 1); return;
    //fragColor = vec4(vec3(bounceCount / MAX_BOUNCE), 1); return;
    //fragColor = vec4(vec3(firstHit.steps / 100.), 1); return;

    col /= MAX_DISPERSE;
        
    fragColor = vec4(col, range(4., 12., firstLen));
	return fragColor; 
 } 



/*

    Crystal Tetrahedron
    -------------------

    Brute force dispersion with a fractal zoom.

    The colours come from performing the whole render multiple times,
    each with a different refractive index and colour tint, then
    combining the results. With so few iterations this usually results
    in a lot of banding, but the addition of some animated blue noise
    gives an acceptable result with just 2 samples.

    If it's running slow, drop MAX_DISPERSE in Buffer A to 2.
	It looks really pretty with 30 bounces and 10 dispersion samples.

	Getting good camera and lighting angles was tricky, so there are a
	bunch of hardcoded matrices and rotations from my tool. 

	Check out the high quality capture on my Instagram
	* https://www.instagram.com/p/CA09i10HgDN/

	Inspired by a couple of things:
	* The shading in this GIF by Tatsuya M https://twitter.com/TatsuyaBot/status/1208415211360665607
	* This animation by ayquo https://www.shadertoy.com/view/3d2BRW

*/


// http://tuxedolabs.blogspot.com/2018/05/bokeh-depth-of-field-in-single-pass.html

vec2 uPixelSize; //The size of a pixel: vec2(1.0/width, 1.0/height)
float uFar = 1.; // Far plane

const float GOLDEN_ANGLE = 2.39996323;
const float MAX_BLUR_SIZE = 10.;
const float RAD_SCALE = 1.; // Smaller = nicer blur, larger = faster

float getBlurSize(float depth, float focusPoint, float focusScale) {
    float coc = clamp((1.0 / focusPoint - 1.0 / depth)*focusScale, -1.0, 1.0);
    return abs(coc) * MAX_BLUR_SIZE;
}

vec3 depthOfField(vec2 texCoord, float focusPoint, float focusScale) {
    vec4 centerTex = texture(BuffA, texCoord);
    float centerDepth = centerTex.a * uFar;
    float centerSize = getBlurSize(centerDepth, focusPoint, focusScale);
    vec3 color = centerTex.rgb;
    
    #ifdef DISABLE_DOF
    	return color;
    #endif

    float tot = 1.0;

    float radius = RAD_SCALE;
    for (float ang = 0.; ang < 10000.; ang += GOLDEN_ANGLE) {
        if (radius >= MAX_BLUR_SIZE) break;

        vec2 tc = texCoord + vec2(cos(ang), sin(ang)) * uPixelSize * radius;
        vec4 sampleTex = texture(BuffA, tc);
        vec3 sampleColor = sampleTex.rgb;
        float sampleDepth = sampleTex.a * uFar;
        float sampleSize = getBlurSize(sampleDepth, focusPoint, focusScale);
        //if (sampleSize < centerSize) break;
        if (sampleDepth > centerDepth) {
            sampleSize = clamp(sampleSize, 0.0, centerSize*2.0);
        }
        float m = smoothstep(radius-0.5, radius+0.5, sampleSize);
        color += mix(color/tot, sampleColor, m);
        tot += 1.0;
        radius += RAD_SCALE/radius;
        
        // modification: exit early when we're in focus
       // if (centerDepth < uFar / 3. && m == 0.) break;
    }
    return color /= tot;
}


// http://filmicworlds.com/blog/filmic-tonemapping-operators/
vec3 tonemap2(vec3 texColor) {
    texColor /= 2.;
   	texColor *= 16.;  // Hardcoded Exposure Adjustment
   	vec3 x = max(vec3(0),texColor-0.004);
   	return (x*(6.2*x+.5))/(x*(6.2*x+1.7)+0.06);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    uPixelSize = vec2(.002) / (RENDERSIZE.xy / RENDERSIZE.y);

    //fragColor = vec4(texture2D(BuffA, uv).aaa, 1.); return;
    
    vec3 col = depthOfField(uv, .65, 1.);

    col = pow(col, vec3(1.25)) * 2.5;
    col = tonemap2(col);

    fragColor = vec4(col, 1);    
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}