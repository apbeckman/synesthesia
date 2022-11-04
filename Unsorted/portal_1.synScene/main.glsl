vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********

// Find t so that mix(a,b,t) = x
float unmix(float a, float b, float x) {
    return (x - a)/(b - a);
}

float rand(in vec2 p) {
    return fract(sin(dot(p,vec2(12.9898,78.233))) * 43758.5453);
}

float noise(in vec2 p) {
    vec2 pi = floor(p);
    vec2 pf = fract(p);

    float r00 = rand(vec2(pi.x    ,pi.y    ));
    float r10 = rand(vec2(pi.x+1.0,pi.y    ));
    float r01 = rand(vec2(pi.x    ,pi.y+1.0));
    float r11 = rand(vec2(pi.x+1.0,pi.y+1.0));

    return mix(mix(r00, r10, pf.x), mix(r01, r11, pf.x), pf.y);
}

float fbm(vec2 uv) {
    vec2 p = uv*256.0;
    
    float v = noise(p);

    p *= 0.3;
    v = mix(v, noise(p), 0.8);

    p *= 0.3;
    v = mix(v, noise(p), 0.8);

    p *= 0.3;
    v = mix(v, noise(p), 0.8);

    return v;
}

float nearInt(float x) {
    return pow(0.5 + 0.5*cos(x*6.28),10.0);
}

float fbmWithBorder(vec2 uv) {
    float r = fbm(uv);
    return mix(1.0, fbm(uv), (1.0-nearInt(uv.x))*(1.0-nearInt(uv.y)) );
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    vec2 uv = fragCoord/RENDERSIZE.xy;

    float heightTerrain = fbmWithBorder(uv);
    float cloudIntensity = heightTerrain;
    
    float starVal = noise(uv*1024.0);
    float starIntensity = unmix(0.95, 1.0, starVal) * step(0.95, starVal);

    fragColor = vec4(heightTerrain, cloudIntensity, starIntensity, 0.0);
	return fragColor; 
 } 



#define TWOPI 6.283185307179586

#define IMAGE_ASPECT_WIDTH_OVER_HEIGHT 1.0

#define FIXED_UP vec3(0.0, 1.0, 0.0)
#define TAN_HALF_FOVY 0.7673269879789604
#define CAM_Z_NEAR 0.1
#define CAM_Z_FAR 50.0

#define MIN_DIST 0.005
#define MAX_DIST 50.0
#define RAY_STEPS 40
#define RAY_STEPS_SHADOW 10
#define POM_QUALITY 100
#define POM_QUALITY_REFL 40

#define TEX_SCALE 40.0
#define BUMP_TEX_DEPTH 0.12
#define BOUNDARY_RADIUS 0.2
#define PLANE_DEPTH 3.0
#define NV_PLANE_N vec3(0.0, 1.0, 0.0)

#define TWIST
#define TWIST_EXPONENT 1.0
#define DROSTE

// =============================================
// Some generic helpers
// =============================================

// Find t so that mix(a,b,t) = x
/*
float unmix(float a, float b, float x) {
    return (x - a)/(b - a);
}
*/
float lensq(vec3 p, vec3 q) {
    vec3 pq = q - p;
    return dot(pq, pq);
}

float hitPlane(vec3 planePoint, vec3 nvPlaneN, vec3 p, vec3 v) {
    return dot(planePoint - p, nvPlaneN) / dot(v, nvPlaneN);
}

mat4 getClipToWorld(float aspectWoverH, vec3 nvCamFw) {
    mat4 clipToEye = mat4(
        aspectWoverH * TAN_HALF_FOVY, 0.0, 0.0, 0.0,
        0.0, TAN_HALF_FOVY, 0.0, 0.0,
        0.0, 0.0,  0.0, (CAM_Z_NEAR - CAM_Z_FAR)/(2.0 * CAM_Z_NEAR * CAM_Z_FAR),
        0.0, 0.0, -1.0, (CAM_Z_NEAR + CAM_Z_FAR)/(2.0 * CAM_Z_NEAR * CAM_Z_FAR)
    );

    vec3 nvCamRt = normalize(cross(nvCamFw, FIXED_UP));
    vec3 nvCamUp = cross(nvCamRt, nvCamFw);
    mat4 eyeToWorld = mat4(
         nvCamRt, 0.0,
         nvCamUp, 0.0,
        -nvCamFw, 0.0,
        0.0, 0.0, 0.0, 1.0
    );

    return eyeToWorld * clipToEye;
}

vec3 nvDirFromClip(mat4 clipToWorld, vec2 clip) {
    vec4 world = clipToWorld * vec4(clip, 1.0, 1.0);
    return normalize(world.xyz / world.w);
}

void computeLighting(
    in float diffuseCoefficient,
    in float specularCoefficient,
    in float specularExponent,
    in vec3 lightColor,
    in vec3 texColor,
    in vec3 nvNormal,
    in vec3 nvFragToLight,
    in vec3 nvFragToCam,
    out vec3 diffuse,
    out vec3 specular
) {
    float valDiffuse = max(0.0, dot(nvNormal, nvFragToLight));
    valDiffuse *= diffuseCoefficient;

    vec3 blinnH = normalize(nvFragToLight + nvFragToCam);
    float valSpecular = pow(max(0.0, dot(nvNormal, blinnH)), specularExponent);
    valSpecular *= specularCoefficient;

    diffuse = valDiffuse * texColor * lightColor;
    specular = valSpecular * lightColor;
}

// =============================================
// "Portal object" rendering
// =============================================

vec3 cpPlane(vec3 planePoint, vec3 nvPlaneN, vec3 p) {
    float t = dot(p - planePoint, nvPlaneN);
    return p - t*nvPlaneN;
}

vec3 cpSeg(vec3 q0, vec3 q1, vec3 p) {
    vec3 vEdge = q1 - q0;
    float t = dot(p - q0, vEdge) / dot(vEdge, vEdge);
    return q0 + clamp(t, 0.0, 1.0)*vEdge;
}

vec3 cpTuple2(vec3 q0, vec3 q1, vec3 p) {
    vec3 q = q0;
    return mix(q, q1, step( lensq(p,q1), lensq(p,q) ));
}

vec3 cpTuple3(vec3 q0, vec3 q1, vec3 q2, vec3 p) {
    vec3 q = cpTuple2(q0,q1, p);
    return mix(q, q2, step( lensq(p,q2), lensq(p,q) ));
}

vec3 cpTuple4(vec3 q0, vec3 q1, vec3 q2, vec3 q3, vec3 p) {
    vec3 q = cpTuple3(q0,q1,q2, p);
    return mix(q, q3, step( lensq(p,q3), lensq(p,q) ));
}

vec3 cpTriBoundary(vec3 q0, vec3 q1, vec3 q2, vec3 p) {
    return cpTuple3(cpSeg(q0,q1, p), cpSeg(q1,q2, p), cpSeg(q2,q0, p), p);
}

vec3 cpQuadBoundary(vec3 q0, vec3 q1, vec3 q2, vec3 q3, vec3 p) {
    return cpTuple4(
        cpSeg(q0,q1, p), cpSeg(q1,q2, p),
        cpSeg(q2,q3, p), cpSeg(q3,q0, p), p
    );
}

float pointInTri(vec3 q0, vec3 q1, vec3 q2, vec3 p) {
    vec3 v01 = cross(q1-q0, p-q0);
    vec3 v12 = cross(q2-q1, p-q1);
    vec3 v20 = cross(q0-q2, p-q2);
    return step(0.0, dot(v01,v12)) * step(0.0, dot(v01,v20));
}

vec3 cpTri(vec3 q0, vec3 q1, vec3 q2, vec3 p) {
    vec3 nvPlaneN = normalize(cross(q1-q0, q2-q0));
    vec3 xp = cpPlane(q0, nvPlaneN, p);
    return mix(cpTriBoundary(q0,q1,q2, p), xp, pointInTri(q0,q1,q2, xp));
}

vec3 cpQuad(vec3 q0, vec3 q1, vec3 q2, vec3 q3, vec3 p) {
    return cpTuple2(cpTri(q0,q1,q2, p), cpTri(q0,q2,q3, p), p);
}

float sdQuadBoundary(vec3 q0, vec3 q1, vec3 q2, vec3 q3, float r, vec3 p) {
    vec3 x = cpQuadBoundary(q0,q1,q2,q3, p);
    return distance(x, p) - r;
}

vec3 normalQuadBoundary(vec3 q0, vec3 q1, vec3 q2, vec3 q3, vec3 p) {
    vec3 x = cpQuadBoundary(q0,q1,q2,q3, p);
    return normalize(p - x);
}

float sdQuad(vec3 q0, vec3 q1, vec3 q2, vec3 q3, float r, vec3 p) {
    vec3 x = cpQuad(q0,q1,q2,q3, p);
    return distance(x, p) - r;
}

void hitObject(
    in vec3 startPos, in vec3 nvRayDir,
    in vec3 q00, in vec3 q10, in vec3 q11, in vec3 q01,
    out float didHit, out vec3 hitPos
){
    didHit = 0.0;
    float travel = 0.0;
    vec3 curPos = startPos;

    for (int k = 0; k < RAY_STEPS; k++) {
        float sdCur = sdQuadBoundary(q00,q10,q11,q01, BOUNDARY_RADIUS, curPos);

        if (sdCur < MIN_DIST) {
            didHit = 1.0;
            break;
        }

        curPos += sdCur * nvRayDir;
        travel += sdCur;
        if (travel > MAX_DIST) {
            break;
        }
    }

    hitPos = curPos;
}

void hitShadow(
    in vec3 startPos, in vec3 nvRayDir,
    in vec3 q00, in vec3 q10, in vec3 q11, in vec3 q01,
    out float lightPercent
){
    lightPercent = 1.0;
    float travel = 0.0;
    vec3 curPos = startPos;

    for (int k = 0; k < RAY_STEPS_SHADOW; k++) {
        float sdCur = sdQuad(q00,q10,q11,q01, BOUNDARY_RADIUS, curPos);

        float curLightPercent = abs(sdCur)/(0.02*travel);
        lightPercent = min(lightPercent, curLightPercent);

        if (sdCur < MIN_DIST) {
            lightPercent = 0.0;
            break;
        }

        curPos += sdCur * nvRayDir;
        travel += sdCur;
        if (travel > MAX_DIST) {
            break;
        }
    }
}

// =============================================
// Terrain and sky rendering
// =============================================

vec4 bumpTex(vec2 uv) {
    float r = textureLod(BuffA, uv, 0.0).r;

    float ang = (r + fract(uv.x) + fract(uv.y)) * TWOPI;
    vec2 q = uv + 0.1*vec2( cos(ang), sin(ang) );

    float height = textureLod(BuffA, q, 0.0).r;
    vec3 color = mix(vec3(0.0), vec3(1.0, 0.2, 0.1), height);

    return vec4(color, height);
}

vec2 skyTex(vec2 uv) {
    vec4 data = textureLod(BuffA, uv, 0.0);

    float starIntensity = data.b;
    
    float ang = (data.g + smoothTimeC*0.1) * TWOPI;
    vec2 q = uv + 0.05*vec2( cos(ang), sin(ang) );
    float cloudIntensity = textureLod(BuffA, q, 0.0).g;
    
    return vec2(cloudIntensity, starIntensity);
}

vec3 bumpTexNormal(vec2 uv) {
    vec2 uvPixel = 1.0 / RENDERSIZE.xy;
    float hSA = bumpTex(uv + vec2(-uvPixel.s, 0.0)).a;
    float hSB = bumpTex(uv + vec2( uvPixel.s, 0.0)).a;
    float hTA = bumpTex(uv + vec2(0.0,-uvPixel.t)).a;
    float hTB = bumpTex(uv + vec2(0.0, uvPixel.t)).a;

    vec2 dhdt = vec2(hSB-hSA, hTB-hTA) / (2.0 * uvPixel);
    vec2 gradh = BUMP_TEX_DEPTH * dhdt;
    return normalize(vec3( -gradh, 1.0 ));
}

vec3 skyColor(vec3 nvDir) {
    float yy = clamp(nvDir.y+0.1, 0.0, 1.0);
    float horiz0 = pow(1.0 - yy, 30.0);
    float horiz1 = pow(1.0 - yy, 5.0);
    
    vec3 sv = nvDir - vec3(0.0, -1.0, 0.0);
    vec2 uvCloud = 0.25*(sv.xz / sv.y) + vec2(0.5);
    vec2 skyTexVal = skyTex(uvCloud);

    float cloudIntensity = pow(skyTexVal.x, 2.0);
    float starIntensity = pow(skyTexVal.y, 2.0);

    vec3 c = vec3(0.0);
    c = mix(c, vec3(0.2, 0.0, 0.5), horiz1);
    c = mix(c, vec3(1.0), horiz0);
    c = mix(c, vec3(0.45, 0.5, 0.48), (1.0-horiz0)*cloudIntensity);
    c = mix(c, vec3(1.0), (1.0-horiz1)*starIntensity);
    return c;
}

void getParallaxMaxOffsets(
    in vec3 tangentS,
    in vec3 tangentT,
    in vec3 nvNormal,
    in vec3 camToFrag,
    in float depthMax,
    out vec2 maxTexOffset,
    out vec3 maxPosOffset
){
    // Scale camToFrag so that its depth into the surface == depthMax
    float camDist = -dot(camToFrag, nvNormal);
    maxPosOffset = (depthMax / camDist) * camToFrag;

    // maxPosOffset = a*tangentS + b*tangentT + N <=> maxTexOffset = [a,b]
    float dss = dot(tangentS, tangentS);
    float dst = dot(tangentS, tangentT);
    float dtt = dot(tangentT, tangentT);
    float dcs = dot(maxPosOffset, tangentS);
    float dct = dot(maxPosOffset, tangentT);
    float invDet = 1.0 / (dss * dtt - dst * dst);
    maxTexOffset = invDet * vec2(dtt*dcs - dst*dct, -dst*dcs + dss*dct);
}

float getParallaxDepthFactor(vec2 uvInitial, vec2 maxTexOffset, int steps) {
    vec2 uvMax = uvInitial + maxTexOffset;
    float dt = 1.0 / float(steps);

    float tOld = 0.0, depthOld = 0.0;
    float tCur = 0.0, depthCur = 0.0;

    for(int i=0; i<=steps; ++i){
        tOld = tCur;
        tCur = float(i)*dt;

        depthOld = depthCur;
        depthCur = 1.0 - bumpTex(mix(uvInitial, uvMax, tCur)).a;

        if(tCur > depthCur){
            tCur = mix(tOld, tCur, unmix(depthOld-tOld, depthCur-tCur, 0.0));
            break;
        }
    }

    return tCur;
}

void terrainAndSky(
    in vec3 startPos, in vec3 nvRayDir, in vec3 lightPos, in int pomSteps,
    in vec3 q00, in vec3 q10, in vec3 q11, in vec3 q01,
    out vec3 hitColor
) {
    float tPlane = hitPlane(
        vec3(0.0, -PLANE_DEPTH, 0.0), NV_PLANE_N,
        startPos, nvRayDir
    );
    float didHitPlane = step(0.0, tPlane);

    if (didHitPlane > 0.5) {

        vec3 hitPos = startPos + tPlane*nvRayDir;
        vec2 hitTex = hitPos.xz / TEX_SCALE;

        vec2 maxTexOffset;
        vec3 maxPosOffset;
        getParallaxMaxOffsets(
            vec3(TEX_SCALE, 0.0, 0.0),
            vec3(0.0, 0.0, TEX_SCALE),
            vec3(0.0, 1.0, 0.0),
            hitPos - startPos,
            BUMP_TEX_DEPTH * TEX_SCALE,
            maxTexOffset,
            maxPosOffset
        );
        float depthPct = getParallaxDepthFactor(hitTex, maxTexOffset, pomSteps);
        vec2 hitTexBump = hitTex + depthPct*maxTexOffset;
        vec3 hitPosBump = hitPos + depthPct*maxPosOffset;

        vec3 bumpColor = bumpTex(hitTexBump).rgb;
        vec3 nvNormal = bumpTexNormal(hitTexBump);
        vec3 nvBumpNormal = normalize(vec3(nvNormal.x, 1.0, nvNormal.y));

        if (distance(startPos, hitPosBump) < MAX_DIST) {

            // Hit terrain: Compute fog, lighting, and shadow

            vec3 nvBumpPosToLight = normalize(lightPos - hitPosBump);
            vec3 vBumpPosToStart = startPos - hitPosBump;
            float dHit = length(vBumpPosToStart);
            vec3 nvBumpPosToStart = vBumpPosToStart / dHit;

            float lightPercent = 1.0;
            hitShadow(
                hitPosBump, nvBumpPosToLight,
                q00,q10,q11,q01,
                lightPercent
            );
            lightPercent *= mix(1.0, 0.0, depthPct);

            vec3 diffuse;
            vec3 specular;
            computeLighting(
                0.8, 0.3, 5.0,
                vec3(1.0), bumpColor,
                nvBumpNormal, nvBumpPosToLight, nvBumpPosToStart,
                diffuse, specular
            );
            vec3 ambient = 0.1 * bumpColor;
            vec3 matColor = ambient + lightPercent*(diffuse + specular);

            float fogMin = 0.0;
            float fogMax = 0.8;
            vec3 fogColor = vec3(0.9, 0.8, 1.0);
            float fogT = mix(fogMin, fogMax, unmix(0.0, MAX_DIST, dHit) );
            float fogAmount = pow(fogT, 1.5);

            hitColor = mix(matColor, fogColor, fogAmount);

        } else {

            // Too far, use skybox
            hitColor = skyColor(nvRayDir);

        }

    } else {

        // Didn't hit plane, use skybox
        hitColor = skyColor(nvRayDir);
    }
}

// =============================================
// Overall scene
// =============================================

float subrectSize() {
    float isMousePressed = clamp(iMouse.z, 0.0, 1.0);
    float valDefault = 0.25;
    float valPressed = mix(0.04, 0.6, iMouse.x/RENDERSIZE.x);
    return mix(valDefault, valPressed, isMousePressed);
}

vec4 scene(vec2 p) {

    // Camera
    // ---------------------------
    vec3 camPos = 4.0 * vec3(cos(smoothTimeC*0.2), 0.0, sin(smoothTimeC*0.2));
    camPos += vec3(0.0, 0.75 + 0.5*cos(smoothTimeC*0.5), 0.0);
    vec3 lookTarget = vec3(0.0);

    vec3 movement = vec3(2.0, 0.0, -smoothTimeC*2.0);
    camPos += movement;
    lookTarget += movement;

    vec3 nvCamFw = normalize(lookTarget - camPos);

    mat4 clipToWorld = getClipToWorld(IMAGE_ASPECT_WIDTH_OVER_HEIGHT, nvCamFw);
    vec3 nvCamDir = nvDirFromClip(clipToWorld, p);

    // Portal geometry
    // ---------------------------
    float a = subrectSize();

    vec3 nv00 = nvDirFromClip(clipToWorld, vec2(-a,-a));
    vec3 nv10 = nvDirFromClip(clipToWorld, vec2( a,-a));
    vec3 nv01 = nvDirFromClip(clipToWorld, vec2(-a, a));
    vec3 nv11 = nvDirFromClip(clipToWorld, vec2( a, a));

    float minY = -PLANE_DEPTH + 2.5*BOUNDARY_RADIUS;
    float tL = hitPlane(vec3(0.0,  minY, 0.0),  NV_PLANE_N, camPos, nv00);
    float tU = hitPlane(vec3(0.0, -minY, 0.0), -NV_PLANE_N, camPos, nv01);
    float tPortal = min(12.0, min(
        mix(MAX_DIST, tL, step(0.0, tL)),
        mix(MAX_DIST, tU, step(0.0, tU))
    ));

    vec3 q00 = camPos + tPortal*nv00;
    vec3 q10 = camPos + tPortal*nv10;
    vec3 q01 = camPos + tPortal*nv01;
    vec3 q11 = camPos + tPortal*nv11;

    vec3 portalVX = normalize(q10 - q00);
    vec3 portalVY = normalize(q01 - q00);
    vec3 portalVZ = cross(portalVX, portalVY);

    q00 += BOUNDARY_RADIUS * ( -portalVX - portalVY );
    q10 += BOUNDARY_RADIUS * (  portalVX - portalVY );
    q01 += BOUNDARY_RADIUS * ( -portalVX + portalVY );
    q11 += BOUNDARY_RADIUS * (  portalVX + portalVY );

    // Light placement
    // ---------------------------
    vec3 lightPos = 0.5*(q01 + q11) + 1.0*portalVY + 5.0*portalVZ;

    // Render scene
    // ---------------------------
    vec3 sceneColor = vec3(0.0);

    float didHitPortal;
    vec3 hitPos;
    hitObject(camPos, nvCamDir, q00,q10,q11,q01, didHitPortal, hitPos);

    if (didHitPortal > 0.5) {

        vec3 n = normalQuadBoundary(q00,q10,q11,q01, hitPos);
        vec3 nvRefl = normalize(reflect( hitPos-camPos, n ));

        vec3 diffuse;
        vec3 specular;
        computeLighting(
            0.2, 0.8, 20.0,
            vec3(1.0), 0.4*vec3(1.0, 0.5, 1.0),
            n, normalize(lightPos - hitPos), normalize(camPos - hitPos),
            diffuse, specular
        );
        vec3 matColor = diffuse + specular;

        vec3 terrainColor;
        terrainAndSky(
            hitPos, nvRefl, lightPos, POM_QUALITY_REFL,
            q00,q10,q11,q01,
            terrainColor
        );

        sceneColor = matColor + 0.8*terrainColor;

    } else {

        terrainAndSky(
            camPos, nvCamDir, lightPos, POM_QUALITY,
            q00,q10,q11,q01,
            sceneColor
        );

    }

    return vec4(clamp(sceneColor, 0.0, 1.0), 1.0);
}

// ------------------------------------------
// Helpers to achieve "Escher effect"
// ------------------------------------------

vec2 cmul(vec2 a, vec2 b) {
    return vec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

vec2 clog(vec2 a) {
    return vec2(0.5*log(dot(a,a)), atan(a.y, a.x));
}

vec2 cexp(vec2 a) {
    return exp(a.x)*vec2(cos(a.y), sin(a.y));
}

vec2 twist(vec2 p) {
#ifdef TWIST
    float a = subrectSize();

    vec2 r = vec2(TWIST_EXPONENT, log(a)/TWOPI);
    return cexp(cmul(r, clog(p)));
#else
    return p;
#endif
}

vec2 droste(vec2 p) {
#ifdef DROSTE
    float a = subrectSize();

    float isMousePressed = clamp(iMouse.z, 0.0, 1.0);
    float apow = mix(
        fract(smoothTime*0.2), 4.0*iMouse.y/RENDERSIZE.y, isMousePressed
    );
    p *= pow(a, apow);

    vec2 log_a = log(abs(p)) / log(a);
    float adjust = min(floor(log_a.x), floor(log_a.y));
    log_a -= vec2(adjust);

    return sign(p) * pow(vec2(a), log_a);
#else
    return p;
#endif
}

// ------------------------------------------
// Put it all together!
// ------------------------------------------

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy;

    vec2 p = (2.0*uv - 1.0);
    p.x *= (RENDERSIZE.x / RENDERSIZE.y) / IMAGE_ASPECT_WIDTH_OVER_HEIGHT;

    vec2 radv = uv - vec2(0.5, 0.5);
    float dCorner = length(radv);
    float vignetteFactor = 1.0 - mix(0.0, 0.3, smoothstep(0.2, 0.707, dCorner));

    fragColor = vignetteFactor * scene(droste(twist(p)));
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