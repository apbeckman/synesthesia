vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********


#define PI 3.14159265359
#define rot(a) mat2(cos(a + PI*0.5*vec4(0,1,3,0)))

#define SAMPLES 6

const float motionBlurAmount = 1.50;
const float nearZ = 2.0;

#define warm
//#define complex

#ifdef warm
const vec3 sunCol = vec3(0.9, 0.55, 0.1)*9.0;
const vec3 subCol = vec3(0.5, 0.8, 0.89)*0.5;
const float sunRadius = 1.50;
const float sunExponent = -40.0;
const float timeScale = 3.0;
#else
const vec3 sunCol = vec3(0.6, 0.2, 0.85)*9.0;
const vec3 subCol = vec3(0.5, 0.8, 0.9)*0.5;
const float sunRadius = 1.5;
const float sunExponent = -30.0;
const float timeScale = 2.0;
#endif
float smoothTime = (smooth_basstime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.75;
float smoothTimeB = (smooth_hightime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.95;
float smoothTimeC = (smooth_midtime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.95;

float time;

vec3 hash33(vec3 p3) {
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy + p3.yxx)*p3.zyx);
}

vec3 hash( in vec3 center ) {
    return center + (hash33(center)-0.5) * 0.5;
}

float B2( vec2 _P ) {
    return mod( 2.0*_P.y + _P.x + 1.0, 4.0 );
}

float B8( vec2 _P ) {
    vec2	P1 = mod( _P, 2.0 );					// (P >> 0) & 1
    vec2	P2 = floor( 0.5 * mod( _P, 4.0 ) );		// (P >> 1) & 1
    vec2	P4 = floor( 0.25 * mod( _P, 8.0 ) );	// (P >> 2) & 1
    return 4.0*(4.0*B2(P1) + B2(P2)) + B2(P4);
}

vec3 voronoi( in vec3 uv, in vec3 no, inout float rough ) {
    
    vec3 center = floor(uv) + 0.5;
    vec3 bestCenterOffset = vec3(0);
    float bestDist = 9e9;
    vec3 bestCenterOffset2 = vec3(0);
    float bestDist2 = 9e9;
    
    for (float x = -0.5 ; x < 1.0 ; x+=1.0)
    for (float y = -0.5 ; y < 1.0 ; y+=1.0)
    for (float z = -0.5 ; z < 1.0 ; z+=1.0) {
		vec3 offset = vec3(x, y, z);
        vec3 newCenter = center + offset;
        vec3 newCenterOffset = hash(newCenter);
        vec3 temp = newCenterOffset - uv;
        float distSq = dot(temp, temp);
        if (distSq < bestDist) {
    		bestCenterOffset2 = bestCenterOffset;
    		bestDist2 = bestDist;
            bestCenterOffset = newCenterOffset;
            bestDist = distSq;
        } else if (distSq < bestDist2) {
            bestCenterOffset2 = newCenterOffset;
            bestDist2 = distSq;
        }
    }
    
    vec3 n1 = normalize(no + hash33(bestCenterOffset)-0.5);
    vec3 n2 = normalize(no + hash33(bestCenterOffset2)-0.5);
    float d = (sqrt(bestDist)-sqrt(bestDist2));
    float aad = 0.02;
    return mix(n1, n2, smoothstep(-aad, +aad, d*2.0));
}

float smin( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); 
}

float cross4D(vec4 p) {
    float aa = length(p.xy)-0.1;
    float bb = length(p.zw)-0.02;
    float cc = length(p.yw)-0.03;
    float de = smin(smin(aa, bb, 0.1), cc, 0.1)-0.03;
    return de;
}

float de4d(vec4 p) {
    p.zw *= rot(smoothTime*0.1);
    p.xz *= rot((0.2*smoothTime)*0.26241);
    p.wy *= rot(smoothTimeB*0.187);
    p.x += smoothTime*0.2687;
   	p.y -= (0.125*smoothTime);
    vec4 inG = (fract(p) - 0.5);
    return cross4D(inG);
}

float de(vec3 p) {
    #ifdef complex
    float d = dot(p, p);
    d = de4d(vec4(p/d, d))*d;
    #else
    vec4 pp = vec4(p, length(p));
    float d = dot(pp, pp);
    d = de4d(pp/d)*d;
    #endif
    return d;
}

vec3 getNormal(vec3 p) {
	vec3 e = vec3(0.0, 0.0001, 0.0);
	return normalize(vec3(
		de(p+e.yxx)-de(p-e.yxx),
		de(p+e.xyx)-de(p-e.xyx),
		de(p+e.xxy)-de(p-e.xxy)));	
}

float DistributionGGX(vec3 N, vec3 H, float roughness) {
    float a = roughness*roughness;
    float a2 = a*a;
    float NdotH = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;
    float num = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;
    return num / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness) {
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;
    float num = NdotV;
    float denom = NdotV * (1.0 - k) + k;
    return num / denom;
}

float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness) {
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2 = GeometrySchlickGGX(NdotV, roughness);
    float ggx1 = GeometrySchlickGGX(NdotL, roughness);
    return ggx1 * ggx2;
}

vec3 fresnelSchlick(float cosTheta, vec3 F0) {
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

vec3 computeLighting(in vec3 normal, in vec3 viewDir,
                     in vec3 albedo, in float metallic, in float roughness,
                     in vec3 lightDir, in vec3 radiance) {
    vec3 result = vec3(0);
    vec3 halfwayDir = normalize(viewDir + lightDir);
    vec3 F0 = vec3(0.004);
    F0 = mix(F0, albedo, metallic);
    float NDF = DistributionGGX(normal, halfwayDir, roughness);
    float G = GeometrySmith(normal, viewDir, lightDir, roughness);
    vec3 F = fresnelSchlick(max(dot(halfwayDir, viewDir), 0.0), F0);
    vec3 kS = F;
    vec3 kD = 1.0 - kS;
    kD *= 1.0 - metallic;
    vec3 numerator = NDF * G * F;
    float denominator = 4.0 * max(dot(normal, viewDir), 0.0) * max(dot(normal, lightDir), 0.0);
    vec3 specular = numerator / max(denominator, 0.0001);
    float NdotL = max(dot(normal, lightDir), 0.0);
    result += (kD * albedo / PI + specular) * radiance * NdotL;
    return result;
}

// fake subsurface scattering
vec3 computeSSS(in vec3 normal, in vec3 viewDir, 
                in vec3 albedo, in float trans, in float index,
                in vec3 lightDir, in vec3 radiance) {
    float add = 1.0 - index;
    add *= add;
    add *= add;
    add *= add;
    add *= add;
    float fr = dot(viewDir, normal)*0.5+0.5;
    float lu = dot(viewDir, lightDir)*-0.5+0.5;
    add *= fr*fr;
    add *= lu;
    return radiance*add*0.750*trans*albedo;
}

// approximation of the error function
float erf( in float x ) {
    //return tanh(1.202760580 * x);
	float sign_x = sign(x);
	float t = 1.0/(1.0 + 0.47047*abs(x));
	float result = 1.0 - t*(0.3480242 + t*(-0.0958798 + t*0.7478556))*exp(-(x*x));
	return result * sign_x;
}

float getIntegral(vec3 start, vec3 dir, float dist) {
    const float a = sunExponent;
	const float b = sunRadius;
    float k = start.x;
    float l = dir.x;
    float m = start.y;
    float n = dir.y;
    float o = start.z;
    float p = dir.z;
    float res = sqrt(PI);
    res *= exp(b+a*(+k*k*(n*n+p*p)
                    -m*m*(-1.0+n*n)
                    -o*o*(-1.0+p*p)
                    -2.0*k*l*o*p
                    -2.0*m*n*(k*l+o*p) ));
    res *= - erf(sqrt(-a)*dot(start, dir)) + erf(sqrt(-a)*(dot(start, dir)+dist));
    res /= 2.0 * sqrt(-a);
    res *= 500.0;
    return res;
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    vec3 acc = vec3(0);
    
    for (int i = 0 ; i < SAMPLES ; i++) {
    	
        float bay = B8(floor(fragCoord))/64.0;
        
        #if SAMPLES > 1
        vec3 rnd = hash33(vec3(fragCoord, FRAMECOUNT*SAMPLES+i));
        bay = rnd.z;
        #endif
        
        time = (smoothTimeC*0.125 + (1.0/120.0)*bay*motionBlurAmount)*timeScale;

        vec2 uv = (fragCoord - RENDERSIZE.xy * 0.5) / RENDERSIZE.y;
		#if SAMPLES > 1
        uv += (rnd.xy-0.5) / RENDERSIZE.y;
        #endif
        
        vec3 from = vec3(0, 0, -1.4);
        vec3 dir = normalize(vec3(uv*0.5, 1.0));
        
        mat2 rotxz = rot(smoothTime*0.125);
		mat2 rotxy = rot(sin(time*0.1));
        if (iMouse.z > 0.5) {
            vec2 delt = iMouse.xy-iMouse.zw;
            rotxz *= rot(-delt.x*0.01);
            rotxy *= rot(delt.y*0.01);
        }
        
        from.zy *= rotxy;
        from.xz *= rotxz;
        dir.zy  *= rotxy;
        dir.xz  *= rotxz;

        float totdist = 0.0;
        totdist += (0.1+pow(bay, 0.2)*0.4)*nearZ;
        totdist += de(from+dir*totdist)*bay;

        float ao = 0.0;
        for (int steps = 0 ; steps < 100 ; steps++) {
            vec3 p = from + totdist * dir;
            float dist = de(p);
            totdist += dist*(0.65+bay*0.1);
            if (dist < 0.0001 || length(p) > 1.4) {
                ao = float(steps)/149.0;
                break;
            }
        }

        vec3 result = vec3(0);
        vec3 p = from + totdist * dir;

        if (length(p) < 1.4) {

            vec3 n = -getNormal(p);

            vec3 sunDir = normalize(p);
            const vec3 subDir = normalize(vec3(2, -7, 3));

            float rough = 0.0;
            vec3 vor = voronoi(p*600.0, n, rough);
            n = normalize(n+vor*0.5);

            vec4 albedo = vec4(0.98, 0.8, 0.5, 0.4);

            result += computeLighting(n, dir, albedo.rgb, 0.9, albedo.a, sunDir, sunCol);
            result += computeLighting(n, dir, albedo.rgb, 0.9, albedo.a, subDir, subCol);
            result += computeSSS(n, dir, albedo.rgb, albedo.a, ao, sunDir, sunCol);
            result += computeSSS(n, dir, albedo.rgb, albedo.a, ao, subDir, subCol);

        } else {

            // background
            float rough = 0.0;
            vec3 vor = voronoi(dir*800.0, vec3(0), rough);
            result = vec3(pow(abs(vor.x), 500.0))*3.0;

        }

        fragColor.rgb = result;
        float sun = getIntegral(from, dir, totdist);
        fragColor.rgb += sun*0.005*sunCol;
        
        acc += fragColor.rgb;
        
    }
    
    acc /= float(SAMPLES);
    fragColor.rgb = clamp(acc, vec3(0), vec3(10));
    

	return fragColor; 
 } 



//#define PI 3.14159265359
#define PHI 1.61803398875

//#define SAMPLES 5
#define BLOOM_RADIUS 35.0

const mat3 ACESInputMat = mat3(
    0.59719, 0.35458, 0.04823,
    0.07600, 0.90834, 0.01566,
    0.02840, 0.13383, 0.83777
);

const mat3 ACESOutputMat = mat3(
     1.60475, -0.53108, -0.07367,
    -0.10208,  1.10813, -0.00605,
    -0.00327, -0.07276,  1.07602
);

vec3 RRTAndODTFit( in vec3 v ) {
    vec3 a = v * (v + 0.0245786) - 0.000090537;
    vec3 b = v * (0.983729 * v + 0.4329510) + 0.238081;
    return a / b;
}

vec3 ACESFitted( in vec3 color ) {
    color = color * ACESInputMat;
    color = RRTAndODTFit(color);
    color = color * ACESOutputMat;
    color = clamp(color, 0.0, 1.0);
    return color;
}

/*vec3 hash33(vec3 p3) {
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy + p3.yxx)*p3.zyx);
}
*/
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    vec3 bloom = vec3(0);
    float totfac = 0.0;
    
    vec3 rnd = hash33(vec3(fragCoord, FRAMECOUNT));
    float offset = rnd.x*2.0*PI;
    
    // bloom
    for (int i = 0 ; i < SAMPLES ; i++) {
        float theta = 2.0*PI*PHI*float(i) + offset;
        float radius = sqrt(float(i)) / sqrt(float(SAMPLES));
        radius *= BLOOM_RADIUS + pow(syn_Level,2.0);
        vec2 offset = vec2(cos(theta), sin(theta))*radius;
        vec2 delta = vec2( 1.0+exp(-abs(offset.y)*0.1) , 0.5);
        offset *= delta;
        vec4 here = textureGrad(BuffA,(fragCoord+offset)/RENDERSIZE.xy, 
                                vec2(0.001, 0)*BLOOM_RADIUS, vec2(0, 0.001)*BLOOM_RADIUS);
        float fact = smoothstep(BLOOM_RADIUS, 0.0, radius);
        bloom += here.rgb*0.05*fact;
        totfac += fact;
    }
    
    bloom /= totfac;
    
    vec2 uv = fragCoord/RENDERSIZE.xy;
    vec2 mo = uv*2.0-1.0;
    mo *= 0.01;
    fragColor.r = textureLod(BuffA, uv-mo*0.1, 0.0).r;
    fragColor.g = textureLod(BuffA, uv-mo*0.6, 0.0).g;
    fragColor.b = textureLod(BuffA, uv-mo*1.0, 0.0).b;
    
    fragColor.rgb += bloom*bloom*100.0;
    vec2 vi = fragCoord / RENDERSIZE.xy * 2.0 - 1.0;
    fragColor.rgb *= (1.0-sqrt(dot(vi,vi)*0.45));
    fragColor.rgb = ACESFitted(fragColor.rgb);
    fragColor.rgb = pow( fragColor.rgb, vec3(1.0/2.2) );
    fragColor.rgb += (rnd-0.5)*0.1;
    
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