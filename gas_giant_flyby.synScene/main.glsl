

			//******** BuffA Code Begins ********


#define PI 3.14159265359
#define rot(a) mat2(cos(a+PI*vec4(0,1.5,0.5,0)))

#define MARCH_EPS 0.1
#define MARCH_MIN_STEP 0.02
#define MARCH_STEPS 300
#define FOV (PI*0.7)

#define BACKGROUND vec3(0.25, 0.37, 0.4)
#define COLOR_SUN vec3(0.8, 0.7, 0.5)
const vec3 sunDir = normalize(vec3(-1, 6, 3));
float time = 0.0;
float yCamera = 0.0;

// Dave_Hoskins hash
#define HASHSCALE1 .1031
#define HASHSCALE3 vec3(.1031, .1030, .0973)
float hash11( in float p ) {
	vec3 p3 = fract(vec3(p) * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}
vec3 hash33( in vec3 p3 ){
	p3 = fract(p3 * HASHSCALE3);
    p3 += dot(p3, p3.yxz+19.19);
    return fract((p3.xxy + p3.yxx)*p3.zyx);
}

// raytrace a plane, used for the rings
float raytracePlane(in vec3 from, in vec3 dir,
                    in vec3 planePos, in vec3 planeNormal ) {
    float originDist = -dot(planePos, planeNormal);
    float t = -(dot(planeNormal, from)+originDist)/dot(planeNormal, dir);
    if ( t > 0.0 )
        return t;
    return -1.0;
}

// iq's smin
float smin( in float d1, in float d2, in float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h);
}

// iq's smax
float smax( in float d1, in float d2, in float k ) {
    float h = clamp( 0.5 - 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) + k*h*(1.0-h);
}

// smooth abs from above
/*float smoothAbs( in float x, in float k ) {
    float h = clamp( 0.5 + x/k, 0.0, 1.0 );
    return mix( -x, x, h ) + k*h*(1.0-h);
}*/

// return center of hex tile
/*vec2 hexCenter( in vec2 p ) {
    #define HEX vec2(1, 1.73205080757)
    vec2 centerA = (floor(p.xy*HEX)+0.5)/HEX;
    vec2 centerB = (floor((p.xy+HEX*0.5)*HEX)+0.5)/HEX-HEX*0.5;
    vec2 a = p.xy-centerA.xy; vec2 b = p.xy-centerB.xy;
    return dot(a,a)<dot(b,b) ? centerA : centerB;
}*/

// noise with smooth derivative
vec4 snoise( in vec3 x, const in float lod ) {
    float dim = 32.0 / exp2(lod);
    x = x * dim;
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
    x = (p+f+0.5) / dim;
    return textureLod(iChannel0, x, lod);
}

// filtered noise
vec4 noise( in vec3 x, in float freq, in float r ) {
    x *= (1.0/32.0) * freq;
    float lod = log2(r*freq+1.0);
    return textureLod(iChannel0, x, lod);
}

// fog color
vec3 getFog( in vec3 dir ) {
    float d = dot(dir, sunDir) * 0.5 + 0.5;
    d *= d; d *= d; d *= d;
    return mix(BACKGROUND, COLOR_SUN, d);
}

// iq's analytical fog
float getFog( in vec3 rgb, in float dist, in vec3 rayOri, in vec3 rayDir ) {
    const float c = 0.3;
    const float b = 1.0;
    float fogAmount = c * exp(-rayOri.z*b) * (1.0-exp( -dist*rayDir.z*b ))/rayDir.z;
    fogAmount = 1.0 - exp(-fogAmount*1.0);
    return min(fogAmount, 1.0);
}

// wrap tunnel around camera
vec3 cameraPath( in float y ) {
    return vec3(12.75 - sin(y*0.1312 - time*0.2)*3.0, y, 1.0 - sin(y*0.1)*1.0);
}

// return thunder for this location
vec4 getThunder( in float y ) {
    float cell = floor(y*0.1)+0.5;
    float rnd = hash11(cell);
    float per = (rnd*4.0+1.5);
    float frac = fract((time + rnd)/per);
    float ti = floor((time + rnd)/per);
    float alpha = smoothstep(0.0, 0.03, frac);
    alpha *= smoothstep(0.2, 0.05, frac);
    // return thunder strength in w
    vec3 pos = vec3(cameraPath(cell/0.1).xy, -3.0 + hash11(ti)*1.0);
    return vec4(pos, alpha);
}

// sample the procedural field inside a sphere
float map( in vec3 p, in float r, out vec4 color ) {
    
    // curve the world a tiny bit around the camera
    float l = abs(yCamera - p.y);
    p.z += l*l*0.00005;
    
    // bounding box optimisation
    vec2 inBB = vec2(abs(p.x), 35.0 - p.z);
    float bb = -inBB.y;
    bb = max(bb, inBB.x-80.0);
    if (bb > 1.0) {
		color = vec4(-1);
        return bb;
    }
    
    vec3 camera = cameraPath(p.y);
    float tunnel = 0.3 - length(p.xz-camera.xz);
    float dist = p.z+100.0;
	
    // keep position and scale
    vec4 pf = vec4(p, 0.01);
    pf.xyz *= pf.w;
    // add some perturbation
    float q = sin(pf.x*14.331)*1.0 + sin(p.y*0.1312-time*0.2)*1.5 + sin(pf.y*9.312)*0.3;
    // keep tornados uv
    vec4 inTorn = pf;
    
    for (float i = 0.99 ; i < 5.0 ; i += 1.0) {
        if (i < 3.0) {
        	pf.xyz += (snoise(pf.xyz*vec3(2.124, 2.2314, 10.347), 2.0).xyz-0.5)*0.06;
        }
        
        pf *= 4.5 - i*0.25;
        pf.x = abs(pf.x) - 0.5 + q*0.07;
        pf.y = abs(fract(pf.y)-0.5);
        vec4 potTorn = pf - vec4(2,0,0,0);
        
        // rounded corner
        float d = length(max(vec2(1.5, -0.5)-pf.xz*vec2(1,-1)*rot(0.15),0.0))-0.1;
        // rounded cylinder
        float m = sin(pf.z*50.0)*0.01 + sin(pf.z*q*2.0)*0.1;
        float dd = length(min(vec2(0.4 + m,0.9)-(vec2(length(potTorn.xy),pf.z)),0.0))-0.1;
        
        d = smin(d, dd, 0.1 + i*0.1) / pf.w;
        
        if (d < dist) inTorn = potTorn;
        dist = smin(d, dist, 0.015);
    }
    
    // add lightning
    vec4 thu = getThunder(p.y);
    vec3 inThun = p-thu.xyz;
    vec2 th = vec2(length(inThun.xz)-1.5 + sin(atan(inThun.x, inThun.z) * 12.0 +
                                               time*12.0)*0.03, inThun.y);
    float d = length(th)-0.15*(thu.w)+0.1;
    if (d < 0.0) {
        color = vec4(10);
        return 0.0;
    }
    
    dist = smax(dist, tunnel, 0.2);
    dist = min(dist, d);
    
    if (dist > r) {
        color = vec4(0);
        return dist;
    }
	
    pf.xyz /= pf.w;
    vec3 pp = vec3(pf.x, p.y, pf.z);
    pp.y -= time*0.5;
    pp.y *= 0.5;
    pp.zx *= rot(2.4412);
    pp.yz *= rot(0.2317);
    float alpha = noise(pp, 0.35, r*2.0).y*0.75;
    pp.xy *= rot(1.3187);
    pp.zy *= rot(0.2317);
    alpha += noise(pp, 1.5, r*2.0).y*0.5;
    pp.yx *= rot(4.5187);
    pp.zy *= rot(3.8317);
    alpha += noise(pp, 8.5, r*2.0).y*0.7;
    
    alpha *= alpha;
    alpha *= alpha;
    alpha = atan(1.0 * alpha) / (0.5 * PI);
    alpha *= alpha;
    alpha *= 8.0;
    alpha *= 1.0 - exp(dist*1.0*pow(inTorn.w, 1.5));
    alpha = min(alpha, 1.0);
    
    float col = pf.z-inTorn.z*10.0;
    col = sin(col*0.6) + sin(col*0.5) + sin(col*0.4);
    col = col * 0.2 + 0.5;
    color.rgb = vec3( sin(col * vec3(3, 4, 5)) * 0.5 + 0.5 );
    color.a = alpha;
    color.a *= smoothstep(r, -r, dist);
    color = clamp(color, 0.0, 1.0);
    
    return dist;
}

// return sun and rings
vec3 getBack( in vec3 dir ) {
    float d = dot(dir, sunDir);
    
    vec3 fog = getFog(dir);
    vec3 base = mix(fog*0.1, fog*0.8, smoothstep(0.8, 0.0, dir.z));
    
    base += exp((d-1.0)*5000.0) * COLOR_SUN * 10.0;
    base += exp((d-1.0)*400.0) * COLOR_SUN * 2.0;
    
    float planeDist = raytracePlane(vec3(0, 0, 1), dir, vec3(0),
                                    normalize(vec3(-9, -1.8, 2)));
    if (planeDist > 0.0) {
        vec3 p = vec3(0, 0, 1) + dir*planeDist;
        
        float l = length(p);
        float s = smoothstep( 0.66, 0.55, abs(l-2.0) );
        s += smoothstep( 0.01, 0.00, abs(l-1.33) );
        s = min(s, 1.0);
        float r = texture(iChannel0, vec3(0, 0, l*0.3)).r*0.5;
        r += texture(iChannel0, vec3(0, 0, l*2.0)).r*0.85;
        vec3 col = vec3( sin(r * 0.9 * vec3(3, 4, 5)) * 0.5 + 0.5 );
        col = mix(col, vec3(1), 0.7);
        float al = texture(iChannel0, vec3(0, 0, l*0.8)).r;
        al = smoothstep(0.2, 0.8, al);
        float m = smoothstep( 0.2, 0.0, abs(l-2.9) )*0.7;
        base = mix(base, col*vec3(1.0, 0.9, 0.6), s*al + m);
    }
    
    base = mix(base, fog, smoothstep(0.1, 0.0, dir.z));
    
    return base;
}

// volumetric rendering at half resolution, in linear space
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    // render at half resolution (checkerboard rendering)
    if (fragCoord.x > RENDERSIZE.x * 0.5) {
        fragColor = vec4(0);
        return;
    } else {
    	fragCoord.x = floor(fragCoord.x) * 2.0 + mod(fragCoord.y, 2.0);
    }
    
    
    vec2 uv = (fragCoord.xy - RENDERSIZE.xy*0.5) / RENDERSIZE.x;
    float fovScale = tan(FOV*0.5);
	vec3 dir = (vec3(uv.x*fovScale, 1.0, uv.y*fovScale));
	dir.y += exp(length(uv)*2.0) * 0.1;
    dir = normalize(dir);
    
	vec2 mouse = vec2(0);
    if (_mouse.z > 0.5) {
        mouse = (_mouse.xy-_mouse.zw)*0.01;
    }
    
    // get random values
    vec3 dither = hash33(vec3(fragCoord, FRAMECOUNT));
    // add motion blur
    time = TIME + dither.y * (1.0/60.0);
    // camera position
    yCamera = time*2.0 - 25.0;
    vec3 from = cameraPath(yCamera);
    vec3 target = cameraPath(yCamera+0.1);
    vec3 targetDir = normalize(target-from);
    dir.yz *= rot(-asin(targetDir.z));
    dir.xy *= rot(atan(targetDir.x, targetDir.y));
    dir.yz *= rot(-mouse.y + 0.1);
    dir.xy *= rot(mouse.x);
    
    // sine of the pixel angle
    float sinPix = sin(FOV/RENDERSIZE.x*0.5)*2.0;
    // accumulate color front to back
    vec4 acc = vec4(0, 0, 0, 1);
    
    vec4 dummy = vec4(0);
    float d = map(from, 0.0, dummy);
    float totdist = max(MARCH_MIN_STEP, d)*dither.x*1.0;
    vec4 fog = vec4(getFog(dir), 0);
    vec3 offset = (hash33(vec3(0, 0, FRAMECOUNT))-0.5)*0.4;
    
    for (int i = min(0, FRAMECOUNT) ; i < MARCH_STEPS ; i++) {
        
        // sample the field here
        vec3 p = from + totdist * dir;
        float r = max(MARCH_MIN_STEP*0.5, sinPix*totdist);
        vec4 color = vec4(0);
        float d = map(p, r, color);
        
        if (color.a < -0.5) {
            totdist = 9e9;
            break;
        } else if (d > r) {
            // away from the surface, skip
            totdist += max(MARCH_MIN_STEP, d*0.9);
        } else {
            
            bool arc = color.a > 5.0;
            
            // blend sample point toward sun/lightning
            vec4 thun = getThunder(p.y);
            vec3 inThun = thun.xyz-p + offset;
            float lenThun = length(inThun);
            float len = smoothstep(4.0, 0.5, lenThun);
            vec3 sampDir = normalize(mix(sunDir, inThun/lenThun, thun.w*len*len));
            
            // get another sample and run lighting
            vec3 samp = p + sampDir*MARCH_EPS;
            vec4 colorD = vec4(0);
            float dD = map(samp, r, colorD);
            vec2 derivV = vec2(color.a-colorD.a, dD-d)/MARCH_EPS;
            float deriv = mix(derivV.x, derivV.y,
                              (1.0-thun.w*len)*smoothstep(4.0, 0.0, abs(d)));
            
            // ad hoc formula for lighting
            deriv = atan(6.0 * deriv) / (0.5 * PI) * 0.5 + 0.5;
            deriv *= deriv;
            deriv = mix(deriv, 1.0, smoothstep(3.5, 0.0, lenThun)*thun.w);
            
            color.rgb *= mix(deriv, 1.0, thun.w*len)*COLOR_SUN;
            color.rgb += exp(-lenThun*7.0+8.0)*thun.w*1000.0 * vec3(0.3, 0.6, 0.8);
            
            // apply fog
            fog.a = getFog( color.rgb, totdist, from, dir );
            color.rgb = mix(color.rgb, fog.rgb, fog.a);
            
            // inside the surface, enable cone-tracing
            color.a = 1.0 - pow(1.0 - color.a, r*20.0);
            // overwrite for lightning
            if (arc) color = vec4(vec3(3), 1);
            color.a *= smoothstep(r, -r, d);
            
            // accumulate color
            acc.rgb += acc.a * (color.a*color.rgb);
        	acc.a *= (1.0 - color.a);
            
            // we can break if the accumulated alpha is 1
            if ( acc.a < 0.01 ) break;
            
            // go forward
            totdist += r*2.0;
        }
    }
    
    // add volumetric rendering on top of the background
    vec3 back = getBack(dir);
    fog.a = getFog( back, totdist, from, dir);
    back.rgb = mix(back.rgb, fog.rgb, fog.a);
    fragColor.rgb = acc.a*back+acc.rgb;
    fragColor.rgb = max(fragColor.rgb, 0.0);
    
    fragColor.a = 1.0;
    
	return fragColor; 
 } 



// Buffer A shader runs the volumetric rendering at half the resolution.
// Image shader reconstructs the frame and apply tonemapping/gamma/dithering.
// The cloud landscape is a repeated distance field deformed with texture based noise.
// Outside the shape, the cone can raymarch and skip empty space very quickly.
// Inside the shape, we raymarch and accumulate opacity along a cone.
// We break when the accumulated opacity is close to 1 or when we hit the sky.
// Filtered value noise (through 3D textures) is used for opacity/coloring.

// Lighting is done by taking the derivative of the opacity/distance field in the
// direction of the light source. During lightning, the direction is blended
// toward the closest lightning bolt.
// https://iquilezles.org/articles/derivative

// Fog is tweaked from iq formula taken from here.
// https://iquilezles.org/articles/fog

// from https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
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

// Dave_Hoskins hash
#define HASHSCALE3 vec3(.1031, .1030, .0973)
vec3 hash33( in vec3 p3 ){
	p3 = fract(p3 * HASHSCALE3);
    p3 += dot(p3, p3.yxz+19.19);
    return fract((p3.xxy + p3.yxx)*p3.zyx);
}

// checkerboard reconstruction and tonemapping
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    // render at half resolution (checkerboard rendering)
    vec2 ch = vec2( step(mod(fragCoord.x + fragCoord.y, 2.0), 0.5), 0 );
    vec4 a = texelFetch(BuffA, ivec2( (fragCoord+ch.xy) * vec2(0.5,1) ), 0);
    vec4 b = texelFetch(BuffA, ivec2( (fragCoord+ch.yx) * vec2(0.5,1) ), 0);
    vec4 c = texelFetch(BuffA, ivec2( (fragCoord-ch.xy) * vec2(0.5,1) ), 0);
    vec4 d = texelFetch(BuffA, ivec2( (fragCoord-ch.yx) * vec2(0.5,1) ), 0);
    fragColor = (a+b+c+d)*0.25;
    
    // exposition
    fragColor.rgb *= 1.4;
    // vignette
    vec2 uv = fragCoord / RENDERSIZE.xy * 2.0 - 1.0;
    fragColor.rgb *= (1.0-dot(uv,uv)*0.1);
    // tonemapping
    fragColor.rgb = ACESFitted(fragColor.rgb);
    // gamma correction
    fragColor.rgb = pow( fragColor.rgb, vec3(1.0/2.2) );
    // dithering
    fragColor.rgb += (hash33(vec3(fragCoord, FRAMECOUNT))-0.5)*0.03;
    
    fragColor.a = 1.0;
    
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