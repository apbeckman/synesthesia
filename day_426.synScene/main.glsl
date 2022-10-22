

			//******** Common Code Begins ********


float r11(float i){ return fract(sin(i*15.126)*115.6);}

float ss( float c, float power, float bias){
    c = clamp(c,-0.,1.);
    //c = smoothstep(0.,1.,c);
    
    c = pow(c,1. + bias);
    
    float a = pow( abs(c), power);
    float b = 1.-pow( abs(c - 1.), power);
    
    return mix(a,b,c);
}

float valueNoise(float i, float p){ return mix(r11(floor(i)),r11(floor(i) + 1.), ss(fract(i), p,0.6));}

float valueNoiseStepped(float i, float p, float steps){ return mix(  floor(r11(floor(i))*steps)/steps, floor(r11(floor(i) + 1.)*steps)/steps, ss(fract(i), p,0.6));}

#define kSegment (floor(pow(mod(TIME,30.)/30.,2.)*2.))

#define pi acos(-1.)

#define R (RENDERSIZE.xy)
#define T(u) texture(iChannel2,(u)/R)
#define T0(u) texture(iChannel0,(u)/R)

//#define sino(a) sin(a)3

//#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))

#define pmod(p,d) mod(p - (d)*0., (d)) - 0.5*(d)

#define pal(a,b,c,d,e) (a + (b)*sin((c)*(d) + (e)))

vec2 r12(float i){float r=r11(i );  return vec2(r,r11(i + r + 2.));}

#define xor(a,b,c) min(max((a),-(b)), max((b),-(a) - c)) 
#define pal(a,b,c,d,e) (a + (b)*sin((c)*(d) + (e)))
#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))

#define pi acos(-1.)

mat3 getOrthogonalBasis(vec3 direction){
    direction = normalize(direction);
    vec3 right = normalize(cross(vec3(0,1,0),direction));
    vec3 up = normalize(cross(direction, right));
    return mat3(right,up,direction);
}
float cyclicNoise(vec3 p, bool turbulent, float time){
    float noise = 0.;
    
    p.yz *= rot(.5);
    p.xz *= rot(.5);
    float amp = 1.;
    float gain = 0.4 + sin(p.x*1.5 + time)*0.;
    const float lacunarity = 1.24;
    const int octaves = 4;
    
     float warp = .1 + sin(time*0.5)*0.05;
    float warpTrk = 1.8 ;
    const float warpTrkGain = 0.51;
    
    vec3 seed = vec3(-2,-2.,0.5);
    mat3 rotMatrix = getOrthogonalBasis(seed);
    
    for(int i = 0; i < octaves; i++){
        
        p -= sin(p.zxy*warpTrk + vec3(0,-time*0.5,0) - .4*warpTrk)*warp; 
        noise += sin(dot(cos(p), sin(p.zxy + vec3(0,time*0.5,0))))*amp;
    
        p *= rotMatrix;
        p *= lacunarity;
        
        warpTrk *= warpTrkGain;
        amp *= gain;
    }
    
    if(turbulent){
        return 1. - abs(noise)*0.5;
    
    }{
        return (noise*0.25 + 0.5);

    }
}



float cyclicNoiseB(vec3 p, bool turbulent, float time){
    float noise = 0.;
    
    p.yz *= rot(1.);
    float amp = 1.;
    float gain = 0.8 + sin(p.z*0.2)*0.2;
    const float lacunarity = 1.6;
    const int octaves = 2;
    
    const float warp =.4;    
    float warpTrk = 1.5 ;
    const float warpTrkGain = .2;
    
    vec3 seed = vec3(-1,-2.,0.5);
    mat3 rotMatrix = getOrthogonalBasis(seed);
    
    for(int i = 0; i < octaves; i++){
        
        p += sin(p.zxy*warpTrk + vec3(0,-time*2.,0) - 2.*warpTrk)*warp; 
        noise += sin(dot(cos(p), sin(p.zxy + vec3(0,time*0.3,0))))*amp;
    
        p *= rotMatrix;
        p *= lacunarity;
        
        warpTrk *= warpTrkGain;
        amp *= gain;
    }
    
    if(turbulent){
        return 1. - abs(noise)*0.5;
    
    }{
        return (noise*0.25 + 0.5);

    }
}

vec3 sdgBox( in vec2 p, in vec2 b )
{
    vec2 w = abs(p)-b;
    vec2 s = vec2(p.x<0.0?-1:1,p.y<0.0?-1:1);
    float g = max(w.x,w.y);
    vec2  q = max(w,0.0);
    float l = length(q);
    return vec3(   (g>0.0)?l  :g,
                s*((g>0.0)?q/l:((w.x>w.y)?vec2(1,0):vec2(0,1))));
}


float sdSq(vec2 p, vec2 s){
    p = abs(p) - s;
    return max(p.x,p.y);
}


float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); }

			//******** BuffA Code Begins ********


float envcnt = 0.;

float getEnv(float t, float speed, float pa, float pb, float jumpAmt, bool cnt){
    //return pow(sin((t - 0.5)*3.14),1.)*0.5 + 0.5;
    t = clamp(t*speed,0.,1.);
    
    envcnt += float(t > 0.99 && cnt);
    //t = smoothstep(0.,1.,t);
    pa += 1.;
    pb += 1.;
    
    float c = cos(t*3.14);
    float a = 1.- ((pow(abs(c),pa)*sign(c))*0.5 + 0.5);
    float b = 1.-((pow(abs(c),pb)*sign(c))*0.5 + 0.5);
    
    a = pow(sin(t*3.14/2.),pa);
    b = 1.-pow(sin((-t + 1.)*3.14/2.),pb);
    
    b *= 1. + (
            smoothstep(0.,1.,t) *smoothstep(0.99,0.7,t)*jumpAmt
        );
    return mix( a, b,t);
}


vec4 renderPassA() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;

    vec2 uv = (U - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    vec3 col = vec3(0.,0.,0.);
    
    
    float segmentMod =  + sin(floor(kSegment))*20.;
    
    float vn = valueNoise(TIME + uv.x*0.01 + segmentMod, 2.)*1.;
    float t = TIME*1.3 - vn + 2. + segmentMod;
        
    float env = getEnv(t, 1., 2., 2., 1., false);
    
    vec2 offs = vec2(sin(t)*0.5,cos(t*0.6)*1.6) ;
    
    vec2 ouv = uv;
    uv -= offs*0.2;
    float nb = cyclicNoise(vec3(uv*15.5 + t*0. + offs*0.3,t*0.2), false, t*0.);
    
    float n = cyclicNoise(vec3(uv*((22.5 ) - dot(ouv,ouv)*4.4),t*1.5 - nb*1.) , true, t*1.);
    n *= 1.;
    
    float r = .16 - n*2. + nb*5.;
    int didx = 3;
    vec2 dfn = vec2(T(U + vec2(1.,0)*r)[didx] - T(U - vec2(1.,0)*r)[didx],T(U + vec2(0.,1)*r)[didx] - T(U - vec2(0.,1)*r)[didx]);

    r = 15. - n*0. + nb*4.;
    didx = 3;
    
    vec2 dfnb = vec2(T(U + vec2(1.,0)*r)[didx] - T(U - vec2(1.,0)*r)[didx],T(U + vec2(0.,1)*r)[didx] - T(U - vec2(0.,1)*r)[didx]);
    
    
    col = mix(col, vec3(0.1+ sin(length(dfn)*44. + TIME + n*3.)*0.5,0.4,0.5),0.4*smoothstep(-0.1 + n*0. + nb*0.0,0. + length(dfnb)*1.,length(dfn)));
    
    //col += 1.-abs(length(dfn.xyx))*4.;
    
    //col = mix(col,0.*vec3(0.,0.,0.1 + 0.2*sin(nb*10.)*0.2),smoothstep(-0.1,.1,-length(dfn)));
    
    vec3 c = vec3(0.4 + sin(length(dfn)*14. + TIME + n*13.)*0.4,0.3 + sin(length(dfn)*24. + TIME*5. +length(dfnb)*7.)*0.1,0.3);
    col = mix(col,3.*c,smoothstep(0.12,.25,abs(abs(n - 0.5) - 0.2) + length(dfn)*0.7));
    
    
    
    //vec3 c = vec3(0.+ sin(length(dfnb)*5. + TIME+ n*10.)*0.4+ sin(dfn.x*10.)*0.3 ,0.1 + n*0.05 + sin(nb*20. + t)*0.05 ,0.1  );
    //c *= 1.9;
    
    
    
    
    
    C = vec4(col,n);
	return C; 
 } 


			//******** BuffB Code Begins ********


vec4 renderPassB() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;

    vec4 fr = texture(BuffA,(U)/RENDERSIZE.xy);
   
    float r = 0.2 
        + sin(fr.w*5.)*.4
        + sin(fr.x*1.)*.71;
    int didx = 1;
    
    vec2 dfn = vec2(T0(U + vec2(1.,0)*r)[didx] - T0(U - vec2(1.,0)*r)[didx],T0(U + vec2(0.,1)*r)[didx] - T0(U - vec2(0.,1)*r)[didx]);
    
    vec2 sc = vec2(0)
        + pow(smoothstep(0.,1.,length(dfn.x)*9.),0.6)*.6;
    
    sc *= dfn*0.5;
    
    
    C.x =texture(BuffA,(U + sc*vec2(0,9))/RENDERSIZE.xy).x;
    
    C.y =texture(BuffA,(U + sc*vec2(0,-4))/RENDERSIZE.xy).y;
    
    C.z =texture(BuffA,(U + sc*vec2(-3,-5))/RENDERSIZE.xy).z;
    
    
	return C; 
 } 


			//******** BuffC Code Begins ********

// Fork of "Day 424" by jeyko. https://shadertoy.com/view/wldfDB
// 2021-02-16 09:06:35

// dont fullscreen

vec4 renderPassC() {
	vec4 C = vec4(0.0);
	vec2 fragCoord = _xy;

    fragCoord -= 0.5*RENDERSIZE.xy;
    //fragCoord *= 0.99;
    fragCoord += 0.5*RENDERSIZE.xy;
    
    float n1d = texelFetch(iChannel1,ivec2(mod(fragCoord + vec2(float(FRAMECOUNT),0.),256.)),0).x;
    
    C -= C;

    vec4 bloom = vec4(0);
    //bloom += texture(BuffB,fragCoord/R,0.);
    //bloom += texture(BuffB,fragCoord/R,1.);
    bloom += texture(BuffB,fragCoord/R,0.);
    //bloom += texture(BuffB,fragCoord/R,3.);
    bloom += texture(BuffB,fragCoord/R,9.);
    bloom += texture(BuffB,fragCoord/R,5.);
    bloom += texture(BuffB,fragCoord/R,6.);
    bloom += texture(BuffB,fragCoord/R,7.);
    bloom += texture(BuffB,fragCoord/R,6.);
    bloom += texture(BuffB,fragCoord/R,9.);
    bloom += texture(BuffB,fragCoord/R,12.);
    
    //bloom = smoothstep(0.4,1.,1.7*bloom/8.)*2.;
    bloom/=6.;
    C = bloom;
	return C; 
 } 


			//******** BuffD Code Begins ********

// Fork of "Day 424" by jeyko. https://shadertoy.com/view/wldfDB
// 2021-02-16 09:06:35

// dont fullscreen

vec4 renderPassD() {
	vec4 C = vec4(0.0);
	vec2 fragCoord = _xy;

    fragCoord -= 0.5*RENDERSIZE.xy;
    //fragCoord *= 0.99;
    fragCoord += 0.5*RENDERSIZE.xy;
    
    float n1d = texelFetch(BuffB,ivec2(mod(fragCoord + vec2(float(FRAMECOUNT),0.),256.)),0).x;
    
    C -= C;

    vec4 bloom = vec4(0);
    //bloom += texture(BuffC,fragCoord/R,0.);
    //bloom += texture(BuffC,fragCoord/R,1.);
    bloom += texture(BuffC,fragCoord/R,5.);
    //bloom += texture(BuffC,fragCoord/R,3.);
    bloom += texture(BuffC,fragCoord/R,4.);
    bloom += texture(BuffC,fragCoord/R,5.);
    bloom += texture(BuffC,fragCoord/R,3.);
    bloom += texture(BuffC,fragCoord/R,7.);
    bloom += texture(BuffC,fragCoord/R,2.);
    bloom += texture(BuffC,fragCoord/R,12.);
    bloom += texture(BuffC,fragCoord/R,9.);
    
    
            
    bloom = smoothstep(0.4,1.5,1.5*bloom/8.)*1.;
    
    
    C = texture(BuffB,fragCoord/R);
    
    C = mix(C,bloom*1.4,pow(smoothstep(0.,1.5,length(bloom)),2.)*0.5);
            
    if(kSegment > 0.5){
        C = 1. - C*0.8;
    }
	return C; 
 } 


// Fork of "day 425" by jeyko. https://shadertoy.com/view/3ltBWj
// 2021-02-17 09:03:04

// Fork of "Day 424" by jeyko. https://shadertoy.com/view/wldfDB
// 2021-02-16 09:06:35

// dont fullscreen, bloom and some other things don't work well there.

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    fragCoord -= 0.5*RENDERSIZE.xy;
    //fragCoord *= 0.99;
    fragCoord += 0.5*RENDERSIZE.xy;
    
    float n1d = texelFetch(image30,ivec2(mod(fragCoord + vec2(float(FRAMECOUNT),0.),256.)),0).x;
    vec3 n  = texelFetch(image30,ivec2(mod(fragCoord  + n1d*200. ,256.)),0).xyz;
    
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    
    fragColor.xyz =texture(BuffD,fragCoord/RENDERSIZE.xy).xyz;
    
    
    fragColor.xyz *= 1. - dot(uv,uv)*0.7;
    fragColor.xyz = pow(max(fragColor.xyz,0.), vec3(1.15,1.3,1.15) + dot(uv,uv)*0.6);
    
    //fragColor.xyz = 1. - fragColor.xyz;
    
    fragColor.xyz = pow(fragColor.xyz, vec3(0.4545 + n*0.15));
    
    
    
    //fragColor = texture(iChannel2,fragCoord/RENDERSIZE.xy);
    
    
    
    fragColor.xyz += smoothstep(1.,0.,length(fragColor))*n*0.04;
    
    fragColor.xyz -= smoothstep(0.,1.,length(fragColor))*n*0.05;
       
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderPassC();
	}
	if(PASSINDEX == 3){
		return renderPassD();
	}
	if(PASSINDEX == 4){
		return renderMainImage();
	}
}