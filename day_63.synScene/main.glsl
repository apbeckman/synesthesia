

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


vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
	vec3 dir = normalize(lookAt - ro);
    vec3 right = normalize(cross(vec3(0,1,0), dir));
    vec3 up = normalize(cross(dir, right));
    return normalize(dir + right*uv.x + up * uv.y);
}
#define pmod0(p,x) mod(p,x) - 0.5*x
#define mx (TIME*0. + 20.*_mouse.x/RENDERSIZE.x)

#define pi acos(-1.)
#define tau (2.*pi)
#define time (TIME*0.125)
#define rot0(x) mat2(cos(x),-sin(x),sin(x),cos(x))
#define pal0(a,b,c,d,e) (a + b*cos(tau*(c*d + e)))

vec3 pA = vec3(0);

vec2 map(vec3 p){
	vec2 d = vec2(10e6);
	
    vec3 q = p;
    

    float sc = 3.4 - sin(TIME*0.5)*2. ;
    float dp = dot(p,p);
    p /= dp;
    p*= sc;
    p=sin(p+vec3(-time*tau*2.2,1.4 - 1.*time*tau,.1 + time*tau*3. + sin(time*tau)*1.5));
    pA = p;
    d.x = 0.;
    d.x = length(p) - 0.7 + length(q)*0.3;
    d.x *= 0.5;
	return d*dp/sc;
}

vec3 glow = vec3(0);
vec2 trace(vec3 ro, vec3 rd,inout vec3 p,inout float t, inout bool hit){
	vec2 d = vec2(10e6);
	t = 0.; hit = false; p = ro;
    
    for(int i = 0; i < 170; i++){
    	d = map(p);
        glow += exp(-d.x*90.);
        if(d.x < 0.0001){
        	hit = true;
            break;
        }
        t += d.x;
        p = ro + rd*t;
    }
    
    
    return d;
}

vec3 getNormal(vec3 p){
    vec2 t = vec2(0.001,0);
	return normalize(map(p).x - vec3(
    	map(p - t.xyy).x,
    	map(p - t.yxy).x,
    	map(p - t.yyx).x
    ));
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    vec3 col = vec3(0);
    col += vec3(0.04,0.1,0.3)*0.02;
    
    col -= vec3(0.01,0.07,0.3)*length(uv)*0.03;
    uv *= 1. + dot(uv,uv)*0.24;
    vec3 ro = vec3(0.6,0.,1.)*0.9;

    
    //ro.x += sin(TIME)*0.1;
    ro.z += 0.5-sin(TIME*0.5)*0.3;
    vec3 lookAt = vec3(0.);
    vec3 rd = getRd(ro, lookAt, uv);
    rd.xy *= rot0(sin(TIME*0.5)*0.6);
    rd.xz *= rot0(sin(TIME*0.75)*0.2);
    vec3 p; float t; bool hit;
    vec2 d = trace(ro, rd, p, t, hit);
    
    if(hit){
        vec3 pAA = pA;
        float modD = 0.1;
        float id = floor(pA.x/modD);
        pA = pmod0(pA, modD);
        //col += pal0(0.,vec3(0.7,0.8 ,0.8)*1., vec3(3.1,9.5,4.1 ), 1.5, id*1.2 + pAA.z*0.9 + pAA.y*0.9 + TIME*0.92)*1.;
        vec3 balls = pal0(0.5,vec3(0.6,0.4 ,0.8)*1., vec3(3.1,9.5,4.1 ), 1.5, id*1.2 + pAA.z*0.9 + pAA.y*0.1 + TIME*0.52 )*1.;
    	balls = mix(balls, 1. - balls,0. + 1.*pow(abs(uv.x)*0.55,2.9)*4.5);
        col += balls;
        col *= step(abs(sin(id*20.))*1., 0.7);
        col -= exp((abs(pA.x) - modD*0.5)*100.);

        col -= exp(-length(p)*15.)*10.;
        //col += 0.1;
        //col += smoothstep(0.01,0., length(pA.x) - modD*0.175);
    	col -= glow*0.025;
    }
    
    col = clamp(col, 0., 1.);

    col *= 1.1;
    col = pow(col, vec3(0.45));
    col *= 1. + 1.*pow(abs(uv.x)*0.55,2.9)*3.5;
    col *= 1. - 1.*pow(abs(uv.y)*1.0,2.9)*0.5;
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


// Fork of "Day 61" by jeyko. https://shadertoy.com/view/WlKXRR
// 2020-02-20 13:48:59

// thanks to mla and Kali!
// They have super nice examples on inversion

// it's really simple, basically
// p /= dot(p,p);
// p = sin(p);
// SDFs
// return distance*dot(p,p);


// and ofc Inigo quilez for pallete!


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    fragColor = texture(BuffA, fragCoord/RENDERSIZE.xy);
    fragColor += texture(BuffA, fragCoord/RENDERSIZE.xy, 6.)*0.3;
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