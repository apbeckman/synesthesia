

			//******** Common Code Begins ********

#define R (RENDERSIZE.xy)
#define T(u) texture(iChannel0,(u)/R)
#define T1(u) texture(iChannel1,(u)/R)
#define T2(u) texture(iChannel2,(u)/R)
#define T3(u) texture(iChannel3,(u)/R)

#define TT(u,T) texture(T,(u)/res)

#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))


#define kTimeCoeff 60./167.85*4.

vec2 getGradient(sampler2D tex, vec2 u, float offset, int channel, vec2 res){
    return vec2(
         TT(u + vec2(1,0)*offset,tex)[channel] - TT(u - vec2(1,0)*offset,tex)[channel],
         TT(u + vec2(0,1)*offset,tex)[channel] - TT(u - vec2(0,1)*offset,tex)[channel]  
    );
}
//#define Neighbordhood() vec4 me = T() 



float sdSegment( in vec2 p, in vec2 a, in vec2 b )
{
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}

vec4 sharpen(sampler2D channel,vec2 fragCoord, vec2 resolution){
    float kernel [9];vec2 offset [9];
    
    vec2 step = vec2(1);
    
    offset[0] = vec2(-step.x, -step.y);
    offset[1] = vec2(0.0, -step.y);
    offset[2] = vec2(step.x, -step.y);
    
    offset[3] = vec2(-step.x, 0.0);
    offset[4] = vec2(0.0, 0.0);
    offset[5] = vec2(step.x, 0.0);
    
    offset[6] = vec2(-step.x, step.y);
    offset[7] = vec2(0.0, step.y);
    offset[8] = vec2(step.x, step.y);
    
    
    kernel[0] = 0.0; kernel[1] = -0.25; kernel[2] = 0.0;
    kernel[3] = -0.25; kernel[4] = 1.0; kernel[5] = -0.25;
    kernel[6] = 0.0; kernel[7] = -0.25; kernel[8] = 0.0;
    
    vec4 sum = texture(channel, (fragCoord)/resolution);
    
    for (int i = 0; i < 9; i++) {
        vec4 color = texture(channel, (fragCoord + offset[i])/resolution);
        sum += color * kernel[i];
    }
    
    sum = clamp(sum,0.,1.);
    
    return sum;
}

vec4 blur(sampler2D channel,vec2 fragCoord, vec2 resolution){
    
    float kernel [9];vec2 offset [9];

     vec2 step = vec2(0.5);
    
    offset[0] = vec2(-step.x, -step.y);
    offset[1] = vec2(0.0, -step.y);
    offset[2] = vec2(step.x, -step.y);
    
    offset[3] = vec2(-step.x, 0.0);
    offset[4] = vec2(0.0, 0.0);
    offset[5] = vec2(step.x, 0.0);
    
    offset[6] = vec2(-step.x, step.y);
    offset[7] = vec2(0.0, step.y);
    offset[8] = vec2(step.x, step.y);
    
    kernel[0] = 1.0; kernel[1] = 1.0; kernel[2] = 1.0;
    kernel[3] = 1.0; kernel[4] = 1.0; kernel[5] = 1.0;
    kernel[6] = 1.0; kernel[7] = 1.0; kernel[8] = 1.0;
    
    vec4 sum = vec4(0);
    
    for (int i = 0; i < 9; i++) {
        vec4 color = texture(channel, (fragCoord + offset[i])/resolution);
        sum += color * kernel[i];
    }
    
    sum /= 9.;
    sum = clamp(sum,0.,1.);
    
    return sum;
}

			//******** BuffA Code Begins ********


vec4 renderPassA() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;

    
    U -= 0.5*R;
    U *= 1. - fract(TIME*kTimeCoeff)*0.002;
    U *= rot(0.1*pow(fract(TIME*kTimeCoeff*0.125/2.),17.1)*0.1);
    
    U += 0.5*R;
    vec2 grad = getGradient( BuffB, U, 3., 0, R);
    vec2 grady = getGradient( BuffB, U, 3., 1, R);
    vec2 gradz = getGradient( BuffB, U, 3., 2, R);
    //vec2 gradw = getGradient( BuffB, U, 3., 3, R);
    
    //grad *= rot(.2);
    
    U -= grad*.2*sin(TIME);
    
    float id = floor(TIME*kTimeCoeff/10.); 
    float md = mod(TIME*kTimeCoeff, 10.);
    if(md < 1.){
        if (id == 0.){
            U -= grady*14.5;
        
        } else {
            U -= gradz*14.5;
    
        }
    }
    
    C = blur(BuffB, U, R);

    
    if(int(FRAMECOUNT)%2100 == 0){
    
        C = 1. - C;
    }
    
    
    if(FRAMECOUNT == 0){
        C = T3(U*0.2);
    }
	return C; 
 } 


			//******** BuffB Code Begins ********


vec4 renderPassB() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;

    C = sharpen(BuffA, U, R);
    
    vec2 uv = (U - 0.5*R)/R.y;   
    
    vec2 muv = (_mouse.xy - 0.5*R)/R.y;   
    if(_mouse.z > 0.)
        C = mix(C,pow(T2(U)*1.,vec4(5.)),smoothstep(0.01,0.,length(uv - muv) - .1));
    
    
    vec4 r = texture(image30,(vec2(FRAMECOUNT%256,floor(float(FRAMECOUNT)/256.)) + 0.5 )/256.).xyzw;
    
    if(FRAMECOUNT % 40 < 2){
        C = mix(C,vec4(0),
                smoothstep(0.01,0.,sdSegment( uv, vec2(r.x,r.y)*2. - 1., vec2(r.z,r.w)*2. - 1. ) - 0.0
            ));

    }
   
	return C; 
 } 



// blur and sharpen from https://www.shadertoy.com/view/MtdXW4

// LOOKS BETTER ON 144 HZ

vec4 renderMainImage() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;

    C = C*0.;
    
    vec2 uv = (U - 0.5*R)/R.y;   
    
    //C += T1(U);
    C = T1(U);
    
    
    float md = kTimeCoeff;
    float fac = fract(TIME*md*0.25);
    
    fac = pow(fac,10.5)*smoothstep(1.,0.96,fac);
    C = mix(C.xxxx, C.yyyy,fac);
    
    //C = mix(C,T1(U).zxzy,dot(uv,uv)*0.1);

    
    float n1d = texelFetch(image30,ivec2(mod(U + vec2(float(FRAMECOUNT),0.),256.)),0).x;
    vec3 n  = texelFetch(image30,ivec2(mod(U  + n1d*200. ,256.)),0).xyz;
    
    C.xyz = pow(max(C.xyz,0.), vec3(1.,1.,1.) + dot(uv,uv)*0.6);
    
    
    
    
    C.xyz += smoothstep(1.,0.,length(C))*n*0.1;
    
    C.xyz -= smoothstep(0.,1.,length(C))*n*0.05;
    
    if(mod(TIME*kTimeCoeff,10.) < 1.){
        C = 1. - C;
    }
    
	return C; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderMainImage();
	}
}