

			//******** Common Code Begins ********

#define R (RENDERSIZE.xy)
#define T(u) texture(BuffA,fract((u)/R))
#define T1(u) texture(BuffB,fract((u)/R))
#define T2(u) texture(image30,fract((u)/R))
//#define T3(u) texture(iChannel3,fract((u)/R))

#define TT(u,T) texture(T,fract((u)/res))

#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))
#define pi acos(-1.)

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

    vec2 uv = (U - 0.5*R)/R.y;   
    
    vec2 OU = U;
    //U -= 0.5*_uvc*(syn_BassLevel*(1.0+syn_Intensity));
    U.xy += _uvc*PI*Zoom*((syn_BassLevel+syn_MidLevel*0.75)*(1.+syn_Intensity));
    
    //U -= 0.5*R;   
    //U *= 1. - fract(TIME*kTimeCoeff)*0.002;
    //U *= rot(0.1*pow(fract(TIME*kTimeCoeff*0.125/2.),17.1)*0.1);
    
    //U *= 1. - dot(uv,uv)*0.001;
    //U += 0.5*R;
    //vec2 grady = getGradient( BuffA, U, 3., 1, R);
    //vec2 gradz = getGradient( BuffA, U, 3., 2, R);
        
    
    float neighrange = 25. - sin(smoothTime*0.00025)*2.;
    neighrange *= 0.65+syn_BassLevel;
    vec2 muv = (_mouse.xy - 0.5*R)/R.y;   
    
    if(_mouse.z > 0.)
        neighrange *= 1. + 1.*smoothstep(0.01,0.,length(uv - muv) - .3);

    
    
    vec2 grad = getGradient( BuffA, U, neighrange*0.4, 0, R);
       
    
    
    
      
    float neighs = 0.;
    
    float iiters = 128.*0.25;            
    float jiters = 128.*0.25;
    
    for(float i = 0.; i < iiters; i++){
        vec2 p = U; 
        
        vec2 offs = vec2(0,1.)*rot(2.*pi*(i + 0.)/iiters);
        
        float samp = 0.;
        float jiters = 32.;
        for(float j = 0.; j < jiters; j++){
            vec2 ioffs = p + offs*mix(neighrange*0.2,neighrange,j/jiters);
            samp += texture(BuffA,fract(ioffs/R),0.5).x/jiters;
        
        }
        
        neighs += samp/iiters*4.;
    }
   
    C = T(U);
    
    float deadness = smoothstep(1.,0.,abs(C.x));
    float aliveness = smoothstep(0.,1.,abs(C.x));
    //deadness = 1. - C.x;
    //aliveness = C.x;
    
    //U += grad*.01*smoothstep(1.8,0.1,abs(neighs - 3.))*aliveness;
    //U += grad*1.9*smoothstep(2.4,1.,abs(neighs - 2.))*aliveness;
    
    
    
    C = T(U);
    
    //deadness = smoothstep(1.,0.,abs(C.x));
    //aliveness = smoothstep(0.,1.,abs(C.x));
    //deadness = 1. - C.x;
    //aliveness = C.x;
    
    
    
    float speed = 0.8;
    C = mix(C,vec4(0),smoothstep(.5,0.,neighs - 1.)*aliveness*speed);
    C = mix(C,vec4(0),smoothstep(-0.125,0.,neighs - 3.)*aliveness*speed);
    
    
    C = mix(C,vec4(1),
            smoothstep(0.5,0.,abs(neighs - 2.5))*
        deadness*
        speed);
    
    
    //vec2 gradw = getGradient( BuffA, U, 3., 3, R);
    
    //grad *= rot(.2);
    
    //U += grad[(FRAMECOUNT/30)%0]*.1*sin(TIME);
    
    
    //if(_mouse.z > 0.)
    
    //C = blur(BuffA, U, R);
    
    
    
    if(int(FRAMECOUNT)%4 >= 0){
        //C = T(OU);
    }
    
    
    if(FRAMECOUNT <= 1){
        C = T2(U*2.2);
    }
	return C; 
 } 


			//******** BuffB Code Begins ********


vec4 renderPassB() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;

    vec4 fr = texture(BuffA,(U)/RENDERSIZE.xy);
   
    float r = 17.4
        + sin(fr.x*2. - 1.)*1.4;
    int didx = 0;
    
    vec2 dfn = vec2(T(U + vec2(1.,0)*r)[didx] - T(U - vec2(1.,0)*r)[didx],T(U + vec2(0.,1)*r)[didx] - T(U - vec2(0.,1)*r)[didx]);
    
    vec2 sc = vec2(0)
        + pow(smoothstep(0.,1.,length(dfn.xy)*4.),.2)*.15*dfn;
    
    
    
    C.x =texture(BuffA,(U + sc*vec2(0,2))/RENDERSIZE.xy).x;
    
    C.y =texture(BuffA,(U + sc*vec2(0,-5))/RENDERSIZE.xy).y;
    
    C.z =texture(BuffA,(U + sc*vec2(5,-5.))/RENDERSIZE.xy).z;
    
	return C; 
 } 


// Fork of "Day 431" by jeyko. https://shadertoy.com/view/ttKBzD
// 2021-02-23 15:15:41

// tried to make a continous Game of Life
// ended up being too wormy to be Game of Life, but cool in its own right! 

vec4 renderMainImage() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;

    C = C*.0;
    
    vec2 uv = (U - 0.5*R)/R.y;   
    
    //C += T1(U);
    C = T1(U);
    
    
    
    
    float n1d = texelFetch(image30,ivec2(mod(U + vec2(float(FRAMECOUNT),0.),256.)),0).x;
    vec3 n  = texelFetch(image30,ivec2(mod(U  + n1d*400. ,256.)),0).xyz;
    
    
    C *= 1. - dot(uv,uv*0.125);
    C = smoothstep(0.,1.,C);
   // C.xyz += _uvc.yyy*_uvc.xxx*C.xyz*_uvc.xxy;
    C.xyz = pow(max(C.xyz,0.), vec3(0.55) + dot(uv,uv)*0.6);
    
    
    
    
    C.xyz += smoothstep(1.,0.,length(C))*n*0.15;
    
    C.xyz -= smoothstep(0.,1.,length(C))*n*0.05;
    
    
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