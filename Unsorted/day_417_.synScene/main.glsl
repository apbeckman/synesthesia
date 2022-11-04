//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** Common Code Begins ********


#define pi acos(-1.)


#define sint(a) (asin(sin(a))*2. - 1.)

#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))

#define pmod(p,d) mod(p - (d)*0.5, (d)) - 0.5*(d)

float r11(float i){ return fract(sin(i*12.126)*12.6);}

#define xor(a,b,c) min(max((a),-(b)), max((b),-(a) - c)) 

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



			//******** BuffA Code Begins ********



mat3 getOrthogonalBasis(vec3 dir){
    vec3 right = normalize(cross(vec3(0,1,0), dir));
    vec3 up = normalize(cross(dir, right));
    return mat3(right,up,dir);
}

float cyclicNoise(vec3 p){
    
    //p.yz *= rot(1.4);
    
    float n = 0.;
    float amp = 1.;
    float gain = 0.5;
    float lac = 1.3 ;
    
    vec3 seed = normalize(vec3(3,-1,2));
    mat3 rotm = getOrthogonalBasis(seed);

    for(int i = 0; i < 5; i++){
        p -= cos(p.zxy*1.5*gain*2. + TIME)*0.1;
        n += (dot(sin(p), cos(p.zxy)))*amp;
    
        amp *= gain;
        p *= lac*rotm;

    }
    return n;
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = ( fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    
    vec2 oouv = uv;
    
    
    uv *= 0.7;
    vec3 col = vec3(0);
    
    float n = cyclicNoise(vec3(uv*5.,TIME*0.5));
    
    n = cyclicNoise(vec3(uv*10. + n*1.15,TIME));
        
        
    float nb = cyclicNoise(vec3(uv*5.,TIME));
        
        
    
    
    float vn = valueNoise(TIME/2.,5.);
    float vnb = valueNoiseStepped(TIME/2. + n*0.4*pow(vn,4.) + 5.,2.,4.);
    float vnc = valueNoiseStepped(TIME/8.,4.,4.);
    
        
    uv.x += (n*2. - 1.)*0.02;
    uv.y += (nb*2. - 1.)*0.4;    
    
    float dfuv = fwidth(uv.x);
    
    float dfn = dFdx(n);
    
    uv.x = pmod(uv.x,0.04 - vn*0.02);
    //uv *= .6;
    
    
    float d = abs(uv.x);
    vec2 oouvb = oouv; 
    oouv *= rot(pi*vnb);
    
    float w = vn*0.3;
    oouv = abs(oouv) - w;
    float db = max(oouv.x, oouv.y);
    
     
    oouvb *= rot(-pi*vnb);
    
    oouvb = abs(oouvb) - w;
    float dc = max(oouvb.x, oouvb.y);
    
    
    dc = abs(dc) - vnc*0.;
    
    
    db = xor(db,abs(dc - 0.1*vnc) - vn*0.1 , +0.0);
    db -= 0.01;
    float lineW =  0.001 + smoothstep(fwidth(db)*15., 0., db)*0.004;
    
    
    //ouv = pmod(ouv,0.2);
    
    col += smoothstep(dfuv*(1. + sin(TIME + uv.x)*0.4 + 0.4),0.,d - lineW);
    
    
    if(floor(mod(TIME/2.,2.)) == 0. ){
        col = 1. - col;
    }
    
    
    col = mix(col,1.- col, smoothstep(fwidth(oouv.x)*1., 0., db));
    
    col = mix(col, 1. - col, 
        step(0.6,valueNoiseStepped(TIME*25.,2.,4.)) *
        step(0.2,valueNoiseStepped(TIME*2.5,2.,2.))
        
        );
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 



			//******** BuffB Code Begins ********


vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    float sc = 0. + valueNoise(TIME*18.,2.)*0.2;
    
    fragColor.x =texture(BuffA,(fragCoord + sc*vec2(0,8))/RENDERSIZE.xy).x;
    
    fragColor.y =texture(BuffA,(fragCoord + sc*vec2(0,-1))/RENDERSIZE.xy).y;
    
    fragColor.z =texture(BuffA,(fragCoord + sc*vec2(0,-4))/RENDERSIZE.xy).z;
    
    
	return fragColor; 
 } 


// influenced by blackle's https://www.shadertoy.com/view/3dXfR2 !!


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    fragCoord -= 0.5*RENDERSIZE.xy;
    fragCoord *= 0.95;
    fragCoord += 0.5*RENDERSIZE.xy;
    
    float n1d = texelFetch(image30,ivec2(mod(fragCoord + vec2(float(FRAMECOUNT),0.),256.)),0).x;
    vec3 n  = texelFetch(image30,ivec2(mod(fragCoord  + n1d*200. ,256.)),0).xyz;
    
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    
    fragColor.xyz =texture(BuffB,fragCoord/RENDERSIZE.xy).xyz;
    
    fragColor.xyz = pow(fragColor.xyz, vec3(1.,1.2,1.2));
    
    //fragColor.xyz = 1. - fragColor.xyz;
    
    //fragColor.xyz *= 1. - dot(uv,uv)*0.5;
    fragColor.xyz = pow(fragColor.xyz, vec3(0.4545 - n*0.15));
    
    
    
    
    
    fragColor.xyz += smoothstep(1.,0.,length(fragColor))*n*0.35;
    
    fragColor.xyz -= smoothstep(0.,1.,length(fragColor))*n*0.15;
    
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
		return renderMainImage();
	}
}