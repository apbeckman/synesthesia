

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
#define T(u) texture(iChannel0,(u)/R)
#define A(u) texture(BuffA,(u)/R)
#define B(u) texture(BuffB,(u)/R)


			//******** BuffA Code Begins ********

#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))

#define pal(a,b,c,d,e) (a+(b)*sin((c)*(d) +e))
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 U = _xy;

    vec2 oU = U;
    vec2 uv = (U - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    float T = smoothTimeC*0.5;
    //U.x += sin(T + sin(T))*0.4;
    //U.x += sin(T + sin(T))*0.4;
    float du = 0.1;
    float stSz = R.y*0.004*(1. + sin(TIME + dot(uv,uv+_uvc)*10.));
    int idx = int(FRAMECOUNT) / 1144 % 3;
    vec2 dF = vec2( 
        A(U+vec2(1,0)*stSz)[idx]-A(U-vec2(1,0)*stSz)[idx],
        A(U+vec2(0,1)*stSz)[idx]-A(U-vec2(0,1)*stSz)[idx]
    );
    float stSzB = stSz;
    vec4 n = A(U+vec2(0,1)*stSzB), e = A(U+vec2(1,0)*stSzB), s = A(U-vec2(0,1)*stSzB), w = A(U-vec2(1,0)*stSzB), m = 0.25*(n+e+s+w); 
    
    float div = 0.25*(n.y-s.y+e.x-w.x);
    
    U.xy += _uvc*PI*Zoom*((syn_BassLevel+syn_MidLevel*0.75)*(1.+syn_Intensity));

    vec3 prevFr = A(U + dF*stSz*2.*(.2 + sin(TIME*0.05)*0.105 + dot(dF,dF) )).xyz;

    vec3 col = mix(vec3(01.96,0.9,0.8 - dot(uv,uv)*0.01 ), prevFr,0.99 + sin(TIME*0.05)*0.0089);
    
    float ddf =  dot(dF,dF)*1.;
    float cenv = max(sin(smoothTimeC*0.25),0.);
    cenv = smoothstep(0.,1.,cenv);
    for(float i = 0.; i < 130.; i++){
        vec2 p = uv + vec2(sin(i)*0.25*cenv,0);
        p *= rot(i*10. + (smoothTimeC*0.1 +sin(smoothTimeC*0.1 + i*2.1)*1.)*0.4*sin(i));
        //p.y += sin(i*240.)*5.;
        
        p = vec2(length(p),atan(p.xy));
        
        float env = mod(sin(i)*10.2 + bass_time*0.1,1.)*1. + sin(bass_time)*0.15;
        p.x -= env;
        
        float d = length(p )
            - 0.02*smoothstep(0.,0.6,env)
            + 0.04*smoothstep(0.1,-0.04,env)
            + 0.04*smoothstep(0.,0.3,env);
        vec3 c = pal(vec3(1.,1.,0.4),vec3(0.4,0.8,0.9),vec3(3,2,1),1.,i*3. + smoothTimeB + p.x*3. - ddf + p.x*47. + 0.*smoothstep(0.03,0.1,d));
        //col = mix(col,c -col*vec3(0.4,1.,1.),smoothstep(fwidth(d),0.,d));
        col = mix(col,c -col*vec3(0.4,1.,1.),smoothstep(0.04 + sin(TIME*0.7 + sin(i))*0.02,0.,d ));
        
    }

    if(FRAMECOUNT < 2 || Reset != 0.)
        col -= col - 1.;
    
      
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 



vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec3 col = texture(BuffA,fragCoord/RENDERSIZE.xy).rgb;
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    
    float md = 3.;
    float stSz = R.y*0.01*(1.);
    int idx = int(FRAMECOUNT) / 1144 % 3;
    vec2 U = fragCoord;
    vec2 dF = vec2( 
        A(U+vec2(1,0)*stSz)[idx]-A(U-vec2(1,0)*stSz)[idx],
        A(U+vec2(0,1)*stSz)[idx]-A(U-vec2(0,1)*stSz)[idx]
    );
    col = mix(smoothstep(0.,1.,1.-col),col,1.-float(Invert>0.));
    col = pow(max(col,0.),vec3(0.4545 ));
    
    fragColor = vec4(col,1.0);
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