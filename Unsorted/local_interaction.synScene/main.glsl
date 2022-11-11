

			//******** BuffA Code Begins ********

#define A(COORD) texture(BuffA,(COORD)/RENDERSIZE.xy)
float growthFactor = normalize(pow((syn_BassLevel*0.5)+(syn_MidLevel*0.35)+(syn_Level*0.15), 2.0));
vec4 image = texture(syn_UserImage,(_xy)/RENDERSIZE.xy);
vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = (_xy);

    U-=.5*RENDERSIZE.xy;
    U *= (.997-growthFactor*Zoom*0.125)-Zoom;
    float a = .003*sin(.125*smoothTimeC-2.*length(U-0.5*RENDERSIZE.xy*Zoom*_uvc)/RENDERSIZE.y)*Spin;
    U *= mat2(cos(a),-sin(a),sin(a),cos(a));
    U+=.5*RENDERSIZE.xy;
    U.xy += Dir_XY.xy;
    Q  =  A(U);
    // Neighborhood :
    vec4 pX  =  A(U + vec2(1,0))+(image*growthFactor)*0.25;
    vec4 pY  =  A(U + vec2(0,1))+(image*growthFactor)*0.25;
    vec4 nX  =  A(U - vec2(1,0))-(image*growthFactor)*0.25;
    vec4 nY  =  A(U - vec2(0,1))-(image*growthFactor)*0.25;
    vec4 m = 0.25*(pX+nX+pY+nY+Test);
    float b = mix(1.,abs(Q.z),.8);
    
    Q.xyz += (1.-b+0.0125*basshits)*(0.25*vec3(pX.z-nX.z,pY.z-nY.z,-pX.x+nX.x-pY.y+nY.y)- Q.xyz)*(Diffusion+growthFactor*0.125);

    
    Q = mix(Q,m,b);
    
    if (length(Q.xy)>0. || Reset > 0.) Q.xy = normalize(Q.xy);
    
    if(FRAMECOUNT <= 1 || Reset > 0.) Q = sin(.01*length(U-(0.5*RENDERSIZE.xy))*vec4(1,2,3,4));
    
    if (_mouse.z>0.&&length(U-_mouse.xy)<(.09+basshits*0.1)*normalize(RENDERSIZE.y)*0.5) Q *= 0.;
    
	return Q; 
 } 
/*bool Media = bool(media);
#ifdef Media
#define A1(COORD) texture(syn_UserImage,(COORD)/RENDERSIZE.xy)
#else*/
#define A1(COORD) texture(BuffA,(COORD)/RENDERSIZE.xy)
//  #endif
float ln (vec3 p, vec3 a, vec3 b) {return length(p-a-(b-a)*min(dot(p-a,b-a),0.)/dot(b-a,b-a));}
vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q  =  A1(U);
    vec4 pX  =  A1(U + vec2(1,0));
    vec4 pY  =  A1(U + vec2(0,1));
    vec4 nX  =  A1(U - vec2(1,0));
    vec4 nY  =  A1(U - vec2(0,1));
    vec3 n = normalize(vec3(pX.z-nX.z,pY.z-nY.z,1));
    vec3 r = reflect(n,vec3(0,0,-1));
    Q = (0.5+0.5*sin(smoothTimeB+atan(Q.x,Q.y)*vec4(3,2,1,4)));
    float d = ln(vec3(.4,.4,6)*RENDERSIZE.xyy,
                 vec3(U,0),vec3(U,0)+r)/RENDERSIZE.y-(0.25*pow(syn_HighLevel*0.4 + syn_MidHighLevel*0.4 + syn_Hits*0.2, 2.));
    Q *= exp(-d*d)*.5+.5*exp(-3.*d*d);
	return Q; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}