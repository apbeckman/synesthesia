

			//******** BuffA Code Begins ********

#define A(COORD) texture(BuffA,(COORD)/RENDERSIZE.xy)

float growthFactor =0.25+2*normalize(pow((syn_BassLevel*0.5)+(syn_MidLevel*0.25)+(syn_Intensity*0.25), 2.0));
//vec4 image = texture(syn_UserImage,(_xy)/RENDERSIZE.xy);
vec4 image = vec4(0.);
vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = (_xy);
  U-= _uvc*Zoom*(2*growthFactor);

  vec3 logoCol = vec3(0.0);
  if (_exists(syn_UserImage)){
    logoCol = _loadUserImageAsMask().rgb;
  }
    U += _uvc*Stretch*PI*2.;


    U-=.5*RENDERSIZE.xy;
    U *= .99375+(growthFactor*0.00625*Zoom);

    float a = .003*sin(.125*smoothTimeC-2.*length(U-0.5*RENDERSIZE.xy)/RENDERSIZE.y)*Spin;
    U *= mat2(cos(a),-sin(a),sin(a),cos(a));
    U+=.5*RENDERSIZE.xy;
    U.xy += Dir_XY.xy*growthFactor*0.5;
    Q  =  A(U);
    image *= MediaImpact;
    image =vec4(logoCol, 0.);

    // Neighborhood :
    vec4 pX  =  A(U + vec2(1,0))+(image*(0.975+(0.25*growthFactor)))*0.375;
    vec4 pY  =  A(U + vec2(0,1))+(image*(0.975+(0.25*growthFactor)))*0.375;
    vec4 nX  =  A(U - vec2(1,0))-(image*(0.975+(0.25*growthFactor)))*0.375;
    vec4 nY  =  A(U - vec2(0,1))-(image*(0.975+(0.25*growthFactor)))*0.375;
    vec4 m = 0.25*(pX+nX+pY+nY+Intensity);
    float b = mix(1.,abs(Q.z),.78);
    
    Q.xyz += (1.-b+0.0125*basshits)*(0.25*vec3(pX.z-nX.z,pY.z-nY.z,-pX.x+nX.x-pY.y+nY.y)- Q.xyz)*(Diffusion+growthFactor*0.125);

    
    Q = mix(Q,m,b);
    
    if (length(Q.xy)>0. || Reset > 0.) Q.xy = normalize(Q.xy);
    
    if(FRAMECOUNT <= 1 || Reset > 0.) Q = sin(.01*length(U-(0.5*RENDERSIZE.xy+PI*0.25*_uvc))*vec4(1,2,3,4));
    
    if (_mouse.z>0.&&length(U-_uvc-_mouse.xy)<(RENDERSIZE.y)*0.125) Q *= 0.;
    
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
    //U-= _uvc*Zoom*(2*growthFactor);

    Q  =  A1(U);
    //U+= _uvc*PI*0.1;

    vec4 pX  =  A1(U + vec2(1,0));
    pX.xy+= _uvc*PI;

    vec4 pY  =  A1(U + vec2(0,1));
    vec4 nX  =  A1(U - vec2(1,0));
    vec4 nY  =  A1(U - vec2(0,1));
    vec3 n = normalize(vec3(pX.z-nX.z,pY.z-nY.z,1));
    vec3 r = reflect(n,vec3(0,0,-2));
    Q = (0.25+0.75*sin(smoothTimeB*0.25+atan(Q.x,Q.y)*vec4(3,2,1,4)));
    float d = ln(vec3(.4,.4,16)*RENDERSIZE.xyy+_uvc.xyy,
                 vec3(U,0),vec3(U,0)+r)/RENDERSIZE.y;
    Q *= exp(-d*d)*.5+.5*exp(-3.*d*d)+(syn_Intensity*pow(syn_HighLevel*0.4 + syn_MidHighLevel*0.4 + syn_Intensity*0.2, 2.))*Flash;
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