

			//******** BuffA Code Begins ********

#define A(COORD) texture(BuffA,(COORD)/RENDERSIZE.xy)

vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    U-=.5*RENDERSIZE.xy;
    U *= .997;
    float a = .003*sin(.5*smoothTime-2.*length(U-0.5*RENDERSIZE.xy)/RENDERSIZE.y);
    U *= mat2(cos(a),-sin(a),sin(a),cos(a));
    U+=.5*RENDERSIZE.xy;
    Q  =  A(U);
    // Neighborhood :
    vec4 pX  =  A(U + vec2(1,0));
    vec4 pY  =  A(U + vec2(0,1));
    vec4 nX  =  A(U - vec2(1,0));
    vec4 nY  =  A(U - vec2(0,1));
    vec4 m = 0.25*(pX+nX+pY+nY);
    float b = mix(1.,abs(Q.z),.8);
    Q.xyz += (1.-b)*(0.25*vec3(pX.z-nX.z,pY.z-nY.z,-pX.x+nX.x-pY.y+nY.y)- Q.xyz);

    
    Q = mix(Q,m,b);
    
    if (length(Q.xy)>0.) Q.xy = normalize(Q.xy);
    
    if(FRAMECOUNT <= 1) Q = sin(.01*length(U-0.5*RENDERSIZE.xy)*vec4(1,2,3,4)+basshits);
    
    if (_mouse.z>0.&&length(U-_mouse.xy)<.1*RENDERSIZE.y) Q *= 0.;
    
	return Q; 
 } 


#define A(COORD) texture(BuffA,(COORD)/RENDERSIZE.xy)
float ln (vec3 p, vec3 a, vec3 b) {return length(p-a-(b-a)*min(dot(p-a,b-a),0.)/dot(b-a,b-a));}
vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q  =  A(U);
    vec4 pX  =  A(U + vec2(1,0));
    vec4 pY  =  A(U + vec2(0,1));
    vec4 nX  =  A(U - vec2(1,0));
    vec4 nY  =  A(U - vec2(0,1));
    vec3 n = normalize(vec3(pX.z-nX.z,pY.z-nY.z,1));
    vec3 r = reflect(n,vec3(0,0,-1));
    Q = (0.5+0.5*sin(smoothTimeC+atan(Q.x,Q.y)*vec4(3,2,1,4)));
    float d = ln(vec3(.4,.4,6)*RENDERSIZE.xyy,
                 vec3(U,0),vec3(U,0)+r)/RENDERSIZE.y;
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