

			//******** BuffA Code Begins ********
#define A(COORD) texture(BuffA,(COORD)/RENDERSIZE.xy)
#define B(COORD) texture(BuffB,(COORD)/RENDERSIZE.xy)
#define C(COORD) texture(BuffC,(COORD)/RENDERSIZE.xy)
#define D(COORD) texture(BuffD,(COORD)/RENDERSIZE.xy)
vec4 Q = vec4(0.0);
vec2 U = _xy;
vec2 mouse = (_mouse.xy*_uvc);
float a = .001*sin(.005*smoothTime)/(1.+length(U-0.5*RENDERSIZE.xy*_uvc)/RENDERSIZE.y);
float b = .225;
vec2 aa = _uvc*growthFactor;
vec4 image = texture(syn_UserImage,(_xy)/RENDERSIZE.xy);

vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    U+=(Zoom*growthFactor*0.00125);

    vec2 X = 0.5*RENDERSIZE.xy*moveXY;
    //X += moveXY;
    if (_mouse.z>0.) X = _mouse.xy;
    U-=X;
    U *= mat2(cos(a),-sin(a),sin(a),cos(a));
    U *= .999;
    U+=X;
    Q  =  D(U);
    // Neighborhood :
    vec4 pX  =  D(U + vec2(1,0)+aa);
    vec4 pY  =  D(U + vec2(0,1)+aa);
    vec4 nX  =  D(U - vec2(1,0)+aa);
    vec4 nY  =  D(U - vec2(0,1)+aa);
    vec4 m = (0.25)*(pX+nX+pY+nY);
    //float b = .3*(1.0-growthFactor*0.12);
    Q += (1.-b)*(0.25*vec4(
       
        (pX.z-nX.z) - 5.*Q.z*(pY.z-nY.z),
       	(pY.z-nY.z) + 5.*Q.z*(pX.z-nX.z),
        -pX.x+nX.x-pY.y+nY.y,1
    )- Q);

    
    Q = mix(Q,m,b);
    
    if (length(Q.xy)>0.) Q.xy = normalize(Q.xy);
    
    if(FRAMECOUNT <= 1 || Reset > 0.) Q = sin(.1*U.xxxx);
    
     
	return Q; 
 } 


			//******** BuffB Code Begins ********


vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    U+=(Zoom*growthFactor*0.00125);

    vec2 X = 0.5*RENDERSIZE.xy*moveXY;
    if (_mouse.z>0.) X = _mouse.xy;
    U-=X-_uvc;
    //float a = .001*sin(.1*TIME)/(1.+length(U-0.5*RENDERSIZE.xy)/RENDERSIZE.y);
    U *= mat2(cos(a),-sin(a),sin(a),cos(a));
    //U *= .999;

    U+=X+_uvc;
    Q  =  A(U);
    // Neighborhood :
    vec4 pX  =  A(U + vec2(1,0))*image;
    vec4 pY  =  A(U + vec2(0,1))*image;
    vec4 nX  =  A(U - vec2(1,0))*image;
    vec4 nY  =  A(U - vec2(0,1))+growthFactor*image;
    vec4 m = 0.25*(pX+nX+pY+nY);
    //float b = .3;
    Q += (1.-b)*(0.25*vec4(
       
        (pX.z-nX.z) - 5.*Q.z*(pY.z-nY.z),
       	(pY.z-nY.z) + 5.*Q.z*(pX.z-nX.z),
        -pX.x+nX.x-pY.y+nY.y,1
    )- Q);

    
    Q = mix(Q,m,b)+0.6*growthFactor*image;
    
    if (length(Q.xy)>0.) Q.xy = normalize(Q.xy);
    
    if(FRAMECOUNT <= 1 || Reset > 0.) Q = sin(.1*U.xxxx);
    
     
	return Q; 
 } 


			//******** BuffC Code Begins ********


vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    U+=(Zoom*growthFactor*0.00125);

    vec2 X = 0.5*RENDERSIZE.xy*moveXY;
    if (_mouse.z>0.) X = _mouse.xy;
    U-=X;
    //float a = .001*sin(.1*TIME)/(1.+length(U-0.5*RENDERSIZE.xy)/RENDERSIZE.y);
    U *= mat2(cos(a),-sin(a),sin(a),cos(a));
    //U *= .999;

    U+=X;
    Q  =  B(U);
    // Neighborhood :
    vec4 pX  =  B(U + vec2(1,0))+growthFactor*image;
    vec4 pY  =  B(U + vec2(0,1))+growthFactor*image;
    vec4 nX  =  B(U - vec2(1,0))+growthFactor*image;
    vec4 nY  =  B(U - vec2(0,1))+growthFactor*image;
    vec4 m = 0.25*(pX+nX+pY+nY)+growthFactor*image;
    //float b = .3;
    Q += (1.-b)*(0.25*vec4(
       
        (pX.z-nX.z) - 5.*Q.z*(pY.z-nY.z),
       	(pY.z-nY.z) + 5.*Q.z*(pX.z-nX.z),
        -pX.x+nX.x-pY.y+nY.y,1
    )- Q);

    
    Q = mix(Q,m,b);
    
    if (length(Q.xy)>0.) Q.xy = normalize(Q.xy);
    
    if(FRAMECOUNT <= 1 || Reset > 0.) Q = sin(.1*U.xxxx);
    
     
	return Q; 
 } 


			//******** BuffD Code Begins ********


vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    vec2 X = 0.5*RENDERSIZE.xy*moveXY;
    if (_mouse.z>0.) X = _mouse.xy;
    U-=X;
    //float a = .001*sin(.1*TIME)/(1.+length(U-0.5*RENDERSIZE.xy)/RENDERSIZE.y);
    U *= mat2(cos(a),-sin(a),sin(a),cos(a));
    //U *= .999;
    U*=(1.0-(Zoom*growthFactor*0.005)-Zoom*0.00125); 
    U+=X;
    Q  =  C(U+_uvc);
    // Neighborhood :
    vec4 pX  =  C(U + vec2(1,0));
    vec4 pY  =  C(U + vec2(0,1));
    vec4 nX  =  C(U - vec2(1,0));
    vec4 nY  =  C(U - vec2(0,1));
    vec4 m = 0.25*(pX+nX+pY+nY);
    float b = .3;
    Q += (1.-b)*(0.25*vec4(
       
        (pX.z-nX.z) - 5.*Q.z*(pY.z-nY.z),
       	(pY.z-nY.z) + 5.*Q.z*(pX.z-nX.z),
        -pX.x+nX.x-pY.y+nY.y,1
    )- Q);

    
    Q = mix(Q,m,b);
    
    if (length(Q.xy)>0.) Q.xy = normalize(Q.xy);
    
    if(FRAMECOUNT <= 1 || Reset > 0.) Q = sin(.1*U.xxxx);
    
     
	return Q; 
 } 


//#define A(COORD) texture(BuffA,(COORD)/RENDERSIZE.xy)
float ln (vec3 p, vec3 a, vec3 b) {return length(p-a-(b-a)*min(dot(p-a,b-a),0.)/dot(b-a,b-a));}
vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    vec2 X = vec2(0.5);
    X+=vec2(cos(smoothTimeB*0.5)*0.125, sin(smoothTimeB*0.5)*0.125);

    if (_mouse.z>0.) X = _mouse.xy/RENDERSIZE.xy;
    Q  =  A(U+_uvc);
    vec4 pX  =  A(U + vec2(1,0)+aa);
    vec4 pY  =  A(U + vec2(0,1)+aa);
    vec4 nX  =  A(U - vec2(1,0)+aa);
    vec4 nY  =  A(U - vec2(0,1)+aa);
    vec3 n = normalize(vec3(pX.z-nX.z,pY.z-nY.z,1));
    vec3 r = reflect(n,vec3(0,0,-1));
    n,r *= 1.0-growthFactor*0.25;

    Q = 0.8+0.1*sin(19.*vec4(1,2,3,4)*(Q.z)+pow(highhits*0.75, 2.));
    float d = ln(vec3(X-.5,10)*RENDERSIZE.xyy,
                 vec3(U,0),vec3(U,0)+r)/RENDERSIZE.y;
    Q *= exp(-.5*d*d);
	return Q; 
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