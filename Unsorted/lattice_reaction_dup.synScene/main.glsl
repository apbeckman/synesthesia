

			//******** BuffA Code Begins ********
float growthFactor = 0.5+pow((syn_BassLevel*0.25)+(syn_MidLevel*0.35)+(syn_Intensity*0.15), 2.0);
vec4 image = texture(syn_UserImage,(_xy)/RENDERSIZE.xy);

#define R RENDERSIZE.xy
#define A(U) texture(BuffC, (U)/R)
vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    U += _uvc*PI*Succ*growthFactor;
    U += MoveXY;
    U += _uvc*Stretch*PI;
    
    for (int i = 0; i < 10; i++) {
        vec4 a = A(U);
        U -= a.x*a.zw;
    }
    Q = A(U);
    //Q *=(1.0 + 0.025* _loadUserImageAsMask());
    Q += Q*0.025* _loadUserImageAsMask()*(1.0+syn_BassLevel);
    
    vec4 
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0)),
        a = A(U+vec2(1,1)),
        b = A(U+vec2(1,-1)),
        c = A(U-vec2(1,1)),
        d = A(U-vec2(1,-1)),
        dQ = 0.125*(n+e+s+w+a+b+c+d)-Q*(1.0+highhits*0.125);
    Q = A(U);
    Q += vec4(0.125,1,1,1)*dQ*(1.0+midhits);
    float x = .1*Q.y*Q.x*(1.-Q.x);
    Q.x = Q.x+x-0.00+(e.z*e.x-w.z*w.x+n.w*n.x-s.w*s.x)*(1.0+growthFactor);
    Q.y = Q.y-x+0.04*Q.y*(1.-Q.y);
    Q.xy = max(Q.xy,0.);
    Q.zw = Q.zw - .025*vec2((e.y)-(w.y),(n.y)-(s.y));
    if (U.x < 1. || U.y < 1. || R.x-U.x<1. || R.y-U.y<1.)Q.x*=0.;
    //if (_mouse.z>0.&&length(_mouse.xy-U)<40.)Q*=0.;
    if (FRAMECOUNT <= 1 || Reset != 0.) Q = vec4(exp(-.1*length(U-0.5*R*0.1)),1,0,0);
  
    
	return Q; 
 } 


			//******** BuffB Code Begins ********
/*
#define R RENDERSIZE.xy
*/
#define AA(U) texture(BuffA, (U)/R)
vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    for (int i = 0; i <12; i++) {
        vec4 a = AA(U);
        U -= a.x*a.zw;
    }
    Q = AA(U);
    vec4 
        n = AA(U+vec2(0,1)),
        e = AA(U+vec2(1,0)),
        s = AA(U-vec2(0,1)),
        w = AA(U-vec2(1,0)),
        a = AA(U+vec2(1,1)),
        b = AA(U+vec2(1,-1)),
        c = AA(U-vec2(1,1)),
        d = AA(U-vec2(1,-1)),
        dQ = 0.125*(n+e+s+w+a+b+c+d)-Q;
    Q = AA(U);
    Q += vec4(0.5,1,1,1)*dQ;
    float x = .1*Q.y*Q.x*(1.-Q.x);
    Q.x = Q.x+x-0.00+(e.z*e.x-w.z*w.x+n.w*n.x-s.w*s.x);
    Q.y = Q.y-x+0.04*Q.y*(1.-Q.y)*(1.0+basshits);
    Q.xy = max(Q.xy,0.);
    Q.zw = Q.zw - .025*vec2(e.y-w.y,n.y-s.y);
    if (U.x < 1. || U.y < 1. || R.x-U.x<1. || R.y-U.y<1.)Q.x*=0.;
    if (_mouse.z>0.&&length(_mouse.xy-U)<40.)Q*=0.2;
    if (FRAMECOUNT <= 1 || Reset != 0.) Q = vec4(exp(-.1*length(U-0.5*R)),1,0,0);
  
    
	return Q; 
 } 


			//******** BuffC Code Begins ********

#define R RENDERSIZE.xy
#define AB(U) texture(BuffB, (U)/R)
vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    U += _uvc*PI*Succ*growthFactor;

    for (int i = 0; i < 15; i++) {
        vec4 a = AB(U);
        U -= a.x*a.zw;
    }
    Q = AB(U);
    vec4 
        n = AB(U+vec2(0,1)),
        e = AB(U+vec2(1,0)),
        s = AB(U-vec2(0,1)),
        w = AB(U-vec2(1,0)),
        a = AB(U+vec2(1,1)),
        b = AB(U+vec2(1,-1)),
        c = AB(U-vec2(1,1)),
        d = AB(U-vec2(1,-1)),
        dQ = 0.125*(n+e+s+w+a+b+c+d)-Q;
    Q = AB(U);
    Q += vec4(0.5,1,1,1)*dQ;
    float x = .1*Q.y*Q.x*(1.-Q.x);
    Q.x = Q.x+x-0.00+(e.z*e.x-w.z*w.x+n.w*n.x-s.w*s.x);
    Q.y = Q.y-x+0.04*Q.y*(1.-Q.y);
    Q.xy = max(Q.xy,0.);
    Q.zw = Q.zw - .025*vec2(e.y-w.y,n.y-s.y);
    if (U.x < 1. || U.y < 1. || R.x-U.x<1. || R.y-U.y<1.)Q.x*=0.;
    if (_mouse.z>0.&&length(_mouse.xy-U)<40.)Q*=0.2;
    if (FRAMECOUNT <= 1 || Reset != 0.) Q = vec4(exp(-.1*length(U-0.5*R)),1,0,0);
  
    
	return Q; 
 } 


			//******** BuffD Code Begins ********

//#define R RENDERSIZE.xy
float A1 (vec2 U) {
	return -10.*texture(BuffA, U/R).x;
}
float ln (vec3 p, vec3 a, vec3 b) {return length(p-a-(b-a)*dot(p-a,b-a)/dot(b-a,b-a));}

vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    U += _uvc*PI*Succ*growthFactor;
    
   vec3 p = vec3(1.5*R,-R.x);
   vec3 d = normalize(p-vec3(U,0));
   vec3 li = vec3(1.3*R,1.3*R.x);
   li += vec3(sin(smoothTimeB), cos(smoothTimeB), 0.);
   p += d*dot(-p,vec3(0,0,1))/dot(d,vec3(0,0,1));
   for (int i = 0; i < 22; i++) {
   	p += .6*d*(p.z-A1(p.xy));
   }
   vec3 q = li;
   q.xy += _uvc*PI*PI;
   vec3 c = normalize(p-li);
   for (int i = 0; i < 30; i++) {
    q += .6*c*(q.z-A1(q.xy));
   }
    U = p.xy;
    float 
        n = A1(U+vec2(0,1)),
        e = A1(U+vec2(1,0)),
        s = A1(U-vec2(0,1)),
        w = A1(U-vec2(1,0));
   float a = abs(A1(U));
   vec3 g = normalize(vec3(e-w,n-s,-1));
   vec3 r = reflect(d,g);
   Q = (0.6+0.124*sin((.2*smoothTimeB)+(a*0.5)*(1.+.1*vec4(1,2,3,4))));
    float o = ln( li, p, p+r );
    float len = length(q-p);
    float h = 0.85+0.5*dot(normalize(p-li),g);
   Q *= h*(exp(-.04*o)*6.+0.3*exp(-.00125*o))*(exp(-.13*len)*(0.5+1.0*syn_HighLevel*0.8+0.2*syn_Intensity));
	return Q; 
 } 


#define R RENDERSIZE.xy
#define D(U) texture(BuffD, (U)/R)
vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    U += _uvc*PI*Succ;
	Q = (4.*(1+highhits*0.125))*D(U);
    float d = 0.1;
    float a = -0.1, si = sin(a), co = cos(a);
    mat2 m = mat2(co,-si,si,co);
    vec2 
        n = m*vec2(0,1),
        e = m*vec2(1,0),
        s = m*vec2(0,-1),
        w = m*vec2(-1,0);
    for (float i = 0.; i < 20.; i++) {
    	Q += d*D(U+i*n)*exp(-.2*i);
    	Q += d*D(U-i*e)*exp(-.2*i);
    	Q += d*D(U+i*s)*exp(-.2*i);
    	Q += d*D(U-i*w)*exp(-.2*i);
    }
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