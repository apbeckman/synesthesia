

			//******** BuffA Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA, (U)/R)
#define B(U) texture(BuffB, (U)/R)
#define C(U) texture(BuffC, (U)/R)
#define D(U) texture(BuffD, (U)/R)

vec4 mediaEdges = _edgeDetectSobel(syn_UserImage, _uv);
vec4 mixedEdges = mix( _loadMedia()*0.7+0.3*mediaEdges, _loadMedia(), 1.0-edge_mix);
bool mediaOn = _exists(syn_UserImage);
vec4 mediaMix(vec4 Q, vec4 image) {
    image = mix(Q, image, int(mediaOn)*0.01*syn_BassLevel);
    return mix(Q, image, 0.0125);
}
float AA (vec2 U) {
	return -10.*texture(BuffA, U/R).x;
}
vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    U += _uvc*stretch;
    for (int i = 0; i < 5; i++) {
        vec4 a = C(U);
        U -= a.x*a.zw;
    }
    Q = C(U);

    vec4 
        n = C(U+vec2(0,1)),
        e = C(U+vec2(1,0)),
        s = C(U-vec2(0,1)),
        w = C(U-vec2(1,0)),
        a = C(U+vec2(1,1)),
        b = C(U+vec2(1,-1)),
        c = C(U-vec2(1,1)),
        d = C(U-vec2(1,-1)),
        dQ = 0.125*(n+e+s+w+a+b+c+d)-Q;
    Q = C(U);
    
    Q += vec4(0.5,1,1,1)*dQ;\

    float x = .1*Q.y*Q.x*(1.-Q.x);
    Q.x = Q.x+x-0.00+(e.z*e.x-w.z*w.x+n.w*n.x-s.w*s.x);
    Q.y = Q.y-x+0.04*Q.y*(1.-Q.y);
    Q.xy = max(Q.xy,0.);
    Q.zw = Q.zw - .025*vec2(e.y-w.y,n.y-s.y);
    if (U.x < 1. || U.y < 1. || R.x-U.x<1. || R.y-U.y<1.)Q.x*=0.;
    if (_mouse.z>0.&&length(_mouse.xy-U)<40.)Q*=0.;
    if (FRAMECOUNT <= 1 || Reset != 0) Q = vec4(exp(-.1*length(U-0.5*R)),1,0,0);
    
	return Q; 
 } 


			//******** BuffB Code Begins ********

#define R RENDERSIZE.xy
vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    for (int i = 0; i < 5; i++) {
        vec4 a = A(U);
        U -= a.x*a.zw;
    }
    Q = A(U);
    vec4 
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0)),
        a = A(U+vec2(1,1)),
        b = A(U+vec2(1,-1)),
        c = A(U-vec2(1,1)),
        d = A(U-vec2(1,-1)),
        dQ = 0.125*(n+e+s+w+a+b+c+d)-Q;
    Q = A(U);
    Q += vec4(0.5,1,1,1)*dQ;
    float x = .1*Q.y*Q.x*(1.-Q.x);
    Q.x = Q.x+x-0.00+(e.z*e.x-w.z*w.x+n.w*n.x-s.w*s.x);
    Q.x += mix((0.0), sin(length(mediaEdges))*PI, 0.3*(media_impact*0.2+0.5*syn_BassLevel*media_impact));
    Q.y = Q.y-x+0.04*Q.y*(1.-Q.y);
    Q.xy = max(Q.xy,0.);
    Q.zw = Q.zw - .025*vec2(e.y-w.y,n.y-s.y);
    if (U.x < 1. || U.y < 1. || R.x-U.x<1. || R.y-U.y<1.)Q.x*=0.;
    if (_mouse.z>0.&&length(_mouse.xy-U)<40.)Q*=0.;
    if (FRAMECOUNT <= 1 || Reset != 0) Q = vec4(exp(-.1*length(U-0.5*R)),1,0,0);
  
    
	return Q; 
 } 


			//******** BuffC Code Begins ********

#define R RENDERSIZE.xy
vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    for (int i = 0; i < 5; i++) {
        vec4 a = B(U);
        U -= a.x*a.zw;
    }
    Q = B(U);
    vec4 
        n = B(U+vec2(0,1)),
        e = B(U+vec2(1,0)),
        s = B(U-vec2(0,1)),
        w = B(U-vec2(1,0)),
        a = B(U+vec2(1,1)),
        b = B(U+vec2(1,-1)),
        c = B(U-vec2(1,1)),
        d = B(U-vec2(1,-1)),
        dQ = 0.125*(n+e+s+w+a+b+c+d)-Q;
    Q = B(U);
    Q += vec4(0.5,1,1,1)*dQ;
    float x = .1*Q.y*Q.x*(1.-Q.x);
    Q.x = Q.x+x-0.00+(e.z*e.x-w.z*w.x+n.w*n.x-s.w*s.x);
    Q.y = Q.y-x+0.04*Q.y*(1.-Q.y);
    Q.xy = max(Q.xy,0.);
    Q.zw = Q.zw - .025*vec2(e.y-w.y,n.y-s.y);
    if (U.x < 1. || U.y < 1. || R.x-U.x<1. || R.y-U.y<1.)Q.x*=0.;
    if (_mouse.z>0.&&length(_mouse.xy-U)<40.)Q*=0.;
    if (FRAMECOUNT <= 1 || Reset != 0) Q = vec4(exp(-.1*length(U-0.5*R)),1,0,0);
  
    
	return Q; 
 } 


			//******** BuffD Code Begins ********


float ln (vec3 p, vec3 a, vec3 b) {return length(p-a-(b-a)*dot(p-a,b-a)/dot(b-a,b-a));}

vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

   vec3 p = vec3(.5*R,-R.x);
   vec3 d = normalize(p-vec3(U,0));
   vec3 li = vec3(1.3*R,.3*R.x);
   p += d*dot(-p,vec3(0,0,1))/dot(d,vec3(0,0,1));
   for (int i = 0; i < 10; i++) {
   	p += .6*d*(p.z-AA(p.xy));
   }
   vec3 q = li;
   vec3 c = normalize(p-li);
   for (int i = 0; i < 40; i++) {
    q += .6*c*(q.z-AA(q.xy));
   }
    U = p.xy;
    float 
        n = AA(U+vec2(0,1)),
        e = AA(U+vec2(1,0)),
        s = AA(U-vec2(0,1)),
        w = AA(U-vec2(1,0));
   float a = AA(U);
   vec3 g = normalize(vec3(e-w,n-s,-1));
   vec3 r = reflect(d,g);
   Q = (0.6+0.4*sin((.2*TIME)+(a)*(1.+.1*vec4(1,2,3,4))));
    float o = ln( li, p, p+r );
    float len = length(q-p);
    float h = 0.5+0.5*dot(normalize(p-li),g);
   Q *= h*(exp(-.04*o)*6.+0.3*exp(-.001*o))*(exp(-.3*len));
	return Q; 
 } 


#define R RENDERSIZE.xy
vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

	   Q = 3.*D(U);
    float d = 0.1;
    float a = -0.1, si = sin(a), co = cos(a);
    mat2 m = mat2(co,-si,si,co);
    vec2 
        n = m*vec2(0,1),
        e = m*vec2(1,0),
        s = m*vec2(0,-1),
        w = m*vec2(-1,0);
    for (float i = 0.; i < 15.; i++) {
    	Q += d*D(U+i*n)*exp(-.2*i);
    	Q += d*D(U-i*e)*exp(-.2*i);
    	Q += d*D(U+i*s)*exp(-.2*i);
    	Q += d*D(U-i*w)*exp(-.2*i);
    }
    Q = mix(Q*1.25, mixedEdges, 0.25*int(mediaOn)*media_impact*media_color_mix) ;
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