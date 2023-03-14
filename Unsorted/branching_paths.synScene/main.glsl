

			//******** BuffA Code Begins ********

#define R RENDERSIZE.xy
#define D(U) texture(BuffD,(U)/R)
vec4 dX (inout vec2 c,inout float m, vec2 u, vec2 r) {
    vec4 n = D(u+r);
    m += n.z;
	if (length(u+r-n.xy)<length(u-c)) c = n.xy;
    return n;
} 
vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = D(U);
    vec2 c = Q.xy;
    float m=0.;
    vec4 
        n = dX(c,m,U,+vec2(0,1)),
        e = dX(c,m,U,+vec2(1,0)),
        s = dX(c,m,U,-vec2(0,1)),
        w = dX(c,m,U,-vec2(1,0));
    dX(c,m,U,vec2(1,1));
    dX(c,m,U,vec2(1,-1));
    dX(c,m,U,vec2(-1,1));
    dX(c,m,U,vec2(-1,-1));
    vec2 g = vec2(e.z-w.z,n.z-s.z);
    Q.xy = c;
    Q.z += (m/8.-Q).z+.05*Q.w - .0001*Q.z;
    Q.w -= 0.001*Q.w;
    Q.zw = max(Q.zw,vec2(2,1)*smoothstep(4.,0.0,length(U-c)));
    Q.xy -= (0.25+syn_BassLevel*0.25)*g;
    if (length(U-_mouse.xy)<.1*R.y&&_mouse.z>0.) Q = vec4(-R,Q.z,0);
    if (FRAMECOUNT <= 1) Q = vec4(clamp(floor(U/2.)*2.,0.5*R-2.,0.5*R+2.),0,0);
    
	return Q; 
 } 


			//******** BuffB Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
vec4 aX (inout vec2 c,inout float m, vec2 u, vec2 r) {
    vec4 n = A(u+r);
    m += n.z;
	if (length(u+r-n.xy)<length(u-c)) c = n.xy;
    return n;
} 
vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = A(U);
    vec2 c = Q.xy;
    float m=0.;
    vec4 
        n = aX(c,m,U,+vec2(0,1)),
        e = aX(c,m,U,+vec2(1,0)),
        s = aX(c,m,U,-vec2(0,1)),
        w = aX(c,m,U,-vec2(1,0));
    aX(c,m,U,vec2(1,1));
    aX(c,m,U,vec2(1,-1));
    aX(c,m,U,vec2(-1,1));
    aX(c,m,U,vec2(-1,-1));
    vec2 g = vec2(e.z-w.z,n.z-s.z);
    Q.xy = c;
    Q.z += (m/8.-Q).z+.05*Q.w - .0001*Q.z;
    Q.w -= 0.001*Q.w;
    Q.zw = max(Q.zw,vec2(2,1)*smoothstep(4.,0.0,length(U-c)));
    Q.xy -= 0.25*g;
    if (length(U-_mouse.xy)<.1*R.y&&_mouse.z>0.) Q = vec4(-R,Q.z,0);
    if (FRAMECOUNT <= 1) Q = vec4(clamp(floor(U/2.)*2.,0.5*R-2.,0.5*R+2.),0,0);
    
	return Q; 
 } 


			//******** BuffC Code Begins ********

#define R RENDERSIZE.xy
#define B(U) texture(BuffB,(U)/R)
vec4 bX (inout vec2 c,inout float m, vec2 u, vec2 r) {
    vec4 n = B(u+r);
    m += n.z;
	if (length(u+r-n.xy)<length(u-c)) c = n.xy;
    return n;
} 
vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = B(U);
    vec2 c = Q.xy;
    float m=0.;
    vec4 
        n = bX(c,m,U,+vec2(0,1)),
        e = bX(c,m,U,+vec2(1,0)),
        s = bX(c,m,U,-vec2(0,1)),
        w = bX(c,m,U,-vec2(1,0));
    bX(c,m,U,vec2(1,1));
    bX(c,m,U,vec2(1,-1));
    bX(c,m,U,vec2(-1,1));
    bX(c,m,U,vec2(-1,-1));
    vec2 g = vec2(e.z-w.z,n.z-s.z);
    Q.xy = c;
    Q.z += (m/8.-Q).z+.05*Q.w - .0001*Q.z;
    Q.w -= 0.001*Q.w;
    Q.zw = max(Q.zw,vec2(2,1)*smoothstep(4.,0.0,length(U-c)));
    Q.xy -= 0.25*g;
    if (length(U-_mouse.xy)<.1*R.y&&_mouse.z>0.) Q = vec4(-R,Q.z,0);
    if (FRAMECOUNT <= 1) Q = vec4(clamp(floor(U/2.)*2.,0.5*R-2.,0.5*R+2.),0,0);
    
	return Q; 
 } 


			//******** BuffD Code Begins ********

#define R RENDERSIZE.xy
#define C(U) texture(BuffC,(U)/R)
vec4 cX (inout vec2 c,inout float m, vec2 u, vec2 r) {
    vec4 n = C(u+r);
    m += n.z;
	if (length(u+r-n.xy)<length(u-c)) c = n.xy;
    return n;
} 
vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = C(U);
    vec2 c = Q.xy;
    float m=0.;
    vec4 
        n = cX(c,m,U,+vec2(0,1)),
        e = cX(c,m,U,+vec2(1,0)),
        s = cX(c,m,U,-vec2(0,1)),
        w = cX(c,m,U,-vec2(1,0));
    cX(c,m,U,vec2(1,1));
    cX(c,m,U,vec2(1,-1));
    cX(c,m,U,vec2(-1,1));
    cX(c,m,U,vec2(-1,-1));
    vec2 g = vec2(e.z-w.z,n.z-s.z);
    Q.xy = c;
    Q.z += (m/8.-Q).z+.05*Q.w - .0001*Q.z;
    Q.w -= 0.001*Q.w;
    Q.zw = max(Q.zw,vec2(2,1)*smoothstep(4.,0.0,length(U-c)));
    Q.xy -= 0.25*g;
    if (length(U-_mouse.xy)<.1*R.y&&_mouse.z>0.) Q = vec4(-R,Q.z,0);
    if (FRAMECOUNT <= 1) Q = vec4(clamp(floor(U/2.)*2.,0.5*R-2.,0.5*R+2.),0,0);
    
	return Q; 
 } 


#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
vec4 X (inout vec2 c, vec2 u, vec2 r) {
    vec4 n = A(u+r);
	if (length(u+r-n.xy)<length(u-c)) c = n.xy;
    return n;
} 
vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    vec4 a = A(U);
    vec2 c = Q.xy;
    vec4 
        n = X(c,U,+vec2(0,1)),
        e = X(c,U,+vec2(1,0)),
        s = X(c,U,-vec2(0,1)),
        w = X(c,U,-vec2(1,0)),
        m = 10.25*(n+e+s+w);
 	vec3 g = normalize(vec3(e.z-w.z,n.z-s.z,.3));
    g = reflect(g,vec3(0,0,1));
 	vec3 b = normalize(vec3(e.w-w.w,n.w-s.w,1));
    float d = dot(g,normalize(vec3(0,1,.5)));
    Q = (exp(-4.*d*d))*m.w*abs(sin(2.+vec4(1+_uvc.y,2,3,4)*(1.+2.*m.z)))*3;

    Q = .0+.2*g.x-1.8*(1.+0.5*(b.x+b.y)) * a.w * sin(1.75+0.75*(g.z)*vec4(1+_uvc.x,2,3,8)*1.75);
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