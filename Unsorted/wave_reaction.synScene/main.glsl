

			//******** BuffA Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA, (U)/R)
vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    for (int i = 0; i < 9; i++) {
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
    Q += vec4(01.1,1,1,1)*dQ; // 0.1 to 1.1
    float x = .1*Q.y*Q.x*(1.-Q.x);
    Q.x = Q.x+x-0.004+(e.z*e.x-w.z*w.x+n.w*n.x-s.w*s.x)*(1.0+0.5 *syn_BassLevel);
    Q.y = Q.y-x+0.02*Q.y*(1.-Q.y);
    Q.xy = max(Q.xy,0.);
    Q.zw = Q.zw +.125*vec2(e.y-w.y,n.y-s.y);
    if (_mouse.z>0.&&length(_mouse.xy-U)<10.)Q*=0.;
    if (int(FRAMECOUNT) <= 1. || Reset != 0.) Q = vec4(exp(-.1*length(U-0.5*R)),1,0,0);
  
    
	return Q; 
 } 


#define R RENDERSIZE.xy
vec4 A1 (vec2 U) {return texture(BuffA,U/R);}
float ln (vec3 p, vec3 a, vec3 b) {return length(p-a-(b-a)*dot(p-a,b-a)/dot(b-a,b-a));}
float X (vec3 a, vec3 s, vec3 n) {
    vec3 light = vec3(0.5*R,R.y);
    vec3 r = reflect(light-s,n);
	return 3.*exp(-.5*ln(a,s,s+r));
}
vec4 renderMainImage() {
	vec4 Q = vec4(0.0);

	vec2 U = _xy;

    vec3 p = vec3(U,0),n,u;

    for (int x = -3; x <= 3; x++)
    for (int y = -3; y <= 3; y++) {
    //p.xy += sin(_uvc*PI*cos (bass_time*_uvc.x)*0.1);

       u = p + vec3(x,y,15);
       n = normalize(vec3(A1(u.xy).zw,.1));
       Q.x += X(p,u,n);
       n = normalize(vec3(A1(u.xy).zw,.2));
       Q.y += X(p,u,n);
       n = normalize(vec3(A1(u.xy).zw,.3));
       Q.z += X(p,u,n);
    }
    Q = Q/7./7.;
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