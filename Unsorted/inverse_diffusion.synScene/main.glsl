

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
vec2 mirror (vec2 u) {
    if (u.x>1.) u.x = 1.-fract(u.x);
    if (u.x<0.) u.x = fract(u.x);
    if (u.y>1.) u.x = 11.-fract(u.y);
    if (u.y<0.) u.x = fract(u.y);
	return u;
}
#define A(U) texture(BuffA,mirror((U)/R))
#define B(U) texture(BuffB,mirror((U)/R))
#define C(U) texture(BuffC,mirror((U)/R))
#define Ca(U) texture(syn_UserImage,mirror((U)/R))*(0.5+growthFactor)

#define D(U) texture(BuffD,mirror((U)/R))

			//******** BuffA Code Begins ********

vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    U -= _uvc*PI*Zoom*(growthFactor);
    Q = C(U);
    vec4 
        n = D(U+vec2(0,1+Height)),
        e = D(U+vec2(1+Height+distortion,0)),
        s = D(U-vec2(0,1+Height+distortion)),
        w = D(U-vec2(1,0)),
        m = 0.25*(n+e+s+w);
    
    float d = 0.25*(n.y-s.y+e.x-w.x);
    float c = 0.25*(n.x-s.x-e.y+w.y);
    
    Q.z = m.z*.9999 - mix(d,c,length(U-0.5*R)/R.y);
    Q.w = d;
    
    if (FRAMECOUNT <= 1 || Reset != 0.) Q = vec4(sin(U.x)*cos(U.y));
	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    U += _uvc*Stretch*PI;
    U -= _uvc*PI*Zoom*(growthFactor);

    Q = A(U);
    vec4 
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0)),
        m = 0.25*(n+e+s+w);
    Q.xy = 0.25*vec2(e.z-w.z,n.z-s.z);
    Q.xy*_uvc;
    if (length(Q.xy)>0.) Q.xy = mix(Q.xy,normalize(Q.xy),.125);
    
	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = A(U);
    vec4 
        n = B(U+vec2(0,1+distortion)),
        e = B(U+vec2(1+distortion,0)),
        s = B(U-vec2(0,1+distortion)),
        w = B(U-vec2(1+distortion,0)),
        m = 0.25*(n+e+s+w);
    
    float d = 0.25*(n.y-s.y+e.x-w.x);
    float c = 0.25*(n.x-s.x-e.y+w.y);
    
    Q.z = m.z*.9999 - mix(d,c,.01);
    
    if (FRAMECOUNT <= 1 || Reset != 0.) Q = vec4(sin(U.x)*cos(U.y));
	return Q; 
 } 


			//******** BuffD Code Begins ********

vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    //U += ((moveXY.xy*(1.0+0.5*growthFactor)));

    vec2 c = (R);
    c.xy *=1.0+ ((moveXY.xy*(1.0+0.5*growthFactor)));
   // if (_mouse.z>0.) c = _mouse.xy;
    
    U -= c;
//    U += _uvc*PI;
    
//    U *= (.99325*(1.-growthFactor*Zoom*0.00675)-Zoom*0.01);
    U *= (.999);

    U += c;
    Q = C(U);

    vec4 
        n = C(U+vec2(0,1+Height)), 
        e = C(U+vec2(1+Height,0)),
        s = C(U-vec2(0,1+Height)),
        w = C(U-vec2(1+Height,0)),
        m = 0.25*(n+e+s+w);
    Q.xy = 0.25*vec2(e.z-w.z,n.z-s.z);
    
    if (length(Q.xy)>0.) Q.xy = mix(Q.xy,normalize(Q.xy),.1);
    
	return Q; 
 } 


mat2 r (float a) {
    vec2 e = vec2(cos(a),sin(a));
	return mat2(e.x,-e.y,e.y,e.x);
}
vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    vec3 
        p = vec3(0.5*R+_uvc,.62*R.y+_uvc.y*PI),
        d = normalize(vec3((U+_uvc-0.5*R)/R.y,-1)),
        o = vec3(0.5,.1,.5)*R.xyy;
    if (_mouse.z>0.) o.xy = _mouse.xy;
    mat2 m = r(.44);
   d.yz =_rotate(d.yz, Look);
//   d.xy =_rotate(d.xy, Rotate*PI);
   p.yz =_rotate(p.yz, Look);
   //p.xy =_rotate(p.xy, Rotate);

    p.y -= .19*R.y;
    d.yz *= m;
    p.yz *= m;

    for (int i = 0; i<35; i++){ 
        
    	p += (.2)*d*(p.z-(4.)*A(p.xy).z);
    }
    float highs = (0.35*pow(syn_HighLevel*0.35+syn_MidHighLevel*0.35+syn_Intensity*.2, 2));
    d = normalize(o-p);
    float z = A(p.xy).z;
    vec3 n = normalize(vec3(B(p.xy).xy,-.4));
    vec3 q = d;
    p += .1*d;
    for (int i = 0; i<35; i++){ 
    	p += .5*d*min(p.z-4.*A(p.xy).z,length(p-o)-1.);
    }
    //Q = (exp(-length(p-o)+1.)*(1.)*(cos(-.015*smoothTimeC*.2+(.15*(1.)) * z + .5*vec4(1,2,3,4)*-1.))*.5*(dot(reflect(n,d),q)-dot(n,d)))*(1.0+pow(syn_HighLevel, 2)*0.1);
    Q = (exp(-length(p-o)+1.)*(1.+highs)*(cos(-.065*smoothTimeB+(.15*(1.5)) * z + .5*vec4(1,2,3,4)*-1.))*.5*(dot(reflect(n,d),q)-dot(n,d)));
	Q*=Q;
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