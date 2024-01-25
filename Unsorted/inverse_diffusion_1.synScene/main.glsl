

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
vec2 mirror (vec2 u) {
    if (u.x>1.) u.x = 1.-fract(u.x);
    if (u.x<0.) u.x = fract(u.x);
    if (u.y>1.) u.x = 1.-fract(u.y);
    if (u.y<0.) u.x = fract(u.y);
	return u;
}
#define A(U) texture(BuffA,mirror((U)/R))
#define B(U) texture(BuffB,mirror((U)/R))
#define C(U) texture(BuffC,mirror((U)/R))
#define D(U) texture(BuffD,mirror((U)/R))

			//******** BuffA Code Begins ********
vec2 Mirror = vec2(x_flip, y_flip);

vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = C(U);
    vec4 
        n = D(U+vec2(0,1)),
        e = D(U+vec2(1,0)),
        s = D(U-vec2(0,1)),
        w = D(U-vec2(1,0)),
        m = 0.25*(n+e+s+w);
    
    float d = 0.25*(n.y-s.y+e.x-w.x);
    float c = 0.25*(n.x-s.x-e.y+w.y);
    
    Q.z = m.z*.999 - mix(d,c,length(U-0.5*R)/R.y);
    Q.w = d;
    if (int(FRAMECOUNT) <= 1 || Reset != 0.) Q = vec4(sin(U.x)*cos(U.y));
	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    U -= _uvc*PI*Zoom*(1.0+low);
    U += vec2(dot(cos(_uvc.x),sin(_uvc.y)))*Pull*PI;
    U += _uvc*Stretch*PI*(1.0+low);
    U += moveXY.xy*(1.0+0.5*low);

    Q = A(U);
    vec4 
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0)),
        m = 0.25*(n+e+s+w);
    Q.xy = 0.25*vec2(e.z-w.z,n.z-s.z);
    vec2 uvc_noise = vec2(_noise(_uvc.x), _noise(_uvc.y));

    if (length(Q.xy)>0.) Q.xy = mix(Q.xy+uvc_noise,normalize(Q.xy),.2-syn_BassLevel*0.05);
    
	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = A(U);
    vec4 
        n = B(U+vec2(0,1)),
        e = B(U+vec2(1,0)),
        s = B(U-vec2(0,1)),
        w = B(U-vec2(1,0)),
        m = 0.25*(n+e+s+w);
    
    float d = 0.25*(n.y-s.y+e.x-w.x);
    float c = 0.25*(n.x-s.x-e.y+w.y);
    
    Q.z = m.z*.999 - mix(d,c,.2);
    
    if (int(FRAMECOUNT) <= 1) Q = vec4(sin(U.x)*cos(U.y));
	return Q; 
 } 


			//******** BuffD Code Begins ********

vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    vec2 c = R;
    
    // if (_mouse.z>0.) c = _mouse.xy;
    c = c+((moveXY.xy*(1.0+0.5*low)));
    c *= 1.0-vec2(_noise(TIME*0.7), _noise(TIME*0.7));
   // c += _noise(c);

    U -= c;
    U -= _uvc*PI*Zoom*(1.0+low);

//    U *= .99;
    U += c;
    Q = C(U);
    vec2 a = vec2(0, Height);
    vec2 b = vec2(Height, 0);
    vec4 
        n = C(U+vec2(0,1)+a),
        e = C(U+vec2(1,0)),
        s = C(U-vec2(0,1)+a),
        w = C(U-vec2(1,0)),
        m = 0.25*(n+e+s+w);
    Q.xy = 0.25*vec2(e.z-w.z,n.z-s.z);
    if (length(Q.xy)>0.) Q.xy = mix(Q.xy,normalize(Q.xy),.2);
    
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
        p = vec3(0.5*R,.62*R.y),
        d = normalize(vec3((U-0.5*R)/R.y,-1)),
        o = vec3(0.5,.1,.5)*R.xyy;
    if (_mouse.z>0.) o.xy = _mouse.xy;
    mat2 m = r(.44);
    d.xy = mix(abs(d.xy), d.xy, 1.0 - Mirror);
    p.y -= .19*R.y;
    
    d.yz *= m;
    p.yz *= m;
    for (int i = 0; i<35; i++){ 
    	p += .2*d*(p.z-4.*A(p.xy).z);
    }
    d = normalize(o-p);
    float z = A(p.xy).z;
    vec3 n = normalize(vec3(B(p.xy).xy,-.25));
    vec3 q = d;
    p += .1*d;
    for (int i = 0; i<35; i++){ 
    	p += .5*d*min(p.z-4.*A(p.xy).z,length(p-o)-1.);
    }
    // Q = (exp(-length(p-o)+1.))*(cos(-.1*TIME+.1*z+.5*vec4(1,2,3,4)))*.5*(dot(reflect(n,d),q)-dot(n,d));
    Q = (exp(-length(p-o)+1.))*(cos(-.075*smoothTimeC+.1*z+.5*vec4(1,2,3,4)))*.5*(dot(reflect(n,d),q)-dot(n,d));
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