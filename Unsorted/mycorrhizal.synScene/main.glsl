

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
#define o vec3(1,0,-1)
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
//#define C(U) texture(iChannel2,(U)/R)
#define D(U) texture(BuffD,(U)/R)
#define Main void mainImage(out vec4 Q, vec2 U)
float ln (vec2 p, vec2 a, vec2 b) {
	return length(p-a-(b-a)*clamp(dot(p-a,b-a)/dot(b-a,b-a),0.,.9));
}
#define norm(u) ((u)/(1e-9+length(u)))
float growthFactor = pow((syn_BassLevel*0.35)+(syn_MidLevel*0.35)+(syn_Level*0.2), 2.0)*syn_Intensity;
			//******** BuffA Code Begins ********

void X (inout vec4 Q, vec2 U, vec2 r) {
    vec4 n = A(U+r);
	if (ln(U,n.xy,n.zw)<ln(U,Q.xy,Q.zw)) Q = n;
}
vec4 renderPassA() {
   	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = A(U);
    for (int x = -2;x <=2; x++)
    for (int y = -2;y <=2; y++)
    X(Q,U,vec2(x,y));
    Q.xy = mix(Q.xy,A(Q.xy).xy,1.);
    Q.zw = mix(Q.zw,A(Q.zw).zw,.01);
    Q.xy += D(Q.xy).xy;
    Q.zw += D(Q.zw).xy;
    
    if (length(Q.xy-Q.zw) > 2.5) {
        vec2 m = 0.5*(Q.xy+Q.zw);
        if (length(U-Q.xy) > length(U-Q.zw)) 
        	Q.xy = m;
        else Q.zw = m;
    }
    if (FRAMECOUNT<=1) {
        Q = vec4(0.51*R,0.49*R);
        vec4 a =vec4(vec2(0.49,.51)*R,vec2(.51,.49)*R);
        if (ln(U,a.xy,a.zw)<ln(U,Q.xy,Q.zw))
            Q = a;
    }

	return Q; 
 } 


			//******** BuffB Code Begins ********

//Mouse
vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

	vec4 C = vec4(0.0);

    vec4 p = texture(BuffB,U/RENDERSIZE.xy);
   	if (_mouse.z>0.) {
      if (p.z>0.) C =  vec4(_mouse.xy,p.xy);
    else C =  vec4(_mouse.xy,_mouse.xy);
   }
    else C = vec4(-RENDERSIZE.xy,-RENDERSIZE.xy);
	return C; 
 } 


			//******** BuffD Code Begins ********

vec4 T(vec2 U) {
	vec4 Q = vec4(0.0);

	U -= .5*D(U).xy;
	U -= .5*D(U).xy;
    return D(U);
}
vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = T(U);
    vec4 
        n = T(U+o.yx),
        e = T(U+o.xy),
        s = T(U+o.yz),
        w = T(U+o.zy),
        m = 0.25*(n+e+s+w);
    Q.xy = m.xy-0.25*vec2(e.z-w.z,n.z-s.z);
	Q.z = Q.z-0.25*(n.y+e.x-s.y-w.x);
    vec4 a = A(U);
    float l = ln(U,a.xy,a.zw);
    float v = exp(-l);
    Q.z -= .001*v;
    Q.xy = mix(Q.xy,-norm(a.xy-a.zw),.1*v*exp(-Q.z*Q.z));
    Q.xy *= .5*(1.0+growthFactor);
    if (_mouse.z>0.) Q.xy += exp(-.1*length(U-_mouse.xy))*(U-_mouse.xy);
    if (U.x<1.||R.x-U.x<1.||U.y<1.||R.y-U.y<1.) Q.z *= 0.;
	if (FRAMECOUNT <= 1) Q = vec4(0,0,0,0);
	return Q; 
 } 


// Fork of "Line Tracking Fluid" by wyatt. https://shadertoy.com/view/tsKXzd
// 2020-12-10 19:40:02

vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    vec4 a = A(U);
    vec4 d = D(U);
    float l = ln(U,a.xy,a.zw);
    Q = max(cos(d.z*vec4(1,2,3,4)),.1)*(exp(-l));
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
		return renderPassD();
	}
	if(PASSINDEX == 3){
		return renderMainImage();
	}
}