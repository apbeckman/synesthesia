

			//******** Common Code Begins ********


// Dave Hoskins
// https://www.shadertoy.com/view/4djSRW
float hash13(vec3 p3)
{
	p3  = fract(p3 * .1031);
    p3 += dot(p3, p3.zyx + 31.32);
    return fract((p3.x + p3.y) * p3.z);
}


// Inigo Quilez
// https://iquilezles.org/articles/distfunctions
float sdBox( vec3 p, vec3 b ) {
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

// rotation matrix
mat2 rot(float a) { return mat2(cos(a),-sin(a),sin(a),cos(a)); }

#define repeat(p,r) (mod(p,r)-r/2.)

			//******** BuffA Code Begins ********

// FLUID EVOLUTION
#define R RENDERSIZE.xy
#define T(U) texture(BuffA,(U)/R)
#define D(U) texture(BuffD,(U)/R)
#define B(U) texture(BuffB,(U)/R)
// Velocity
vec2 v (vec4 b) {
	return vec2(b.x-b.y,b.z-b.w);
}
// Pressure
float p (vec4 b) {
	return 0.25*(b.x+b.y+b.z+b.w);
}
// TRANSLATE COORD BY Velocity THEN LOOKUP STATE
vec4 A(vec2 U) {
    U-=.5*v(T(U));
    U-=.5*v(T(U));
	return T(U);
}
vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    // THIS PIXEL
    Q = A(U);
    // NEIGHBORHOOD
    vec4 
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0));
    // GRADIENT of PRESSURE
    float px = 0.25*(p(e)-p(w));
    float py = 0.25*(p(n)-p(s)); 
    
    		// boundary Energy exchange in :   
    Q += 0.25*(n.w + e.y + s.z + w.x)
        	// boundary Energy exchange out :
        	-p(Q)
        	// dV/dt = dP/dx,  dEnergy In dTime = dEnergy in dSpace
        	-vec4(px,-px,py,-py);
    
    // get value from picture buffer
    float z = .8-length(B(U).xyz);
    // some kind of viscolsity thing 
    Q = mix(mix(Q,0.25*(n+e+s+w),.01),vec4(p(Q)),.01*(1.-z));
    // gravity polarizes energy! pretty cool imo
    Q.zw -= 0.001*z*vec2(1,-1);
    // Init with no velocity and some pressure
    if (FRAMECOUNT <= 1||(_mouse.z>0.&&length(U-_mouse.xy)<R.y/5.)) Q = vec4(.2);
    // At boundarys turn all kinetic energy into potential energy
    if(U.x<3.||R.x-U.x<3.||U.y<3.||R.y-U.y<3.)Q = vec4(p(Q));
	return Q; 
 } 


			//******** BuffB Code Begins ********

// LOOK UP PICTURE IN LOCATION FROM BUFFER D
//#define R RENDERSIZE.xy
//#define T(U) texture(iChannel0,(U)/R)
//#define D(U) texture(BuffD,(U)/R)
vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = texture(image7,D(U).xy/R);
	return Q; 
 } 


			//******** BuffC Code Begins ********

// Lighting on Buffer B

//#define T(U) texture(BuffB,(U)/R)

vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

   Q =  1.2-2.2*B(U);
    Q.xyz = Q.xyz+.5*normalize(Q.xyz);
   float
       n = length(B(U+vec2(0,1))),
       e = length(B(U+vec2(1,0))),
       s = length(B(U-vec2(0,1))),
       w = length(B(U-vec2(1,0)));
    vec3 no = normalize(vec3(e-w,n-s,1));
    float d = dot(reflect(no,vec3(0,0,1)),normalize(vec3(1)));
    Q *= 8.*exp(-3.*d*d);
	return Q; 
 } 


			//******** BuffD Code Begins ********

// TRANSLATE LOCATION FIELD WITH v(A(coord)), INIT WITH FragCoord
/*
#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
*/
#define d(U) texture(BuffD,(U)/R)
#define C(U) texture(BuffC,(U)/R)
/*
vec2 v (vec4 b) {
	return vec2(b.x-b.y,b.z-b.w);
}
*/
vec4 X(vec2 U) {
    U-=.5*v(A(U));
    U-=.5*v(A(U));
	return d(U);
}
vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = X(U);
    
    vec4 
        q = A(U),
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0)),
        N = X(U+vec2(0,1)),
        E = X(U+vec2(1,0)),
        S = X(U-vec2(0,1)),
        W = X(U-vec2(1,0));
    Q += 0.25*((n.w-q.z)*(N-Q) + (e.y-q.x)*(E-Q) + (s.z-q.w)*(S-Q) + (w.x-q.y)*(W-Q));
    
    if (FRAMECOUNT <= 1||(_mouse.z>0.&&length(U-_mouse.xy)<R.y/5.)) Q = vec4(U,0,0);
	return Q; 
 } 


// LENS FLAIR EFFECT
//#define R RENDERSIZE.xy
//#define T(U) texture(BuffC,(U)/R)
vec4 F (vec2 U,vec2 r) {
	vec4 t = C(U+r);
    return exp(-.01*dot(r,r))*(exp(2.*t)-1.);
}
vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

   
   Q = vec4(0);
    for (float i = 0.; i < 7.; i+=1.1) {
    	Q += F(U,+vec2(-i,i));
    	Q += F(U,+vec2(i,i));
    	Q += F(U,-vec2(-i,i));
    	Q += F(U,-vec2(i,i));
    }
    Q = C(U)*0.15+ 1e-5*Q;
    Q = atan(Q);
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