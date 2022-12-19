

			//******** Common Code Begins ********

// Wavelengths of Each Color
#define _A (10000./600.)
#define _B (10000./500.)
#define _C (10000./400.)
// Depth of Bubbles
#define D 1.
// Iterations
#define I 8.


// Window of Samples
#define W 90.
// Number of Samples
#define N 80.

// Brightness Decay
#define Z mix(0.,.999,step(.01,mod(T/.2,1.)+0.05*syn_Hits))
// Adjust for Brightness
#define F .06



vec2 R;
float T;
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
#define C(U) texture(BuffC,(U)/R)
#define Main void mainImage(out vec4 Q, in vec2 U)
#define ei(a) mat2(cos(a),-sin(a),sin(a),cos(a))
#define pi 3.14159265359
#define _sin(a) sin(mod(a,2.*pi))
#define _cos(a) cos(mod(a,2.*pi))
float map (vec2 u) {
    u = 4.*(u-.5*R)/R.y;
    float d = 0.;
    float s = 1.;
    for (float i = 1.; i < I; i++)
    {
        
        u *= (1.7+.5*sin(1.74*floor(.2*T)))*ei(1.-floor(.2*T));
        u = abs(u)-2.;
        
        float l = 1.-length(u);
        if (l>0.)
            d += exp2(-i)*sqrt(l);
        
        
    }
    
        
    return D*d;
}
// Dave H
vec2 hash22(vec2 p)
{
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+33.33);
    return fract((p3.xx+p3.yz)*p3.zy);

}


			//******** BuffA Code Begins ********

vec4 renderPassA() 
 
{
    vec4 Q = vec4(0.0);
	vec2 U = _xy;

    R = RENDERSIZE.xy;
    T = smoothTimeC;
    Q = A(U);
    vec3 eye = vec3(U,R.y);
    for (float i = 0.; i<N; i++)
    {
        vec2 v = U+W*(hash22(vec2(i,FRAMECOUNT))-.5);
        vec3 u = vec3(v,map(v));
        float l = _A*(length(eye-u));
        Q.xy += F*vec2(_cos(l),_sin(l))/N;
    }
    Q *= Z;
    if (FRAMECOUNT < 1) Q = vec4(0);
	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() 
{
  	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    R = RENDERSIZE.xy;
    T = smoothTimeC;
    Q = B(U);
    vec3 eye = vec3(U,R.y);
    for (float i = 0.; i<N; i++)
    {
        vec2 v = U+W*(hash22(vec2(i,FRAMECOUNT))-.5);
        vec3 u = vec3(v,map(v));
        float l = _B*length(eye-u);
        Q.xy += F*vec2(_cos(l),_sin(l))/N;
    }
    Q *= Z;
    if (FRAMECOUNT < 1) Q = vec4(0);
	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() 
{
  	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    R = RENDERSIZE.xy;
    T = smoothTimeC;
    Q = C(U);
    vec3 eye = vec3(U,R.y);
    for (float i = 0.; i<N; i++)
    {
        vec2 v = U+W*(hash22(vec2(i,FRAMECOUNT))-.5);
        vec3 u = vec3(v,map(v));
        float l = _C*length(eye-u);
        Q.xy += F*vec2(_cos(l),_sin(l))/N;
    }
    Q *= Z;
    if (FRAMECOUNT < 1) Q = vec4(0);
	return Q; 
 } 



vec4 renderMainImage() 
{
  	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    R = RENDERSIZE.xy;
    Q.x = length(A(U).xy);
    Q.y = length(B(U).xy);
    Q.z = length(C(U).xy);
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
		return renderMainImage();
	}
}