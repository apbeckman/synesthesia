

			//******** Common Code Begins ********

// Wavelength
#define H 8.
//Brightness
#define F 70.
// Iterations
#define I 20.


// Window of Samples
#define W .2*R.y
// Number of Samples
#define N 15.




vec2 R;
float T;
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
#define C(U) texture(BuffC,(U)/R)
//#define Main void mainImage(out vec4 Q, in vec2 U)
#define ei(a) mat2(cos(a),-sin(a),sin(a),cos(a))
#define pi 3.14159265359
#define _sin(a) sin(mod(a,2.*pi))
#define _cos(a) cos(mod(a,2.*pi))
vec4 map (vec2 u) {
    u = 1.1*(u-.5*R)/R.y;
    vec4 c = vec4(0);
    float j = 0.;
    for (float i = 1.; i < I; i++)
    {
        float d = length(u)-.2;
        c += .02/i/i*exp(-1e2*d*d)*(.5+.5*
            sin(2.+i+vec4(1,2,3,4))
        );
        u.x = abs(u.x);
        u -= vec2(.6,0);
        u = 2.4*vec2(u.x*u.x-u.y*u.y,2.*u.x*u.y);
        if (length(u)>4.)break;
        
    }
    return F*c;
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
    vec4 Q = vec4(0.);
    vec2 U = _xy;

    R = RENDERSIZE.xy;
    T = TIME;
    Q = A(U);
    vec3 eye = vec3(U,R.y);
    for (float i = 0.; i<N; i++)
    {
        vec2 v = U+H*(hash22(vec2(i,smoothTime))-.5);
        vec4 m = map(v);
        vec3 u = vec3(v,m.x);
        float l = W*(length(eye-u));
        Q.xy += vec2(_cos(l),_sin(l))/N;
    }
    if (FRAMECOUNT <= 1) Q = vec4(0);
	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() 
{
    vec4 Q = vec4(0.);
    vec2 U = _xy;

    R = RENDERSIZE.xy;
    T = TIME;
    Q = B(U);
    vec3 eye = vec3(U,R.y);
    for (float i = 0.; i<N; i++)
    {
        vec2 v = U+H*(hash22(vec2(i,smoothTime))-.5);
        vec4 m = map(v);
        vec3 u = vec3(v,m.y);
        float l = W*(length(eye-u));
        Q.xy += vec2(_cos(l),_sin(l))/N;
    }
    if (FRAMECOUNT <= 1) Q = vec4(0);
	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC()
{
    vec4 Q = vec4(0.);
    vec2 U = _xy;
    R = RENDERSIZE.xy;
    T = TIME;
    Q = C(U);
    vec3 eye = vec3(U,R.y);
    for (float i = 0.; i<N; i++)
    {
        vec2 v = U+H*(hash22(vec2(i,smoothTime))-.5);
        vec4 m = map(v);
        vec3 u = vec3(v,m.z);
        float l = W*(length(eye-u));
        Q.xy += vec2(_cos(l),_sin(l))/N;
    }
    if ((FRAMECOUNT) <= 1) Q = vec4(0);
	return Q; 
 } 


// Fork of "Microscopy 101" by wyatt. https://shadertoy.com/view/NlKSzG
// 2022-01-04 06:32:20


vec4 renderMainImage() 
{
    vec4 Q = vec4(0.);
    vec2 U = _xy;

    R = RENDERSIZE.xy;T = TIME;
    Q.x = length(A(U).xy)+sin(TIME);
    Q.y = length(B(U).xy);
    Q.z = length(C(U).xy);
    Q /= float(FRAMECOUNT);
    //=Q = map(U);
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