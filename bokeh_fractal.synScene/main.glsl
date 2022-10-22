

			//******** Common Code Begins ********

vec2 R; float T;
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
#define Main void mainImage(out vec4 Q, in vec2 U) {R=RENDERSIZE.xy,T=TIME;

#define D 9. 
#define N exp2(D)

#define ei(a) mat2(cos(a),-sin(a),sin(a),cos(a))

void shape (in float i, out vec3  u, out vec4 c) 
{
   	vec4 Q = vec4(0.0);
	vec2 U = _xy;


    u = vec3(0);
    c = vec4(0);
    float r = .1;
    for (float j = 0.; j < D; j++) {
         float a = mod(i,2.)*2.-1.;
         c += .2*sin(1.3*a+vec4(1,2,3,4));
         u.xy  *= ei(-a+.001*T);
         u.xz *= ei(.2+.1*sin(.31*T));
         u.y += r;
         u.x -= a*r;
         r *= .5;
         i = floor(i/2.);
    }
    u.xz  *= ei(.174*T);
    u *= 2.5;
}

			//******** BuffA Code Begins ********

vec4 renderPassA(){
	vec4 Q = vec4(0.0);
	vec2 U = _xy;


    vec3 u;
    float i = floor(U.x);
    shape(i,u,Q);
    
    Q.xyz = u;

	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB(){
	vec4 Q = vec4(0.0);
	vec2 U = _xy;


    vec3 u;
    float i = floor(U.x);
    shape(i,u,Q);
    
    //Q.xyz = u;

	return Q; 
 } 


vec4 renderMainImage(){
   	vec4 Q = vec4(0.0);
	vec2 U = _xy;


    //Q = vec4(0.);
    vec3 u;
    vec4 c;
    U = (U-.5*R)/R.y;
    for (float i = 0.; i < N; i++)
    {
        u = A(vec2(i+.5,.5)).xyz;
        c = B(vec2(i+.5,.5));
        float d = .2*abs(u.z);
        float l = (length(u.xy-U)-d);
        Q += c*(smoothstep(1e-2,0.,l)+.5*exp(-1e4*l*l))/(1.+400.*d);
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
		return renderMainImage();
	}
}