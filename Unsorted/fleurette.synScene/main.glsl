

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
//#define C(U) texture(iChannel2,(U)/R)
//#define D(U) texture(iChannel3,(U)/R)
// /#define Main void mainImage(out vec4 Q, in vec2 U)


vec2 hash23(vec3 p3) {  // Dave H
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+33.33);
    return fract((p3.xx+p3.yz)*p3.zy);
}

			//******** BuffA Code Begins ********

vec4 renderPassA() {
  	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = vec4(0);
    U = 4.*(U-.5*R)/R.y;
    for (float i = 0.; i < 9.; i++) {
        U = vec2(U.x*U.x-U.y*U.y,2.*U.x*U.y)-vec2(-.6-(cos(smoothTime*0.085)*0.25-0.5),sin(.085*smoothTime));
        U /= .5+0.3*dot(U,U);
        Q += .05*length(U)*vec4(.5+sin(2.*U.x+vec3(1,2,3)),1);
    }

	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
   	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = .99*B(U);
    for (int i = 0; i < 30; i++) {
        vec2 h = 20.*(hash23(vec3(U+R,i+30*FRAMECOUNT))*2.-1.);
        vec4 c = A(U+h),
             n = A(U+h+vec2(0,1)),
             e = A(U+h+vec2(1,0)),
             s = A(U+h-vec2(0,1)),
             w = A(U+h-vec2(1,0));

        vec2 g = 2.*R.x*vec2(e.w-w.w,n.w-s.w);

        Q += exp(-length(h-g))*c;
   }
	return Q; 
 } 


vec4 renderMainImage() {
   	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = .05*B(U);

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