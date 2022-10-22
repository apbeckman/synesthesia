

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
#define Main void mainImage(out vec4 Q, in vec2 U)
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
#define ei(a) mat2(cos(a),sin(a),-sin(a),cos(a))
float segment (vec2 p, vec2 a, vec2 b) {
    return length(p-a-(b-a)*clamp(dot(p-a,b-a)/dot(b-a,b-a),0.,1.));
}

			//******** BuffA Code Begins ********

vec4 renderPassA() {
   	vec4 Q = vec4(0.0);
    vec2 U = _xy;

    Q = A(U);
    if (FRAMECOUNT <= 1) Q = vec4(length(U-.5*R));
    
    if (_mouse.z>0.) {
        vec4 b = B(U);
        Q = vec4(min(Q.x,segment(U,b.xy,b.zw)));
    }

	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    Q = B(U);
    if (_mouse.z>0.) {
        if (Q.x<R.x) {
            Q = vec4(Q.zw,_mouse.xy);
        }
        else Q = _mouse.xyxy;
    }
    else Q = vec4(1e9);
    if (Q.x<.5*R.x) Q.x = R.x-Q.x;
    if (Q.z<.5*R.x) Q.z = R.x-Q.z;
	return Q; 
 } 


float map (vec3 p) {
    vec2 u = vec2(length(p.xz),.5*p.y);
    u = u*R.y*0.5+0.5*R;
    return A(u).x/R.x-.001;
    
   
}

vec3 norm (vec3 p) {
    vec3 e = vec3(1e-3,0,0);
    return normalize(vec3(
        map(p+e.xyz)-map(p-e.xyz),
        map(p+e.zxy)-map(p-e.zxy),
        map(p+e.yzx)-map(p-e.yzx)
    ));
}

vec4 renderMainImage() {
 	vec4 Q = vec4(0.0);
    vec2 U = _xy;

    vec3 d = normalize(vec3(2.*(U-.5*R)/R.y,1));
    vec3 p = vec3(0,0,-3);
    p.yz *= ei(.75*sin(.5*TIME));
    d.yz *= ei(.75*sin(.5*TIME));
    Q = vec4(0);
    for (float x = 0.; x < 150.; x++) {
        float l = map(p);
        p += d*l;
        if (l < .001) Q.xyz = .5-.5*norm(p);
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