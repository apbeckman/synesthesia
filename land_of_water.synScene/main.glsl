

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA, (U)/R)
#define B(U) texture(BuffB, (U)/R)
#define C(U) texture(BuffC, (U)/R)
#define D(U) texture(BuffD, (U)/R)

#define N 3.

#define PRECIPITATION 1.*(1.0+syn_Intensity)//*(1.+basshits)
#define EVAPORATION .0001*(1.0+syn_BassLevel)

			//******** BuffA Code Begins ********

// Calculate forces and pressure
vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = B(U);
    vec4 
        n = B(U+vec2(0,1)),
        e = B(U+vec2(1,0)),
        s = B(U-vec2(0,1)),
        w = B(U-vec2(1,0)),
        c = D(U),
        nc = D(U+vec2(0,1)),
        ec = D(U+vec2(1,0)),
        sc = D(U-vec2(0,1)),
        wc = D(U-vec2(1,0));
    Q.xy -= 1./(1.+2.*sqrt(c.y))*(
        // slope force
        0.25*(vec2(ec.x-wc.x,nc.x-sc.x)*.5+vec2(ec.y-wc.y,nc.y-sc.y))+
        // pressure force
        (0.25+basshits)*vec2(e.z-w.z,n.z-s.z)+
        // magnus force
        0.25*Q.w*vec2(n.w-s.w,e.w-w.w));
    Q.xy *= min(1.,c.y);
    // divergence
    Q.z  = 0.25*(s.y-n.y+w.x-e.x+n.z+e.z+s.z+w.z)+basshits;
    // curl
    Q.w = 0.25*(s.x-n.x+w.y-e.y)+basshits;
    if (length(Q.xy) > .8) Q.xy = .8*normalize(Q.xy);
    
    //Boundary conditions
    if (FRAMECOUNT<=1) Q = vec4(0);
    
	return Q; 
 } 


			//******** BuffB Code Begins ********

// Advect along velocity and curl feild
vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    for (float i = 0.; i< N;i++) {
        Q = A(U);
        float co = cos(Q.w/N), si = sin(Q.w/N);
        U -= Q.xy*mat2(co,-si,si,co)/N;
    }
    Q = A(U);
	return Q; 
 } 


			//******** BuffC Code Begins ********

// Advect along velocity and curl feild
vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    for (float i = 0.; i< N;i++) {
        Q = A(U);
        float co = cos(Q.w/N), si = sin(Q.w/N);
        U -= Q.xy*mat2(co,-si,si,co)/N;
    }
    Q = D(U);
	return Q; 
 } 


			//******** BuffD Code Begins ********

vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = C(U);
    // neighborhood
    vec4 
        a = A(U),
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0)),
    	nc = C(U+vec2(0,1)),
        ec = C(U+vec2(1,0)),
        sc = C(U-vec2(0,1)),
        wc = C(U-vec2(1,0));
    Q = mix(Q,0.25*(nc+ec+sc+wc),vec4(0.,.1,0,0));
    // divergence 
    Q += 0.25*(s.y*sc-n.y*nc+w.x*wc-e.x*ec);
    
    // x : height y : water, z : sediment 
    float m = 0.25*(D(U+vec2(0,1)).x+D(U+vec2(1,0)).x+D(U-vec2(1,0)).x+D(U-vec2(0,1)).x);
    float me = D(U).x;
    float l = m-me;
   	Q.x = me + 0.03*l*l*l;
    if (_mouse.z>0.) Q.x += .25/(1.+.0005*dot(U-_mouse.xy,U-_mouse.xy));
    float x = .01*(Q.y*length(a.xy)*(1.-Q.z)+.01*Q.z);
    Q.z += x;
    Q.x -= x;
    Q.y = Q.y*(1.-EVAPORATION) + PRECIPITATION/R.x;
    Q = max(Q,0.);
    // boundary conditions
    if (FRAMECOUNT<5)Q = vec4(10.+1.-U.y/R.y+.5*B(U).x+10.*exp(-length(U-0.5*R)/R.y),.3,0,0);
    if (U.x<2.||U.y<2.||R.y-U.y<2.||R.x-U.x<2.) Q*=0.;

	return Q; 
 } 


// Controls in Common

vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    vec4 
        n = (C(U+vec2(0,1))),
        e = (C(U+vec2(1,0))),
        s = (C(U-vec2(0,1))),
        w = (C(U-vec2(1,0)));
    vec4 c = C(U);
    vec3 no = normalize(vec3(e.x-w.x+e.y-w.y,n.x-s.x+n.y-s.y,.5));
    Q = abs(sin(1.+3.*c.z+sqrt(c.y)*(1.+.5*vec4(1,2,3,4))));
    float a = .5;
    no.zy *= mat2(cos(a),-sin(a),sin(a),cos(a));
    no.zx *= mat2(cos(a),-sin(a),sin(a),cos(a));
    Q*=max(0.,.2+.8*no.z);
	
    
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