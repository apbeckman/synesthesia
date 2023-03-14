

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
#define C(U) texture(BuffC,(U)/R)
#define D(U) texture(BuffD,(U)/R)
#define Main void mainImage(out vec4 Q, in vec2 U)

            #define r 1.2
#define N 15.
#define S vec4(4,7,1,1)
#define Gaussian(i) 0.3989422804/S*exp(-.5*(i)*(i)/S/S)


			//******** BuffA Code Begins ********

vec4 renderPassA() {
  	vec4 Q = vec4(0.0);
	vec2 U = _xy;


    if (int(FRAMECOUNT)%2<=1) {
        Q = vec4(0);
        for (int x = -1; x <= 1; x++)
        for (int y = -1; y <= 1; y++)
        {
            vec2 u = vec2(x,y);
            vec4 a = A(U+u);
            vec2 w1 = clamp(U+u+a.xy-0.5*r,U - 0.5,U + 0.5),
                 w2 = clamp(U+u+a.xy+0.5*r,U - 0.5,U + 0.5);
            float m = (w2.x-w1.x)*(w2.y-w1.y)/(r*r);
            Q.xyz += m*a.w*a.xyz;
            Q.w += m*a.w;
        }
        if (Q.w>0.)
            Q.xyz/=Q.w;
        if ((FRAMECOUNT) <= 1) 
        {
            Q = vec4(0,0,.1,0);
            if (length(U-vec2(0.5)*R)<.3*R.y)Q.w = .8;
        }
        if (_mouse.z>0.&&length(U-_mouse.xy)<20.) Q.xw = vec2(.25,.3);
        if (U.x<1.||U.y<1.||R.x-U.x<1.||R.y-U.y<1.) Q.xy *= 0.;
    } else {
    	Q = A(U);vec4 q = Q;
    for (int x = -1; x<=1; x++)
	for (int y = -1; y<=1; y++)
    if (x!=0||y!=0)
    {
        vec2 u = vec2(x,y);
        vec4 a = A(U+u), b = B(U+u), d = D(U+u);
        u = (u)/dot(u,u);
        Q.xy -= q.w*0.125*(d.y+d.x+a.w*(a.w-.8))*u;
    }
    if (Q.w < 1e-3) Q.z *= 0.;
    Q.y -= 1e-4*step(1e-3,Q.w);
    Q.xy *= 0.999;
    }
	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
  	vec4 Q = vec4(0.0);
	vec2 U = _xy;


    if (int(FRAMECOUNT)%2<1) {
        Q = vec4(0);
        for (int x = -1; x <= 1; x++)
        for (int y = -1; y <= 1; y++)
        {
            vec2 u = vec2(x,y);
            vec4 a = B(U+u);
            vec2 w1 = clamp(U+u+a.xy-0.5*r,U - 0.5,U + 0.5),
                 w2 = clamp(U+u+a.xy+0.5*r,U - 0.5,U + 0.5);
            float m = (w2.x-w1.x)*(w2.y-w1.y)/(r*r);
            Q.xyz += m*a.w*a.xyz;
            Q.w += m*a.w;
        }
        if (Q.w>0.)
            Q.xyz/=Q.w;
        if ((FRAMECOUNT) <= 1) 
        {
            Q = vec4(0,0,.1,0);
            if (length(U-vec2(0.5)*R)<.3*R.y)Q.w = .8;
        }
        if (_mouse.z>0.&&length(U-_mouse.xy)<20.) Q.xw = vec2(.25,.3);
        if (U.x<1.||U.y<1.||R.x-U.x<1.||R.y-U.y<1.) Q.xy *= 0.;
    } else {
    	Q = B(U);vec4 q = Q;
    for (int x = -1; x<=1; x++)
	for (int y = -1; y<=1; y++)
    if (x!=0||y!=0)
    {
        vec2 u = vec2(x,y);
        vec4 a = B(U+u), b = A(U+u), d = D(U+u);
        u = (u)/dot(u,u);
        Q.xy -= q.w*0.125*(d.x-.01*b.w+d.y+a.w*(a.w-1.))*u;
        Q.xy -= q.w*0.125*(abs(q.z)*a.z)*vec2(-u.y,u.x);
        Q.z  -= q.w*0.125*a.w*(dot(vec2(-u.y,u.x),a.xy-q.xy));
    }
    if (Q.w < 1e-3) Q.z *= 0.;
    Q.y -= 1e-3*step(1e-3,Q.w);
    Q.xy *= 0.999;
    }
	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

	Q = vec4(0);
    for (float i = -N; i <= N; i++) {
    	vec4 a = A(U+vec2(i,0)),
             b = B(U+vec2(i,0));
        Q += Gaussian(i)*vec4(a.w,b.w,0,0);
    }
	return Q; 
 } 


			//******** BuffD Code Begins ********

vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;


	Q = vec4(0);
    for (float i = -N; i <= N; i++) {
        Q += 2.*Gaussian(i)*C(U+vec2(0,i));
    }
	return Q; 
 } 


// Fork of "4-Substance" by wyatt. https://shadertoy.com/view/3lffzM
// 2020-08-03 02:14:45

// Fork of "Multi-Substance" by wyatt. https://shadertoy.com/view/WtffRM
// 2020-08-01 02:57:11

vec4 renderMainImage() {
   	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    vec4
        n = A(U+vec2(0,1))+B(U+vec2(0,1))+C(U+vec2(0,1))+D(U+vec2(0,1)),
        e = A(U+vec2(1,0))+B(U+vec2(1,0))+C(U+vec2(1,0))+D(U+vec2(1,0)),
        s = A(U-vec2(0,1))+B(U-vec2(0,1))+C(U-vec2(0,1))+D(U-vec2(0,1)),
        w = A(U-vec2(1,0))+B(U-vec2(1,0))+C(U-vec2(1,0))+D(U-vec2(1,0));
    vec3 norm = 
        normalize(vec3(e.z*e.w-w.z*w.w,n.z*n.w-s.z*s.w,10)),
        ref = reflect(vec3(0,0,-1),norm);
   
	vec4 a = A(U), b = B(U), c = C(U), d = D(U);
    Q = vec4(a.w,b.w,c.w,1)+d.w;
    Q.xyz *= 1.8*mat3(.3,.5,.5,.3,.7,.2,.5,.6,.3);
	Q *= 1.+0.5*norm.x;
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