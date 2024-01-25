

			//******** BuffA Code Begins ********
bool fc = (FRAMECOUNT <= 1);
#define R RENDERSIZE.xy
vec4 A (vec2 U) {return texture(BuffA,U/R);}
vec4 B (vec2 U) {return texture(BuffB,U/R);}
vec4 C (vec2 U) {return texture(BuffC,U/R);}
vec4 D (vec2 U) {return texture(BuffD,U/R);}

float ln (vec3 p, vec3 q, vec3 b) {return length(p-q-(b-q)*dot(p-q,b-q)/dot(b-q,b-q));}
float a = 1.;

float o = a*.7071;

vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    Q = D(U);
    float g = TIME+3.*length(U), c = cos(g), s = sin(g);
    mat2 m = mat2(c,-s,s,c);
    vec4 n = D(U+vec2(0,a)*m)+D(U+vec2(a,0)*m)+D(U+vec2(0,-a)*m)+D(U+vec2(-a,0)*m)+D(U+vec2(-o,o)*m)+D(U+vec2(o,-o)*m)+D(U+vec2(-o,-o)*m)+D(U+vec2(o,o)*m);
    n *= .125;
    vec4 dx = n-Q;
    Q += dx*vec4(1.,.3,1.,1);
    
    float x = .3*Q.x*Q.y*(1.-Q.y);
    
    Q.y = Q.y+x-0.025-.1*Q.z;
    Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;
    Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;
    
    Q = clamp(Q,0.,1.);
    
    if (_mouse.z> 0. && length(U-_mouse.xy) < 40.) Q.y=1.;
    
if (fc) {
        
        Q = vec4(1,0,0,0);
        if (length(U-0.5*R)<10.) Q.y += 1.;
    }
	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    Q = A(U);
    float a = TIME+3.*length(U);
    float c = cos(a);
    float s = sin(a);
    mat2 m = mat2(c,-s,s,c);
    vec4 n = A(U+vec2(0,a)*m)+A(U+vec2(a,0)*m)+A(U+vec2(0,-a)*m)+A(U+vec2(-a,0)*m)+A(U+vec2(-o,o)*m)+A(U+vec2(o,-o)*m)+A(U+vec2(-o,-o)*m)+A(U+vec2(o,o)*m);
    n *= .125;
    vec4 dx = n-Q;
    Q += dx*vec4(1.,.3,1.,1);
    
    float x = .3*Q.x*Q.y*(1.-Q.y);
    
    Q.y = Q.y+x-0.025-.1*Q.z;
    Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;
    Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;
    
    Q = clamp(Q,0.,1.);
    
    if (_mouse.z> 0. && length(U-_mouse.xy) < 40.) Q.y=1.;
    
    if (fc) {
        
        Q = vec4(1,0,0,0);
        if (length(U-0.5*R)<10.) Q.y += 1.;
    }
	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    a += 1.0;
    o = a*.7071;

    Q = B(U);
    float g = TIME+3.*length(U), c = cos(g), s = sin(g);
    mat2 m = mat2(c,-s,s,c);
    vec4 n = B(U+vec2(0,a)*m)+B(U+vec2(a,0)*m)+B(U+vec2(0,-a)*m)+B(U+vec2(-a,0)*m)+B(U+vec2(-o,o)*m)+B(U+vec2(o,-o)*m)+B(U+vec2(-o,-o)*m)+B(U+vec2(o,o)*m);
    n *= .125;
    vec4 dx = n-Q;
    Q += dx*vec4(1.,.3,1.,1);
    
    float x = .3*Q.x*Q.y*(1.-Q.y);
    
    Q.y = Q.y+x-0.025-.1*Q.z;
    Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;
    Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;
    
    Q = clamp(Q,0.,1.);
    
    if (_mouse.z> 0. && length(U-_mouse.xy) < 40.) Q.y=1.;
    
    if (fc) {
        
        Q = vec4(1,0,0,0);
        if (length(U-0.5*R)<10.) Q.y += 1.;
    }
	return Q; 
 } 


			//******** BuffD Code Begins ********
vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = C(U);
    float g = TIME+3.*length(U), c = cos(g), s = sin(g);
    mat2 m = mat2(c,-s,s,c);
    a += 1.0;

    o = a*.7071;

    vec4 n = C(U+vec2(0,a)*m)+C(U+vec2(a,0)*m)+C(U+vec2(0,-a)*m)+C(U+vec2(-a,0)*m)+C(U+vec2(-o,o)*m)+C(U+vec2(o,-o)*m)+C(U+vec2(-o,-o)*m)+C(U+vec2(o,o)*m);
    n *= .125;
    vec4 dx = n-Q;
    Q += dx*vec4(1.,.3,1.,1);
    
    float x = .3*Q.x*Q.y*(1.-Q.y);
    
    Q.y = Q.y+x-0.025-.1*Q.z;
    Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;
    Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;
    
    Q = clamp(Q,0.,1.);
    
    if (_mouse.z> 0. && length(U-_mouse.xy) < 40.) Q.y=1.;
    
    if (fc) {
        
        Q = vec4(1,0,0,0);
        if (length(U-0.55*R)<10.) Q.y += 1.;
        if (length(U-0.45*R)<10.) Q.y += 1.;

    }
	return Q; 
 } 


vec4 renderMainImage() {

	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    vec4
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0));
    vec2 h = vec2(e.w-w.w,n.w-s.w);
    vec2 g = vec2(e.y-w.y,n.y-s.y);
    vec4 b = A(U);
   	vec4 dx = 0.25*(n+e+s+w) - b;
	Q = (b.w+b.y)*normalize(abs(sin(4.+.5*(b.z-b.x)*vec4(3,2,1,4))));
    
    vec3 no = normalize(vec3(g-h,-1));
    vec3 r = reflect(vec3(0,0,1),no);
    vec3 l = vec3(2.*R.x,2.*R.y,10.*R.y);
    vec3 u = vec3(U,0);
    vec3 lu = l-u;
    float o = length(r-lu*dot(r,lu)/dot(lu,lu));
    Q = (Q)*(exp(-o)+2.*exp(-7.5*o)+1e2*exp(-30.*o));
   // Q *= 1.+texture(iChannel1,r);
     Q *= 1.+texture(prism,r.xz);

    no = normalize(vec3(g-h,1));
    r = refract(vec3(0,0,1),no,.2);
    l = vec3(3.*R.x,3.*R.y,10.*R.y);
    u = vec3(U,0);
    lu = l-u;
    o = length(r-lu*dot(r,lu)/dot(lu,lu));
    Q = 2.*dx+(Q*0.9+0.1)*(0.4+exp(-o)+10.*exp(-9.5*o));
   // Q *= 1.+0.2*sin(Q+3.*texture(iChannel1,r).x*vec4(1,2,3,4));
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