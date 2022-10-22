

			//******** Common Code Begins ********

#define N 7.
vec2 R;
#define R3D vec3(R/N,N*N)
vec2 d2 (vec3 U) {
    U = mod(U,R3D);
    U.z = floor(U.z);
	return U.xy+vec2(mod(U.z,N),floor(U.z/N))*R/N;
}
vec3 d3 (vec2 u) {
    vec2 o = floor(u/R*N);
	return vec3(mod(u,R/N),o.x+o.y*N);
}
mat2 e (float a) {
	float c = cos(a), s = sin(a);
    return mat2(c,-s,s,c);
}
#define A(U) texelFetch(iChannel0,ivec2(d2(U)),0)
#define B(U) texelFetch(iChannel1,ivec2(d2(U)),0)
#define Sampler vec4 T(vec3 U) {return mix(texture(iChannel0,d2(vec3(U.xy,floor(U.z)))/R),texture(iChannel0,d2(vec3(U.xy, ceil(U.z)))/R),fract(U.z));}
#define Main void mainImage( out vec4 Q, in vec2 u )
#define _3D  R = RENDERSIZE.xy; vec3 U = d3(u);
#define Me Q = A(U);
#define Them vec4 M = 1./6.*(A(U+vec3(1,0,0))+A(U+vec3(0,1,0))+A(U+vec3(0,0,1))+A(U-vec3(1,0,0))+A(U-vec3(0,1,0))+A(U-vec3(0,0,1)));
#define Init  if (FRAMECOUNT < 1) Q = vec4(exp(-.3*length(U-0.5*R3D)),0,1,0);

#define Diffuse Q += vec4(.2,0,1,0)*(M-Q);
#define Grow float g = .13*Q.z*Q.x*(1.-Q.x);
#define React Q.x += g-.006; Q.z += .006*(1.-Q.z)-g;


			//******** BuffA Code Begins ********

renderPassA
{
    _3D
    Me
    Them
    
    Diffuse
    
    Grow
    
    React
    
    Q = max(Q,0.);
        
    Init  
   

	return ; 
 } 


			//******** BuffB Code Begins ********

renderPassB
{
    _3D
    Me
    Them
    
    Diffuse
    
    Grow
    
    React
    
    Q = max(Q,0.);
        
    Init  
   

	return ; 
 } 


			//******** BuffC Code Begins ********

renderPassC
{
    _3D
    Me
    Them
    
    Diffuse
    
    Grow
    
    React
    
    Q = max(Q,0.);
        
    Init  
   

	return ; 
 } 


			//******** BuffD Code Begins ********

renderPassD
{
    _3D
    Me
    Them
    
    Diffuse
    
    Grow
    
    React
    
    Q = max(Q,0.);
        
    Init  
   

	return ; 
 } 


Sampler
Main
{
    R = RENDERSIZE.xy;
    vec3 mi = 0.5*vec3(R/N,N*N);
    vec3 p = vec3(0,0,-R.x/N);
    vec3 d = normalize(vec3((u-0.5*R)/R.y,1));
    if (_mouse.z>0.) {
 		p.zx *= e(6.2*_mouse.x/R.x);
		d.zx *= e(6.2*_mouse.x/R.x);
        p.yz *= e(6.2*_mouse.y/R.y);
		d.yz *= e(6.2*_mouse.y/R.y);
    } else {
		p.yz *= e(.2*TIME);
		d.yz *= e(.2*TIME);
	}
    Q = vec4(0);
    for (int i = 0; i < 140; i++) {
        vec3 o = abs(p)-mi;
        float m = length(max(o,0.));
        if (m<.01)
        { 	
            vec4 a = T(p+mi);
            float aa = a.x;
            Q += .02*(1.-exp(-aa))*abs(.75+0.25*sin(.1*TIME+6.*a.x+vec4(1,2,3,4)));
            p += d*exp(-10.*aa);
            p = mod(p+mi,R3D)-mi;
        } else p += d*m;
        
 	}
	return ; 
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