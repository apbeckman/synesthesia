

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
#define C(U) texture(BuffC,(U)/R)
#define D(U) texture(image47,(U)/R)


// oneshade:
//https://www.shadertoy.com/view/7sKSRh
float erf(in float x) {
    x *= .8;
    //return sign(x) * sqrt(1.0 - exp(-1.239192 * x * x));
    return sign(x) * sqrt(1.0 - exp2(-1.787776 * x * x)); // likely faster version by @spalmer
}
float erfstep (float a, float b, float x) {
    return .5*(erf(b-x)-erf(a-x));
}
// Dave_Hoskins https://www.shadertoy.com/view/4djSRW
vec3 hash31(float p)
{
   vec3 p3 = fract(vec3(p) * vec3(.1031, .1030, .0973));
   p3 += dot(p3, p3.yzx+33.33);
   return fract((p3.xxy+p3.yzz)*p3.zyx); 
}


			//******** BuffA Code Begins ********

vec4 renderPassA()
{
    vec4 Q = vec4(0.0);
    vec2 U = _xy;
    U -= _uvc*Zoom;
    vec4 dQ = Q = vec4(0);
    for (float x = -4.; x<=4.;x++)
    for (float y = -4.; y<=4.;y++)
    {
        vec2 u = vec2(x,y);
        vec4 a = B(U+u);
        vec2 v = u+a.xy;
        float w = erfstep(-.5,.5,v.x)*
                  erfstep(-.5,.5,v.y);
        dQ.xyz += w*a.w*a.xyz;
        dQ.w   += w*a.w;
    }
    if (dQ.w>0.)
    {
        dQ.xyz/=dQ.w;
        Q = dQ;
    }
    Q.xy *= .995;
    Q.w += 4e-4;
    Q.w *= .999;
    Q.w = min(Q.w,2.);
    if (U.x < 1.||U.y<1.||R.x-U.x<1.||R.y-U.y<1.)Q.w = 0.;
    if (_mouse.z>0.&&length(U-_mouse.xy)<10.)Q.w = 1.;
    if (int(FRAMECOUNT) <= 1) {Q = vec4(0,0,0,.5+.01*sin(3.1*U.x)+.01*sin(3.1*U.y));}
	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    Q = A(U);
    vec4 dQ = vec4(0);
    for (float x = -1.; x<=1.;x++)
    for (float y = -1.; y<=1.;y++)
    if(x!=0.||y!=0.)
    {
        vec2 u = vec2(x,y);
        vec4 a = A(U+u);
        vec4 b = C(U+u);
        float f = sqrt(a.w)+b.x;
        dQ.xy -= .5*f*u/dot(u,u);
    }
    Q += dQ;
	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    vec4 a = B(U);
    Q = C(U);
    Q.x -= 2e-2*min(length(a.xy)*sqrt(a.w),1.);
    Q.x += 1e-3;
    
    if (FRAMECOUNT < 100)
        Q.x = D(U).x+20.*(1.-length(U-.5*R)/R.x);

	return Q; 
 } 


// Fork of "Liquid Disco" by wyatt. https://shadertoy.com/view/Dt2XWK
// 2023-02-20 06:14:33

vec4 renderMainImage() 
{
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    vec4 f = A(U);
    Q = f.wwww*(sin(-1.1+6.2*f.z+vec4(1,2,3,4)));
    
    vec4 n = C(U+vec2(0,1));
    vec4 e = C(U+vec2(1,0));
    vec4 s = C(U-vec2(0,1));
    vec4 w = C(U-vec2(1,0));
    Q.xyz = .5+.5*normalize(vec3(e.x-w.x,n.x-s.x,1)).zzz;
    Q *= .5+.5*sin(.1*C(U).x+vec4(1,2,3,4));
    Q.w = 1.;

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