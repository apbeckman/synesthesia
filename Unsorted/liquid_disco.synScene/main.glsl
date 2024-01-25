

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
#define Main void mainImage(out vec4 Q, in vec2 U)


// oneshade:
//https://www.shadertoy.com/view/7sKSRh
float erf(in float x) {
    x *= xMultiplier;
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
    vec4 Q = vec4(0.);
    vec2 U = _xy;
    vec4 dQ = Q = vec4(0);
    U -= _uvc*Zoom*(2.0+low);
    U += Drift;
    U += _uvc*Stretch;
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
        Q *= 1.0- _edgeDetectSobel(syn_UserImage);

    Q.xy *= .9;
    Q.w += 1e-14;
    if (_mouse.z>0.&&length(U-_mouse.xy)<10.)Q.w = 1.;
    if (FRAMECOUNT <= 1 || Reset != 0.0) {Q = vec4(0,0,0,.5+.01*sin(3.1*U.x)+.01*sin(3.1*U.y));}
	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB()
{
    vec4 Q = vec4(0.);
    vec2 U = _xy;
        
    Q = A(U);
    vec4 dQ = vec4(0);
    for (float x = -1.; x<=1.;x++)
    for (float y = -1.; y<=1.;y++)
    if(x!=0.||y!=0.)
    {
        vec2 u = vec2(x,y);
        vec4 a = A(U+u);
        vec4 b = B(U+u);
        float f = a.w*(a.w-2.);
        dQ.xy -= (.125+low*0.125+syn_Level*0.25)*f*u/dot(u,u);
    }
    Q += dQ;
	return Q; 
 } 


vec4 renderMainImage()
{
    vec4 Q = vec4(0.);
    vec2 U = _xy;

    vec4 f = A(U);
    //Q = f.wwww*abs(sin(TIME+6.2*f.z+vec4(1,2,3,4))); //original
    Q = f.wwww*abs(sin(f.z+7*vec4(1,2,3,4)));
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