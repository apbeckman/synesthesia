

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
//#define Main void mainImage(out vec4 Q, in vec2 U)
#define box for(int x=-1;x<=1;x++)for(int y=-1;y<=1;y++)
float dt = .2;
float _a = .05;
float _b = .05;

#define A(U) texture(BuffB,(U)/R)
#define B(U) texture(BuffD,(U)/R)
#define A1(U) texture(BuffA,(U)/R)
// /#define A(U) texture(BuffA,(U)/R)

#define B2(U) texture(BuffD,(U)/R)

#define A3(U) texture(BuffC,(U)/R)
#define B3(U) texture(BuffB,(U)/R)

			//******** BuffA Code Begins ********

// Forces 🦜
vec4 renderPassA()
{
	vec4 Q = vec4(0.0);
    vec2 U = _xy;

    Q = A(U);
    vec4 b = B(U);
    Q.w = b.w;
    vec4 T = Q;
    box if(abs(x)!=abs(y))
    {
        vec2 u = vec2(x,y)*(1.0+basshits);
        vec4 a = A(U+u);
        float f = dt*.25*(a.w*(a.w-1.+a.z)*(1.0+basshits));
        Q.xy -= f*u.xy;
        Q.z  += .125*.25*(a.z-T.z)-Q.w*f*dot(a.xy-T.xy+u.xy,u.xy); 
    }
    float M = dt*Q.z*Q.w*b.x*b.y*(1.0+basshits);
    Q.z += 10.*M;
    Q.y -= .1/R.y;
    if (_mouse.z>0.&&length(U-_mouse.xy)<50.)
        Q.x = .5;
    if (FRAMECOUNT <= 1)Q = vec4(0,0,.1,.1);
    if (U.x < 4.||R.x-U.x<4.) Q.xyz *= 0.;
    if (U.y < 4.||R.y-U.y<4.) Q.xyz *= 0.;
    Q = clamp(Q,vec4(-.5,-.5,0,0),vec4(.5,.5,2,Q.w));
    if (length(Q.xy)>.5-basshits*0.45) Q.xy = (.5)*normalize(Q.xy);
    
	return Q; 
 } 


			//******** BuffB Code Begins ********

// Advect 🐿
vec4 renderPassB()
{
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    Q = A1(U);
    vec4 dQ = vec4(0);
    box if(abs(x)!=abs(y))
    {
        vec2 u = vec2(x,y);
        vec4 q = A1(U+u);
        vec2 a = Q.xy,
             b = q.xy+u;
       float ab = dot(u,b-a);
       float i = dot(u,(0.5*u-a))/ab;
       float j = .5+.5*max(1.-5.*Q.w*q.w,0.);
       float k = .5+.5*max(1.-5.*Q.w*q.w,0.);
       float wa = 0.25*Q.w*min(i,j)/j;
       float wb = 0.25*q.w*max(k+i-1.,0.)/k;
        dQ.xyz += Q.xyz*wa+q.xyz*wb;
        dQ.w += wa+wb;
        
    }
    if (dQ.w>0.)dQ.xyz/=dQ.w;
    Q = dQ;
    
    if (U.x < 4.||R.x-U.x<4.) Q.xy *= 0.;
    if (U.y < 4.||R.y-U.y<4.) Q.xy *= 0.;
	return Q*(1.0+basshits*0.001); 
 } 


			//******** BuffC Code Begins ********

// Advect 🐿
vec4 renderPassC()
{
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    Q = A(U);
    vec4 Qb = B(U);
    vec4 dQ = vec4(0);
    box if(abs(x)!=abs(y))
    {
        vec2 u = vec2(x,y);
        vec4 q = A(U+u),
             qb = B(U+u);
        vec2 a = Q.xy,
             b = q.xy+u;
       float ab = dot(u,b-a);
       float i = dot(u,(0.5*u-a))/ab;
       float j = .5;//+.2*max(1.-5.*Q.w*q.w,0.);
       float k = .5;//+.2*max(1.-5.*Q.w*q.w,0.);
       float wa = 0.25*Q.w*min(i,j)/j;
       float wb = 0.25*q.w*max(k+i-1.,0.)/k;
        dQ.xyz += Qb.xyz*wa+qb.xyz*wb;
        dQ.w += wa+wb;
        
    }
    if (dQ.w>0.)dQ.xyz/=dQ.w;
    Q = dQ;
	return Q; 
 } 


			//******** BuffD Code Begins ********

// Reaction 🔥
vec4 renderPassD() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    Q = A3(U);
    vec4 b = B3(U);
    Q.w = b.w;
    
    float M = dt*b.z*b.w*Q.x*Q.y;
    Q.x += -_a*M;
    Q.y += -_b*M;
    Q.z += M;
    Q.w += (1.-_a-_b)*M*(1.0+basshits);
    
    
    if (FRAMECOUNT<=1) {
        float l = length(U-vec2(.5,.1)*R);
        Q = vec4(step(l,.1*R.y),
                 1,
                 0,.5);
    }
	return Q; 
 } 


// Fork of "Colors!" by wyatt. https://shadertoy.com/view/7dX3z7
// 2021-06-02 16:58:24

// Fork of "Transport Dynamics II" by wyatt. https://shadertoy.com/view/sdl3RN
// 2021-03-14 01:41:52

// Display 🌵
vec4 renderMainImage()
{
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    vec4 f = A1(U), c = A3(U);
    Q = c*min(f.w,1.3);
    Q = ((.5-.5*cos(5.*c.x+c.y*vec4(1.,2,3,4))))*atan(f.wwww);
    Q += f.z*f.w*vec4(1,.5,0,1);
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