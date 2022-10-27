

			//******** Common Code Begins ********

vec2 R;
vec4 M;
int I;
#define N 9.
#define R3D vec3(R/N,N*N)
#define e(a) mat2(cos(a),-sin(a),sin(a),cos(a))
#define d2(U) ((U).xy+vec2(mod(floor((U).z),N),floor(floor((U).z)/N))*R/N)
#define d3(u) vec3(mod(u,R/N),floor(u/R*N).x+floor(u/R*N).y*N)
#define _3D  vec3 U = d3(u)
#define Sampler vec4 T(vec3 U) {return mix(texture(BuffA,d2(vec3(U.xy,floor(U.z)))/R),texture(BuffA,d2(vec3(U.xy, ceil(U.z)))/R),fract(U.z));}
#define A(U) texture(cha,d2(mod(U,R3D))/R)
#define B(U) texture(chb,d2(mod(U,R3D))/R)
#define Main void mainImage (out vec4 Q, in vec2 U)
float signe (float x) {return atan(100.*x);}
void prog (vec3 U, out vec4 a, out vec4 b, sampler2D cha, sampler2D chb) {
	
    a = vec4(0); b = vec4(0);
    for (int x = -1; x <= 1; x++)
    for (int y = -1; y <= 1; y++)
    for (int z = -1; z <= 1; z++)
    {
        vec3 u = vec3(x,y,z);
    	vec4 aa = A(U+u), bb = B(U+u);
        aa.xyz += bb.xyz;
        #define q 1.125
		vec3 w1 = clamp(aa.xyz-0.5*q,U - 0.5,U + 0.5),
             w2 = clamp(aa.xyz+0.5*q,U - 0.5,U + 0.5);
        float m = (w2.x-w1.x)*(w2.y-w1.y)*(w2.z-w1.z)/(q*q*q);
        aa.xyz = 0.5*(w1+w2);
        a.xyz += aa.xyz*aa.w*m;
        b.xyz += bb.xyz*aa.w*m;
        a.w += aa.w*m;
    }
    if (a.w>0.) {
        a.xyz/=a.w;
        b.xyz/=a.w;
    }
}
void prog2 (vec3 U, out vec4 a, out vec4 b, sampler2D cha, sampler2D chb) {
	
    a = A(U); b = B(U);
    vec3 f = vec3(0); float m = 0.;
    for (int x = -1; x <= 1; x++)
    for (int y = -1; y <= 1; y++)
    for (int z = -1; z <= 1; z++)
    {
        vec3 u = vec3(x,y,z);
        float l = length(u);
        if (l>0.) {
    		vec4 aa = A(U+u), bb = B(U+u);
            f += 1e-2*(aa.w*(1.-.2*aa.w))*u/l;
            m += aa.w;
        }
    }
    if (m>0.) b.xyz += f/m;
    
    
    // Boundaries:
   	b.xyz -= 1e-3*signe(a.w)*(a.xyz-0.5*R3D)*sin(1e-5*float(I));

    
    if (FRAMECOUNT<=1.||U.x<1.||R3D.x-U.x<1.||R3D.y-U.y<1.||R3D.x-U.x<1.||U.z<1.||R3D.z-U.z<1.) {
    	a = vec4(U,0);
        b = vec4(0);
        if (length(U-0.5*R3D) < 0.2*R3D.y) a.w = 26.;
    }
}

			//******** BuffA Code Begins ********

vec4 renderPassA() {
    vec4 Q = vec4 (0.);
	R = RENDERSIZE.xy;
    M = _mouse;
    I = int(FRAMECOUNT);
   	vec4 a, b;
    vec2 U = _xy;
   	prog (d3(U),a,b,BuffC,BuffD);
    
    Q = a;
	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
    vec4 Q = vec4 (0.);
    vec2 U = _xy;

	R = RENDERSIZE.xy;
    M = _mouse;
    I = int(FRAMECOUNT);
   	vec4 a, b;
    
   	prog (d3(U),a,b,BuffC,BuffD);
    
    Q = b;
	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
    vec4 Q = vec4 (0.);
    vec2 U = _xy;

	R = RENDERSIZE.xy;
    M = _mouse;
    I = int(FRAMECOUNT);
   	vec4 a, b;
    
   	prog2 (d3(U),a,b,BuffA,BuffB);
    
    Q = a;
	return Q; 
 } 


			//******** BuffD Code Begins ********

vec4 renderPassD() {
    vec4 Q = vec4 (0.);
    vec2 U = _xy;

	R = RENDERSIZE.xy;
    M = _mouse;
    I = int(FRAMECOUNT);
   	vec4 a, b;
    
   	prog2 (d3(U),a,b,BuffA,BuffB);
    
    Q = b;
	return Q; 
 } 


// Fork of "Michael's Cool Algorithm" by wyatt. https://shadertoy.com/view/ttsyzs
// 2020-06-25 01:31:08

Sampler
vec4 renderMainImage()
{
    vec4 Q = vec4 (0.);

    vec2 U = _xy;
    R = RENDERSIZE.xy;
    vec3 mi = 0.5*vec3(R/N,N*N);
    vec3 p = vec3(0,0,-R.x/N);
    vec3 d = normalize(vec3((U-0.5*R)/R.y,1));
    if (_mouse.z>0.) {
 		p.zx *= e(6.2*_mouse.x/R.x);
		d.zx *= e(6.2*_mouse.x/R.x);
        p.yz *= e(6.2*_mouse.y/R.y);
		d.yz *= e(6.2*_mouse.y/R.y);
    } else {
		p.xz *= e(.1*smoothTime);
		d.xz *= e(.1*smoothTime);
	}
    Q = vec4(0);
    for (int i = 0; i < 100; i++) {
        vec3 o = abs(p)-mi;
        float m = length(max(o,0.));
        if (m<.1)
        { 	
            vec4 a = T(p+mi).wwww;
            float aa = length(a);
            Q += 2e-2*(1.-exp(-aa))*(a)*(0.5+0.5*sin(a+vec4(1,2,3,4)));
            p += d*(.01+exp(-.1*aa*aa));
           //p = mod(p+mi,R3D)-mi;
        } else p += d*m;
        
 	}
    Q = .8*atan(log(1.+Q));
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