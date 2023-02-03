

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
#define C(U) texture(BuffC,(U)/R)
#define D(U) texture(BuffD,(U)/R)
#define Main void mainImage( out vec4 Q, in vec2 U )
#define NeighborhoodD vec4 n = D(U+vec2(0,1)), e = D(U+vec2(1,0)), s = D(U-vec2(0,1)), w = D(U-vec2(1,0)), m = 0.25*(n+e+s+w); 
#define NeighborhoodC vec4 n = C(U+vec2(0,1)), e = C(U+vec2(1,0)), s = C(U-vec2(0,1)), w = C(U-vec2(1,0)), m = 0.25*(n+e+s+w); 
#define NeighborhoodA vec4 n = A(U+vec2(0,1)), e = A(U+vec2(1,0)), s = A(U-vec2(0,1)), w = A(U-vec2(1,0)), m = 0.25*(n+e+s+w); 
#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))
#define grad 0.125*vec2(e.x-w.x,n.x-s.x)
#define div 0.125*(n.y-s.y+e.x-w.x)

#define N 32.
#define For for (float i = -(N); i<=(N); i++)
#define S vec4(2,2,7.+12,8)
#define Gaussian(i) (0.39989422804*GaussianMod)/S*exp(Scale*0.5+(-.5)*(i)*(i)/S/S)
#define Init if (FRAMECOUNT <= 1 || Reset > 0.) 
#define Mouse if (_mouse.z>0.&&abs(length(U-_mouse.xy))<30.) 
vec4 image = vec4(0.);

			//******** BuffA Code Begins ********

vec4 renderPassA(){
    vec4 Q = vec4(0.);
    vec2 U = _xy;
    U += moveXY;

    U-=0.5*R;
    float a = 0.02*sin(.01225*smoothTimeC)*exp(-3.*length(U)/R.y);
    U *= (1.-(Zoom)*(.0024*(1.0+syn_BassLevel*syn_MidLevel*syn_Intensity))*exp(exp(-length(U)/(R.y+R.y))))*rot(a*0.1);
    U+=0.5*R;
    NeighborhoodD
    
    Q = A(U+.25*rot(Flip)*D(U).xy)-div;
    
    Q = sin(Q);
    Q *= 1.0 + (0.35*_loadUserImage()-0.4*_loadUserImageAsMask()*(1.0+0.5*syn_Intensity));

    Mouse Q.x = 1.;
    Init Q = sin((.1+syn_BassLevel)*U.xxxx)*cos(U.y);
	return Q; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB(){
    vec4 Q = vec4(0.);
    vec2 U = _xy;

    Q = vec4(0);
    For Q += Gaussian(i) * A(U+vec2(i,i)).x;
	return Q; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC(){
    vec4 Q = vec4(0.);
    vec2 U = _xy;
    U += Stretch*_uvc*PI;

    Q = vec4(0);
    For Q += Gaussian(i) * B(U+_uvc+vec2(-i,i));
	return Q; 
 } 


			//******** BuffD Code Begins ********

vec4 renderPassD(){
    vec4 Q = vec4(0.);
    vec2 U = _xy;
    NeighborhoodC
    Q.xy = 
        10.*vec2(e.x-w.x,n.x-s.x)+
        6.*vec2(e.y-w.y,n.y-s.y)+
        7.*vec2(e.z-w.z,n.z-s.z)+
        1.*vec2(e.w-w.w,n.w-s.w);
	return Q; 
 } 


vec4 renderMainImage(){
    vec4 Q = vec4(0.);
    vec2 U = _xy;
    NeighborhoodA
    vec3 no = normalize(vec3(grad*10.5*_uvc.x,.5));
    vec3 re = reflect(no,vec3(0-_uvc.x/PI,0+_uvc.y/PI,1));
    re.xz *= rot(.03*smoothTimeB);
    re.zy *= rot(.0125*smoothTimeB)*0.4;
    Q = (0.5*(1.0+syn_HighHits)-01.7*A(U).xxxx)*(0.5+0.7*re.xzyz);
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