

			//******** Common Code Begins ********


// Dave Hoskins
// https://www.shadertoy.com/view/4djSRW
float hash13(vec3 p3)
{
	p3  = fract(p3 * .1031);
    p3 += dot(p3, p3.zyx + 31.32);
    return fract((p3.x + p3.y) * p3.z);
}
float dist(vec2 p0, vec2 pf){return sqrt((pf.x-p0.x)*(pf.x-p0.x)+(pf.y-p0.y)*(pf.y-p0.y));}
float d2 = dist(RENDERSIZE.xy*0.5,_xy.xy)*0.00125;


// Inigo Quilez
// https://iquilezles.org/articles/distfunctions
float sdBox( vec3 p, vec3 b ) {
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

// rotation matrix
mat2 rot(float a) { return mat2(cos(a),-sin(a),sin(a),cos(a)); }

#define repeat(p,r) (mod(p,r)-r/2.)

			//******** BuffA Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
#define C(U) texture(BuffC,(U)/R)
#define D(U) texture(BuffD,(U)/R)
vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    U += Drift;
    U += Stretch*_uvc*PI+Zoom*_uvc*PI*0.5;
    U += _uvc*Zoom*PI*0.5;
    vec2 u = .5+.1*vec2(sin(.25*twist_time),cos(.25*twist_time));
    if (_mouse.z>0.) u = _mouse.xy/R;
    U -= u*R; 
    U += abs(d2)*_uvc*Fisheye*PI* (1.0 +  low);

    float r = length(U)/R.y,
        a = .01*sin(.05*twist_time)/(1.+5.*r);

    U *= (.999)*mat2(cos(a),-sin(a),sin(a),cos(a));
    U += u*R;
    Q = A(U);
    vec4
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0));
    Q = mix(Q,0.25*(n+e+s+w),.1);
    Q = mix(Q,D(U),.0005);
    if (FRAMECOUNT <= 1) Q = vec4(.0005)*sin(U.y/R.y*60.+U.x/R.x);
    if (_mouse.z>0.&&length(U-_mouse.xy)<2.) Q=vec4(.03);
    
	return Q; 
 } 


			//******** BuffB Code Begins ********
/*
#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(iChannel1,(U)/R)
#define C(U) texture(iChannel2,(U)/R)
#define D(U) texture(BuffD,(U)/R)
*/
vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = A(U);
    vec4
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0));
    Q.xy = vec2(-e.x+w.x,-n.x+s.x);

    if (length(Q.xy)>0.) Q.xy = mix(Q.xy,normalize(Q.xy),.023);
	return Q; 
 } 


			//******** BuffC Code Begins ********
/*
#define R RENDERSIZE.xy
#define A(U) texture(BuffB,(U)/R)
#define B(U) texture(iChannel1,(U)/R)
#define C(U) texture(iChannel2,(U)/R)
#define D(U) texture(iChannel3,(U)/R)
*/
vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = B(U);
    vec4
        n = B(U+vec2(0,1)),
        e = B(U+vec2(1,0)),
        s = B(U-vec2(0,1)),
        w = B(U-vec2(1,0));
        
    Q = .25*(n+e+s+w);
    
	return Q; 
 } 


			//******** BuffD Code Begins ********
/*
#define R RENDERSIZE.xy
#define A(U) texture(BuffC,(U)/R)
#define B(U) texture(iChannel1,(U)/R)
#define C(U) texture(iChannel2,(U)/R)
#define D(U) texture(iChannel3,(U)/R)
*/
vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = C(U);
    vec4
        n = C(U+vec2(0,1)),
        e = C(U+vec2(1,0)),
        s = C(U-vec2(0,1)),
        w = C(U-vec2(1,0));
    Q = vec4(.25)*(n.y-s.y+e.x-w.x);
	return Q; 
 } 

/*
#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(iChannel1,(U)/R)
#define C(U) texture(BuffC,(U)/R)
#define D(U) texture(BuffD,(U)/R)
*/
vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    vec4
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0));
    vec3 no = normalize(vec3(-e.x+w.x,-n.x+s.x,.000001));
    vec3 re = reflect(normalize(vec3(0,0,1)),no);
    
    Q = vec4(0.7+0.5*no.x)*(0.25+.25*(sin(120.*D(U).x*vec4(1,2,3,4)))*(1.0+syn_HighLevel*(1.0*syn_Intensity)));
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