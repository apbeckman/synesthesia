

			//******** Common Code Begins ********



// Scene control:
//
// $0003	BA	1..0	- camera distance 0..3
// $000C	DC	3..2	- scenery 	00:forest		01:city
//									10:terrain		11:landscape
//	 								-0:circle		-1:triangle
//
// $0010	 E	 4		- 0:glowing		1:solid
// $0020	 F	 5		- large voronoi range / end phase of city destruction
// $0040	 G	 6		- city destruction
// $0080	 H	 7		- black sphere (only flag H = black screen)
// $0100	 I	 8		- camera angle
// $0200	 J	 9		- voronoi destruction
// $0400	 K	10		- camera distance + 0.5
// $0800	 L  11		- water (with C=0 D=1)
// $1000	 M  12		- darkness
// $2000	 N	13		- sharp OUT transition
// $4000	 O	14		- sharp IN transition
// $8000	 P	15		- shape particles
//
//


const ivec3 SCENES[] = ivec3[](
#define SCENE(n,d,b)	ivec3((d),(b),(n))

	//intro
	SCENE(2, 16, 0x050A | 0x0000),
	SCENE(3, 16, 0x040B | 0x0000),

	//PART 1 - forest+rocks
	SCENE(4, 12, 0x0013 | 0x4000),
	SCENE(5, 11, 0x0051 | 0x2000),
	SCENE(6, 1, 0x0410 | 0x6000),
	SCENE(7, 8, 0x0010 | 0x6000),
	SCENE(8, 12, 0x0110 | 0x4000),
	SCENE(9, 12, 0x001F | 0x0000),
	SCENE(10, 12, 0x001C | 0x6000),

	//the ball!
	SCENE(11, 8, 0x0080 | 0x0000), //d
	SCENE(13, 20, 0x0418 | 0x2000),

	//PART 2 - city
	SCENE(14, 16, 0x0407 | 0x4000),
	SCENE(15, 16, 0x0014 | 0x0000),
	SCENE(16, 20, 0x0124 | 0x6000),
	SCENE(17, 11, 0x0044 | 0x6000),

	//the ball!
	SCENE(18, 13, 0x0080 | 0x0000), //d        
	SCENE(20, 20, 0x0818 | 0x4000),

	//PART 3 - destruction
	SCENE(80, 16, 0x0041 | 0x2000),
	SCENE(21, 16, 0x0200 | 0x6000),
	SCENE(33, 16, 0x0074 | 0x0000),
	SCENE(22, 16, 0x040c | 0x7000),
	SCENE(24, 19, 0x063c | 0x3000),

	SCENE(26, 7, 0x0080 | 0x0000), //d

								   //outtro - water
	SCENE(28, 12, 0x1908 | 0x0000),
	SCENE(27, 16, 0x0808 | 0x1000),
	SCENE(97, 10, 0x0C08 | 0xB000),

	SCENE(96, 8, 0x0080 | 0x0000), //d

								   //the ball! 
//	SCENE(30, 12, 0x0018 | 0x4000),
//	SCENE(31, 16, 0x0098 | 0x0000),
	SCENE(32, 24, 0x0098 | 0x6000),
	SCENE(0, 255, 0x0080 | 0x0000) //d	
#undef SCENE
	);



// ================================ Helper functions ================================


float sat(float x)
{
	return clamp(x, 0., 1.);
}


// ================================ Noises ================================

// 3D noise - method by iq/shane
float noise(vec3 p)
{
	vec3 ip=floor(p);
	p-=ip;
	vec3 s=vec3(7, 157, 113);
	vec4 h=vec4(0, s.yz, s.y+s.z)+dot(ip, s);
	p=p*p*(3.-2.*p);
	h=mix(fract(43758.5*sin(h)), fract(43758.5*sin(h+s.x)), p.x);
	h.xy=mix(h.xz, h.yw, p.y);
	return mix(h.x, h.y, p.z);
}


// ================================ Voronoi ================================

vec2 vcore(vec2 m, vec2 p, vec2 s)
{
	vec2 c = floor(2.5*p+s);		// 1./.4   r
	c += fract(43758.5*sin(c+17.*c.yx));

	float v = length(.4*c-p);	// r
	return v<m.x ? vec2(v, m.x) : v<m.y ? vec2(m.x, v) : m;
}



// ================================ Patterns ================================


float lattice(vec3 p)
{
	p=abs(p);
	p=max(p, p.yzx);
	p=min(p, p.yzx);
	p=min(p, p.yzx);
	return p.x;
}



// ================================ SDF merge functions ================================


void dmin(inout vec3 d, float x, float y, float z)
{
	if( x < d.x ) d = vec3(x, y, z);
}



// ================================ Domain operations ================================


// rotation
void pR(inout vec2 p, float a)
{
	a *= 6.283;
	p = cos(a)*p+sin(a)*vec2(p.y, -p.x);
}


// 3D repetition
vec3 rep(vec3 p, float r)
{
	return (fract(p/r-.5)-.5)*r;
}

// diffuse reflection hash - method by fizzer
vec3 hashHs(vec3 n, inout float seed)
{
	vec2 uv = (seed=32.+seed*fract(seed))+vec2(78.233, 10.873);
	uv = fract(.1031*uv);
	uv *= 19.19+uv;
	uv = fract(2.*uv*uv);

	float u = 2.*uv.x-1.;

	vec3 v = vec3(sqrt(1.-u*u), 0., u);
	pR(v.xy, uv.y);
	return normalize(n+v);
}


// ================================ Complex SDFs ================================

float vines(vec3 p, float s)
{
	p.y=abs(p.y);
	pR(p.xz, .1*p.y); p=abs(p); p.xz -= .06*s;
	pR(p.xz, -.16*p.y); p=abs(p); p.xz -= .05*s;
	pR(p.xz, .4*p.y);
	return length(abs(p.xz) - .04*(s*.5+.5));
}


			//******** BuffA Code Begins ********

#define R RENDERSIZE.xy
#define D(U) texture(BuffD,(U)/R)
vec4 X (inout vec2 c,inout float m, vec2 u, vec2 r) {
    vec4 n = D(u+r);
    m += n.z;
	if (length(u+r-n.xy)<length(u-c)) c = n.xy;
    return n;
} 
vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = D(U);
    vec2 c = Q.xy;
    float m=0.;
    vec4 
        n = X(c,m,U,+vec2(0,1)),
        e = X(c,m,U,+vec2(1,0)),
        s = X(c,m,U,-vec2(0,1)),
        w = X(c,m,U,-vec2(1,0));
    X(c,m,U,vec2(1,1));
    X(c,m,U,vec2(1,-1));
    X(c,m,U,vec2(-1,1));
    X(c,m,U,vec2(-1,-1));
    vec2 g = vec2(e.z-w.z,n.z-s.z);
    Q.xy = c;
    Q.z += (m/8.-Q).z+.05*Q.w - .00001*Q.z;
    Q.w -= 0.0001*Q.w;
    Q.zw = max(Q.zw,vec2(2,1)*smoothstep(4.,0.0,length(U-c)));
    Q.xy -= 0.25*g;
    if (length(U-_mouse.xy)<.1*R.y&&_mouse.z>0.) Q = vec4(-R,Q.z,0);
    if (FRAMECOUNT <= 1) Q = vec4(clamp(floor(U/2.)*2.,0.5*R-2.,0.5*R+2.),0,0);
    
	return Q; 
 } 


			//******** BuffB Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
/*
vec4 X (inout vec2 c,inout float m, vec2 u, vec2 r) {
    vec4 n = A(u+r);
    m += n.z;
	if (length(u+r-n.xy)<length(u-c)) c = n.xy;
    return n;
} 
*/
vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = A(U);
    vec2 c = Q.xy;
    float m=0.;
    vec4 
        n = X(c,m,U,+vec2(0,1)),
        e = X(c,m,U,+vec2(1,0)),
        s = X(c,m,U,-vec2(0,1)),
        w = X(c,m,U,-vec2(1,0));
    X(c,m,U,vec2(1,1));
    X(c,m,U,vec2(1,-1));
    X(c,m,U,vec2(-1,1));
    X(c,m,U,vec2(-1,-1));
    vec2 g = vec2(e.z-w.z,n.z-s.z);
    Q.xy = c;
    Q.z += (m/8.-Q).z+.05*Q.w - .0001*Q.z;
    Q.w -= 0.001*Q.w;
    Q.zw = max(Q.zw,vec2(2,1)*smoothstep(4.,0.0,length(U-c)));
    Q.xy -= 0.25*g;
    if (length(U-_mouse.xy)<.1*R.y&&_mouse.z>0.) Q = vec4(-R,Q.z,0);
    if (FRAMECOUNT <= 1) Q = vec4(clamp(floor(U/2.)*2.,0.5*R-2.,0.5*R+2.),0,0);
    
	return Q; 
 } 


			//******** BuffC Code Begins ********

#define R RENDERSIZE.xy
#define B(U) texture(BuffB,(U)/R)
/*
vec4 X (inout vec2 c,inout float m, vec2 u, vec2 r) {
    vec4 n = B(u+r);
    m += n.z;
	if (length(u+r-n.xy)<length(u-c)) c = n.xy;
    return n;
} 
*/
vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = B(U);
    vec2 c = Q.xy;
    float m=0.;
    vec4 
        n = X(c,m,U,+vec2(0,1)),
        e = X(c,m,U,+vec2(1,0)),
        s = X(c,m,U,-vec2(0,1)),
        w = X(c,m,U,-vec2(1,0));
    X(c,m,U,vec2(1,1));
    X(c,m,U,vec2(1,-1));
    X(c,m,U,vec2(-1,1));
    X(c,m,U,vec2(-1,-1));
    vec2 g = vec2(e.z-w.z,n.z-s.z);
    Q.xy = c;
    Q.z += (m/8.-Q).z+.05*Q.w - .0001*Q.z;
    Q.w -= 0.001*Q.w;
    Q.zw = max(Q.zw,vec2(2,1)*smoothstep(4.,0.0,length(U-c)));
    Q.xy -= 0.25*g;
    if (length(U-_mouse.xy)<.1*R.y&&_mouse.z>0.) Q = vec4(-R,Q.z,0);
    if (FRAMECOUNT <= 1) Q = vec4(clamp(floor(U/2.)*2.,0.5*R-2.,0.5*R+2.),0,0);
    
	return Q; 
 } 


			//******** BuffD Code Begins ********

#define R RENDERSIZE.xy
#define C(U) texture(BuffC,(U)/R)
/*
vec4 X (inout vec2 c,inout float m, vec2 u, vec2 r) {
    vec4 n = C(u+r);
    m += n.z;
	if (length(u+r-n.xy)<length(u-c)) c = n.xy;
    return n;
}
*/ 
vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    Q = C(U);
    vec2 c = Q.xy;
    float m=0.;
    vec4 
        n = X(c,m,U,+vec2(0,1)),
        e = X(c,m,U,+vec2(1,0)),
        s = X(c,m,U,-vec2(0,1)),
        w = X(c,m,U,-vec2(1,0));
    X(c,m,U,vec2(1,1));
    X(c,m,U,vec2(1,-1));
    X(c,m,U,vec2(-1,1));
    X(c,m,U,vec2(-1,-1));
    vec2 g = vec2(e.z-w.z,n.z-s.z);
    Q.xy = c;
    Q.z += (m/8.-Q).z+.05*Q.w - .0001*Q.z+basshits;
    Q.w -= 0.001*Q.w;
    Q.zw = max(Q.zw,vec2(2,1)*smoothstep(4.,0.0,length(U-c)));
    Q.xy -= 0.25*g;
    if (length(U-_mouse.xy)<.1*R.y&&_mouse.z>0.) Q = vec4(-R,Q.z,0);
    if (FRAMECOUNT <= 1) Q = vec4(clamp(floor(U/2.)*2.,0.5*R-2.,0.5*R+2.),0,0);
    
	return Q; 
 } 


#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
vec4 X (inout vec2 c, vec2 u, vec2 r) {
    vec4 n = A(u+r);
	if (length(u+r-n.xy)<length(u-c)) c = n.xy;
    return n;
} 
vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    vec4 a = A(U);
    vec2 c = Q.xy;
    vec4 
        n = X(c,U,+vec2(0,1)),
        e = X(c,U,+vec2(1,0)),
        s = X(c,U,-vec2(0,1)),
        w = X(c,U,-vec2(1,0)),
        m = 0.125*(n+e+s+w);
 	vec3 g = normalize(vec3(e.z-w.z,n.z-s.z,.3));
    g = reflect(g,vec3(0,0,1));
 	vec3 b = normalize(vec3(e.w-w.w,n.w-s.w,1));
    float d = dot(g,normalize(vec3(0,1,.5)));
    Q = (exp(-4.*d*d))*m.w*abs(sin(2.+vec4(1,2,3,4)*(1.+2.*m.z)));

    Q = .8+.2*g.x-.8*(1.+0.5*(b.x+b.y)) * a.w * sin(2.+0.5*(g.z)*vec4(1,2,3,4));
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