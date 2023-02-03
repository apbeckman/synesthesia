

			//******** Common Code Begins ********

#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
#define C(U) texture(BuffC,(U)/R)
#define D(U) texture(BuffD,(U)/R)
#define Main void mainImage (out vec4 Q, vec2 U) 
#define Init if (FRAMECOUNT < 1) 
#define Border if (U.x<1.||R.x-U.x<1.||U.y<1.||R.y-U.y<1.)
#define T(U) A((U)-dt*A(U).xy)
#define NeighborhoodT vec4 n = T(U+vec2(0,1)), e = T(U+vec2(1,0)), s = T(U-vec2(0,1)), w = T(U-vec2(1,0)), m = 0.25*(n+e+s+w);
#define Neighborhood vec4 n = A(U+vec2(0,1)), e = A(U+vec2(1,0)), s = A(U-vec2(0,1)), w = A(U-vec2(1,0)), m = 0.25*(n+e+s+w);
#define grd 0.25*vec2(e.z-w.z,n.z-s.z)
#define div 0.25*(e.x-w.x+n.y-s.y)
#define dt (int(FRAMECOUNT)<500?1.:0.1)
#define I 4
//Dave H :
vec3 hash33(vec3 p3)
{
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy + p3.yxx)*p3.zyx);
}
float pie (vec2 p, vec2 a, vec2 b) {
	vec2 m = 0.5*(a+b); // midpoint
    if (length(a-b)<1e-3) return 1e3; // ignore self
	return abs(dot(p-m,b-m)/dot(b-m,b-m)); // pojection
} 

			//******** BuffA Code Begins ********

vec4 renderPassA() {
    vec4 Q = vec4 (0.);
    vec2 U = _xy;
	// Fluid:
if (int(FRAMECOUNT)%2==0) {
    Q = T(U);
    NeighborhoodT;
    Q.xy -= dt*grd;
    Init Q = vec4(0);
    Border Q.xy *= 0.;
} else {
    Q = A(U);
    Neighborhood;
    Q.z -= dt*div;
}
	vec2 mo = 0.5*R;
    if (_mouse.z>0.) mo = _mouse.xy;
    Q = mix(Q,vec4(0.5*sin(dt*2.*TIME)*
            vec2(sin(dt*TIME),0),
            Q.z,
            1.
        ),
     exp(-.2*length(U-mo)));
	return Q; 
 } 


			//******** BuffB Code Begins ********

// Translate Coordinate by gradient of pressure
// Strength of translation determined by wavelength
vec4 renderPassB() {
    vec4 Q = vec4 (0.);
    vec2 U = _xy;

    if (int(FRAMECOUNT)%I==0) {
        Neighborhood;
        Q.xyz = hash33(vec3(U,int(FRAMECOUNT)));
        Q.xy = U+Q.xy*2.-1.;
        Q.xy += 3e3*grd*(1.+Q.z);
    } else Q = B(U);
	return Q; 
 } 


			//******** BuffC Code Begins ********

// Hierarchical Sort
void X (inout vec4 Q, inout vec4 r, vec2 U, vec2 u) {
        vec4 b,c;vec2 i; float l;
        if (int(FRAMECOUNT)%I==0) {
            i = U+u; 
            b = B(i);
            l = length(b.xy-U);
            if (l<r.x) {
                Q.xy = i;
                r.x = l;
            } else if (l<r.y) {
                Q.zw = i;
                r.y = l;
            }
        } else {
            c = C(U+u);
            b = B(c.xy);
            vec4 bb = B(Q.xy);
            i = c.xy;
            l = length(b.xy-U);
            if (l<r.x) {
                Q.xy = i;
                r.x = l;
            } else if (l<r.y) {
                Q.zw = i;
                r.y = l;
            }
            b = B(c.zw);
            bb = B(Q.xy);
            i = c.zw;
            l = length(b.xy-U);
            if (l<r.x) {
                Q.xy = i;
                r.x = l;
            } else if (l<r.y) {
                Q.zw = i;
                r.y = l;
            }
        }
    	
}
vec4 renderPassC() {
    vec4 Q = vec4 (0.);
    vec2 U = _xy;

    Q = vec4(1e9);
    vec4 r = vec4(1e9);
    if (int(FRAMECOUNT)%I<0) {
    	Q = C(U);
        r = vec4(
        	length(U-B(Q.xy).xy),
        	length(U-B(Q.zw).xy),
            0,0
        );
    }
    float k = exp2(float(I-(int(FRAMECOUNT)%I)));
    X(Q,r,U,vec2(0,k));
    X(Q,r,U,vec2(k,0));
    X(Q,r,U,vec2(0,-k));
    X(Q,r,U,vec2(-k,0));
    
    X(Q,r,U,vec2(k,k));
    X(Q,r,U,vec2(-k,k));
    X(Q,r,U,vec2(k,-k));
    X(Q,r,U,vec2(-k,-k));
	return Q; 
 } 


			//******** BuffD Code Begins ********

// Draw Photons to the screen
vec4 P (vec2 U, vec3 p) {
	return 0.0033/(1.+dot(U-p.xy,U-p.xy))*max(cos(p.z*6.2+vec4(1,2,3,4)),0.);
}
vec4 renderPassD() {
    vec4 Q = vec4 (0.);
    vec2 U = _xy;
    Q = D(U);
    if (int(FRAMECOUNT)%I==0) Q *= 0.9;
    for (int x = -3; x <=3; x++) {
        for (int y = -3; y<=3; y++) {
            vec2 u = U+vec2(x,y);
        	vec4 c = C(u);
            vec4 b = B(c.xy);
            Q += P(U,b.xyz);
            b = B(c.zw);
            Q += P(U,b.xyz);
            b = B(u);
            Q += P(U,b.xyz);
        }
    }
	return Q; 
 } 


vec4 renderMainImage() {
    vec4 Q = vec4 (0.);
    vec2 U = _xy;

	if (int(FRAMECOUNT)%I<I-1) discard;
    Q = A(U);
    Q *= sqrt(Q);
    
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