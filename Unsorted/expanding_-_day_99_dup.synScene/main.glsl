vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** Common Code Begins ********

#define PI acos(-1.)

#define T true
#define F false

#define init vec2 s = vec2(0)
#define ret return s
#define bpm 140.

#define tick (60./bpm)

#define beat tick
#define bar  beat*4.
#define hbar  bar/2.
#define msr  bar*4.
#define hbeat tick/2.
#define qbeat tick/4.





void pR(inout vec2 p, float a) {
    p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

// Reflect space at a plane
float pReflect(inout vec3 p, vec3 planeNormal, float offset) {
    float t = dot(p, planeNormal)+offset;
    if (t < 0.) {
        p = p - (2.*t)*planeNormal;
    }
    return sign(t);
}

int Type=5;
vec3 nc,pab,pbc,pca;
void initIcosahedron() {
    float cospin=cos(PI/float(Type)), scospin=sqrt(0.75-cospin*cospin);
	nc=vec3(-0.5+test2,-cospin,scospin);
	pab=vec3(0.,0.,1.);
	pbc=vec3(scospin,0.,0.5);
	pca=vec3(0.,scospin,cospin);
	pbc=normalize(pbc+vec3(0., 0.+Fold*0.5,Fold));	pca=normalize(pca+vec3(0., 0+Fold*0.5,Fold));
}

vec3 bToC(vec3 A, vec3 B, vec3 C, vec3 barycentric) {
	return barycentric.x * A + barycentric.y * B + barycentric.z * C;
}




vec3 pIcosahedron(inout vec3 p, int subdivisions) {
    p = abs(p);
	pReflect(p, nc, 0.);
    p.xy = abs(p.xy);

	pReflect(p, nc, 0.);
    p.xy = abs(p.xy);

	pReflect(p, nc, 0.);
    float smoothSin = (sin(smoothTimeC*0.125)+sin(smoothTimeC*0.125+1.)-1.)*0.25+0.5;
    if (subdivisions > 0) {

        vec3 A = pbc;
       	vec3 C = reflect(A, normalize(cross(pab, pca)));
        vec3 B = reflect(C, normalize(cross(pbc, pca)));
       
        vec3 n;


        float d = .5;
        
        vec3 p1 = bToC(A, B, C, vec3(1.-d, .0, d));
        vec3 p2 = bToC(A, B, C, vec3(1.-d, d, .0));
        n = normalize(cross(p1, p2));
        pReflect(p, n, 0.);
        
        if (subdivisions > 1) {

            A = reflect(A, n);
            B = p1-smoothSin*0.5*Sin;
            C = p2-smoothSin*0.4*Sin;

            p1 = bToC(A, B, C, vec3(.5, .0, .5));
            p2 = bToC(A, B, C, vec3(.5, .5, .0));
            n = normalize(cross(p1, p2));
            pReflect(p, n, 0.);
            

			p2 = bToC(A, B, C, vec3(.0, .5, .5));
            p1 = bToC(A, B, C, vec3(.5, .5, .0));
            n = normalize(cross(p1, p2));
            pReflect(p, n, 0.);
        }
    }
    
    return pca;
}






			//******** BuffA Code Begins ********


#define pal(a,b,c,d,e) ((a) + (b)*sin((c)*(d) + (e)))
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
float kick;
vec4 r11(float x){return texture(image30,vec2(x)/512.);}
float noise( in vec2 x )
{
    vec2 p = floor(x);
    vec2 f = fract(x);
	f = f*f*(3.0-2.0*f);
	
	vec2 uv = (p.xy) + f.xy;
	return textureLod( image30, (uv+ 0.5)/512.0, 0.0 ).x;
}

float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
	f = f*f*(3.0-2.0*f);
	
	vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
	vec2 rg = textureLod( image30, (uv+ 0.5)/512.0, 0.0 ).yx;
	return mix( rg.x, rg.y, pow(f.z, 0.05) );
}
vec3 glow = vec3(0);

vec3 ro;
#define pmod(p,x) mod(p,x) - 0.5*x
vec2 map(vec3 p, float t){
	vec2 d = vec2(10e6);

    
    vec3 q = p;
    
    
   	vec3 n = pIcosahedron(p, int(1.));
    n = pIcosahedron(p, int(1.));
    
    float dP = length(p.yx) - 0.01;
    
    p.x -= 0.12;
    
    n = pIcosahedron(p, int(2.));
    
	
    
    p.z = log(p.z + 0.) - (smoothTime*0.1);
    
    p.z = sin(0.5*p.z*6.28) - 0.5;
    
    float modD = 1.;
    
    float id = floor(p.z/modD);
    
    //p.z = pmod(p.z, modD);
    
    //p.z -= 0.1;
    
    p.x *= log(p.y)*1.;
    //p.xz *= rot(-1.4);
    p.x -= 0.1;
    p.xz *= rot(-1.4);
    
    //p.x -= 0.4;
    //p.y -= 0.4;
    
    //d.x =  length(p) - 0.44;
    
    
    
    float dc = length(p.z) - 0.06;
    
    
    
    //dc = max(dc, -abs(p.x) + 0.09);
    dc = max(dc, -abs(p.x) + 0.7);
    //dc = max(dc, -abs(p.y) + 0.14);
    //dc = max(dc, -abs(p.y) + 0.44);
    d.x =  min(d.x, dc);
    
    
    d.x = abs(d.x) + 0.0001;
    
    //vec3 c = pal(0.1, 0.5, vec3(0.9,0.2,0.1), 0.6,4.7 + t*6. + TIME + p.z*0.2);
    //vec3 c = pal(0.4, 0.5, vec3(0.9,0.2,0.05), 1.4, t*2. + TIME + p.z*0.2);
    vec3 c = pal(0.4, 0.5, vec3(1.4,0.9,0.3), 1.4, t*2. + smoothTimeB*1.25 + p.z*0.2);
    c = max(c, 0.);
    vec3 ga = 0.0007/(0.1 + d.x*d.x*d.x*d.x*d.x*20000000000.) * c;
    
    vec3 gb = 0.0004/(0.1 + d.x*d.x*d.x*d.x*d.x*2000.) * pal(0.6,0.1, vec3(0.7,0.2,0.5),0.67,0.2);
    
    //glow += mix(ga,gb, smoothstep(0.,1.,length(ga)*400.));
    glow += ga;
    
    d.x =  min(d.x, abs(dP) + 0.01);
    float at = pow(abs(sin(q.z - (smoothTimeB*0.00125) + n.y)), 100.)*smoothstep(0.,1.,length(q.z)*0.4);
    //glow  += 0.001/(0.006 + dP*dP*dP*100.) * pal(0.5,0.5, vec3(0.7,0.2,0.5),0.67,0.2)*at;
    
    d.x *= 0.6;
    
    return d;
}
float it = 0.;
float dith = 0.;
vec2 march(vec3 ro, vec3 rd,inout vec3 p, inout float t, inout bool hit){
	vec2 d = vec2(10e6);
    
    t = 0.; hit = false; p = ro;
    
    for(it = 0.; it < 90. + min(float(FRAMECOUNT), 0.) ; it++){
    	d = map(p,t);
     	d.x += 0.01*pow(smoothstep(1.,0.,t*0.2), 2.5);
        d.x *= dith;
        t += d.x;
    	p = ro + rd*t;
    }
    
    
	return d;
}

vec3 getRd(vec3 ro,vec3 lookAt,vec2 uv){
	vec3 dir = normalize(lookAt - ro);
    vec3 right = normalize(cross(vec3(0,1,0), dir));
    vec3 up = normalize(cross(dir, right));
	return normalize(dir + (right*uv.x + up *uv.y)*0.8);
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    initIcosahedron();
    
    
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    dith = mix(0.9,1.,texture(image30, 4000.*(uv + 0.5 + smoothTimeC*0.15)/256.).x);
    uv *= 1. + dot(uv,uv)*0.2;
    vec3 col = vec3(0);

    
    vec3 ro = vec3(0);
    
    
    float sp = 0.3;
    
    float tt = ((smoothTime*0.35))*sp ;
    ro = vec3(cos(tt+ iMouse.x/RENDERSIZE.x*2.),sin(tt+ iMouse.y/RENDERSIZE.x*2.)*1.6,sin(tt+ iMouse.x/RENDERSIZE.x*2.))*(4.1);
    
    vec3 lookAt = vec3(0);
    
    vec3 rd = getRd(ro, lookAt, uv);
    
    bool hit; float t; vec3 p;
    vec2 d = march(ro, rd, p, t, hit);
    
    
    
    
    
    glow = smoothstep(0.,1.,glow*1.);
    //glow = smoothstep(0.,1.,glow*1.);
    col += glow;
    
    
    //col -= smoothstep(0.,1.,it*0.0009)*0.8;
    //col -= pow(smoothstep(0.,1.,it*0.0045), 10.)*2.4;
    
    
    //col = mix(col, vec3(0.,0.1,0.2)*0.2,pow(smoothstep(0.,1.,t*0.2), 2.));
    //col = mix(col, vec3(0.1,0.,0.1)*0.002,pow(smoothstep(0.,1.,t*0.1), 1.));
    
    
    col = max(col, 0.);
    col *= 9.;
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


// All geodesic reflection functions and are from tdhooper! shadertoy.com/view/4tG3zW
// They are in common tab! 

// radiual blur in this buffer

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    uvn *=_uvc;
	fragColor = vec4(0);
    // Radial blur
    float steps = 16.;
    float scale = 0.00 + pow(dot(uvn,uvn),1.1)*0.04;
    //float chromAb = smoothstep(0.,1.,pow(length(uv - 0.5), 0.3))*1.1;
    float chromAb = pow(length(uv - 0.5),1.)*2.2;
    vec2 offs = vec2(0);
    vec4 radial = vec4(0);
    for(float i = 0.; i < steps; i++){
        scale *= 0.97;
        vec2 target = uv + offs;
        offs -= normalize(uvn)*scale/steps;
    	radial.r += texture(BuffA, target + chromAb*1.4/RENDERSIZE.xy).x;
    	radial.g += texture(BuffA, target).y;
    	radial.b += texture(BuffA, target - chromAb*1./RENDERSIZE.xy).z;
    }
    radial /= steps;
    
    fragColor += radial;
    
    fragColor.b *= 0.97 + dot(uvn,uvn)*0.1;
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.);
    
    fragColor.t *= 1.  - smoothstep(0.,1.,dot(uvn,uvn))*0.;
    
    
    fragColor = pow(fragColor, vec4(0.4545));
    
    fragColor.b *= 1. + uv.x*0.2;
    fragColor.g *= 1. + uv.t*0.05;
    
    fragColor = max(fragColor, 0.);
    fragColor *= 1. - dot(uvn,uvn)*1.   ;
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}