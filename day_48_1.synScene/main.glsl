

			//******** Common Code Begins ********


//#define time 
    #define T (TIME + 2.5)
    
    #define SPEED 2.

vec3 glow = vec3(0);
    
#define TIME (TIME + 10.)
#define mx (TIME*0.9 + 20.*_mouse.x/RENDERSIZE.x)
#define my (20.*_mouse.y/RENDERSIZE.x)
#define pal(a,b,c,d,e) (a + b*sin(c*d + e))
    
#define dmin(a, b) a.x < b.x ? a : b
#define pmod(a,x) mod(a,x) - x*0.5

vec3 path(float z){
    z *= 0.4;
	return vec3(
    	sin(z + cos(z)*0.6),
    	cos(z + sin(z*0.8)*0.5),
    	0.
    )*1.;
}
float opSmoothUnion(float d1, float d2, float k) {
    float h = clamp(0.5 + 0.5 * (d2 - d1)/k,0.,1.);
    return mix(d2, d1, h) - k*h*(1. - h);  
}
// Hex code from BigWings! He has a tutorial on them.
float HexDist(vec2 p) {
	p = abs(p);
    float c = dot(p, normalize(vec2(1,1.73)));
    c = max(c, p.x);
    return c;
}


vec4 HexCoords(vec2 uv) {
	vec2 r = vec2(1, 1.73);
    vec2 h = r*.5;
    vec2 a = mod(uv, r)-h;
    vec2 b = mod(uv-h, r)-h;
    vec2 gv = dot(a, a) < dot(b,b) ? a : b;
    float x = atan(gv.x, gv.y);
    float y = .5-HexDist(gv);
    vec2 id = uv-gv;
    return vec4(x, y, id.x,id.y);
}



float modu;
float moduB;
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
vec2 map(vec3 p){
	vec2 d = vec2(10e6);
	p -= path(p.z);
    vec2 n = normalize(p.xy);
    #define modDist 1.
    #define tunnW 0.7
    #define pipeW 0.04
    
    vec3 g = p;
    p.z = pmod(p.z, modDist);
    
    vec3 o = p;
    //p.xy *= rot(0.4 + p.z*(0.1 + sin(TIME*0.1) )+ TIME*0.3);
    vec2 pC = vec2(atan(p.y,p.x), length(p.xy));
    
    vec3 q = vec3(pC, p.z);

    
    // hex
    //vec4 hc = HexCoords(vec2(pC.x, z)*20.);
    vec4 hc = HexCoords(vec2(pC.x, p.z*2.)*1.);
    
    float dHex = hc.y - 0.1 + sin(o.z)*0.06;
    dHex = max(dHex, -length(o.xy*1.) + tunnW*0.6);
    //dHex = max(dHex, length(o.xy*1.) - tunnW*0.97);
    //d = dmin(d, vec2(dHex, 2.));    
    
    q.y -= tunnW*1.;
    
    q.x = pmod(q.x*4.,1.);
    q.z = pmod(q.z,0.2);
        
    // pipes
    float dPipes = min(length(q.xy) - pipeW,length(q.zy) - pipeW );
    
    // dots
    #define pi acos(-1.)
    q = vec3(pC, g.z);
    q.z = pmod(q.z, 0.14);
    
    q.x = pmod(q.x*3.*pi + 0.5, 1.);
    
    q.y -= tunnW*0.9;
    q.y *= 15.;
    q.z *= 8.;
    float dDots = length(q) - 0.4;
    //dDots *= 0.7;
    
    
    // mod
    float mm = sin(TIME*0.5 + g.z*0.5 + p.z);
    modu = (mm/sqrt(0.02 + mm*mm ))*0.5 + 0.5;
    float mmB = sin(TIME*1.25 + g.z*0.25 + p.z*0.8 + p.y);
    moduB = (mmB/sqrt(0.01 + mmB*mmB ))*0.5 + 0.5;
    //moduB = 1. - moduB;
    moduB *= 0.16;
    //modu = 1.;
    
    // thing
    dPipes = mix(dPipes,dDots,moduB*0.7  + 0.);
    
    float dThing = mix(
    	dHex, dPipes, modu*1. + 0.3
    );
    d = dmin(d, vec2(dThing, 2.));
    
    // tunnel
    float dTunn = -length(o.xy*1.) + tunnW;
    dTunn = max(length(o.xy*1.) - tunnW - 0.02, dTunn);
    dTunn = max(dTunn, -dThing);
    d = dmin(d, vec2(dTunn, 10.));
    
    
    
    d.x *= 0.14;
	return d;
}

vec2 march(vec3 ro,vec3 rd,inout vec3 p,inout float t,inout bool hit){
	hit = false;
    p = ro;
    t = 0.;
    vec2 d;
    for(int i = 0; i < 200; i++){
    	d = map(p);
        if (d.y < 10.){
            glow += mix(
                mix(
                    exp(-d.x*5.)*pal(1.3,0.7,vec3(1.8+modu*0.5,0.4,0.8), 3.9 +modu*0.2 + sin(p.z)*0.5,2. + t*0.1)*2.,
                    exp(-d.x*0.01)*pal(2.8,1.4,vec3(5.8,2.4+modu*0.2,0.8+ sin(p.y + TIME)*.09), 6.9 +modu*0.2 + sin(p.z)*.5,0.5 + t*0.5)*2.,
                    moduB
                ),exp(-d.x*20.)*pal(0.8,0.7,vec3(1.8+modu*0.5,0.4,0.8), 3.9 +modu*0.2 + sin(p.z)*0.5,1.9 - t*0.1)*2.
                , 1. - modu
            );
        	
        
        } else {
        	glow += exp(-d.x*2.)*pal(0.2,1.9,vec3(1.8,0.4,0.8), 5.6 ,2.- t*0.1)*0.4;
        }
        
        if(d.x < 0.001){
        	hit = true;
            break;
        }
        t += d.x;;
    	p = ro + rd*t;
    }
	return d;
}

vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
	vec3 dir = normalize(lookAt - ro);
	vec3 right = cross(vec3(0,1,0), dir);
	vec3 up = cross(dir, right);
	return normalize(dir + right*uv.x + up*uv.y);
}
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    uv *= 1. - dot(uv,uv)*0.2;
    
    vec3 col = vec3(0);

    vec3 ro = vec3(0);
    ro.z += mx;
    ro += path(ro.z);
    vec3 lookAt = ro + vec3(0,0,4);
    lookAt += path(lookAt.z);
    
    vec3 rd = getRd(ro, lookAt, uv);
    float t; bool hit;
    vec3 p;
    
    vec2 d = march(ro, rd, p, t, hit);
    
    if(hit){
        //vec3 N = getNormal
		col += 0.04;
    }
        
    col *= 0.4;
    col = pow(col,vec3(0.45));
    
    col += glow*0.011;
    uv *= 0.8;
    col *= 1. - dot(uv,uv);
    //col *= 1. - t*0.2;    
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}