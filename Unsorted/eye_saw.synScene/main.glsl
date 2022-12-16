

			//******** Common Code Begins ********

// Created by SHAU - 2020
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
//-----------------------------------------------------

#define PI 3.141592
#define R RENDERSIZE.xy
#define T mod(TIME, 640.0)
#define CAM_POS vec2(0.5, 0.5)
#define LOOK_AT vec2(1.5, 0.5)
#define ANIM vec2(2.5, 0.5)

#define UI0 1597334673U
#define UI1 3812015801U
#define UI2 uvec2(UI0, UI1)
#define UI3 uvec3(UI0, UI1, 2798796415U)
#define UIF (1.0 / float(0xffffffffU))

//Dave Hoskins - Hash without sin
//https://www.shadertoy.com/view/XdGfRR
float hash11(float p)
{
	uvec2 n = uint(int(p)) * UI2;
	uint q = (n.x ^ n.y) * UI0;
	return float(q) * UIF;
}

float hash12(vec2 p) 
{
	uvec2 q = uvec2(ivec2(p)) * UI2;
	uint n = (q.x ^ q.y) * UI0;
	return float(n) * UIF;
}

vec3 hash32(vec2 q)
{
	uvec3 n = uvec3(ivec3(q.xyx)) * UI3;
	n = (n.x ^ n.y ^ n.z) * UI3;
	return vec3(n) * UIF;
}

//Fabrice - compact rotation
mat2 rot(float x) {return mat2(cos(x), sin(x), -sin(x), cos(x));}


//standard path deform - check out any tunnel by Shane or Aiekick for some good examples
vec3 path(float t) 
{
    float a = sin(t * PI / 24. + 2.7);
    float b = cos(t * PI / 24.);
    return vec3(a * 3., b * a, 0.);
} 

			//******** BuffA Code Begins ********

// Created by SHAU - 2020
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
//-----------------------------------------------------

//CAMERA

vec4 renderPassA() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;

    float AT = mod(T, 40.0),
          scene = 1.0,
          fade = min(AT*0.5, 1.0) - clamp(AT - 19.0, 0.0, 1.0) +
                 clamp(AT - 20.0, 0.0, 1.0) - clamp(AT*0.5 - 19.0, 0.0, 1.0);
    
    vec3 lookAt = vec3(2.0, 0.0, T*3.0 + sin(T*0.2)*20.0),
         camPos = vec3(0.0, 0.0, T*3.0 - 5.0);
    lookAt += path(lookAt.z);
    camPos += path(camPos.z);
    
    if (AT>20.0)
    {
        scene = 0.0; 
        lookAt = vec3(0.0, 0.0, 0.0),
        camPos = vec3(0.0, sin(T*0.23)*1.8, -5.0 - sin(T*0.1));
        camPos.xz *= rot(T*0.2);
    }
    
    if (U==ANIM)
    {
        C = vec4(scene, fade, 0.0, 0.0);
    }
    else if (U==CAM_POS)
    {
        C = vec4(camPos, 1.0);    
    }
    else if (U==LOOK_AT)
    {
        C = vec4(lookAt, 1.0);    
    }
    
	return C; 
 } 


			//******** BuffB Code Begins ********

// Created by SHAU - 2020
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
//-----------------------------------------------------

//RENDERING

#define EPS 0.005
#define FAR 100.0
#define ZERO (min(int(FRAMECOUNT),0))

#define S(a,b,v) smoothstep(a,b,v)
#define TORUS vec2(5.0, 3.0)
#define TS vec2(0.0625)

//SHADERTOY - TEXT
const vec2[9] text = vec2[](vec2(9.0, 10.0),  //Y
                            vec2(15.0, 11.0), //O
                            vec2(4.0, 10.0),  //T
                            vec2(2.0, 10.0),  //R
                            vec2(5.0, 11.0),  //E
                            vec2(4.0, 11.0),  //D
                            vec2(1.0, 11.0),  //A
                            vec2(8.0, 11.0),  //H
                            vec2(3.0, 10.0)); //S

//extract letter from texture
vec4 letter(vec2 uv, vec2 id, float s)
{
    vec4 l = vec4(0.0);
    uv *= s; //scale
    uv += TS*id - vec2(0.5*s) + TS*0.5; //center
    
    if (uv.x>(id.x*TS.x) && uv.x<((id.x + 1.0)*TS.x) && 
        uv.y>(id.y*TS.y) && uv.y<((id.y + 1.0)*TS.y))
    {
        l = texture(image49, uv, -1.0);
    }
    return l;
}

//noise IQ - Shane
float n3D(vec3 p) 
{    
	const vec3 s = vec3(7, 157, 113);
	vec3 ip = floor(p); 
    p -= ip; 
    vec4 h = vec4(0., s.yz, s.y + s.z) + dot(ip, s);
    p = p * p * (3. - 2. * p);
    h = mix(fract(sin(h) * 43758.5453), fract(sin(h + s.x) * 43758.5453), p.x);
    h.xy = mix(h.xz, h.yw, p.y);
    return mix(h.x, h.y, p.z);
}

//IQ - Usefull little functions
//https://iquilezles.org/articles/functions
float expImpulse(float x, float k)
{
    float h = k*x;
    return h*exp(1.0-h);
}

//Shane - Perspex Web Lattice - one of my favourite shaders
//https://www.shadertoy.com/view/Mld3Rn
//Standard hue rotation formula... compacted down a bit.
vec3 rotHue(vec3 p, float a)
{
    vec2 cs = sin(vec2(1.570796, 0) + a);

    mat3 hr = mat3(0.299,  0.587,  0.114,  0.299,  0.587,  0.114,  0.299,  0.587,  0.114) +
        	  mat3(0.701, -0.587, -0.114, -0.299,  0.413, -0.114, -0.300, -0.588,  0.886) * cs.x +
        	  mat3(0.168,  0.330, -0.497, -0.328,  0.035,  0.292,  1.250, -1.050, -0.203) * cs.y;
							 
    return clamp(p*hr, 0., 1.);
}

//IQ - SDF primitives
//https://iquilezles.org/articles/distfunctions
float sdTorus(vec3 p, vec2 t) 
{
    vec2 q = vec2(length(p.xz) - t.x, p.y);
    return length(q) - t.y;
}

vec3 tunnel(vec3 p)
{
    p -= path(p.z); 
	return vec3(3.0 - length(p.xy), p.xy);    
}

void ring(inout vec3 pc, vec3 p, vec3 col)
{
    float z = fract(p.z*0.04 + -T*0.3),
          atn1 = (1.0 / (1.0 + z*z*100000.0)) + (1.0 / (1.0 + (1.0 - z)*(1.0 - z)*100000.0)),
          atn2 = (1.0 / (1.0 + z*z*6000.0)) + (1.0 / (1.0 + (1.0 - z)*(1.0 - z)*6000.0));
    pc += pc*atn2*16.0;
    pc += col * length(pc)*atn2*0.2;
    pc += mix(col, vec3(1), 0.4)*atn1 * n3D(p.zxy*2.3+T*0.2);
}

vec4 renderTunnel(vec3 ro, vec3 rd, vec3 la, vec3 col)
{
    vec3 pc = vec3(0);
    
    float t = 0.0, dof = 0.0;
    vec2 id = vec2(0.0);;
    for (int i=ZERO; i<100; i++) {
        vec3 ns = tunnel(ro + rd*t);
        if (abs(ns.x)<EPS) {
            id = ns.yz;
            break;
        }
        t += ns.x;
        
        if (t>FAR)
        {
            t = -1.0;
            break;
        }
    }
    
    if (t>0.0)
    {
        vec3 p = ro + rd*t;
        dof = 0.0;
		//partition surface
        float a = (atan(id.x, id.y)/6.2831853) + 0.5; //0->1

        vec2 c = vec2(fract(a*60.0), //cell positiom
                      fract(p.z*10.0)),
        cid = vec2(floor(a*60.0), //cell id
                   floor(p.z*10.0));
        vec3 h3 = hash32(cid+T*0.4);

        //RANDOM CHARACTERS
        if (h3.z>0.3)
        {
            pc = col * letter(vec2(1.0 - c.x, c.y), floor(h3.xy*16.0), 0.1).x;
            pc *= hash11(cid.y) * n3D(p*0.5) * n3D(p*5.0+T);
            pc = pc*0.3 + pc*S(0.5, 1.0, n3D(p.zxy*-0.4));
        }
        //SHADERTOY
        vec2 cz = vec2(c.x,
                       fract(fract(p.z)*9.0));
        vec2 czid = vec2(cid.x,
                       floor(fract(p.z)*9.0));
        if (mod(p.z, 10.0)>9.0 && hash12(czid)>0.3)
        {
            pc += col * 2.0 * letter(vec2(1.0 - cz.x, cz.y), text[int(czid.y)], 0.1).x;
        }
        
        //lines
        float l = length(c.x-0.5); //center line
        l = 0.3 / (1.0 + l*l*800.);    
        l *= hash11(cid.x+T*-0.3) * n3D(p*0.5) * n3D(p*5.0+T); //add noise
        l = l*0.3 + l*S(0.5, 1.0, n3D(p.yxz*0.5)); //more noise   
        float ls = l*S(1.0, 0.8, fract(p.z))*S(0.0, 0.2, fract(p.z)); //gaps in lines
        l = mix(ls, l, step(0.8, hash11(cid.x+T*-0.3))); //more noise        
        pc += col * l * S(0.0, 1.0, n3D((p-32.3)*0.7 + T*0.3)) * 2.0;
        
        ring(pc, p, col);
        p.z += 1.4;
        ring(pc, p, col);        
    }

    return vec4(pc, dof);
}                      

vec4 renderTorus(vec3 ro, vec3 rd, vec3 la, vec3 col)
{
    vec3 pc = vec3(0);
    
    float t = 0.0, dof = 0.0;
    for (int i=ZERO; i<100; i++) {
        float ns = abs(sdTorus(ro + rd*t, TORUS));
        if (abs(ns)<EPS) break;
        t += ns;
        if (t>FAR)
        {
            t = -1.0;
            break;
        }
    } 

    if (t>0.0)
    {
        vec3 p = ro + rd*t;
        dof = length(la-p) * 0.03;
        //partition surface
        float aX = (atan(p.x, p.z)/6.2831853) + 0.5, //0->1
            aY = (atan(length(p.xz)-TORUS.x, p.y)/6.2831853) + 0.5, //0->1 from BigWings
            aYR = aY*40.0+T*0.3, //split and rotation
            aYID = mod(floor(aYR), 40.0), //id
            h1 = hash11(aYID)*2.0 - 1.2, //biased hash
            h1r = h1*-1.0, //reversed
            aXR = aX*60.0 + T*(0.5*sign(h1)+h1) + //split and rotation 
            expImpulse(mod(T*(0.6+abs(h1)), 2.0), 2.0) * sign(h1);

        vec2 c = vec2(fract(aXR), fract(aYR)), //cell co-ordinate
            cid = vec2(mod(floor(aXR), 60.0), aYID);  //cell id
        vec3 h3 = hash32(cid); //hash on cell id

        //RANDOM CHARACTERS
        if (h3.z>0.3)
        {
            pc = col * letter(vec2(1.0 - c.x, c.y), floor(h3.xy*16.0), 0.1).x;
            pc *= hash11(cid.y) * n3D(p*0.5) * n3D(p*5.0+T);
            //pc = pc*0.3 + pc*S(0.5, 1.0, n3D(p.zxy*-0.4));
        }

        //SHADERTOY
        if (cid.y==9.0 && cid.x<9.0)
        {
            pc = 2.0 * col * letter(vec2(1.0 - c.x, c.y), text[int(cid.x)], 0.1).x;                
        }
        else if (cid.y==15.0 && cid.x>22.0 && cid.x<32.0)
        {
            pc = 2.0 * col * letter(vec2(1.0 - c.x, c.y), text[int(cid.x) - 23], 0.1).x;                
        }
        else if (cid.y==30.0 && cid.x>39.0 && cid.x<49.0)
        {
            pc = 2.0 * col * letter(vec2(1.0 - c.x, c.y), text[int(cid.x) - 40], 0.1).x;                
        }

        //GLOW LINES
        //with wrapping
        vec2 scl = vec2(30.0, 1.0), //scale
             bp = vec2(fract(T*(0.1 + abs(h1r)*0.1))*scl.x * sign(h1r), 0.5);  //center

        //disk
        float bd = h1r>0.0 ? -scl.x : 2.0*scl.x;
        float b = length(bp - vec2(aX, c.y)*scl);
        b = min(b, length(vec2(bp.x + scl.x, bp.y) - vec2(aX, c.y)*scl));
        b = min(b, length(vec2(bp.x + bd, bp.y) - vec2(aX, c.y)*scl));

        //line
        float l = length(c.y - 0.5);
        float lb1 = length(bp.x/scl.x - aX);
        float lb2 = length((bp.x + scl.x)/scl.x - aX);

        pc += col / (1.0 + b*b*1000.0);
        pc += vec3(1) / (1.0 + b*b*6000.0);            

        float clip1 = h1r>0.0 ? step(aX, bp.x/scl.x) : step(bp.x/scl.x, aX);
        float clip2 = h1r>0.0 ? step(aX, (bp.x + scl.x)/scl.x) : step((bp.x + scl.x)/scl.x, aX);
        pc += col * S(0.02, 0.0, l) / (1.0 + lb1*lb1*4000.0) * clip1;
        pc += col * S(0.02, 0.0, l) / (1.0 + lb2*lb2*4000.0) * clip2;
        pc += vec3(1) * S(0.02, 0.0, l) / (1.0 + lb1*lb1*12000.0) * clip1;
        pc += vec3(1) * S(0.02, 0.0, l) / (1.0 + lb2*lb2*12000.0) * clip2;
    }

    return vec4(pc, dof);
}

vec3 camera(vec2 U, vec3 ro, vec3 la, float fl) 
{
    vec2 uv = (U - R*.5) / R.y;
    vec3 fwd = normalize(la-ro),
         rgt = normalize(vec3(fwd.z, 0., -fwd.x));
    return normalize(fwd + fl*uv.x*rgt + fl*uv.y*cross(fwd, rgt));
}

vec4 renderPassB() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;
     
    vec4 cam = texture(BuffA, CAM_POS/R),
         lookAt = texture(BuffA, LOOK_AT/R),
         anim = texture(BuffA, ANIM/R);
    vec3 col = rotHue(vec3(0,1,0), T*0.2) * 2.0;
    
    vec3 rd = camera(U, cam.xyz, lookAt.xyz, 1.8);

    //debug text lookup
    //vec2 uv = U/R;
    //vec2 lid = vec2(3.0, 10.0); //S
    //pc += letter(uv, lid, 0.06);
    
    if (anim.x==0.0)
    {
        C = renderTorus(cam.xyz, rd, lookAt.xyz, col);
    }
    else if (anim.x==1.0)
    {
        C = renderTunnel(cam.xyz, rd, lookAt.xyz, col);
    }
	return C; 
 } 


// Created by SHAU - 2020
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
//-----------------------------------------------------

/*
   Inspired by some of the animations on Isaac Arthur's YouTube channel 
*/

const float GA =2.399; 

// simplified version of Dave Hoskins blur from Virgill
vec3 dof(sampler2D tex, vec2 uv, float rad) {
	vec3 acc = vec3(0);
    vec2 pixel = vec2(.002*R.y/R.x, .002), angle = vec2(0, rad);;
    rad = 1.;
	for (int j = 0; j < 80; j++) {  
        rad += 1. / rad;
	    angle *= rot(GA);
        vec4 col=texture(tex,uv+pixel*(rad-1.)*angle);
		acc+=col.xyz;
	}
	return acc/80.;
}

vec4 renderMainImage() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;
     
    vec4 anim = texture(BuffA, ANIM/R);
    vec3 pc = vec4(dof(BuffB, U/R, texture(BuffB, U/R).w), 1.).xyz;
    
    pc = pow(pc*2.0, vec3(1.0/2.3)); //gamma correction
    pc *= anim.y; //fade

    C = vec4(pc, 1.0);
	return C; 
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