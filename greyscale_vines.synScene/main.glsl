//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


// Code by Flopine
// Thanks to wsmind, leon, XT95, lsdlive, lamogui, Coyhot, Alkama and YX for teaching me
// Thanks LJ for giving me the love of shadercoding :3

// Thanks to the Cookie Collective, which build a cozy and safe environment for me 
// and other to sprout :)  https://twitter.com/CookieDemoparty


// fix image for fanzine
//#define dt 19.6
/*
float smoothTime = (smooth_basstime * 0.75 + TIME * 0.1 + syn_Time * 0.25);
float smoothTimeB = (smooth_hightime * 0.75 + TIME * 0.125 + syn_Time * 0.25);
float smoothTimeC = (smooth_midtime * 0.75 + TIME * 0.125 + syn_Time * 0.125);
*/
#define dt (smoothTime*1.5)
#define PI 3.141592
#define ITER 96.

float hash11 (float x)
{return fract(sin(x)*1245.5);}

float hash21 (vec2 x)
{return fract(sin(dot(x,vec2(12.45,23.5)))*1245.4);}

mat2 rot (float a)
{return mat2(cos(a),sin(a),-sin(a),cos(a));}

void moda (inout vec2 p, float rep)
{
    float per = 2.*PI/rep;
    float a = atan(p.y,p.x);
    float l = length(p);
    a = mod(a,per)-per*0.5;
    p = vec2(cos(a),sin(a))*l;
}

float smin( float a, float b, float k )
{
    float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
} 

float cyl (vec2 p, float r)
{return length(p)-r;}

float sdHexPrism( vec3 p, vec2 h )
{
    const vec3 k = vec3(-0.8660254, 0.5, 0.57735);
    p = abs(p);
    p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
    vec2 d = vec2(
       length(p.xy-vec2(clamp(p.x,-k.z*h.x,k.z*h.x), h.x))*sign(p.y-h.x),
       p.z-h.y );
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float twist (vec3 p, float w)
{
    p.xz *= rot(p.y*0.8+(moreTwisties*(smoothTimeC)));
    moda(p.xz,7.);
    p.x -= w*2.;
    return cyl(p.xz, w);
}

float twists (vec3 p, float w)
{
    float per = 7.;
    float id = floor(p.y/per);
	if (mod(id,2.) == 0.) p.xz *= rot(PI/6.);
    p.y = mod(p.y, per)-per*0.5;
    
    p.y += sin(p.x+smoothTimeC)*(0.5+distortion);
    p.x += cos(p.z+smoothTimeC*0.43)*(0.4+distortion);
    p.xy *= rot(PI/2.);
    moda(p.yz, 6.);
    return twist(p, w);
}

float room (vec3 p)
{return -sdHexPrism(p.xzy, vec2((15.+wallsOut),1e10));}

float SDF (vec3 p)
{
    return min(smin(twist(p,0.5),twists(p,0.2),0.8),max(-twists(p, 1.),room(p)));
}

vec3 get_cam (vec3 ro, vec3 target, vec2 uv, float fov)
{
    vec3 forward = normalize(target - ro);
    vec3 left = normalize(cross(vec3(0.,1.,0.), forward));
    vec3 up = normalize(cross (forward, left));
    return normalize(forward*fov+ left*uv.x + up*uv.y);
}

float raymarch (vec3 ro, vec3 rd, vec2 uv)
{
    vec3 p = ro;
    float dither = hash21(uv);
    for (float i=0.; i<ITER; i++)
    {
        float d = SDF(p);
        if (d<0.0001)
        {
            float shad = i/ITER;
            return 1.-shad;
        }
        d *= 0.7+dither*0.1;
        p += d*rd;
    }
    return 0.;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (2.*fragCoord-RENDERSIZE.xy)/RENDERSIZE.y;
 
    vec3 ro = vec3(3.,-4. +dt,-10.-Height),
        p = ro,
        tar = vec3(0., dt, 0.),
        rd = get_cam(ro, tar, uv, 1.),
        col = vec3(0.);
        rd.yz = _rotate(rd.yz, lookXY.y*PI);
        rd.xy = _rotate(rd.xy, -1.0*lookXY.x*PI);
        rd.xz = _rotate(rd.xz, lookZ*PI);
    
    // anaglyph technique from leon/ponk 
    float red  = raymarch(ro-vec3(.01,0.,0.),rd,uv);
    float cyan = raymarch(ro+vec3(.01,0.,0.),rd,uv);
    
    col = vec3(red,vec2(cyan));

    fragColor = vec4(pow(col,vec3(2.2)),1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}