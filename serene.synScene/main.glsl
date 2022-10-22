vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


// Code by Flopine

// Thanks to wsmind, leon, XT95, lsdlive, lamogui and Coyhot for teaching me
// Thanks LJ for giving me the love of shadercoding :3

// Cookie Collective rulz


#define ITER 64.
#define PI 3.141592
#define time TIME
#define BPM 25./2.
#define tempo BPM/60.

vec3 palette (float t, vec3 a, vec3 b, vec3 c, vec3 d)
{return a+b*cos(2.*PI*(c*t+d));}

float random (vec2 st)
{return fract(sin(dot(st.xy, vec2(12.2544, 35.1571)))*5418.548416);}

vec2 moda (vec2 p, float per)
{
    float a = atan(p.y, p.x);
    float l = length(p);
    a = mod(a-per/2., per)-per/2.;
    return vec2(cos(a),sin(a))*l;
}

vec2 mo(vec2 p, vec2 d)
{
    p = abs(p)-d;
    if (p.y > p.x) p.xy = p.yx;
    return p;
}

float stmin(float a, float b, float k, float n)
{
    float st = k/n;
    float u = b-k;
    return min(min(a,b), 0.5 * (u+a+abs(mod(u-a+st,2.*st)-st)));
}

float smin( float a, float b, float k )
{
    float h = max( k-abs(a-b), 0.0 );
    return min( a, b ) - h*h*0.25/k;
}

mat2 rot(float a)
{return mat2(cos(a),sin(a),-sin(a),cos(a));}

vec2 path (float t)
{
    float a = sin(t*0.2 + 1.5), b = sin(t*0.2);
    return vec2(a, a*b);
}

float sphe (vec3 p, float r)
{return length(p)-r;}

float od (vec3 p, float d)
{return dot(p, normalize(sign(p)))-d;}

float cyl (vec2 p, float r)
{return length(p)-r;}

float box( vec3 p, vec3 b )
{
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float sc (vec3 p, float d)
{
    p = abs(p);
    p = max(p.xyz, p.yzx);
    return min(p.x, min(p.y,p.z)) - d;
}

float prim1 (vec3 p)
{
    float c = cyl(p.xz, 0.1);
    float per = 2.;
    p.y += (script_time*0.75+TIME*0.2+(syn_MidTime+syn_MidHighTime)*0.125)*tempo;
    p.y = mod (p.y-per/2., per)-per/2.;
    return smin(sphe(p, 0.3), c, 0.5);
  
}

float prim2 (vec3 p)
{
    float s = sphe(p,1.);
    float o = od(p,.9);
    p.xz *= rot(p.y*0.7);
    p.xz = moda(p.xz, 2.*PI/5.);
    p.x -= 1.;
    return smin(prim1(p), max(-o,s), 0.5);
}

float prim3 (vec3 p)
{
    p.xy = mo(p.xy, vec2(2.));
    return prim2(p);
}

float prim4 (vec3 p)
{
    p.yz *= rot(p.x*0.1);
    p.xy = moda(p.xy, 2.*PI/4.);
    p.x -= 4.;
    return prim3(p);
}

float prim5 (vec3 p)
{
    float per = 8.;
    p.xy *= rot(p.z*0.2);
    p.z = mod(p.z-per/2., per) -per/2.;
    return prim4(p);
} 

float prim6 (vec3 p)
{
    p.z -= script_time+script_bass_time*0.35;
    p.xz *= rot((script_time+script_high_time*0.375)*tempo);
    p.xy *= rot((script_time+script_high_time*0.375)*tempo);
    p *= 1.2;
    return stmin(od(p, 1.), sphe(p,1.), 0.5, 4.);
}



float g = 0.;
float SDF(vec3 p)
{
    float d = smin(prim5(p), prim6(p), 0.5);
    g+=0.1/(0.1+d*d); 
    return d;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy;
    uv -= 0.5;
    uv /= vec2(RENDERSIZE.y / RENDERSIZE.x, 1);

    vec3 ro = vec3(0.001,0.001,-10. + (script_time+script_bass_time*0.35)); vec3 p = ro;
    vec3 rd = normalize(vec3(uv,1.));

    float shad = 0.;
    float dither = random(uv);

    for (float i=0.; i<ITER; i++)
    {
        float d = SDF(p);
        if (d<0.001)
        {
            shad = i/ITER;
            break;
        }
        d *= 0.9 + dither*0.1;
        p += d*rd * 0.7;
    }

    float t = length(ro-p);

    vec3 pal = palette
        (length(uv),
         vec3(0.5),
         vec3(0.5),
         vec3(0.5),
         vec3(0.,0.3,0.7));

    vec3 c = vec3(shad) * pal;
    c = mix(c, vec3(0.,0.,0.2), 1.-exp(-0.001 *t *t));
    c += g* 0.04 * (1.-length(uv));

    fragColor = vec4(c,1.);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}