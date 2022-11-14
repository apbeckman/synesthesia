

			//******** Common Code Begins ********

// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0
// Unported License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ 
// or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
// =========================================================================================================

#define GLOW_SAMPLES 80
#define GLOW_DISTANCE 0.2
#define GLOW_POW 1.
#define GLOW_OPACITY 1.86

#define sat(a) clamp(a, 0., 1.)
#define PI 3.14159265
#define TAU (PI*2.0)

mat2 r2d(float a) { float c = cos(a), s = sin(a); return mat2(c, -s, s, c); }
float hash11(float seed)
{
    return mod(sin(seed*123.456789)*123.456,1.);
}

vec3 getCam(vec3 rd, vec2 uv)
{
    float fov = 3.;
    vec3 r = normalize(cross(rd, vec3(0.,1.,0.)));
    vec3 u = normalize(cross(rd, r));
    return normalize(rd+fov*(r*uv.x+u*uv.y));
}

vec2 _min(vec2 a, vec2 b)
{
    if (a.x < b.x)
        return a;
    return b;
}

float _cucube(vec3 p, vec3 s, vec3 th)
{
    vec3 l = abs(p)-s;
    float cube = max(max(l.x, l.y), l.z);
    l = abs(l)-th;
    float x = max(l.y, l.z);
    float y = max(l.x, l.z);
    float z = max(l.x, l.y);
    
    return max(min(min(x, y), z), cube);
}

float _cube(vec3 p, vec3 s)
{
    vec3 l = abs(p)-s;
    return max(l.x, max(l.y, l.z));
}

			//******** BuffA Code Begins ********

// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0
// Unported License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ 
// or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
// =========================================================================================================

float _seed;
float rand()
{
    _seed++;
    return hash11(_seed);
}
vec2 map(vec3 p)
{
    vec3 p2 = p;
    vec2 acc = vec2(10000.,-1.);
    vec2 rep = vec2(.25);
    vec2 id = floor((p.xz+rep*.5)/rep);

    p.xz = mod(p.xz+rep*.5,rep)-rep*.5;
    float sz = .1;
    float h = .5*texture(image17, id*.01+TIME*.0008).x*length(id*.15);
    vec3 op = p;
    p.y = abs(p.y);
    p.y-= h;
    float bar = length(p)-sz;
    bar = min(bar, max(length(p.xz)-sz, abs(op.y)-h));
    bar = max(bar, op.y);
    acc = _min(acc, vec2(bar, 0.));
    
    vec3 rep2 = vec3(1.);
    p2.y += TIME;
    p2.xz *= r2d(sin(p2.y*3.-TIME)*.05);
    vec3 id2 = floor((p2+rep2*.5)/rep2);
    p2 = mod(p2+rep2*.5,rep2)-rep2*.5;
    float bubble = length(p2)+abs(sin(id2.y+id2.x+id2.z*10.));
    acc = _min(acc, vec2(bubble, 3.));
    
    
    return acc;
}

vec3 getNorm(vec3 p, float d)
{
    vec2 e = vec2(0.01, 0.);
    return normalize(vec3(d)-vec3(map(p-e.xyy).x, map(p-e.yxy).x, map(p-e.yyx).x));
}
vec3 accCol;
vec3 trace(vec3 ro, vec3 rd, int steps)
{
    accCol = vec3(0.);
    vec3 p = ro;
    for (int i = 0; i < steps && distance(p, ro) < 10.; ++i)
    {
        vec2 res = map(p);
        if (res.x < 0.01)
            return vec3(res.x, distance(p, ro), res.y);
        p+=rd*res.x*.35;
        vec3 rgb = mix(vec3(.1,.3,.5),vec3(0.078,0.522,0.471), sat(sin(2.*distance(p, ro))));
        accCol += rgb*0.01*sat(.3+sat(sin(p.x*2.)+sin(p.z)+sin(length(p.xz)+TIME)));//*(1.-sat(res.x/0.5))*.01;
    }
    return vec3(-1.);
}

vec3 rdr(vec2 uv)
{
    vec3 col = vec3(0.);
    float d = 2.;
    float a = TIME*.1;
    uv *= r2d(sin(TIME*.2)*.2);
    vec2 offr = (vec2(rand(), rand())-.5)*.05;
    vec3 ro = vec3(sin(a)*d+offr.x,-.8+offr.y+.2*sin(TIME*.1),cos(a)*d);
    vec3 ta = vec3(0.,0.,0.);
    vec3 rd = normalize(ta-ro);
    
    rd = getCam(rd, uv);
    vec3 res = trace(ro, rd, 128);
    float depth = 100.;
    if (res.y > 0.)
    {
        depth = res.y;
        vec3 p = ro+rd*res.y;
        vec3 n = getNorm(p, res.x);
        col = n*.5+.5;
        vec2 rep = vec2(.25);
        vec2 id = floor((p.xz+rep*.5)/rep);
        float h = texture(image17, id*.01+TIME*.0008).x;
        vec3 rgb = mix(vec3(0.102,1.000,0.698), vec3(0.486,0.435,0.094), pow(h,3.));
        col = rgb*(1.-pow(sat(-dot(rd, n)), .25));
        col += rgb * pow(sat(-p.y*1.5),2.);
        col *= pow(h*1.5,2.);
        col = mix(col*.35, 1.5*col.yxz, sat(sin(length(id)*.7-TIME*2.)));
    }
    col = mix(col, accCol, 1.-sat(exp(-depth*0.4)));
    return col;
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-.5*RENDERSIZE.xy)/RENDERSIZE.xx;
    _seed = TIME+texture(image17, uv).x;

    vec3 col = rdr(uv);
    
    vec2 off = vec2(1., -1.)/(RENDERSIZE.x*1.5);

    if (true) // Not so cheap antialiasing
    {
        //col = vec3(1.,0.,0.);
        vec3 acc = col;
        acc += rdr(uv+off.xx);
        acc += rdr(uv+off.xy);
        acc += rdr(uv+off.yy);
        acc += rdr(uv+off.yx);
        col = acc/5.;
        
    }
    col *= 1.9/(col+1.);
    col = sat(col);

    col = mix(col, texture(BuffA, fragCoord/RENDERSIZE.xy).xyz, .9);
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0
// Unported License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ 
// or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
// =========================================================================================================

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy;

    const int steps = GLOW_SAMPLES;
    vec3 col = vec3(0.);
    
    for (int i = 0; i< steps; ++i)
    {
        float f = float(i)/float(steps);
        f = (f -.5)*2.;
        float factor = GLOW_DISTANCE;
        col += texture(BuffA, uv+vec2(f*factor,0.)).xyz/float(steps);
    }
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


// This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0
// Unported License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ 
// or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
// =========================================================================================================

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy;
    const int steps = GLOW_SAMPLES;
    vec3 col = vec3(0.);
    
    for (int i = 0; i< steps; ++i)
    {
        float f = float(i)/float(steps);
        f = (f -.5)*2.;
        float factor = GLOW_DISTANCE;
        col += texture(BuffB, uv+vec2(0.,f*factor)).xyz/float(steps);
    }
    
    vec3 rgb = texture(BuffA, uv).xyz+GLOW_OPACITY*pow(col, vec3(GLOW_POW));
    rgb = pow(rgb*1.2, vec3(2.2));
    vec2 cuv = (fragCoord-.5*RENDERSIZE.xy)/RENDERSIZE.xx;
    
    fragColor = vec4(rgb,1.0);
	return fragColor; 
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