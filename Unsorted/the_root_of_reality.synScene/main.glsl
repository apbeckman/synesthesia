

//                  = The Root of Reality =         
//               by Maximilian Knape ·∑>| 2023            
// -----------------------------------------------------------
// This work is licensed under a Creative Commons Attribution-
//        NonCommercial-ShareAlike 3.0 Unported License

#define MAX_STEPS 500
#define MAX_DIST 350.
#define MIN_DIST 20.
#define STEP_FAC 1.1

#define SURF_DIST .03
#define SURF_MUL 500.
#define SURF_EXP 12.

#define PP_CONT 1.4
#define PP_VIGN 2.0

#define TIME TIME*1.3
#define PI 3.14159265358979
#define TAU 6.28318530717958
#define S(x,y,t) smoothstep(x,y,t)

float smin( float a, float b, float k ) //iq
{
    float h = max(k-abs(a-b),0.);
    return min(a, b) - h*h*.25/k;
}

vec2 Map(vec3 p) 
{    
    float d = MAX_DIST, col = 0.85;
    vec3 pos = p;

    if (length(p - 100. + vec3(pow(sin(TIME*3./4.5), 3.)*20.)) < 50.)
        for(float i = 0.; i < 5.; i++)
        {
            float t = (TIME + 2.*pow(sin(TIME/23.),3.))*3. - i*.8;
            p = pos - 100. + vec3(pow(sin(t/4.5), 3.)*20. + i);
            float s = length(p + sin(t / vec3(4.1,6.3,2.2)) * vec3(11,19,8));
            s -= (4. / (i*i*.1 + 1.)); 
            d = smin(s *.6, d, 3. - i*.5);
        }
    
    pos = pos + vec3(TIME*10. + pow(sin(TIME/15.), 3.));
    p = pos / 12.;
    
    float n1 = length(dot(sin(p/2.5), cos(p*1.8))) + .2 + 1./(TIME/2.+.01);
    col = mix(col, 1.0, step(n1, d));
    d = min(n1, d);
    
    p =  pos / 72. + pow(sin(TIME/124.), 3.)*10.;
    float n2 = length(dot(sin(p/14.5), cos(p))) - (.3 + sin(TIME/13.)*.1);
    col = mix(col, 2.0, step(n2, d));
    d = min(n1, d);
   
    
    return vec2(d, col);
}

vec3 Palette(int index)
{
    switch (index)
    {
        case 0: return vec3(1., 1., 1.);
        case 1: return vec3(.85, .15, .3);
        case 2: return vec3(.3, .65, 1.);
        case 3: return vec3(1., .8, .3);
    }
    return vec3(0.);
}

vec3 RTM(vec3 ro, vec3 rd) 
{
    int steps;
    float sum = 0.;
	float s = 1.;
    float d = MIN_DIST;
    const float a = 1. / float(MAX_STEPS); 
    vec3 p = vec3(0), col = vec3(1);
    
    for(int i = 0; i < MAX_STEPS; i++) 
    {    
        float sd = (SURF_DIST * (pow(d/MAX_DIST, SURF_EXP)*SURF_MUL + 1.));
        if (s < sd || d > MAX_DIST) break;
        
        steps = i;
        p = ro + rd*d;
        
        vec2 map = Map(p);
        col = mix(col, Palette(int(floor(map.y))), .02 * (1.-sum));
        
        s = max(abs(map.x), 2. * sd);
        d += s * STEP_FAC * (1.1 - fract(map.y));
        
        sum += a;
    }
    
    col *= sum - pow(1. - (length(ro - p) / MAX_DIST), SURF_EXP);
    
    return col;
}

mat2 Rot(in float a)
{
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

vec3 R(in vec2 uv, in vec3 p, in vec3 l, in float z)
{
    vec3 f = normalize(l-p),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f,r),
        c = p+f*z,
        i = c + uv.x*r + uv.y*u,
        d = normalize(i-p);
    return d;
}

vec4 PP(vec3 col, vec2 uv)
{
    col = mix(col, S(vec3(0.), vec3(1.), col), PP_CONT);    
    col *= S(PP_VIGN,-PP_VIGN/5., dot(uv,uv)); 
    col = pow(col, vec3(.9)); //vec3(.4545)
    
    return vec4(col, 1.);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-.5 * RENDERSIZE.xy) / RENDERSIZE.y;
	vec2 m = _mouse.xy / RENDERSIZE.xy;
    if (length(m) <= 0.) m = vec2(.5);

    vec3 ro = vec3(-1);
    ro.yz *= Rot(-PI/2. + PI*.5);
    ro.xz *= Rot(-m.x * PI*2. - PI);
    vec3 rd = R(uv, ro, vec3(0), .8);
    
    vec3 col = RTM(ro, rd);
    
    float dr = dot(rd, vec3(.57735));
    col += pow(max(dr, 0.), 100.+sin(TIME*.9)*20.)*Palette(3)*.65;
    col *= pow(dr*.5+.5, .5+sin(TIME/5.)*.2)*S(0., 1., TIME);
    col -= S(1.-1e-3, 1., dr)*.2*sin(TIME/5.);
    col += S(1.-5e-5, 1., dr)*.3;
    
    col += .2 * S(0., 10., TIME) * ((S(-1., 1., dr) - S(0.99, 1., dr)) *
           S(0.3, 3., length(sin(10.*normalize(vec2(-abs(rd.z-rd.x), abs(rd.y-rd.x)) * Rot(TIME/31.))))));  

    fragColor = PP(col, uv);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}