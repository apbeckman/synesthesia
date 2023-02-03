

/*
* License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
* Created by bal-khan
*/

vec2	march(vec3 pos, vec3 dir);
vec3	camera(vec2 uv);
void	rotate(inout vec2 v, float angle);
vec3	calcNormal( in vec3 pos, float e, vec3 dir);
float	loop_circle(vec3 p);
float	circle(vec3 p, float phase);
float	sdTorus( vec3 p, vec2 t, float phase );
float	mylength(vec2 p);
float	mylength(vec3 p);
float	nrand( vec2 n );

float 	t;			// time
vec3	ret_col;	// torus color
vec3	h; 			// light amount

#define I_MAX		200.
#define E			0.0001
#define FAR			200.
//#define PI			3.14159
#define TAU			PI*2.

/*
* Leon's mod polar from : https://www.shadertoy.com/view/XsByWd
*/

vec2 modA (vec2 p, float count) {
    float an = TAU/count;
    float a = atan(p.y,p.x)+an*.5;
    a = mod(a, an)-an*.5;
    return vec2(cos(a),sin(a))*length(p);
}

/*
* end mod polar
*/

vec4 renderMainImage() {
	vec4 c_out = vec4(0.0);
	vec2 f = _xy;

    t  = TIME*.125;
    vec3	col = vec3(0., 0., 0.);
	vec2 R = RENDERSIZE.xy,
          uv  = vec2(f-R/2.) / R.y;
	vec3	dir = camera(uv);
    rotate(dir.xz, -LookXY.x*PI);
    rotate(dir.yz, -LookXY.y*PI);
    vec3	pos = vec3(.0, .0, 0.0);

    h*=0.;
    vec2	inter = (march(pos, dir));
    ret_col = vec3(.490, .482, .470);
    col.xyz = ret_col*((1.-inter.x*.005)+inter.y*.005);
    col += h*.005;

    c_out =  vec4(col,1.0);
    return c_out;
}


float	scene(vec3 p)
{
    float	var;
    float	mind = 1e5;

    p.z-=bass_time*2.;

    p.y-= -4.;
    p.z -= -3.;
    p.x -= 5.;
    p.y -= -9.;
    rotate(p.xz, 1.5);
    vec3 sp = p;
    vec3 op = p;
    
    p = op;
    
    float dd = .030625;
	op.zy = fract(op.zy*dd)-.5;
    op.zy /=dd;

    rotate(op.yz, -smoothTimeC*.05 + op.x*.15);
    float num = .5;
    vec3 rp = op;
    op.x = fract(op.x*num)-.5;
    op.x /= num;
    p = op;
    
    p = p.xyz-vec3(.0, -0., 0.);
    p.z = abs(p.z)-10.;
    float mada = max(max(abs(p.x)-.9, abs(p.y)-.8), abs(p.z)-.62 );
    float light_wave = length(p.zy-vec2(-4.0,-8.))-1.*sin(rp.x*.25+TIME);
    light_wave = abs(light_wave)+.05;
    float light_wave2 = mylength( (fract((sp-vec3(.0+smoothTimeB,.0,0.) ).xyz*.03125)-.5)/.03125 )-16.;
    mada = min(mada, light_wave);
    light_wave2 = abs(light_wave2)+.05;

    h += vec3(.1, .4, .1)*vec3(1.)*.25/max(.001, .051 + 1.*light_wave*light_wave );
    h += vec3(.51, .4, .1)*vec3(1.)*.25/max(.001, .051 + 1.*light_wave2*light_wave2 );
    mind = min(mind, mada);
    
    p = op;
    float ten0 = mylength(p.xy)-.2;
    p = op;
    p.z = fract(p.z*3.)-.5;
    p.z /= 3.;

    rotate(p.xy, (op.z*1.)*1.+smoothTimeC*.5 );

    float ming = mylength(vec3(abs(p.y)-.5, p.xz))-.105;

    float caps = 1e5;
    rotate(p.xz, p.y*5.+TIME*-4.);
    p.xz = modA(p.xz, 2.);
    p.x -= .06125;
    
    caps = max(length(p.xz)-.01, -p.y-0.9*.5);
    caps = max(caps, p.y-.9*.5);
    caps = min(caps, max(max(length(vec2(p.z, (fract(p.y*16.)-.5)/16. ) )-.01, p.x-.01), -p.x-.05) );
    caps = max(caps, +p.y-.5);
    caps = max(caps, -p.y-.5);
    ming = min(ming, caps);
    ming = max(ming, op.z-10.);
    ming = max(ming, -op.z-10.);
    ten0 = ming;
    ten0 = min(ten0, mada);
    p = op;
    p.yx = modA(p.yx, 5.);
    p.y -= 3.;
    p.x = max(abs(fract(p.z)-.5 ), max(abs(p.x), abs(p.y)))-.5125;
    float ten1 = 1e5+length(p.xy)-.02;
    float ten2 = min(ten0, ten1);

    mind = min(mind, ten2);

    h -= (vec3(.05,.05,-0.505))*vec3(1.)*.0125/max(.01, .01+mind*mind);

    return (mind)*.75;
}

vec2	march(vec3 pos, vec3 dir)
{
    vec2	dist = vec2(0.0, 0.0);
    vec3	p = vec3(0.0, 0.0, 0.0);
    vec2	s = vec2(0.0, 0.0);

	    for (float i = -1.; i < I_MAX; ++i)
	    {
	    	p = pos + dir * dist.y;
	        dist.x = scene(p);
	        dist.y += dist.x*.2; // makes artefacts disappear
	        if (dist.x < E || dist.y > FAR)
            {
                break;
            }
	        s.x++;
    }
    s.y = dist.y;
    return (s);
}

// Utilities

float	mylength(vec2 p)
{
	float	ret;

    ret = max(abs(p.x), abs(p.y));

    return ret;
}

float	mylength(vec3 p)
{
	float	ret;

    ret = max(max(abs(p.x), abs(p.y)), abs(p.z));

    return ret;
}

void rotate(inout vec2 v, float angle)
{
	v = vec2(cos(angle)*v.x+sin(angle)*v.y,-sin(angle)*v.x+cos(angle)*v.y);
}

vec2	rot(vec2 p, vec2 ang)
{
	float	c = cos(ang.x);
    float	s = sin(ang.y);
    mat2	m = mat2(c, -s, s, c);
    
    return (p * m);
}

vec3	camera(vec2 uv)
{
    float		fov = 1.;
	vec3		forw  = vec3(0.0, 0.0, -1.0);
	vec3    	right = vec3(1.0, 0.0, 0.0);
	vec3    	up    = vec3(0.0, 1.0, 0.0);

    return (normalize((uv.x) * right + (uv.y) * up + fov * forw));
}

vec3 calcNormal( in vec3 pos, float e, vec3 dir)
{
    vec4 c_out = vec4(0.);
    vec3 eps = vec3(e,0.0,0.0);

    return normalize(vec3(
           march(pos+eps.xyy, dir).y - march(pos-eps.xyy, dir).y,
           march(pos+eps.yxy, dir).y - march(pos-eps.yxy, dir).y,
           march(pos+eps.yyx, dir).y - march(pos-eps.yyx, dir).y ));
	return c_out.rgb; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}