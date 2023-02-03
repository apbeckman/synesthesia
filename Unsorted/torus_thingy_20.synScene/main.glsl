

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
#define FAR			110.
#define PI			3.14159
#define TAU			PI*2.

vec4 renderMainImage() {
	vec4 c_out = vec4(0.0);
	vec2 f = _xy;

    t  = TIME*.125;
    vec3	col = vec3(0., 0., 0.);
	vec2 R = RENDERSIZE.xy,
          uv  = vec2(f-R/2.) / R.y;
	vec3	dir = camera(uv);
    vec3	pos = vec3(.0, .0, 60.0);

    h*=0.;
    vec2	inter = (march(pos, dir));
    if (inter.y >= FAR)
        ret_col = vec3(.90, .82, .70);
    col.xyz = ret_col*(1.-inter.x*.005);
    col += h*.005;
    c_out =  vec4(col,1.0);
    return c_out;
}

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

float  plane(vec3 p, vec3 origin, vec3 normal){ 
   return dot(p - origin,normal);   
}

float	scene(vec3 p)
{  
    p.z += -26.;
    float	var;
    float	mind = 1e5;
    float	cage = 1e5;
    rotate(p.xz, 1.57-.35*TIME );
    rotate(p.yz, 1.57-.25*TIME );
    vec3 op = p;
    var  = atan(op.x,op.y);
    
    p = op;
    vec3 pp = p;
    p.xy = modA(p.xy, 100.);
    p.x -= 10.+.0*abs(5.*sin(TIME*.125));
    p.zx = modA(p.zx, 100.);
    p.z -= 4.+2.*sin(atan(pp.x,pp.y)*2.+atan(p.z)*8.-TIME*0.);
    
    mind = length(p.zx)-.0506125;
    mind = min(mind,
        length(p.zy)-.0506125
        );
    mind = min(mind,
               length(p)-.10010506125
               );
    mind = min(mind,
               length(p.z+.5*sin(TIME*.5)+.45)-.01001
               );

    ret_col = vec3(.90, .82, .70);
    float ball = 1e5;
    h += vec3(.5, .5, .9)*1./(ball*ball*1.+.25);
    h -= vec3(-.50,.1250,1.)*vec3(1.)*.0125/(.01+mind*mind);
    h -= vec3(.05,.05,1.)*vec3(1.)*.0125/(.01+mind*mind);
    return (mind);
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

    ret = max( abs(p.x)+.5*abs(p.y), abs(p.y)+.5*abs(p.x) );
    
    return ret;
}

float	mylength(vec3 p)
{
	float	ret;

    ret = max( abs(p.x)+.5*abs(p.y), abs(p.y)+.5*abs(p.x) );
    ret = max(abs(p.z)+.5*abs(p.x), ret);
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
  	vec3 c_out = vec3(0.0);
	vec2 f = _xy;

    vec3 eps = vec3(e,0.0,0.0);

    return normalize(vec3(
           march(pos+eps.xyy, dir).y - march(pos-eps.xyy, dir).y,
           march(pos+eps.yxy, dir).y - march(pos-eps.yxy, dir).y,
           march(pos+eps.yyx, dir).y - march(pos-eps.yyx, dir).y ));
	return c_out; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}