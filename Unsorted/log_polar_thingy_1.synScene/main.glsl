

/*
* License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
* Created by bal-khan
*/



#define I_MAX		400
#define E			0.0001
#define FAR			20.

#define	FUDGE		1.

vec4	march(vec3 pos, vec3 dir);
vec3	camera(vec2 uv);
vec3	calcNormal(in vec3 pos, float e, vec3 dir);
void	rotate(inout vec2 v, float angle);
float	mylength(vec2 p);
float	mylength(vec3 p);

#define PI			3.14159
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
#define SCALE 2.0/PI

float sdCy( vec3 p, vec2 h )
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

vec3 ilogspherical(in vec3 p)
{
	float erho = exp(p.x);
	float sintheta = sin(p.y);
	return vec3(
		erho * sintheta * cos(p.z),
		erho * sintheta * sin(p.z),
		erho * cos(p.y)
	);
}

void pR(inout vec2 p, float a) {
	p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

float side = .0, shorten = 20.50;
float layer(in vec3 p, in float twost)
{
	float ret = 1e5;
	pR(p.yz, twost);
    vec3 op = p;
	p.xyz = abs(p.xyz) - 2.750;
    float a = .0+p.x*.5;
    /*
    p.xy= modA(p.xy, 3.);p.x-=5.;
    p.zy= modA(p.zy, 3.);p.z-=2.;
    p.yx= modA(p.yx, 3.);p.y-=4.;
    */
    ret = length(vec2(max(p.y, p.x), max(p.z, min(p.x,p.y) ) ) )-.25;

    float sf = .1525+-1.5100251010/(length(p)*length(p)*0.015+10.01)-.0;
    ret = min(ret,
              length(
              fract(sf*vec2(p.z, p.x))/sf-.5
              )-.10105
              );
    ret = min(ret,
              length(
              fract(sf*vec2(p.y, p.x))/sf-.5
              )-.10105
              );
    ret = min(ret,
              length(
              fract(sf*vec2(p.z, p.y))/sf-.5
              )-.105105
              );
    ret = max(ret,
              -(length(fract(sf*p)/sf-.5)-.5 )
              );
    op.xy = modA(op.xy, 50.);
    op.zx = modA(op.zx, 50.);
    float shell_cubes = (mylength(op-vec3(.0,-.0,1.0505)+-.00)-0.0251250905125995);
    ret = min(ret, shell_cubes );
    ret = min(ret, (length(p+-4.50705)-3.50905995) );

	return ret;
}

float sdf(in vec3 pin)
{
    float dens = .25;
    float idens = 1./dens;
    float twist = 1.5, stepZoom = 1.;
	float r = length(pin);
	vec3 p = vec3(log(r), acos(pin.z / length(pin)), .0+TIME*.1*-1.0+atan(pin.y, pin.x));

	// Apply rho-translation, which yields zooming
	p.x -= TIME*1.2;

	// find the scaling factor for the current tile
	float xstep = floor(p.x*dens) + (TIME*1.2)*dens;
	
	// Turn tiled coordinates into single-tile coordinates
	p.x = mod(p.x, idens);

	// Apply inverse log-spherical map (box tile -> shell tile)
	p = ilogspherical(p);

	float ret = 1e5;
	ret = min(ret, layer(p/stepZoom, (xstep+1.0)*twist)*stepZoom);

	// Compensate for scaling applied so far
	ret = ret * exp(xstep*idens) / shorten;

	return ret;
}

vec4 renderMainImage() {
	vec4 c_out = vec4(0.0);
	vec2 f = _xy;

    vec3	col = vec3(0., 0., 0.);
    vec2	R = RENDERSIZE.xy;
	vec2	uv  = vec2((f.xy-.5*R.xy)/R.y);
	vec3	dir = normalize(vec3(uv*vec2(1.,-1.), 1.));
    vec3	pos = vec3(-.0, -.2105017501050*.450*.0+.033*.0, -2.0*.0-1.52525045);
    

    vec4	inter = march(pos, dir);
    if (inter.y <= E*1.)
    {
        vec3	v = pos+(inter.w-E*0.)*dir;
        vec3	n = calcNormal(v, E*1., dir);
        vec3	ev = normalize(v - pos);
		vec3	ref_ev = reflect(ev, n);
        vec3	light_pos   = vec3(-100.0, 60.0, -50.0);
		vec3	light_color = vec3(.1, .4, .7);
        vec3	vl = normalize(light_pos - v);
		float	diffuse  = max(0.0, 1.-dot(vl, n))+.0*max(0.0, dot(vl, n));
		float	specular = pow(max(0.0, dot(vl, ref_ev)), 3.);
        col.xyz += light_color * (specular*1.0)+ diffuse * vec3(.51,.515, .53);
    }
    else
        col *= .0;
    c_out =  vec4(col, 1.);
}

vec4	march(vec3 pos, vec3 dir)
{
    vec2	dist = vec2(0.0, 0.0);
    vec3	p = vec3(0.0, 0.0, 0.0);
    vec4	ret = vec4(0.0, 0.0, 0.0, 0.0);

    for (int i = -1; i < I_MAX; ++i)
    {
    	p = pos + dir * dist.y;
        dist.x = sdf(p);
        dist.y += dist.x*FUDGE;
        if ( dist.x < E || dist.y > FAR)
            break;
    }
    ret.w = dist.y;
    ret.y = dist.x;
    return (ret);
}

// Utilities

float	mylength(vec3 p)
{
	float	ret = 1e5;

    p = abs(p);
    ret = max(p.x, max(p.y, p.z));
    
    return ret;
}

float	mylength(vec2 p)
{
	float	ret = 1e5;

    p = abs(p);
    
    ret = max(p.x, p.y);

    return ret;
}

float smin( float a, float b, float k )
{
	float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
}


void rotate(inout vec2 v, float angle)
{
	v = vec2(cos(angle)*v.x+sin(angle)*v.y,-sin(angle)*v.x+cos(angle)*v.y);
}

vec3 calcNormal( in vec3 pos, float e, vec3 dir)
{
    vec3 eps = vec3(e,0.0,0.0);

	return normalize(vec3(
           sdf(pos+eps.xyy) - sdf(pos-eps.xyy),
           sdf(pos+eps.yxy) - sdf(pos-eps.yxy),
           sdf(pos+eps.yyx) - sdf(pos-eps.yyx) ));
}

vec3	camera(vec2 uv)
{
    float		fov = 1.;
	vec3		forw  = vec3(0.0, 0.0, -1.0);
	vec3    	right = vec3(1.0, 0.0, 0.0);
	vec3    	up    = vec3(0.0, 1.0, 0.0);

    return (normalize((uv.x) * right + (uv.y) * up + fov * forw));
	return c_out; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}