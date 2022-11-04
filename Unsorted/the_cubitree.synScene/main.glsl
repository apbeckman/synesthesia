

			//******** BuffA Code Begins ********

// Created by anatole duprat - XT95/2014
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.


//Maths
//const float PI = 3.14159265;

mat3 rotate( in vec3 v, in float angle)
{
	float c = cos(radians(angle));
	float s = sin(radians(angle));
	
	return mat3(c + (1.0 - c) * v.x * v.x, (1.0 - c) * v.x * v.y - s * v.z, (1.0 - c) * v.x * v.z + s * v.y,
		(1.0 - c) * v.x * v.y + s * v.z, c + (1.0 - c) * v.y * v.y, (1.0 - c) * v.y * v.z - s * v.x,
		(1.0 - c) * v.x * v.z - s * v.y, (1.0 - c) * v.y * v.z + s * v.x, c + (1.0 - c) * v.z * v.z
		);
}

mat3 lookat( in vec3 fw, in vec3 up )
{
	fw = normalize(fw);
	vec3 rt = normalize( cross(fw, normalize(up)) );
	return mat3( rt, cross(rt, fw), fw );
}


//Raymarching 
float map( in vec3 p );

float box( in vec3 p, in vec3 data )
{
    return max(max(abs(p.x)-data.x,abs(p.y)-data.y),abs(p.z)-data.z);
}

float sphere( in vec3 p, in float size)
{
	return length(p)-size;
}

vec4 raymarche( in vec3 org, in vec3 dir, in vec2 nfplane )
{
	float d = 1.0, g = 0.0, t = 0.0;
	vec3 p = org+dir*nfplane.x;
	
	for(int i=0; i<42; i++)
	{
		if( d > 0.001 && t < nfplane.y )
		{
			d = map(p);
			t += d;
			p += d * dir;
			g += 1./42.;
		}
	}
	
	return vec4(p,g);
}

vec3 normal( in vec3 p )
{
	vec3 eps = vec3(0.01, 0.0, 0.0);
	return normalize( vec3(
		map(p+eps.xyy)-map(p-eps.xyy),
		map(p+eps.yxy)-map(p-eps.yxy),
		map(p+eps.yyx)-map(p-eps.yyx)
	) );
}

float ambiantOcclusion( in vec3 p, in vec3 n, in float d )
{
    float dlt = 0.0;
    float oc = 0.0;
    
    for(int i=1; i<=6; i++)
    {
		dlt = d*float(i)/6.;
		oc += (dlt - map(p+n*dlt))/exp(dlt);
    }
    oc /= 6.;
    
    return clamp(pow(1.-oc,d), 0.0, 1.0);
}




//Geometry
float ill = 0.;
#define impulsTime (smoothTimeC+sin(smoothTimeC+PI))
float map( in vec3 p )
{
	float d = p.y;
	vec3 pp = p;
	ill = 0.;
	
	//mirrors
	p = abs(p);
	p = rotate(vec3(-1.,0.,0.),40.)*p;
	p = abs(p);
	p = rotate(vec3(0.,1.,0.),45.)*p;
	p = abs(p);
	
	//make a branch of cubes
	for(int i=0; i<15; i++)
	{
		p -= vec3(.25);
		p = rotate( normalize( vec3(.5, .25, 1.0 ) ), 20.+pp.x+pp.y+pp.z )*p;
		
		
		float size = cos(float(i)/20.*PI*2.-impulsTime);
		float dbox = box( p, vec3( (1.1-float(i)/20.)*.25 + pow(size*.4+.4,10.) ) );
	
		if( dbox < d)
		{
			d = dbox;
			ill = pow(size*.5+.5, 10.);
		}
	
	}
	//add another one iteration with a sphere
	p -= vec3(.25);
	p = rotate( normalize( vec3(.5, .25, 1.0 ) ), 20.+pp.x+pp.y+pp.z )*p;
	d = min(d, sphere(p,.25) );
	
	return d;
}

//Shading
vec3 ldir = normalize( vec3(.267,.358,.90) );
vec3 sky( in vec3 dir )
{
	vec3 col = mix( vec3(40., 34., 30.), vec3(18., 28., 44.), min( abs(dir.y)*2.+.5, 1. ) )/255.*.5;
	col *= (1. + vec3(1.,.7,.3)/sqrt(length(dir-ldir))*4.); //sun
	
	return col;
}
vec3 shade( in vec4 p, in vec3 n, in vec3 org, in vec3 dir )
{		
	//direct lighting
	vec3 col = vec3(.1);
	col += pow(sky(vec3(1.,0.,0.))*max( dot(n, ldir), 0.)*2., vec3(2.));
	
	//illumination of the tree
	col += mix( vec3(1.,.3,.1), vec3(.1, .7, .1), length(p.xyz)/8.)*ill*p.w*5.;
	
	//ao
	col *= pow( ambiantOcclusion(p.xyz,n,1.) , 1.5 );
	
	//fog/sky
	col = mix(col, sky(dir), vec3(1.)*min( pow( distance(p.xyz,org)/20., 2. ), 1. ) );
	
	return col;
}

vec2 hash2( in float n )
{
    return fract(sin(vec2(n,n+1.0))*vec2(43758.5453123,22578.1459123));
}
//Main
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	//screen coords

    vec2 o = hash2( float(FRAMECOUNT) ) - 0.5;
    vec2 v = (-RENDERSIZE.xy + 2.0*(fragCoord+o))/ RENDERSIZE.y;
	vec2 q = fragCoord.xy/RENDERSIZE.xy;
	
	//camera ray
	float ctime = (TIME+140.)*.025;
	vec3 org = vec3( cos(ctime)*10., 2.+cos(ctime), sin(ctime)*10. );
	vec3 dir = normalize( vec3(v.x, v.y, 1.5) );
	dir = lookat( -org + vec3(0., 2., 0.), vec3(0., 1., 0.) ) * dir;
	
	//classic raymarching by distance field
	vec4 p = raymarche(org, dir, vec2(4., 20.) );
	vec3 n = normal(p.xyz);
	vec3 col = shade(p, n, org, dir);
	
    fragColor = mix(vec4(col, 1.), texture(BuffA,q), .8);
	return fragColor; 
 } 


// Created by anatole duprat - XT95/2014
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.


// cheap bloom based on mipmap !
vec3 bloom( sampler2D tex, vec2 p) {
    vec3 col = vec3(0.);
    for (int i=1; i<9; i++)
        col += textureLod(tex, p, float(i)).rgb / float(9-i);
    
    return max(col-.6, vec3(0.));
}


//Main
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


	vec2 uv = fragCoord.xy/RENDERSIZE.xy;
    vec3 col = texture(BuffA, uv).rgb;
    
	//post process
    vec2 q = uv;
    col = pow( col*2., vec3(1.75) );
	col *= sqrt( 32.0*q.x*q.y*(1.0-q.x)*(1.0-q.y) ); //from iq
    col += bloom(BuffA, uv)*2.;
	
	//fragColor = vec4( col*min(TIME/5., 1.), 1. );
    
    
    fragColor = vec4(col, 1.);
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