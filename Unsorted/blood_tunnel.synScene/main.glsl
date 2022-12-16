

// Created by Stephane Cuillerdier - Aiekick/2015 (twitter:@aiekick)
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Tuned via XShade (http://www.funparadigm.com/xshade/)
float hash(float n)
{
	return fract(sin(n) * 43728.1453);
}

vec3 noise(vec3 x)
{
	vec3 p = floor(x);
	vec3 f = fract(x);

	f = f * f * (3.0 - 2.0 * f);
	float n = p.x + p.y * 55.0 + p.z * 101.0 ;

	vec3 n3D;
	n3D = vec3(mix(
	mix(
		mix(hash(n), hash(n + 1.0), f.x),
		mix(hash(n+55.0), hash(n + 56.0), f.x),
		f.y),
	mix(
		mix(hash(n+101.0), hash(n + 102.0), f.x),
		mix(hash(n+156.0), hash(n + 157.0), f.x),
		f.y),
	f.z));
	return n3D;
	
}

vec2 path(float t){	return vec2(cos(t*0.2)*0.5,0.5* sin(t*0.2)) * 2.;}

vec2 df(vec3 p)
{
	p.xy -= path(p.z);
	float hex = (max(abs(p.x) + p.y, -p.y) + max(abs(p.x) - p.y, p.y));
	float tex = texture(image47, (abs(p.xz) + abs(p.yz))*0.28).r * atan(p.x,p.y)/3.14159 * 1.;
	p*=3.;
	float maze = cos(3.14 * (p.x + p.y * sign(sin(1e3 * length(floor(p.xy))))));
	return vec2(0.8 * hex - maze * 0.12 + tex - 4.8 , 2.);
}

vec3 nor( in vec3 pos, float prec )
{
	vec3 eps = vec3( prec, 0., 0. );
	vec3 nor = vec3(
	    df(pos+eps.xyy).x - df(pos-eps.xyy).x,
	    df(pos+eps.yxy).x - df(pos-eps.yxy).x,
	    df(pos+eps.yyx).x - df(pos-eps.yyx).x );
	return normalize(nor);
}

// return color from temperature 
//http://www.physics.sfasu.edu/astro/color/blackbody.html
//http://www.vendian.org/mncharity/dir3/blackbody/
//http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color.html
vec3 blackbody(float Temp)
{
	vec3 col = vec3(255.);
    col.x = 56100000. * pow(Temp,(-3. / 2.)) + 148.;
   	col.y = 100.04 * log(Temp) - 623.6;
   	if (Temp > 6500.) col.y = 35200000. * pow(Temp,(-3. / 2.)) + 184.;
   	col.z = 194.18 * log(Temp) - 1448.6;
   	col = clamp(col, 0., 255.)/255.;
    if (Temp < 1000.) col *= Temp/1000.;
   	return col;
}

// get density of the df at surfPoint
// ratio between constant step and df value
float SubDensity(vec3 surfPoint, float prec, float ms) 
{
	vec3 n = vec3(0);
	float s = 0.;
    const int iter = 10;
	for (int i=0;i<iter;i++)
	{
		n = nor(surfPoint,prec); 
		surfPoint = surfPoint - n * ms; 
		s += df(surfPoint).x;
	}
	return 1.-s/(ms*float(iter)); // s < 0. => inside df
}

float SubDensity(vec3 p, float s) 
{
	vec3 n = nor(p,s); 							// precise normale at surf point
	return df(p - n * s).x;						// ratio between df step and constant step
}

// from shane sahders
// Tri-Planar blending function. Based on an old Nvidia writeup:
// GPU Gems 3 - Ryan Geiss: http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
vec3 tex3D( sampler2D tex, in vec3 p, in vec3 n )
{
    n = max((abs(n) - .2)*7., .001);
    n /= (n.x + n.y + n.z );  
	p = (texture(tex, p.yz)*n.x + texture(tex, p.zx)*n.y + texture(tex, p.xy)*n.z).xyz;
    return p*p;
}

vec4 shade(vec3 ro, vec3 rd, float d, vec3 lp, float li)
{
	vec3 p = ro + rd * d;											// surface point
	float sb = SubDensity(p, 0.01, 0.076);							// deep subdensity (10 iterations)
	vec3 bb = blackbody(40.*sb+0.);								// bb
	vec3 ld = normalize(lp-p); 										// light dir
	vec3 n = nor(p, 0.08);	// normal at surface point
    
    // derived from bumpmap func from shane
    const vec2 e = vec2(0.1, 0);
    mat3 m = mat3( tex3D(image47, noise(e.xyy), n), tex3D(image47, noise(e.yxy), n), tex3D(image47, e.yyx, n));
   	vec3 g = vec3(1,0,0) * m * 20.;
    g -= n * dot(n, g);
    n =  normalize( n + g );
    
	vec3 refl = reflect(rd,n);										// reflected ray dir at surf point 
	float amb = 0.1242; 											// ambiance factor
	float diff = clamp( dot( n, ld ), 0.0, 1.0 ); 					// diffuse
	float fre = pow( clamp( 1. + dot(n,rd),0.0,1.0), 4. ); 			// fresnel
	float spe = pow(clamp( dot( refl, ld ), 0.0, 1.0 ),16.);		// specular
	float sss = 1. - SubDensity(p, 7.8); 							// one step sub density of df
	return vec4(
        (diff + fre + bb.x * sss) * amb * li + spe, 
        (diff + fre + bb * sb * 0.608 + sss * 0.352) * amb * li + spe * 0.612 	
    );
}

vec3 cam(vec2 uv, vec3 ro, vec3 cv, float t)
{
	vec3 cu = normalize(vec3(0,1,0));
  	vec3 z = normalize(cv-ro);
    vec3 x = normalize(cross(cu,z));
  	vec3 y = cross(z,x);
  	return normalize(z + uv.x*x + uv.y*y);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    float t = -smoothTime;

    vec2 si = RENDERSIZE.xy;
    vec2 g = fragCoord;
    
    vec2 uv = (g+g-si)/si.y;
    uv.xy += _uvc*PI*FOV;
    vec3 col = vec3(0.);
    
	vec3 ro = vec3(path(t),t);
  	vec3 cv = vec3(path(t),t+.1);
	vec3 rd = cam(uv, ro, cv, t);
    rd.yz = _rotate(rd.yz, LookXY.y*PI);
    rd.xz = _rotate(rd.xz, LookXY.x*PI);
    rd.xy = _rotate(rd.xy, -Rotate*PI);

	float md = 20.,s = 1.,d = 0.;
	
	const float iter = 250.;
    for(float i=0.;i<iter;i++)
    {      
        if (s<.025*log(d*d/s/500.)||d>md) break;
        s = df(ro+rd*d).x;
		d += s * 0.3;
    }
    
    if (_mouse.z > 0.)
        if (_mouse.x < si.x * 0.5)
			d = min(d,d*d*d/md); // weird but cool also
       	else
    		d = min(d*.6,d*d*d/md); // another flame version
    
    fragColor.rgb = mix(
        shade(ro, rd, d, ro, 1.2).yzw, 
        vec3(.2,0,0), 
        1.-exp(-0.01*d*d)
    );
    
    fragColor.a = 1.;
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}