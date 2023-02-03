

// Created by Stephane Cuillerdier - Aiekick/2017 (twitter:@aiekick)
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Tuned via XShade (http://www.funparadigm.com/xshade/)

const vec3 ld = vec3(0.,1., .5);

float t = 0.;

vec2 m;

float fullAtan(vec2 p)
{
    return step(0.0,-p.x)*3.1415926535 + sign(p.x) * atan(p.x, sign(p.x) * p.y);
}

float fractus(vec2 p, vec2 v)
{
	vec2 z = p;
    vec2 c = v;
	float k = 1., h = 1.0;    
    float t = (sin(TIME * .5)*.5+.5) * 5.;
    for (float i=0.;i<5.;i++)
    {
        if (i > t) break;
        h *= 4.*k;
		k = dot(z,z);
        if(k > 4.) break;
		z = vec2(z.x * z.x - z.y * z.y, 2. * z.x * z.y) + c;
    }
    float d = sqrt(k/h)*log(k);
    
    // next iteration
    if (k > 4.)
    	z = vec2(z.x * z.x - z.y * z.y, 2. * z.x * z.y) + c;
    h *= 4.*k;
    k = dot(z,z);
    float d1 = sqrt(k/h)*log(k);
    
    // df blending
	return mix(d, d1, fract(t));
}

vec4 dfFractus(vec3 p)
{
	float a = fullAtan(p.xz); // axis y
    
    vec2 c;
    c.x = mix(0.2, -0.5, sin(a * 2.));
    c.y = mix(0.5, 0.0, sin(a * 3.));
    
    float path = length(p.xz) - 3.;
    
    vec2 rev = vec2(path, p.y);
    float aa = a + TIME;
    rev *= mat2(cos(aa),-sin(aa),sin(aa),cos(aa)); // rot near axis y
	
	return vec4(fractus(rev, c) - 0.0, rev, 0);
}

vec4 df(vec3 p)
{
	return dfFractus(p); 
}

vec3 nor( vec3 p, float prec )
{
    vec2 e = vec2( prec, 0. );
    vec3 n = vec3(
		df(p+e.xyy).x - df(p-e.xyy).x,
		df(p+e.yxy).x - df(p-e.yxy).x,
		df(p+e.yyx).x - df(p-e.yyx).x );
    return normalize(n);
}

// from iq code
float softshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax )
{
	float res = 1.0;
    float t = mint;
    for( int i=0; i<1; i++ )
    {
		float h = df( ro + rd*t ).x;
        res = min( res, 8.0*h/t );
        t += h*.25;
        if( h<0.001 || t>tmax ) break;
    }
    return clamp( res, 0., 1. );
}

// from iq code
float calcAO( in vec3 pos, in vec3 nor )
{
	float occ = 0.0;
    float sca = 1.0;
    for( int i=0; i<10; i++ )
    {
        float hr = 0.01 + 0.12*float(i)/4.0;
        vec3 aopos =  nor * hr + pos;
        float dd = df( aopos ).x;
        occ += -(dd-hr)*sca;
        sca *= 0.95;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );    
}

//--------------------------------------------------------------------------
// Grab all sky information for a given ray from camera
// from Dave Hoskins // https://www.shadertoy.com/view/Xsf3zX
vec3 GetSky(in vec3 rd, in vec3 sunDir, in vec3 sunCol)
{
	float sunAmount = max( dot( rd, sunDir), 0.0 );
	float v = pow(1.0-max(rd.y,0.0),6.);
	vec3  sky = mix(vec3(.1, .2, .3), vec3(.32, .32, .32), v);
	sky = sky + sunCol * sunAmount * sunAmount * .25;
	sky = sky + sunCol * min(pow(sunAmount, 800.0)*1.5, .3);
	return clamp(sky, 0.0, 1.0);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    m = _mouse.xy/RENDERSIZE.xy;
    if (m.x == 0.) m.x = .1;
    
	fragColor = vec4(1);
	
	vec2 g = fragCoord.xy;
	vec2 si = RENDERSIZE.xy;
	vec2 uv = (2.*g-si)/min(si.x, si.y);
	
    float a = 1.57;
	vec3 rayOrg = vec3(cos(a),1.5,sin(a)) * 5.;
	vec3 camUp = vec3(0,1,0);
	vec3 camOrg = vec3(0,-1.5,0);
	
	float fov = .5;
	vec3 axisZ = normalize(camOrg - rayOrg);
	vec3 axisX = normalize(cross(camUp, axisZ));
	vec3 axisY = normalize(cross(axisZ, axisX));
	vec3 rayDir = normalize(axisZ + fov * uv.x * axisX + fov * uv.y * axisY);
	
	float s = 0.;
    float d = 0.;
	float dMax = 20.;
	float count = 0.;
	for (float i=0.; i<500.; i++)
	{
		if (d*d/s>1e6 || d>dMax) break;
        s = df(rayOrg + rayDir * d).x;
		d += s * 0.2;
		count+=1.;
	}
	
    vec3 sky = GetSky(rayDir, ld, vec3(1.5));
    
   	/*vec4 mat = df(rayOrg + rayDir * d);
   	vec3 tex = texture(image6, mat.yz * 0.8).rgb;
    d += dot(tex, vec3(0.05));*/
	
	if (d<dMax)
	{
        
		vec3 p = rayOrg + rayDir * d;
        vec3 n = nor(p, 0.001);
		vec4 mat = df(p);
        
        vec3 tex = texture(image6, mat.yz).rgb; // wood
        
		// iq lighting
		float occ = calcAO( p, n );
        float amb = clamp( 0.5+0.5*n.y, 0.0, 1.0 );
        float dif = clamp( dot( n, ld ), 0.0, 1.0 ) * (df(p+n*1.16).x);
        float spe = pow(clamp( dot( rayDir, ld ), 0.0, 1.0 ),16.0);
        float sss = df(p - n*0.001).x/0.01;
	
        dif *= softshadow( p, ld, 0.1, 2. );
        
        vec3 brdf = vec3(-0.5) * tex;
        brdf += 0.5*dif*vec3(1.00,0.90,0.60);
        brdf += 0.4*spe*vec3(1.00,0.90,0.60)*dif;
        brdf += 0.3*amb*vec3(0.50,0.70,1.00)*occ;
        brdf += 0.5*(1.-sss)*vec3(0.2,0.7,0.2);
        fragColor.rgb *= brdf;

        fragColor.rgb = mix( fragColor.rgb, sky, 1.0-exp( -0.01*d*d*count/150. ) ); 
	}
	else
	{
		fragColor.rgb = sky;
	}
    
    fragColor.rgb = sqrt(fragColor.rgb * fragColor.rgb * 0.8);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}