

// Created by Stephane Cuillerdier - Aiekick/2015
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.


float random (in float x) {
	return fract(sin(x) * 1e4);
}

// Based on Morgan McGuire @morgan3d
// https://www.shadertoy.com/view/4dS3Wd
float n3D (in vec3 p) {
	const vec3 step = vec3(110.0, 241.0, 171.0);

	vec3 i = floor(p);
	vec3 f = fract(p);

    // For performance, compute the base input to a
    // 1D random from the integer part of the
    // argument and the incremental change to the
    // 1D based on the 3D -> 1D wrapping
	float n = dot(i, step);

	vec3 u = f * f * (3.0 - 2.0 * f);

	return mix( mix(mix(random(n + dot(step, vec3(0,0,0))),
                        random(n + dot(step, vec3(1,0,0))),
                        u.x),
                    mix(random(n + dot(step, vec3(0,1,0))),
                        random(n + dot(step, vec3(1,1,0))),
                        u.x),
                u.y),
                mix(mix(random(n + dot(step, vec3(0,0,1))),
                        random(n + dot(step, vec3(1,0,1))),
                        u.x),
                    mix(random(n + dot(step, vec3(0,1,1))),
                        random(n + dot(step, vec3(1,1,1))),
                        u.x),
                u.y),
                u.z);
}

float Noise = 0.065*n3D(vec3(smoothTime*0.025));

float dstepf = 0.0;
#define uScreenSize RENDERSIZE.xy
#define uTime TIME

const vec2 RMPrec = vec2(0.3, 0.01); // ray marching tolerance precision // low, high
const vec2 DPrec = vec2(1e-3, 30.); // ray marching distance precision // low, high
    
// light
const vec3 LCol = vec3(0.8,0.5,0.2);
const vec3 LPos = vec3(-0.6, 0.7, -0.5);
const vec3 LAmb = vec3( 0. );
const vec3 LDif = vec3( 1. , 0.5, 0. );
const vec3 LSpe = vec3( 0.8 );

// material
const vec3 MCol = vec3(0.);
const vec3 MAmb = vec3( 0. );
const vec3 MDif = vec3( 1. , 0.5, 0. );
const vec3 MSpe = vec3( 0.6, 0.6, 0.6 );
const float MShi =30.;
    
#define mPi PI
#define m2Pi PI*2.0

float flashLevel = mix(0.0, pow(syn_HighHits, 2.0), flashing);

vec2 s,g,uv,m;

vec2 uvs(vec3 p) // uv sphere
{
	p = normalize(p);
	vec2 sp;
	sp.x = atan(p.z, p.x) / (m2Pi+1.27);
	sp.y = asin(p.y) / (mPi);
	return sp;
}

float smin( float a, float b, float k )
{
	float h = clamp( .5+.5*(b-a)/k, 0., 1. );
	return mix( b, a, h ) - k*h*(1.-h);
}

float sdCyl( vec3 p, vec2 h )
{
  	vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  	return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

// return color from temperature 
//http://www.physics.sfasu.edu/astro/color/blackbody.html
//http://www.vendian.org/mncharity/dir3/blackbody/
//http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color.html
vec3 getHotColor(float Temp) // blackbody temperature color
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

vec2 getTemp(vec3 p)
{
	p*=2.;
	float r = fract(p.x+p.z);
	return vec2(dot(p,p)*(1000.)*r,r);
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
	sky += _loadUserImage(sunDir.xz*vec2(0.89,1.0)+vec2(0.175,0.0)).rgb*syn_HighHits;
	return clamp(sky, 0.0, 1.0);
}

vec3 colorPaletteChooser(float colReg, float var){
  vec3 paletteCol = vec3(0.5);
  if (colReg < 1.0){
    paletteCol = _palette(var, vec3(0.500, 0.500, 0.520), vec3(0.500, 0.500, 0.500), vec3(0.780, 0.765, 0.750), vec3(0.360, 0.570, 0.680));
  } else if (colReg < 2.0){
    paletteCol = _palette(var, vec3(0.500, 0.580, 0.500), vec3(0.190, 0.460, 0.500), vec3(0.760, 0.740, 0.580), vec3(1.000, 0.300, 0.490));
  } else if (colReg < 3.0){
    paletteCol = _palette(var,vec3(0.129, 0.359, 0.072), vec3(0.933, 0.561, 0.616), vec3(0.334, 1.013, 0.882), vec3(0.597, 0.050, 0.885));
  } else if (colReg < 4.0){
    paletteCol = _palette(var,  vec3(-0.060, -0.340, 0.100), vec3(0.940, 0.840, 0.713), vec3(0.600, 0.735, 0.7191), vec3(0.500, 0.260, 0.335));
  } else if (colReg < 5.0){
    paletteCol = _palette(1.0-var,vec3(1.040, 0.180, 0.260), vec3(0.053, 0.775, 0.330), vec3(0.142, 0.523, 0.800), vec3(0.242, 0.887, 0.000));
  }
  return paletteCol;
}

vec3 getCol(vec3 p){
	float temp = getTemp(p).x;
	vec3 col = getHotColor(temp);
	col = mix(col, colorPaletteChooser(palette_select, clamp(temp,0.0,1.0))*sin(p.z*0.2+smoothTimeB*0.1), color_mod);
	return col;
}

vec3 opRep( vec3 p, vec3 c )
{
	p.xy += _uvc;
    vec3 q = mod(p,c)-0.5*c;
    return q;
}

float sdSphere( vec3 p, float s )
{
  return length(p*PI*_uvc.xyy)-s;
}

float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return mix(length(q)-t.y, length(p)-t.x*01.25, 1.0);
}

float sdCylinder( vec3 p, vec3 c )
{
	c.xy *=1.0- _uvc*p.xy;
  return length(p.xz-c.xy)-c.z;
}

vec2 map(vec3 p)
{
	p = mix((p), floor(p)+_pixelate(fract(p), voxel_size), voxelate);
	
	vec2 res = vec2(0.);
    
	//float t = sin(smoothTimeC*0.2)*.5+.5;
	float t = sin(smoothTimeC*0.2)*.5+.5;
	t*=2.;

	float rugo = cos(0.5+2.*p.x+smoothTimeC*0.75)*sin(1.5*p.z)*sin(3.*p.y+syn_BeatTime*beat_twitching)*cos(1.3);
    
	vec3 sci = vec3(1.0-0.5*syn_Intensity* syn_Presence,1.25,0.9);// scale in
	vec3 pob = vec3(0.,5.,0.);// pos bottom
	vec3 pocy = vec3(0.,-5.,0.);//pos cyl
        
	vec3 col = getCol(p);

	float diamhole = 2.1*syn_Presence*(0.5+0.5*syn_Intensity)*Hole;
	float spi = length(p*sci) - 4.5 - rugo;//in
	float spb = length(p+pob) - 4.5 + dot(col, vec3(0.01));//bottom
	float spo = spi - 1. + flashLevel*0.8 + 2.0*astral_plane;//out

	float cyl = sdCyl(p+pocy, vec2(diamhole,4.));//top hole
    
	float disp = dot(col, vec3(0.03));
    
	float spicyl = smin(spi,cyl,0.6) + disp;
	float cavern = smin(max(-spicyl, spo ), spb, 3.5);
  
	vec3 dotPos = opRep(p+vec3(0.0, smoothTimeC*0.125+smoothTime*0.125, 0.0), vec3(5.0, 4.0, 3.0));

	cavern = smin(cavern, sdTorus(dotPos+0.25, vec2(01.1,0.001))+step(5.0, p.y)*1000.0+droppers*10000.0, 2.0);
	cavern = smin(cavern, sdCylinder(dotPos+0.25, vec3(0.1,0.001, 0.05))+step(5.0, p.y)*1000.0+cylinders*100000.0, 0.5+syn_BPMSin4*0.5);

  dstepf += 0.01;
    
	return vec2(cavern, 1.0);
}


vec3 nor(vec3 p, float prec)
{
	vec2 e = vec2(prec, 0.);
    
	vec3 n;
    
	n.x = map(p+e.xyy).x - map(p-e.xyy).x; 
	n.y = map(p+e.yxy).x - map(p-e.yxy).x; 
	n.z = map(p+e.yyx).x - map(p-e.yyx).x;  
    
	return normalize(n); 
}

vec3 ads( vec3 p, vec3 n )
{
	vec3 ldif = normalize( LPos - p);
	vec3 vv = normalize( vec3(0.) - p );
	vec3 refl = reflect( vec3(0.) - ldif, n );
    
	vec3 amb = MAmb*LAmb;
	vec3 dif = max(0., dot(ldif, n.xyz)) * MDif * LDif;
	vec3 spe = vec3( 0. );
	if( dot(ldif, vv) > 0.)
		spe = pow(max(0., dot(vv,refl)),MShi)*MSpe*LSpe;
    
	vec3 col = amb*1.2 + dif*1.5 + spe*0.8;
    
	return col;
}

vec4 scn(vec4 col, vec3 ro, vec3 rd)
{
	float s = DPrec.x;
	float d = 0.;
	vec3 p = ro+rd*d;
	vec4 c = col;
    
	float b = 0.35;
    
	vec4 glow = vec4(0.0);

	//MAPPING
	for(int i=0;i<200;i++)
	{
		if(s<DPrec.x||s>DPrec.y) break;
		vec2 dat = map(p);
		s = dat.x;
		// glow += d*texture(syn_UserImage, p.xz*0.2);
		d += s*(s>DPrec.x?RMPrec.x:RMPrec.y);
		p = ro+rd*d;
	}	
    
	float lightIntensity = syn_MidHighLevel* syn_HighLevel*syn_Intensity*.25+0.25;

	if (s<DPrec.x)
	{
		vec2 r = getTemp(p)*(1.0+0.4*syn_HighLevel*syn_Intensity);
	
		vec3 n = nor(p, r.y); 
      	
		c.rgb = getCol(p)*(1.0+0.3*syn_HighLevel*syn_Intensity) + dot(n,rd) + ads(p,n) * lightIntensity;
	}
	else
	{
		vec3 dir = -normalize(vec3(2.,10.,0.));
		vec3 col = vec3(lightIntensity);
		c.rgb = GetSky(rd, dir, col);
	}
	c+=glow*0.01;
	return c;
}

vec3 cam(vec2 uv, vec3 ro, vec3 cu, vec3 org, float persp)
{
	uv.xy += _uvc*PI*uv*PI*Mirror;
	uv.xy *= 1.0 - Mirror*0.5*FOV*0.75;
	vec3 rorg = normalize(org-ro);
	vec3 u =  normalize(cross(cu, rorg));
	vec3 v =  normalize(cross(rorg, u));
	vec3 rd = normalize(rorg + u*uv.x + v*uv.y);
	rd.xz = _rotate(rd.xz, PI*LookXY.x+_uvc.x*PI*Flip*FOV);
	rd.xy = _rotate(rd.xy, PI*Spin);
	rd.yz = _rotate(rd.yz, PI*LookXY.y);
	return rd;
}

vec4 Image(in vec2 fragCoord )
{
	s = uScreenSize;
	g = fragCoord;
	uv = (2.*g-s)/s.y;
	uv += _uvc*PI*FOV;
	//float t = uTime*.2;
	float t = bass_time;
	float ts = sin(t)*.5+.5;
    
	float axz = -t/4; // angle XZ
	float axy = 2.6 + 0.42*ts; // angle XY // inf 3.02 // sup 2.60
	axz = mix(axz, camera_look.x*PI, manual_camera); // angle XZ
	axy = mix(axy, camera_look.y, manual_camera); // angle XY // inf 3.02 // sup 2.60
	float cd = 3.;// cam dist to scene origine

	//axy = 2.6; // on bloque la camera an haut pour mise au point
    
	float ap = 1.; // angle de perspective
	vec3 cu = vec3(0.,1.,0.); // cam up 
	vec3 org = vec3(0., 0.8, 0.); // scn org
	vec3 ro = vec3(cos(axz),sin(axy),sin(axz))*cd; // cam org
    
	vec3 rd = cam(uv, ro, cu, org, ap);
    
	vec4 c = vec4(0.,0.,0.,1.); // col
    
	c = scn(c, ro, rd);//scene
	
  c += dstepf;
    
	return c;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  fragColor = Image(fragCoord);

  fragColor = _gamma(fragColor, 0.3);

	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}