

// Created by Stephane Cuillerdier - @Aiekick/2021
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Tuned via Noodlesplate (https://github.com/aiekick/NoodlesPlate)

const vec3 uColor0 = vec3(0.2,0.8,0.3);
const vec3 uColor1 = vec3(0.8,0.2,0.37);

mat3 rotx(float a){return mat3(1.,0.,0.,0.,cos(a),-sin(a),0.,sin(a),cos(a));}
mat3 roty(float a){return mat3(cos(a),0.,sin(a),0.,1.,0.,-sin(a),0.,cos(a));}
mat3 rotz(float a){return mat3(cos(a),-sin(a),0.,sin(a),cos(a),0.,0.,0.,1.);}

mat3 m1;
mat3 m2;

vec3 path(vec3 p)
{
	p.xy += cos(TIME * 0.1);
	p *= rotz(p.z * 0.2);
    p += sin(p.zxy * 0.5) * 0.5;
	p *= rotz(p.z * 0.2);
    p = sin(p.zxy * 0.1) * 2.;
    return p;
}

float pattern(vec3 p)
{
	p = abs(fract(p*.3) - 0.5);
	return length(max(abs(p.x), abs(p.y)) - p.z);
}

float sdf(vec3 p)
{
	p += path(p);
	p *= rotz(p.z * 0.045);
	vec2 q = mod(p.xz, 3.0) - 1.5;
	vec3 qm1 = p*m1, qm2 = p*m2;
	float d0 = min(pattern(qm1), pattern(qm2));
    float d1 = min(pattern(qm1*3.), pattern(qm2*2.));
   	float dist0 = (1.-clamp(d0,0.,1.));
	float dist1 = (1.-clamp(d1,0.,1.))*d0;
	float sp = 1. - cos(p.z * 1.5) * 0.1 - abs(p.y) + dist0*0.75 + dist1*2.25;
	p.xz = mod(p.xz, 1.5) - 0.75;
	float ct = (cos(TIME * 0.5) * 0.5 + 0.5) * 2.0;
	float st = (sin(TIME * 0.5) * 0.5 + 0.5) * 2.0;
	float bo = length(max(abs(p) - vec3(ct,1000.0,st), 0.0)) - 0.1;
	
	return max(bo,sp);
}

float mat(vec3 p)
{
	p += path(p);
	p *= rotz(p.z * 0.045);
	vec2 q = mod(p.xz, 3.0) - 1.5;
	vec3 qm1 = p*m1, qm2 = p*m2;
	float d0 = min(pattern(qm1), pattern(qm2));
    float d1 = min(pattern(qm1*3.), pattern(qm2*2.));
   	float dist0 = (1.-clamp(d0,0.,1.));
	float dist1 = (1.-clamp(d1,0.,1.))*d0;
	float sp = 1. - cos(p.z * 1.5) * 0.1 - abs(p.y) + dist0*0.75 + dist1*2.25;
	p.xz = mod(p.xz, 1.5) - 0.75;
	float ct = (cos(TIME * 0.5) * 0.5 + 0.5) * 2.0;
	float st = (sin(TIME * 0.5) * 0.5 + 0.5) * 2.0;
	float bo = length(max(abs(p) - vec3(ct,1000.0,st), 0.0)) - 0.1;
	
	if (bo > sp)
	{
		return 1.0;
	}
	return 0.0;
}

vec3 nor(vec3 p, float epsilon)
{
	vec3 eps = vec3(epsilon, 0., 0.);
	vec3 nor = vec3(
		sdf(p + eps.xyy) - sdf(p - eps.xyy),
		sdf(p + eps.yxy) - sdf(p - eps.yxy),
		sdf(p + eps.yyx) - sdf(p - eps.yyx));
	return normalize(nor);
}

float getSha( in vec3 ro, in vec3 rd, in float hn) // iq code
{
    float res = 1.0;
    float t = 0.0005;
	float h = 1.0;
    for( int i=0; i<40; i++ )
    {
        h = sdf(ro + rd*t);
        res = min( res, hn*h/t );
		t += clamp( h, 0.02, 2.0 );
    }
    return clamp(res,0.0,1.0);
}


float getAO( in vec3 p, in vec3 nor ) // iq code
{
	float occ = 0.0;
    float sca = 1.0;
    for( int i=0; i<5; i++ )
    {
        float hr = 0.01 + 0.12 * float(i)/4.0;
        vec3 aopos =  nor * hr + p;
        float dd = sdf( aopos );
        occ += -(dd-hr)*sca;
        sca *= 0.95;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );    
}

vec3 tex3D( sampler2D tex, in vec3 p, in vec3 n ) // from shane
{
    n = max((abs(n) - .2)*7., .001);
    n /= (n.x + n.y + n.z );  
    p = (texture(tex, p.yz)*n.x + texture(tex, p.zx)*n.y + texture(tex, p.xy)*n.z).xyz;
    return p*p;
}


vec3 doBumpMap( sampler2D tx, in vec3 p, in vec3 n, float bf) // from shane
{
    const vec2 e = vec2(0.001, 0);
    mat3 m = mat3( tex3D(tx, p - e.xyy, n), tex3D(tx, p - e.yxy, n), tex3D(tx, p - e.yyx, n));
    vec3 g = vec3(0.299, 0.587, 0.114)*m;
    g = (g - dot(tex3D(tx,  p , n), vec3(0.299, 0.587, 0.114)) )/e.x; g -= n*dot(n, g);
    return normalize( n + g*bf );
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	fragColor = vec4(0);
	
	vec2 uv = (fragCoord * 2.0 - RENDERSIZE.xy)/RENDERSIZE.y;
	
	float ca = 3.0;
    float ce = 8.0;
    float cd = 20.0;
	float maxd = 1000.0;
	
	ca = _mouse.x / RENDERSIZE.x * 6.28318;
	ce = _mouse.y / RENDERSIZE.y * 30.0;
	
    vec3 ro = vec3(0,0, TIME * 6.); 
	ro -= path(ro);
	vec3 cv = ro + vec3(0,0,4); 
	cv -= path(cv);
	vec3 lp = ro + vec3(0,0,3);
    lp -= path(lp);
	
	vec3 cu = normalize(vec3(0,1,0));
  	vec3 z = normalize(cv-ro);
    vec3 x = normalize(cross(cu,z));
  	vec3 y = cross(z,x);
  	vec3 rd = normalize(uv.x * x + uv.y * y + z);
	
	mat3 mx = rotx(-7.);
	mat3 my = roty(-5.);
	mat3 mz = rotz(-3.);
	
    m1 = mx * my * mz;
    m2 = m1 * m1;
	
	float d = 0.0, s = 1.0;
	for (float i = 0.0; i < 500.0; i++)
	{
		if (log(d/1e7) > 0.0 || d > maxd) break;
		s = sdf(ro + rd * d);
		d += s * 1.0;
	}
	
	float fog = 1.0-exp( -0.001*d*d );
		
	if (d < maxd)
	{
		vec3 p = ro + rd * d;
		vec3 n = nor(p, 0.01);
		
		vec3 ld = normalize(lp-p); 										
		
		float mt = mat(p);
		
		vec4 col = vec4(0);									
			
		if (mt < 0.5) // top
		{
			n = doBumpMap(image7, -p*0.5, n, 0.018);
			vec3 refl = reflect(rd,n);
			float diff = clamp( dot( n, ld ), 0.0, 1.0 );
			float fre = pow( clamp( 1. + dot(n,rd),0.0,1.0), 16. );
			float spe = pow(clamp( dot( refl, ld ), 0.0, 1.0 ),25.);
			float sha = 0.5 + 0.5 * getSha(p, n, 2.0);
			float ao = getAO(p, n);
			fragColor = vec4(
				(diff + fre) * 0.8 + diff * ao,
				(diff + fre) * uColor1 + spe * sha
			) * 1.0;
		}
		else // wall
		{
			n = doBumpMap(image3, -p*0.5, n, 0.1);
			vec3 refl = reflect(rd,n);
			float diff = clamp( dot( n, ld ), 0.0, 1.0 );
			float spe = pow(clamp( dot( refl, ld ), 0.0, 1.0 ),25.);
			float sha = 0.5 + 0.5 * getSha(p, n, 2.0);
			float ao = getAO(p, n);
    
			fragColor = vec4(uColor0 * ao + spe * sha, 1);
		}

		fragColor = fragColor.zyww + fragColor.x*0.1;
	}
	
	fragColor = mix( fragColor, vec4(0), fog);
       
	fragColor = mix(fragColor, fragColor.grba, sin(fog*3.));
   	fragColor = sqrt(fragColor*fragColor*fragColor*1.5);
	
	// vigneting from iq
    vec2 q = fragCoord/RENDERSIZE.xy;
    fragColor.rgb *= 0.5 + 0.5*pow( 16.0*q.x*q.y*(1.0-q.x)*(1.0-q.y), 1.0 );

	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}