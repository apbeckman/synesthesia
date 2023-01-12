vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


//**************************************************************************************************
// There´s something under the surface
// Just testing enhanced sphere tracing: lgdv.cs.fau.de/get/2234
// Implemented by chronos in www.shadertoy.com/view/4lj3zK

// Soundtrack "System malfunction" done with 64klang

//**************************************************************************************************

float sminPoly( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}


float vmax(vec3 v) {
	return max(max(v.x, v.y), v.z);
}

float fBox(vec3 p, vec3 b) {
	vec3 d = abs(p) - b;
	return length(max(d, vec3(0))) + vmax(min(d, vec3(0)));
}

void pR(inout vec2 p, float a) {
	p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

vec2 pModGrid2(inout vec2 p, vec2 size) {
    p+=its;
	vec2 c = floor((p + size*0.5)/size);
	p = mod(p + size*01.5, size) - size*0.5;
	p *= mod(c,vec2(2))*2. - vec2(1);
	p -= size/2.;
	if (p.x > p.y) p.xy = p.yx;
	return floor(c/2.);
}

vec2 ku=vec2(0.);

float sdPlane(vec3 p) 
{
return p.y+(0.002*sin(p.x*110.))+(0.002*sin(p.z*112.))-0.3-0.1*sin(smoothTime*0.5+2.*p.z);
}

float fField2(vec3 p) 
{
	vec2 q = pModGrid2(p.xz,vec2(0.8));
	pR(p.xz,smoothTimeC*0.125);
	vec2 q2 = pModGrid2(p.xz,vec2(0.7));
	pR(p.xy,smoothTimeC*0.0755);
	float box = fBox(p-vec3(0),vec3(0.6))-0.08;
	ku=q; return box;
}

float map(vec3 p)
{
    return sminPoly(fField2(p),sdPlane(p),-0.01);
}

float shadowsoft( vec3 ro, vec3 rd, float mint, float maxt, float k )
{
	float t = mint;
	float res = 1.0;
    for ( int i = 0; i < 28; ++i )
    {
        float h = map( ro + rd * t );
        if ( h < 0.001 ) return 0.0;
		res = min( res, k * h / t );
        t += h;
		if ( t > maxt )
			break;
    }
    return res;
}


vec3 calcNormal(vec3 pos)
{
    float eps=0.0001;
	float d=map(pos);
	return normalize(vec3(map(pos+vec3(eps,0,0))-d,map(pos+vec3(0,eps,0))-d,map(pos+vec3(0,0,eps))-d));
}

float castRay(in vec3 ray_origin, in vec3 ray_direction, in bool inside) {
    float relaxation = 1.5;// range [1.0, 2.0]
    float distance_min =  0.01;
    float distance_max = 12.0;
    float precis   = 0.0005;
    float distance = distance_min;
    float previous_radius = 0.0;
    float stepLength = 0.0;
    float function_sign = 1.0;
    if(map(ray_origin) < 0.0) function_sign = -1.0;
    
	for(int i = 0; i < 60; i++ ) {
        float result; 
        if (inside==false) result = map(ray_origin + ray_direction * distance);
        else result = -map(ray_origin + ray_direction * distance);
        float signed_radius = function_sign * result;
        float radius = abs(signed_radius);
        
        bool sorFail = relaxation > 1.0 && (radius + previous_radius) < stepLength;
        if(sorFail) {
        	stepLength -= relaxation * stepLength; // revert last step
            relaxation = 1.0;
        } else {
        	stepLength = signed_radius * relaxation;   
        }
        previous_radius = radius;
        if(!sorFail && radius < precis || distance > distance_max ) break;
        distance += stepLength;
    }
    return  distance;

}


//***************************************************************************************************
// main
//***************************************************************************************************

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv,p;
	uv.xy = gl_FragCoord.xy /RENDERSIZE.xy; p = uv * 2.0 - 1.0;
   
//  cam    
	p.x *= RENDERSIZE.x /RENDERSIZE.y;
	float theta = sin(smoothTime*0.05) * 1.28;
    float x = 3.0 * cos(theta); 
    float z = 3.0 * sin(theta); 
	vec3 ro = vec3(x*2.2, 5.0+2.*sin((smoothTime+15.)*0.125), z*1.14);		

	vec3 cw = normalize(vec3(0., 0.25, 0.) - ro);
    vec3 cp = vec3(0.0, 1.0, 0.0);
    vec3 cu = normalize(cross(cw, cp));
    vec3 cv = normalize(cross(cu, cw));
	vec3 rd = normalize(p.x * cu + p.y * cv + 7.5 * cw);

// 	render:
    vec3 col= vec3(0.0);
    float t = castRay(ro,rd,false);
	vec3 pos = ro + rd *t;
	vec3 nor = calcNormal(pos);
	vec3 ligvec = vec3(-0.5, 0.5, 0.5);
	vec3 lig = normalize(ligvec);	
	float dif = clamp(dot(lig, nor), 0.0, 1.0);
	float spec = pow(clamp(dot(reflect(rd, nor), lig), 0.0, 1.0), 16.0);
	col = vec3(0.2*dif+1.*spec);
    
    float sh = shadowsoft(pos,lig,0.01,0.2,1.2); 
    col *= clamp(sh, 0.0, 1.0);


//  refraction
	vec3 rd2 = refract(rd,nor,0.77);  
    float t2 = castRay(pos,rd2,true);
	vec3 pos2 = pos + rd2* t2;
    vec3 nor2 = calcNormal(pos2);
	float dif2 = clamp(dot(lig, nor2), 0.0, 1.0);
    col.r += (1.0-t2*0.25)+(1.0-t*0.15)*dif2;
	rd2 = refract(rd,nor,0.83);
    t2 = castRay(pos,rd2,true);
	pos2 = pos + rd2* t2;
    nor2 = calcNormal(pos2);
	dif2 = clamp(dot(lig, nor2), 0.0, 1.0);
    col.b += (1.0-t2*0.25)+(1.0-t*0.15)*dif2;
	rd2 = refract(rd,nor,0.8);
    t2 = castRay(pos,rd2,true);
	pos2 = pos + rd2* t2;
    nor2 = calcNormal(pos2);
	dif2 = clamp(dot(lig, nor2), 0.0, 1.0);
    float spec2 = pow(clamp(dot(reflect(rd2, nor2), lig), 0.0, 1.0), 16.0);
    col.g += (1.0-t2*0.25)+(1.0-t*0.15)*dif2;
    col +=0.3*(spec2+highhits*0.5);
    col *= clamp(sh, 0.8, 1.0);

	vec3 ro3 = pos2+rd; 
	vec3 rd3 = rd2+0.05;
    float t3 = castRay(ro3,rd3,false);
	vec3 pos3 = ro3 + rd3* t3;
    vec3 nor3 = calcNormal(pos3);
	float dif3 = clamp(dot(lig, -nor3), 0.0, 1.0);
    col-= 0.2*(1.-dif3);
    col = mix(col, vec3(0.4,0.5,0.5), ku.y*0.2*t2);
    fragColor = vec4(col, 1.0);

	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}