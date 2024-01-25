

// Based on "video heightfield" "by simesgreen: https://www.shadertoy.com/view/Xss3zr#

//const int _Steps = 128;
const vec3 lightDir = vec3(0.577, 0.577, 0.577);

// transforms
vec3 rotateX(vec3 p, float a)
{
    float sa = sin(a);
    float ca = cos(a);
    vec3 r;
    r.x = p.x;
    r.y = ca*p.y - sa*p.z;
    r.z = sa*p.y + ca*p.z;
    return r;
}

vec3 rotateY(vec3 p, float a)
{
    float sa = sin(a);
    float ca = cos(a);
    vec3 r;
    r.x = ca*p.x + sa*p.z;
    r.y = p.y;
    r.z = -sa*p.x + ca*p.z;
    return r;
}
float maxheight = MAXHEIGHT *1.0 + sign(MAXHEIGHT)*basshits*0.3;
bool
intersectBox(vec3 ro, vec3 rd, vec3 boxmin, vec3 boxmax, out float tnear, out float tfar)
{
	// compute intersection of ray with all six bbox planes
	vec3 invR = 1.0 / rd;
	vec3 tbot = invR * (boxmin - ro);
	vec3 ttop = invR * (boxmax - ro);
	// re-order intersections to find smallest and largest on each axis
	vec3 tmin = min (ttop, tbot);
	vec3 tmax = max (ttop, tbot);
	// find the largest tmin and the smallest tmax
	vec2 t0 = max (tmin.xx, tmin.yz);
	tnear = max (t0.x, t0.y);
	t0 = min (tmax.xx, tmax.yz);
	tfar = min (t0.x, t0.y);
	// check for hit
	bool hit;
	if ((tnear > tfar)) 
		hit = false;
	else
		hit = true;
	return hit;
}

float luminance(sampler2D tex, vec2 uv)
{
	vec3 c = texture(tex, uv).xyz;
	return dot(c, vec3(0.33, 0.33, 0.33));
}

vec2 gradient(sampler2D tex, vec2 uv, vec2 texelSize)
{
	float h = luminance(tex, uv);
	float hx = luminance(tex, uv + texelSize*vec2(1.0, 0.0));	
	float hy = luminance(tex, uv + texelSize*vec2(0.0, 1.0));
	return vec2(hx - h, hy - h);
}

vec2 worldToTex(vec3 p)
{
	vec2 uv = p.xz*.5+.5;
	uv.y = 1.0 - uv.y;
	return uv;
}

float heightField(vec3 p)
{
	//return sin(p.x*4.0)*sin(p.z*4.0);
	//return luminance(syn_UserImage, p.xz*0.5+0.5)*2.0-1.0;
	//return luminance(syn_UserImage, worldToTex(p))*0.5;
	vec4 mediaEdges = _edgeDetectSobel(syn_UserImage, _uv);
	return luminance(syn_UserImage, worldToTex(p))*maxheight*OVERDRIVE; //FACTOR IS SCALING in Z
}
bool traceHeightField(vec3 ro, vec3 rayStep, out vec3 hitPos)
{
	vec3 p = ro;
	bool hit = false;
	float pH = 0.0;
	vec3 pP = p;
	for(int i=0; i< 256; i++) {
		
		if (i > int(SLICES)) { break; }
		
		float h = heightField(p);
		if ((p.y < h) && !hit) {
			hit = true;
			//hitPos = p;
			// interpolate based on height
            hitPos = mix(pP, p, (pH - pP.y) / ((p.y - pP.y) - (h - pH)));
		}
		pH = h;
		pP = p;
		p += rayStep;
	}
	return hit;
}

vec3 background(vec3 rd)
{
     return mix(vec3(1.0, 1.0, 1.0), vec3(0.0, 0.5, 1.0), abs(rd.y));
}

#define TWOPI 6.28318530718

vec4 renderMainImage() { 
 	vec4 out_FragColor = vec4(0.0);

    vec2 pixel = (_xy.xy / RENDERSIZE.xy)*2.0-1.0;

    // compute ray origin and direction
    float asp = 1.0; //RENDERSIZE.x / RENDERSIZE.y;
    vec3 rd = normalize(vec3(asp*pixel.x, pixel.y, -2.0));
    vec3 ro = vec3(-1.*MoveXY.x, -1.*MoveXY.y, MOVEZ);
		
	vec2 mouse = vec2(RotateXY.x,RotateXY.y);


	// rotate view
    float ax = (mouse.y + .25) * TWOPI;

    rd = rotateX(rd, ax);
    ro = rotateX(ro, ax);
		
    float ay = (mouse.x + .5) * TWOPI;
    rd = rotateY(rd, ay);
    ro = rotateY(ro, ay);
	
	// intersect with bounding box
    bool hit;	
	vec3 boxMin = vec3(-1.0, -.001, -1.0);
	vec3 boxMax = vec3(1.0, maxheight, 1.0); //FACTOR IS SCALING Z - 2nd term in vec3
	float tnear, tfar;
	hit = intersectBox(ro, rd, boxMin*OVERDRIVE, boxMax*OVERDRIVE, tnear, tfar);

	tnear -= 0.000001;
	vec3 pnear = ro + rd*tnear;
    vec3 pfar = ro + rd*tfar;
	
    float stepSize = length(pfar - pnear) / float(SLICES);
	
    vec4 col = vec4(0,0,0,0);
    if(hit)
    {
    	// intersect with heightfield
		ro = pnear;
		vec3 hitPos;
		hit = traceHeightField(ro, rd*stepSize, hitPos);
		if (hit) {
			
			vec2 uv = worldToTex(hitPos);
			col = texture(syn_UserImage, uv);
				col += highhits*0.5*flashy;

		}
    }

    out_FragColor=vec4(col);

return out_FragColor; 
 } 
vec4 renderMain(){

	if(PASSINDEX == 0){
		return renderMainImage();
	}
}
