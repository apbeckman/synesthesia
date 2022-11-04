

			//******** BuffA Code Begins ********

#define MARCHING_STEP 128

/////
// SDF Operation function
/////

vec3 opRep( in vec3 p, in vec3 c)
{
    vec3 q = mod(p+0.5*c,c)-0.5*c;
    return q;
}

/////
// Scene and primitive SDF function
/////

float sphereSDF(vec3 samplePoint) {
    return length(samplePoint) - 1.0;
}



float sdPlane( vec3 p )
{
    return p.y;
}

#define Scale 2.
#define iteration 15
#define Power (7.+sin(smoothTime/7.)*5.)
#define Bailout 5.
vec3 CSize = vec3(1.3);

float DE(vec3 p) {
    p = p.xzy;
	float scale = 1.;
	for( int i=0; i < iteration;i++ )
	{
		p = 2.0*clamp(p, -CSize, CSize) - p;
		//float r2 = dot(p,p);
        float r2 = dot(p,p+sin(p.x*.1)); //Alternate fractal
		float k = max((2.)/(r2), .027);
		p     *= k;
		scale *= k;
	}
	float l = length(p.xy);
	float rxy = l - 4.0;
	float n = l * p.z;
	rxy = max(rxy, -(n) / 4.);
	return (rxy) / abs(scale);
}

float sceneSDF(vec3 samplePoint) {

  
    float res = DE(samplePoint);
    //res += sdPlane(-0.5, vec4(0.,1.,0.,1.));
    return res;
   
}


/////
// Ray function
/////

vec3 getCameraRayDir(vec2 uv, vec3 camPos, vec3 camTarget)
{
    // Calculate camera's "orthonormal basis", i.e. its transform matrix components
    vec3 camForward = normalize(camTarget - camPos);
    vec3 camRight = normalize(cross(vec3(0.0, 1.0, 0.0), camForward));
    vec3 camUp = normalize(cross(camForward, camRight));
     
    float fPersp = 2.0;
    vec3 vDir = normalize(uv.x * camRight + uv.y * camUp + camForward * fPersp);
 
    return vDir;
}

vec3 rayDir(float fov, vec2 size, vec2 fragCoord)
{
    vec2 xy = fragCoord - size/2.0;
    float z = size.y * 0.5 / tan(radians(fov)/ 2.0);
    return normalize(vec3(xy,-z));
}

vec2 normalizeScreenCoords(vec2 screenCoord)
{
    vec2 result = 2.0 * (screenCoord/RENDERSIZE.xy - 0.5);
    result.x *= RENDERSIZE.x/RENDERSIZE.y;
    return result;
}

/////
// Marching function
/////

float march(vec3 pos, vec3 direction, float start, float end, inout int i)
{
    float depth = start;
    for(i = 0; i < MARCHING_STEP; i++)
    {
        float dist =  sceneSDF(pos + direction * depth);
        if(dist < 0.005f)
        {
            //return depth;
            break;
        }
        depth += dist;
        if(depth >= end)
            return end;
    }
    return depth;
}


/////
// Main function
/////

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    
    vec2 uv = normalizeScreenCoords(fragCoord);
    vec3 pos = vec3(smoothTime * 10.,52.0,52.);
    vec3 at = vec3(pos.x + 0.1, pos.y, pos.z);
    
    int i = 0;
    
    vec3 dir = getCameraRayDir(uv, pos, at);
    
    float dist = march(pos, dir, 0.f,5000.f, i);
    vec3 col = vec3(dist);
    
    if((dist - 5000.f) > 0.001f)
    {
        col = vec3(0.0529, 0.0808, 0.0922);
    }
    else
    {
        col = vec3(dist*0.4); 
        col = vec3(0.75 + sin(smoothTimeB/10.), 0.615, 0.053 + cos(smoothTimeB/10.)) * float(i)/float(MARCHING_STEP);
    }
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec2 computeUV( vec2 uv, float k, float kcube ){
    
    vec2 t = uv - .5;
    float r2 = t.x * t.x + t.y * t.y;
	float f = 0.;
    
    if( kcube == 0.0){
        f = 1. + r2 * k;
    }else{
        f = 1. + r2 * ( k + kcube * sqrt( r2 ) );
    }
    
    vec2 nUv = f * t + .5;
    nUv.y = 1. - nUv.y;
 
    return nUv;
    
}



vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    float lensType = 2.0;
	vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    float k = 0.2 * cos( lensType);
    float kcube = .02 * cos( lensType );
    
    float offset = .1 * cos( lensType * .5 );
    
    float red = texture( BuffA, computeUV( uv, k + offset, kcube ) ).r; 
    float green = texture( BuffA, computeUV( uv, k, kcube ) ).g; 
    float blue = texture( BuffA, computeUV( uv, k - offset, kcube ) ).b; 
    
    fragColor = vec4( red, green,blue, 1. );

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