vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


// ##############################
// BEGIN	IQ methods
// ##############################
float length8( vec2 p )
{
	p = p*p; p = p*p; p = p*p;
	return pow( p.x + p.y, 1.0/8.0 );
}
float sdPlane( vec3 p, vec4 n )
{
	// n must be normalized
	return dot(p,n.xyz) + n.w;
}
float sdTorus88( vec3 p, vec2 t )
{
	vec2 q = vec2(length8(p.xy)-t.x,p.z);
	return length8(q)-t.y;
}

// Union
float opU( float d1, float d2 )
{
	return min(d1,d2);
}
// Substraction
float opS( float d1, float d2 )
{
	return max(-d1,d2);
}
// Repetition
vec3 opRep( vec3 p, vec3 c )
{
	return mod(p,c)-0.5*c;
}
// ##############################
// END		IQ methods
// ##############################

// ##############################
// BEGIN	Camera helpers
// ##############################
float iCamPosX = 0.0;
float iCamPosY = 0.0;
float iCamPosZ = 0.0;
float iCamRotX = 0.0;
float iCamRotY = 0.0;
float iCamRotZ = 0.0;

vec3 calcCameraPos()
{
	return vec3(iCamPosX, iCamPosY, iCamPosZ);
}
void rotateAxis(inout vec2 p, float a)
{
	p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}
vec3 calcCameraRayDir(float fov, vec2 fragCoord, vec2 resolution)
{
	float tanFov = tan(fov / 2.0 * PI / 180.0) / resolution.x;
	vec2 p = tanFov * (fragCoord * 2.0 - resolution.xy);
	vec3 rayDir = normalize(vec3(p.x, p.y, 1.0));
	rotateAxis(rayDir.yz, iCamRotX);
	rotateAxis(rayDir.xz, iCamRotY);
	rotateAxis(rayDir.xy, iCamRotZ);
	return rayDir;
}
// ##############################
// END		Camera helpers
// ##############################

const vec3 repSpacing = vec3( 50.0, 50.0, 1.5 );

//	Calculates distance to nearest object given a point
float distFunc( vec3 point )
{
	vec3 point2 = opRep( point, repSpacing );
	//point2 += vec3(cos(smoothTime*0.25)*0.5+0.35, sin(smoothTime*0.25)*0.5+0.35, 0.);
	rotateAxis( point2.xy, (smoothTimeC*0.5+floor(point.z/repSpacing.z))/4.0 );

	//rotateAxis( point2.xy, (-smoothTimeC*0.25+floor(point.z/repSpacing.z))/2.0 );
	return sdTorus88( point2-vec3( 0, 0, 0 ), vec2( 2.0, 0.5 ) );
}


vec3 getNormal( in vec3 pos )
{
	// IQ
	vec2 e = vec2( 1.0,-1.0 ) * 0.001;
	return normalize( e.xyy*distFunc( pos + e.xyy ) +
					  e.yyx*distFunc( pos + e.yyx ) +
					  e.yxy*distFunc( pos + e.yxy ) +
					  e.xxx*distFunc( pos + e.xxx ) );
}

bool isEdge(const vec3 point)
{
	float d = 0.0125;
	//get points a little bit to each side of the point
	vec3 right = point + vec3(d, 0.0, 0.0);
	vec3 left = point + vec3(-d, 0.0, 0.0);
	vec3 up = point + vec3(0.0, d, 0.0);
	vec3 down = point + vec3(0.0, -d, 0.0);
	vec3 behind = point + vec3(0.0, 0.0, d);
	vec3 before = point + vec3(0.0, 0.0, -d);

	vec3 normRight = getNormal(right);
	vec3 normLeft = getNormal(left);
	vec3 normUp = getNormal(up);
	vec3 normDown = getNormal(down);
	vec3 normBehind = getNormal(behind);
	vec3 normBefore = getNormal(before);

	vec3 normal = getNormal(point);

	const float limit = 0.99;

	// float gradient1 = abs(dot(normal, normRight - normLeft));
	// float gradient2 = abs(dot(normal, normUp - normDown));
	// float gradient3 = abs(dot(normal, normBehind - normBefore));

	if(abs(dot(normal, normRight)) < limit) {
		return true;
	}
	if(abs(dot(normal, normLeft)) < limit) {
		return true;
	}
	if(abs(dot(normal, normUp)) < limit) {
		return true;
	}
	if(abs(dot(normal, normDown)) < limit) {
		return true;
	}
	if(abs(dot(normal, normBehind)) < limit) {
		return true;
	}
	if(abs(dot(normal, normBefore)) < limit) {
		return true;
	}

	return false;
}


//vec3 lightColor = clamp(vec3( afFrequencies[0]+afFrequencies[1], afFrequencies[2]+afFrequencies[3], afFrequencies[4]+afFrequencies[5] ), 0.0, 1.0);


const float lightAttenuation = 0.02;
vec3 getShadedColor( vec3 hitPosition, vec3 normal, vec3 cameraPosition )
{
	//	light relative to camera position
	vec3 lightPosition = vec3(sin(smoothTimeB*0.125), 1.0, cos(smoothTimeB*0.125));
	lightPosition += cameraPosition;

	//	Specular highlight factor
	float materialShininess = 16.0;
	vec3 materialSpecularColor = vec3( 1.0 );

	//	Output color
	vec3 outputColor = vec3( 0.0 );

	//	Calculate eye vector and its reflection
	//vec3 ev = normalize( hitPosition - cameraPosition );
	//vec3 ref_ev = reflect( ev, normal );
	vec3 surfaceToLight = normalize(lightPosition - hitPosition);
	vec3 surfaceToCamera = normalize(cameraPosition - hitPosition);

	//	surface color
	vec3 surfaceColor = vec3( 1.0 );

	//	edge detection
	 if(isEdge(hitPosition))
	 {
	 	surfaceColor = vec3( 0.3 );
	 }

	//	ambient component
    vec3 lightColor = vec3(abs(0.5+cos(smoothTimeB*0.125+10.)*0.125+0.9125), abs(sin(smoothTimeB*0.125)*0.125+0.925), 0.2525+abs(sin(smoothTimeB*0.125)*0.125+0.925))*0.65;
	vec3 ambientColor = surfaceColor * lightColor * 0.0; // ambient factor

	//	diffuse component
	float diffuseCoefficient = max(0.0, dot(normal, surfaceToLight));
	vec3 diffuseColor = diffuseCoefficient * surfaceColor * lightColor;

	//	specular component
	float specularCoefficient = 0.0;
	if(diffuseCoefficient > 0.0) {
		specularCoefficient = pow(max(0.0, dot(surfaceToCamera, reflect(-surfaceToLight, normal))), materialShininess);
	}
	vec3 specularColor = specularCoefficient * materialSpecularColor * lightColor;

	//	light attenuation (falloff based on distance, fog)
	float distanceToLight = length(lightPosition - hitPosition);
	float attenuation = 1.0 / (1.0 + lightAttenuation * pow(distanceToLight, 2.));

	outputColor = ambientColor + attenuation*(diffuseColor + specularColor);

	//	gamma correction
	//vec3 gamma = vec3(1.0/2.2);
	//outputColor = vec3(pow(outputColor, gamma));

	//	return shading result
	return outputColor;
}

const float epsilon = 0.0001;
const int maxSteps = 256;
const float maxT = 125.0;
float trace(vec3 ro, vec3 rd, out vec3 point, out bool objectHit)
{
	float t = 0.0;
	point = ro;

	for(int steps = 0; steps < maxSteps; ++steps)
	{
		//check how far the point is from the nearest surface
		float dist = distFunc(point);
		//if we are very close
		if(epsilon > dist)
		{
			objectHit = true;
			break;
		}
		//not so close -> we can step at least dist without hitting anything
		t += dist;
		// return immediately if maximum t is reached
		if(t > maxT)
		{
			objectHit = false;
			return maxT;
		}
		//calculate new point
		point = ro + t * rd;
	}

	return t;
}

const int reflectionBounces = 1;
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	//	Set up Camera
	vec3 camP = calcCameraPos(); // Camera position
	//camP += vec3( repSpacing.x/2.0+sin( -smoothTime*0.125 )*0.25, repSpacing.y/2.0+cos( -smoothTime*0.125 )*0.25, 2.5*smoothTime );
	camP += vec3( repSpacing.x/2.0+sin( smoothTime*0.25 )*(0.5), repSpacing.y/2.0+cos( smoothTime*0.25 )*(0.5), 2.0*smoothTime );

	vec3 camDir = calcCameraRayDir( 90.0, fragCoord.xy, RENDERSIZE.xy ); // Camera view direction

	//	Set up ray
	vec3 point;		// Set in trace()
	bool objectHit;	// Set in trace()

	//	Initialize color
	vec3 color = vec3(0.0);

	float t = trace(camP, camDir, point, objectHit);
	if(objectHit)
	{
		//	Lighting calculations
		vec3 normal = getNormal(point);
		color = getShadedColor( point, normal, camP );

		//	Reflections
		for(int i = 0; i < reflectionBounces; i++)
		{
			vec3 pointRef;	// Set in trace()
			camDir = reflect(camDir, normal);
			trace(point + camDir*0.001, camDir, pointRef, objectHit);
			if(objectHit)
			{
				// Get color of reflection
				color += 0.3 * getShadedColor( pointRef, getNormal(pointRef), point );
			}
			point = pointRef;
		}
	}

	//	fog
	//vec3 fogColor = vec3( 0.1*(1.0-sin(smoothTimeB)*0.25), 0.1+0.1*(1.0-sin(smoothTimeB)*0.25), 0.1+0.1*(1.0-sin(smoothTimeB)*0.25) );
	vec3 fogColor = vec3(abs(0.45+cos(smoothTimeB*0.125+10.)*0.125+0.9125), abs(sin(smoothTimeB*0.125)*0.125+0.925), 0.525)*1.5*(1.0+pow(normalize(highhits), 2.0)*0.125);
	float FogDensity = 0.025*(1.0+pow(syn_HighLevel*0.875+syn_Hits*0.125, 2.0)*0.5);
	float fogFactor = 1.0 /exp(t * FogDensity);
	fogFactor = clamp( fogFactor, 0.0, 1.0 );
	color = mix(fogColor, color, fogFactor);

	fragColor = vec4(color, 0.0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}