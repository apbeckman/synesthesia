//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


//this is based on https://www.shadertoy.com/view/4tcGDr

const int MAX_MARCHING_STEPS = 1000;
const float MIN_DIST = 0.0;
const float MAX_DIST = 10000.0;
const float EPSILON = 0.001;

const float size = 1000.0;


const float scale = EPSILON*10.0; //to prevent rendering artifacts

/**
 * Signed distance function describing the scene.
 * 
 * Absolute value of the return value indicates the distance to the surface.
 * Sign indicates whether the point is inside or outside the surface,
 * negative indicating inside.
 */


float planet_surface(vec3 p,float i){
    return length(sin((p*(i+1.0)/size))*size) - size;
}

float sceneSDF(vec3 p,float anim) {
    p /= scale;
    float result = 0.0;
    for(float i = 1.0; i < 81.0; i *= 3.0){
    	result += planet_surface(p*i,i)/(i*i);
    }
    //float result = sceneSDF1(p/1000.0+sceneSDF1(p/1000.0));
    return -result*scale/2.0;
}

float sceneSDF(vec3 p){
	return sceneSDF(p,1.0);
}

vec3 surface_color(vec3 p)
{
    p /= scale;
    p /= 200.0;
    return sin(vec3(sceneSDF(p/5.0,0.0),sceneSDF(p*3.0,0.0),sceneSDF(p*2.0,0.0)))/3.0+vec3(.3);
}


/**
 * Return the shortest distance from the eyepoint to the scene surface along
 * the marching direction. If no part of the surface is found between start and end,
 * return end.
 * 
 * eye: the eye point, acting as the origin of the ray
 * marchingDirection: the normalized direction to march in
 * start: the starting distance away from the eye
 * end: the max distance away from the ey to march before giving up
 */
float shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end) {
    float depth = start;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        float dist = sceneSDF(eye + depth * marchingDirection);
        if (dist < EPSILON) {
			return depth;
        }
        depth += dist;
        if (depth >= end) {
            return end;
        }
    }
    return end;
}
            

/**
 * Return the normalized direction to march in from the eye point for a single pixel.
 * 
 * fieldOfView: vertical field of view in degrees
 * size: resolution of the output image
 * fragCoord: the x,y coordinate of the pixel in the output image
 */
vec3 rayDirection(float fieldOfView, vec2 size, vec2 fragCoord) {
    vec2 xy = fragCoord - size / 2.0;
    float z = size.y / tan(radians(fieldOfView) / 2.0);
    return normalize(vec3(xy, -z));
}

/**
 * Using the gradient of the SDF, estimate the normal on the surface at point p.
 */
vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
        sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
        sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
        sceneSDF(vec3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

/**
 * Lighting contribution of a single point light source via Phong illumination.
 * 
 * The vec3 returned is the RGB color of the light's contribution.
 *
 * k_a: Ambient color
 * k_d: Diffuse color
 * k_s: Specular color
 * alpha: Shininess coefficient
 * p: position of point being lit
 * eye: the position of the camera
 * lightPos: the position of the light
 * lightIntensity: color/intensity of the light
 *
 * See https://en.wikipedia.org/wiki/Phong_reflection_model#Description
 */
vec3 phongContribForLight(vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye,
                          vec3 lightPos, vec3 lightIntensity) {
    lightPos = eye;
    vec3 N = estimateNormal(p);
    vec3 L = normalize(lightPos - p);
    vec3 V = normalize(eye - p);
    vec3 R = normalize(reflect(-L, N));
    
    float dotLN = dot(L, N);
    float dotRV = dot(R, V);
    
    if (dotLN < 0.0) {
        // Light not visible from this point on the surface
        return vec3(0.0, 0.0, 0.0);
    } 
    
    if (dotRV < 0.0) {
        // Light reflection in opposite direction as viewer, apply only diffuse
        // component
        return lightIntensity * (k_d * dotLN);
    }
    return lightIntensity * (k_d * dotLN + k_s * pow(dotRV, alpha));
}

/**
 * Lighting via Phong illumination.
 * 
 * The vec3 returned is the RGB color of that point after lighting is applied.
 * k_a: Ambient color
 * k_d: Diffuse color
 * k_s: Specular color
 * alpha: Shininess coefficient
 * p: position of point being lit
 * eye: the position of the camera
 *
 * See https://en.wikipedia.org/wiki/Phong_reflection_model#Description
 */
vec3 phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye) {
    const vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
    vec3 color = ambientLight * k_a;
    
    vec3 light1Pos = vec3(4.0,
                          2.0,
                          4.0);
    vec3 light1Intensity = vec3(0.8);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light1Pos,
                                  light1Intensity);   
    return color;
}

/**
 * Return a transform matrix that will transform a ray from view space
 * to world coordinates, given the eye point, the camera target, and an up vector.
 *
 * This assumes that the center of the camera is aligned with the negative z axis in
 * view space when calculating the ray marching direction. See rayDirection.
 */
mat3 viewMatrix(vec3 eye, vec3 center, vec3 up) {
    // Based on gluLookAt man page
    vec3 f = normalize(center - eye);
    vec3 s = normalize(cross(f, up));
    vec3 u = cross(s, f);
    return mat3(s, u, -f);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec3 viewDir = rayDirection(45.0, RENDERSIZE.xy, fragCoord);
    float speed = 10.0;
	vec3 eye = vec3(TIME*5.0,0,0);
    
    mat3 viewToWorld = -viewMatrix(eye, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
    
    vec3 worldDir = viewToWorld * viewDir;
    
    float dist = shortestDistanceToSurface(eye, worldDir, MIN_DIST, MAX_DIST);
    
    if (dist > MAX_DIST - EPSILON) {
        // Didn't hit anything
        fragColor = vec4(0.0, 0.0, 0.0, 0.0);
		return fragColor;
    }
    
    // The closest point on the surface to the eyepoint along the view ray
    vec3 p = eye + dist * worldDir;
    
    // Use the surface normal as the ambient color of the material
    vec3 K_a = surface_color((p));
    vec3 K_d = K_a;
    vec3 K_s = vec3(1.0, 1.0, 1.0);
    float shininess = 100.0;
    
    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, eye);
    
    fragColor = vec4(color, 1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}