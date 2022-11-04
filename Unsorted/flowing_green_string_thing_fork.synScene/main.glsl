

			//******** Common Code Begins ********

#define MAX_FLOAT 3.402823466e+38
#define MIN_FLOAT 1.175494351e-38
#define MAX_DOUBLE 1.7976931348623158e+308
#define MIN_DOUBLE 2.2250738585072014e-308

// Indicate to 'repeat' function that we don't wish to
#define NEVER 1000000.0

/**
 * Common vectors
 */
const vec3 ORIGIN = vec3(0,0,0);
const vec3 X = vec3(1,0,0);
const vec3 Y = vec3(0,1,0);
const vec3 Z = vec3(0,0,1);

/**
 * Common color values
 */
const vec3 BLACK = vec3(0,0,0);
const vec3 WHITE = vec3(1,1,1);
const vec3 RED   = vec3(1,0,0);
const vec3 GREEN = vec3(0,1,0);
const vec3 BLUE  = vec3(0,0,1);
const vec3 YELLOW  = vec3(1,1,0);
const vec3 CYAN    = vec3(0,1,1);
const vec3 MAGENTA = vec3(1,0,1);

/**
 * For the given 2d screen position, figure out the ray vector
 */
vec3 calculateRay(vec3 res, vec2 screenPos, 
                  vec3 eye, vec3 look_at, vec3 up) {
	vec2 screen_pos = screenPos.xy / res.xy;
    float aspect = res.y / res.x;
    screen_pos -= 0.5;
    screen_pos.y *= aspect;
    vec3 look_center = normalize(look_at - eye);
    vec3 look_right = cross(up, look_center);
    vec3 look_up = cross(look_center, look_right);
        
	vec3 newRay = normalize(look_center + screen_pos.x * look_right + screen_pos.y * look_up);
    return newRay;
}



/*
 * Signed distance functions for object primitives
 */
float sphere(vec3 where, vec3 center, float radius) {
  return length(where - center) - radius;
}

//float torus_around_x(vec3 where, float major, float minor) {
    

float round_box( vec3 where, vec3 sizes, float roundness ) {
	return length(max(abs(where)-sizes,0.0))-roundness;
}

/**
 * centred modulo
 */
float cmod(float x, float r) {
    return mod(x + 0.5 *r, r) - 0.5 *r;
}

vec3 repeat(vec3 where, vec3 repetition) {

    return mod(where, repetition);
}
vec3 repeat_x(vec3 where, float r) {

    where.x = mod(where.x, r);
    return where;
}


#define PI 3.141592653589793
vec3 radial_symmetry_xz(vec3 where, float count) {
    float ang = mod(atan(where.x, where.z) + PI, 2.0 *PI /count);
    float r = length(where.xz);
    return vec3(r *cos(ang), where.y, r * sin(ang));
}

// polynomial smooth min (k = 0.1);
float blend( float a, float b, float k )
{
    float h = max( k-abs(a-b), 0.0 );
    return min( a, b ) - h*h*0.25/k;
}


int hash(int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}


			//******** BuffA Code Begins ********


/**
 * Ray marching parameters for this scene
 */
#define MAX_STEPS  100
#define MAX_DIST   50.0
#define EPSILON    0.001
#define STEP_RATIO 0.01

/**
 * object ids
 */
#define ID_FLOOR 10
#define ID_BIZZO 1
#define ID_BIZZO2  2
#define ID_CAP   4

#define COS_30    0.8660254037844387
#define HEX_EVEN  1.7320508075688772

vec2 dist_bizzo(vec3 where) {
    float d = length(where) -1.0;
    float ratio = 2.0 + 0.1 * sin(smoothTime * 0.13);
    
    float a2 = 1.0 + 1.0 * sin(smoothTime *0.49) + 0.5 * where.y;
    float s2 = sin(a2);
    float c2 = cos(a2);
    
    float a3 = 1.0 * sin(TIME *0.29) + 0.15 * 0.5 * where.x ;
    float s3 = sin(a3);
    float c3 = cos(a3);
    
    float a = 1.0 + sin(TIME * 0.39) + 0.15 * cos(where.z);
    float s = sin(a);
    float c = cos(a);
    
	float depth = 0.0;
    
    vec3 displacement = vec3(
        -6.2 +  c,
        6.5 + 0.9 * s,
        -6.5 + 0.9
	);
    for (int  i = 0; i <6; i++) {
        
        where = vec3(
            c2 * where.x - s2 * where.z,
            where.y,
            s2 * where.x + c2 * where.z
        );
        where = vec3(
            c3 * where.x - s3 * where.y,
            s3 * where.x + c3 * where.y,
            where.z
        );

        where = -abs(where);

        where.xz = -where.zx;
        where.y = -5.0/(1.0 + pow(length(where.yz),2.0));
        where = where * ratio + displacement;
        float j = 0.5*(length(where) -1.0) / ratio;
//        d = blend(d, j, 0.2);
        // blend function body
        float k = 0.2;
    	float h = max( k - abs(d - j), 0.0 );
        float together = min(d, j);
    	d= together - h*h*0.25/k;
        float progress = clamp((j-d)/k, 0.0, 1.0);
        depth+=progress;
	
    }
    return vec2(d, depth);
}

float dist_bizzo2(vec3 where) {
    where.y += 22.0;
    

	where.y += cos(sin(TIME + length(where.xz)) * 5.0);
    where.y += 0.5 * sin(TIME*0.73 + sin(where.z) * 3.0);
    return where.y;
}

/**
 * find the closest object in the scene and return its distance and id
 */
vec3 measure(vec3 where) {
    vec3 closest = vec3(100000.0, 0.0, 0.0);
	float depth = 0.0;
    float dist_floor = where.y + 30.0;
    if (dist_floor <= closest.x) {
        closest = vec3(dist_floor, ID_FLOOR, 0.0);
    }

    vec2 dist_bizzo = dist_bizzo(where);
    if (dist_bizzo.x <= closest.x) {
        closest = vec3(dist_bizzo.x, ID_BIZZO, dist_bizzo.y);
    }

    float dist_bizzo2 = dist_bizzo2(where);
    if (dist_bizzo2 <= closest.x) {
        closest = vec3(dist_bizzo2, ID_BIZZO2, 0.0);
    }

    return closest;
}

/**
 * Figure out coloring for where we hit
 */
const vec4 floor_color0 =  vec4(0.018,0.018,0.022,0.0);
const vec4 bizzo_color0 = vec4(0.1, 0.4, 0.32, 0.2) * 0.5;
const vec4 bizzo2_color0 = vec4(0.25,0.23,0.21,0.2);

const vec4 guts_color0 = vec4(0.5,0.2,0.1,0.0);
const vec4 rod_color0 = vec4(1,0,0,1.0);
const vec4 bone_color0 = vec4(0.6,0.57,0.50,0.2);
const vec4 sky = vec4(0,0,0,0);

vec4 paint(vec3 hit, vec3 where) {

    int who = int(hit.y);
    float ambient = 0.0;
    if (who == ID_FLOOR) {
        return 0.1 * bone_color0;
    }
    if (who == ID_BIZZO) {
        return bizzo_color0 + vec4(hit.z/5.0) * vec4(1,1,0.5,1);
    }
    if (who == ID_BIZZO2) {
        return bizzo2_color0;
    }
    if (who == ID_CAP) {
        return bone_color0;
    }
    return sky;
}

// end of model stuff

vec3 calc_surface_normal(vec3 hit);
float calcSoftshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax );

/**
 * main entrypoint
 */
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec3 eye = vec3(sin(TIME/10.0)*12.0,
                    7.0,
                    cos(TIME/10.0)*12.0);
    vec3 look_at = vec3(0,0,0);
    vec3 up = Y;
    vec3 ray = calculateRay(vec3(RENDERSIZE, 1.0), fragCoord, eye, look_at, up);
    
    vec3 where = eye;
    float total_dist = 0.0;
    vec3 current;
    int who = 0;
    for(int steps = 0; steps < MAX_STEPS; steps++) {
        current = measure(where);
        float current_dist = current.x;
        if (current_dist < EPSILON) {
            who = int(current.y);
            break;
        }
        total_dist += current_dist * STEP_RATIO;
        if (total_dist > MAX_DIST) {
            break;
        }
        where = eye + total_dist * ray;
    }

    vec3 fog_color = vec3(0,0,0);
    if (who == 0){
        fragColor = vec4(fog_color, 1.0);
        return fragColor;
    }
	vec3 hit = where;
    
    vec4 the_paint = paint(current, where);
    vec3 to_light = normalize(vec3(-5,15,-1));
    float shadow = calcSoftshadow(hit, to_light, 0.0, total_dist);
    vec3 surface_normal = calc_surface_normal(hit);
    float dotty = dot(to_light, surface_normal);
    float light_amount = max(0.0, dotty);
    float light_fade = 1.0;
    float ambient = the_paint.w;
    float lighting = ambient + (1.0-ambient) * 
        (shadow*0.5 * (1.0 + light_amount * light_fade));

	vec3 coloring = light_fade *(the_paint.xyz * lighting)
        + fog_color * (1.0-light_fade);
    vec3 reflected = surface_normal * 2.0 * dotty - to_light;
    vec3 toEye = normalize(-ray);
	float specular = shadow * pow(max(0.0, dot(toEye, reflected)), 32.0);
	coloring += vec3(specular, specular, specular);
    fragColor = vec4(coloring, total_dist * 0.2);
    return fragColor;
}

#define NORMAL_DELTA 0.001

vec3 calc_surface_normal(vec3 hit) {
	return normalize(vec3(
            measure(hit+vec3(NORMAL_DELTA, 0.0, 0.0)).x - measure(hit-vec3(NORMAL_DELTA, 0.0, 0.0)).x,
            measure(hit+vec3(0.0, NORMAL_DELTA, 0.0)).x - measure(hit-vec3(0.0, NORMAL_DELTA, 0.0)).x,
            measure(hit+vec3(0.0, 0.0, NORMAL_DELTA)).x - measure(hit-vec3(0.0, 0.0, NORMAL_DELTA)).x
    ));
}


float calcSoftshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax )
{
    //vec4 fragColor = vec4(0.0);
	float res = 1.0;
    float t = mint;
    for( int i=0; i<16; i++ )
    {
		float h = measure( ro + rd*t ).x;
        res = min( res, 8.0*h/t );
        t += clamp( h, 0.02, 0.10 );
        if( res<0.005 || t>tmax ) break;
    }
    return clamp( res, 0.0, 1.0 );
	//return float(fragColor); 
 } 



/**
 * Ray marching parameters for this scene
 */
#define MAX_STEPS1  300000
#define MAX_DIST1   20.0
//#define EPSILON    0.001
//#define STEP_RATIO1 0.20
#define STEP_RATIO1 0.20

#define FRACTAL_ITERATIONS 7
/**
 * object ids
 */
#define ID_FLOOR 10
#define ID_BIZZO 1
#define ID_BIZZO2  2
#define ID_CAP   4

#define COS_30    0.8660254037844387
#define HEX_EVEN  1.7320508075688772

vec2 dist_bizzo1(vec3 where) {
    float d = length(where) -1.0;
//    float ratio = 1.75 + 0.1 * sin(TIME * 0.13);
    float ratio = 1.75 + 0.1 * sin(smoothTime * 0.13);
    
    //float a2 = 1.0 + 1.0 * sin(TIME *0.49) + 0.25 * sin(TIME*0.01)* where.y;
    float a2 = 1.0 + 1.0 * -sin(smoothTime *0.079) + 0.25 * sin(smoothTime*0.01)* where.y;

    float s2 = sin(a2);
    float c2 = cos(a2);

//    float a3 = 1.0 * sin(TIME *0.29) + 0.015 * 0.25 * sin(TIME*0.03)*where.x ;
    float a3 = 1.0 * sin(smoothTime *0.129) + 0.015 * 0.25 * -sin(smoothTime*0.03)*where.x ;
    
    float s3 = sin(a3);
    float c3 = cos(a3);
    
 //   float a = 1.0 + sin(TIME * 0.39) + 0.15 * cos(where.z);
    float a = 1.0 + sin(smoothTimeC * 0.19) + 0.15 * cos(where.z);

    float s = sin(a);
    float c = cos(a);
    
	float depth = 0.0;
    
    vec3 displacement = vec3(
        -2.2 +  c,
        2.5 + 0.9 * s,
        -1.5 + 0.9 * sin(smoothTime * 0.13)
	);
    for (int  i = 0; i <FRACTAL_ITERATIONS; i++) {
        
        where = vec3(
            c2 * where.x - s2 * where.z,
            where.y,
            s2 * where.x + c2 * where.z
        );
        where = vec3(
            c3 * where.x - s3 * where.y,
            s3 * where.x + c3 * where.y,
            where.z
        );

        where = -abs(where);

        where.xz = -where.zx;
        where.y = -5.0/(1.0 + pow(length(where.yz),2.0));
        where = where * ratio + displacement;
        float j = 0.5*(length(where) -1.0) / ratio;
//        d = blend(d, j, 0.2);
        // blend function body
        float k = 0.2;
    	float h = max( k - abs(d - j), 0.0 );
        float together = min(d, j);
    	d= together - h*h*0.25/k;
        float progress = clamp((j-d)/k, 0.0, 1.0);
        depth+=progress;
	
    }
    return vec2(d, depth);
}

float dist_bizzo3(vec3 where) {
    where.y += 22.0;
    

	where.y += cos(sin(smoothTime + length(where.xz)) * 5.0);
    where.y += 0.5 * -sin(smoothTime*0.73 + sin(where.z) * 3.0);
    return where.y;
}

/**
 * find the closest object in the scene and return its distance and id
 */
vec3 measure1(vec3 where) {
    vec3 closest = vec3(100000.0, 0.0, 0.0);
	float depth = 0.0;
    /*
    float dist_floor = where.y + 30.0;
    if (dist_floor <= closest.x) {
        closest = vec3(dist_floor, ID_FLOOR, 0.0);
    }
*/
    vec2 dist_bizzo = dist_bizzo1(where);
    if (dist_bizzo.x <= closest.x) {
        closest = vec3(dist_bizzo.x, ID_BIZZO, dist_bizzo.y);
    }
/*
    float dist_bizzo2 = dist_bizzo2(where);
    if (dist_bizzo2 <= closest.x) {
        closest = vec3(dist_bizzo2, ID_BIZZO2, 0.0);
    }
*/
    return closest;
}

/**
 * Figure out coloring for where we hit
 */
const vec4 floor_color1 =  vec4(0.018,0.018,0.022,0.0);
const vec4 bizzo_color1 = vec4(0.1, 0.4, 0.32, 0.2) * 0.5;
const vec4 bizzo2_color1 = vec4(0.25,0.23,0.21,0.2);

const vec4 guts_color1 = vec4(0.5,0.2,0.1,0.0);
const vec4 rod_color1 = vec4(1,0,0,1.0);
const vec4 bone_color1 = vec4(0.6,0.57,0.50,0.2);
const vec4 sky_top = vec4(0.8, 0.74, 0.74, 1.0);
const vec4 sky_bottom = vec4(BLACK,1.0); //vec4(0.4, 0.37, 0.37, 1.0);

vec4 paint1(vec3 hit, vec3 where) {

    int who = int(hit.y);
    float ambient = 0.0;
    if (who == ID_FLOOR) {
        return 0.1 * bone_color1;
    }
    if (who == ID_BIZZO) {
        return bizzo_color1 + vec4(hit.z/7.0) * vec4(1,1,0.5,1);
    }
    if (who == ID_BIZZO2) {
        return bizzo2_color1;
    }
    if (who == ID_CAP) {
        return bone_color1;
    }
    return sky_top;
}

// end of model stuff

vec3 calc_surface_normal1(vec3 hit);
float calcSoftshadow1( in vec3 ro, in vec3 rd, in float mint, in float tmax );

/**
 * main entrypoint
 */
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec3 eye = vec3(
        sin(smoothTime/10.0)*12.0,
        7.0,
        cos(smoothTime/10.0)*12.0
    );
    vec3 look_at = vec3(0,0,0);
    vec3 up = Y;
    vec3 ray = calculateRay(vec3(RENDERSIZE, 0.), fragCoord, eye, look_at, up);
    
    vec3 where = eye;
    float total_dist = 0.0;
    vec3 current;
    int who = 0;
    for(int steps = 0; steps < MAX_STEPS1; steps++) {
        current = measure1(where);
        float current_dist = current.x;
        if (current_dist < EPSILON) {
            who = int(current.y);
            break;
        }
        total_dist += current_dist * STEP_RATIO1;
        if (total_dist > MAX_DIST1) {
            break;
        }
        where = eye + total_dist * ray;
    }

    vec3 fog_color = vec3(0,0,0);
    if (who == 0){
        float pos = 0.5 * (ray.y + 1.0);
        fragColor = mix(sky_bottom, sky_top, pos);
        return fragColor;
    }
	vec3 hit = where;
    
    vec4 the_paint = paint1(current, where);
    vec3 to_light = normalize(vec3(-5,15,-1));
    float shadow = calcSoftshadow1(hit, to_light, 0.0, total_dist);
    vec3 surface_normal = calc_surface_normal1(hit);
    float dotty = dot(to_light, surface_normal);
    float light_amount = max(0.0, dotty);
    float light_fade = 1.0;
    float ambient = the_paint.w;
    float lighting = ambient + (1.0-ambient) * 
        (shadow*0.5 * (1.0 + light_amount * light_fade));

	vec3 coloring = light_fade *(the_paint.xyz * lighting)
        + fog_color * (1.0-light_fade);
    vec3 reflected = surface_normal * 2.0 * dotty - to_light;
    vec3 toEye = normalize(-ray);
	float specular = shadow * pow(max(0.0, dot(toEye, reflected)), 32.0);
	coloring += vec3(specular, specular, specular);
    fragColor = vec4(coloring, total_dist * 0.2);
    return fragColor;
}

#define NORMAL_DELTA 0.001

vec3 calc_surface_normal1(vec3 hit) {
	return normalize(vec3(
            measure1(hit+vec3(NORMAL_DELTA, 0.0, 0.0)).x - measure1(hit-vec3(NORMAL_DELTA, 0.0, 0.0)).x,
            measure1(hit+vec3(0.0, NORMAL_DELTA, 0.0)).x - measure1(hit-vec3(0.0, NORMAL_DELTA, 0.0)).x,
            measure1(hit+vec3(0.0, 0.0, NORMAL_DELTA)).x - measure1(hit-vec3(0.0, 0.0, NORMAL_DELTA)).x
    ));
}


float calcSoftshadow1( in vec3 ro, in vec3 rd, in float mint, in float tmax )
{
	float res = 1.0;
    float t = mint;
    for( int i=0; i<16; i++ )
    {
		float h = measure1( ro + rd*t ).x;
        res = min( res, 8.0*h/t );
        t += clamp( h, 0.02, 0.10 );
        if( res<0.005 || t>tmax ) break;
    }
    return clamp( res, 0.0, 1.0 );
	//return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}